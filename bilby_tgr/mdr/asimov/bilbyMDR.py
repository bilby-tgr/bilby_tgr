"""BilbyMDR Pipeline specification."""

import os
import re
import subprocess

import time

from asimov import config
from asimov.pipeline import PipelineException, PipelineLogger

from asimov.pipelines import Bilby
from asimov.pipeline import PostPipeline
from asimov import utils

import htcondor


class BilbyMDR(Bilby):
    """
    The BilbyMDR Pipeline.

    Parameters
    ----------
    production : :class:`asimov.Production`
       The production object.
    category : str, optional
        The category of the job.
        Defaults to "C01_offline".
    """
    config_template = os.path.join(os.path.dirname(__file__), "bilbymdr.ini")
    name = "BilbyMDR"
    STATUS = {"wait", "stuck", "stopped", "running", "finished"}

    def __init__(self, production, category=None):
        super(Bilby, self).__init__(production, category)
        self.logger.info("Using the bilbyMDR pipeline")

        if not production.pipeline.lower() == "bilbymdr":
            raise PipelineException

    def build_dag(self, psds=None, user=None, clobber_psd=False, dryrun=False):
        """
        Construct a DAG file in order to submit a production to the
        condor scheduler using bilby_pipe.

        Parameters
        ----------
        production : str
           The production name.
        psds : dict, optional
           The PSDs which should be used for this DAG. If no PSDs are
           provided the PSD files specified in the ini file will be used
           instead.
        user : str
           The user accounting tag which should be used to run the job.
        dryrun: bool
           If set to true the commands will not be run, but will be printed to standard output. Defaults to False.

        Raises
        ------
        PipelineException
           Raised if the construction of the DAG fails.
        """

        cwd = os.getcwd()

        self.logger.info(f"Working in {cwd}")

        if self.production.event.repository:
            ini = self.production.event.repository.find_prods(
                self.production.name, self.category
            )[0]
            ini = os.path.join(cwd, ini)
        else:
            ini = f"{self.production.name}.ini"

        if self.production.rundir:
            rundir = self.production.rundir
        else:
            rundir = os.path.join(
                os.path.expanduser("~"),
                self.production.event.name,
                self.production.name,
            )
            self.production.rundir = rundir

        if "job label" in self.production.meta:
            job_label = self.production.meta["job label"]
        else:
            job_label = self.production.name

        command = [
            os.path.join(config.get("pipelines", "environment"), "bin", "bilby_pipe"),
            ini,
            "--label",
            job_label,
            "--outdir",
            f"{os.path.abspath(self.production.rundir)}",
            "--accounting",
            f"{self.production.meta['scheduler']['accounting group']}",
        ]

        if dryrun:
            print(" ".join(command))
        else:
            self.logger.info(" ".join(command))
            pipe = subprocess.Popen(
                command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
            )
            out, err = pipe.communicate()
            self.logger.info(out)

            if err or "DAG generation complete, to submit jobs" not in str(out):
                self.production.status = "stuck"
                self.logger.error(err)
                raise PipelineException(
                    f"DAG file could not be created.\n{command}\n{out}\n\n{err}",
                    production=self.production.name,
                )
            else:
                time.sleep(10)
                return PipelineLogger(message=out, production=self.production.name)

    def submit_dag(self, dryrun=False):
        """
        Submit a DAG file to the condor cluster.

        Parameters
        ----------
        dryrun : bool
           If set to true the DAG will not be submitted,
           but all commands will be printed to standard
           output instead. Defaults to False.

        Returns
        -------
        int
           The cluster ID assigned to the running DAG file.
        PipelineLogger
           The pipeline logger message.

        Raises
        ------
        PipelineException
           This will be raised if the pipeline fails to submit the job.

        Notes
        -----
        This overloads the default submission routine, as bilby seems to store
        its DAG files in a different location
        """

        cwd = os.getcwd()
        self.logger.info(f"Working in {cwd}")

        self.before_submit()

        try:
            # to do: Check that this is the correct name of the output DAG file for billby (it
            # probably isn't)
            if "job label" in self.production.meta:
                job_label = self.production.meta["job label"]
            else:
                job_label = self.production.name
            dag_filename = f"dag_{job_label}.submit"
            command = [
                # "ssh", f"{config.get('scheduler', 'server')}",
                "condor_submit_dag",
                "-batch-name",
                f"bilbyMDR/{self.production.event.name}/{self.production.name}",
                os.path.join(self.production.rundir, "submit", dag_filename),
            ]

            if dryrun:
                print(" ".join(command))
            else:

                # with set_directory(self.production.rundir):
                self.logger.info(f"Working in {os.getcwd()}")

                dagman = subprocess.Popen(
                    command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
                )

                self.logger.info(" ".join(command))

                stdout, stderr = dagman.communicate()

                if "submitted to cluster" in str(stdout):
                    cluster = re.search(
                        r"submitted to cluster ([\d]+)", str(stdout)
                    ).groups()[0]
                    self.logger.info(
                        f"Submitted successfully. Running with job id {int(cluster)}"
                    )
                    self.production.status = "running"
                    self.production.job_id = int(cluster)
                    return cluster, PipelineLogger(stdout)
                else:
                    self.logger.error("Could not submit the job to the cluster")
                    self.logger.info(stdout)
                    self.logger.error(stderr)

                    raise PipelineException(
                        "The DAG file could not be submitted.",
                    )

        except FileNotFoundError as error:
            self.logger.exception(error)
            raise PipelineException(
                "It looks like condor isn't installed on this system.\n"
                f"""I wanted to run {" ".join(command)}."""
            ) from error

    def after_completion(self):
        post_pipeline = PESummaryPipeline(production=self.production)
        self.logger.info("Job has completed. Running PE Summary.")
        cluster = post_pipeline.submit_dag()
        self.production.meta["job id"] = int(cluster)
        self.production.status = "processing"
        self.production.event.update_data()


class PESummaryPipeline(PostPipeline):
    """
    A postprocessing pipeline add-in using PESummary customized for MDR analysis.
    """

    name = "PESummary"

    def submit_dag(self, dryrun=False):
        """
        Run PESummary on the results of this job.
        """

        psds = {ifo: os.path.abspath(psd) for ifo, psd in self.production.psds.items()}

        if "calibration" in self.production.meta["data"]:
            calibration = [
                os.path.abspath(os.path.join(self.production.repository.directory, cal))
                if not cal[0] == "/"
                else cal
                for cal in self.production.meta["data"]["calibration"].values()
            ]
        else:
            calibration = None

        configfile = self.production.event.repository.find_prods(
            self.production.name, self.category
        )[0]

        # if gr result present, use it to make comparison pages (turned off for now)
        # compare_with_gr = self.production.meta['gr pe info']['available']
        compare_with_gr = False

        if compare_with_gr:
            labels = " ".join([self.production.name, 'gr'])
            samples = " ".join([self.production.pipeline.samples(absolute=True)[0],
                                self.production.meta['gr pe info']['result file path']])
            configs = " ".join([
                os.path.join(
                    self.production.event.repository.directory, self.category, configfile),
                self.production.meta['gr pe info']['config file path']])
            approximants = " ".join([self.production.meta["waveform"]["approximant"],
                                    self.production.meta['gr pe info']["approximant"]])
        else:
            labels = self.production.name
            samples = self.production.pipeline.samples(absolute=True)[0]
            configs = os.path.join(
                self.production.event.repository.directory, self.category, configfile)
            approximants = self.production.meta["waveform"]["approximant"]

        command = [
            "--webdir",
            os.path.join(
                config.get("project", "root"),
                config.get("general", "webroot"),
                self.production.event.name,
                self.production.name,
                "pesummary",
            ),
            "--labels",
            labels,
            "--gw",
            "--approximant",
            approximants,
            "--f_low",
            str(min(self.production.meta["quality"]["minimum frequency"].values())),
            "--f_ref",
            str(self.production.meta["waveform"]["reference frequency"]),
        ]

        if "cosmology" in self.meta:
            command += [
                "--cosmology",
                self.meta["cosmology"],
            ]
        if "redshift" in self.meta:
            command += ["--redshift_method", self.meta["redshift"]]
        if "skymap samples" in self.meta:
            command += [
                "--nsamples_for_skymap",
                str(
                    self.meta["skymap samples"]
                ),  # config.get('pesummary', 'skymap_samples'),
            ]

        if "evolve spins" in self.meta:
            if "forwards" in self.meta["evolve spins"]:
                command += ["--evolve_spins_fowards", "True"]
            if "backwards" in self.meta["evolve spins"]:
                command += ["--evolve_spins_backwards", "precession_averaged"]

        if "nrsur" in self.production.meta["waveform"]["approximant"].lower():
            command += ["--NRSur_fits"]

        if "calculate" in self.meta:
            if "precessing snr" in self.meta["calculate"]:
                command += ["--calculate_precessing_snr"]

        if "multiprocess" in self.meta:
            command += ["--multi_process", str(self.meta["multiprocess"])]

        if "regenerate" in self.meta:
            command += ["--regenerate", " ".join(self.meta["regenerate posteriors"])]

        # Config file
        command += [
            "--config",
            configs,
        ]
        # Samples
        command += ["--samples"]
        command += [samples]
        # Calibration information
        if calibration:
            command += ["--calibration"]
            command += calibration
        # PSDs
        command += ["--psd"]
        for key, value in psds.items():
            command += [f"{key}:{value}"]

        # adding mdr parameter to corner plot
        command += ["--add_to_corner"]
        if 'mass_graviton_eff' in self.production.meta['priors']:
            command += ["mass_graviton"]
        else:
            command += ['A_alpha']

        with utils.set_directory(self.production.rundir):
            with open(f"{self.production.name}_pesummary.sh", "w") as bash_file:
                bash_file.write(
                    f"{config.get('pesummary', 'executable')} " + " ".join(command)
                )

        self.logger.info(
            f"PE summary command: {config.get('pesummary', 'executable')} {' '.join(command)}"
        )

        if dryrun:
            print("PESUMMARY COMMAND")
            print("-----------------")
            print(command)

        submit_description = {
            "executable": config.get("pesummary", "executable"),
            "arguments": " ".join(command),
            "output": f"{self.production.rundir}/pesummary.out",
            "error": f"{self.production.rundir}/pesummary.err",
            "log": f"{self.production.rundir}/pesummary.log",
            "request_cpus": self.meta["multiprocess"],
            "environment": "HDF5_USE_FILE_LOCKING=FAlSE OMP_NUM_THREADS=1 OMP_PROC_BIND=false",
            "getenv": "CONDA_EXE,USER,LAL*,PATH",
            "batch_name": f"PESummary/{self.production.event.name}/{self.production.name}",
            "request_memory": "8192MB",
            # "should_transfer_files": "YES",
            "request_disk": "8192MB",
            "+flock_local": "True",
            "+DESIRED_Sites": htcondor.classad.quote("nogrid"),
        }

        if "accounting group" in self.meta:
            submit_description["accounting_group_user"] = config.get("condor", "user")
            submit_description["accounting_group"] = self.meta["accounting group"]
        else:
            self.logger.warning(
                "This PESummary Job does not supply any accounting"
                " information, which may prevent it running on"
                " some clusters."
            )

        if dryrun:
            print("SUBMIT DESCRIPTION")
            print("------------------")
            print(submit_description)

        if not dryrun:
            hostname_job = htcondor.Submit(submit_description)

            with utils.set_directory(self.production.rundir):
                with open("pesummary.sub", "w") as subfile:
                    subfile.write(hostname_job.__str__())

            try:
                # There should really be a specified submit node, and if there is, use it.
                schedulers = htcondor.Collector().locate(
                    htcondor.DaemonTypes.Schedd, config.get("condor", "scheduler")
                )
                schedd = htcondor.Schedd(schedulers)
            except:  # NoQA
                # If you can't find a specified scheduler, use the first one you find
                schedd = htcondor.Schedd()
            with schedd.transaction() as txn:
                cluster_id = hostname_job.queue(txn)

        else:
            cluster_id = 0

        return cluster_id
