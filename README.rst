|pipeline status| |coverage report|

=========
Bilby_TGR
=========

This is the Bilby Testing General Relativity (TGR) package. It is used to
develop and share analysis done by the LVC TGR group building on the
`Bilby <https://git.ligo.org/lscsoft/bilby>`_ stochastic sampling package.


Installation
------------

To get started developing and using :code:`bilby_tgr`, first clone the repository

.. code-block:: console

   $ git clone git@git.ligo.org:gregory.ashton/bilby_tgr.git
   
To install, enter the cloned directory and run

.. code-block:: console

   $ pip install .
   
To check that you have installed it correctly, open an python prompt and run

.. code-block:: python

   >>> import bilby_tgr
   >>> bilby_tgr.source.non_gr_d_alpha_2_binary_black_hole

If this returns the function, then you have it installed! You can now add new
functions to the sources and access them in the same way.

Running using bilby_pipe
------------------------

Once you have installed :code:`bilby_tgr`, you can use the :code:`bilby_pipe`
package to run stoachastic sampling. For help getting installed and setup with
:code:`bilby_pipe` itself, see `the documentation <https://git.ligo.org/lscsoft/bilby_pipe>`_.
Here, we give an example ini file. Notice that the :code:`frequency-domain-source-model`
is pointing to the :code:`bilby_tgr.source.non_gr_d_alpha_2_binary_black_hole`
function. You can replace this with any new function you care to define in the
:code:`bilby_tgr` package.

.. code-block:: console

   trigger-time = 1186741861.5
   detectors = [H1, L1, V1]
   channel-dict = {H1:DCH-CLEAN_STRAIN_C02, L1:DCH-CLEAN_STRAIN_C02, V1:Hrec_hoft_V1O2Repro2A_16384Hz}
   prior_file = 4s.prior
   time-marginalization=False
   distance-marginalization=True
   phase-marginalization=True
   create-plots=True
   local-generation = True
   psd-dict = {H1:BayesWave_median_PSD_H1.dat, L1:BayesWave_median_PSD_L1.dat, V1:BayesWave_median_PSD_V1.dat}

   label = GW170814
   outdir = d_alpha_2
   accounting = ligo.dev.o3.cbc.pe.lalinference

   duration = 4
   coherence-test = False

   maximum-frequency=1024
   minimum-frequency=20
   sampling-frequency=2048
   reference-frequency = 20
   waveform-approximant = IMRPhenomPV2
   frequency-domain-source-model = bilby_tgr.source.non_gr_d_alpha_2_binary_black_hole

   calibration-model=CubicSpline
   spline-calibration-envelope-dict = {H1:GWTC1_GW170814_H_CalEnv.txt, L1:GWTC1_GW170814_L_CalEnv.txt, V1:GWTC1_GW170814_V_CalEnv.txt}
   spline_calibration-nodes = 10

   deltaT = 0.2
   sampler = dynesty
   sampler-kwargs = {nlive: 1000, sample: rwalk, walks=50, nact=5}
   n-parallel = 4

   transfer-files = False

The other thing you need is a prior file (4s.prior in the above ini).
This will be a standard CBC prior, plus any new parameters.

.. code-block:: python

    chirp_mass = Uniform(name="chirp_mass", minimum=12.299703, maximum=45, unit='$M_{\\odot}$')
    mass_ratio = Uniform(name="mass_ratio", minimum=0.125, maximum=1)
    mass_1 = Constraint(name="mass_1", minimum=1.001398, maximum=1000)
    mass_2 = Constraint(name="mass_2", minimum=1.001398, maximum=1000)
    a_1 = Uniform(name="a_1", minimum=0, maximum=0.88)
    a_2 = Uniform(name="a_2", minimum=0, maximum=0.88)
    tilt_1 = Sine(name="tilt_1")
    tilt_2 = Sine(name="tilt_2")
    phi_12 = Uniform(name="phi_12", minimum=0, maximum=2 * np.pi, boundary="periodic")
    phi_jl = Uniform(name="phi_jl", minimum=0, maximum=2 * np.pi, boundary="periodic")
    luminosity_distance = bilby.gw.prior.UniformSourceFrame(name="luminosity_distance", minimum=1e2, maximum=5e3, unit="Mpc")
    dec = Cosine(name="dec")
    ra = Uniform(name="ra", minimum=0, maximum=2 * np.pi, boundary="periodic")
    theta_jn = Sine(name="theta_jn")
    psi = Uniform(name="psi", minimum=0, maximum=np.pi, boundary="periodic")
    phase = Uniform(name="phase", minimum=0, maximum=2 * np.pi, boundary="periodic")
    d_alpha_2 = Uniform(minimum=-10, maximum=10, latex_label="$\\delta \\alpha_2$")
   
   

.. |pipeline status| image:: https://git.ligo.org/gregory.ashton/bilby_tgr/badges/master/pipeline.svg
   :target: https://git.ligo.org/gregory.ashton/bilby_tgr/commits/master
.. |coverage report| image:: https://gregory.ashton.docs.ligo.org/bilby_tgr/coverage_badge.svg
   :target: https://gregory.ashton.docs.ligo.org/bilby_tgr/htmlcov/
