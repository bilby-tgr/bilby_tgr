import unittest
from copy import copy
from unittest.mock import patch

import astropy.units as u
import bilby
import lal
import numpy as np
from gwpy.frequencyseries import FrequencySeries
from lalsimulation.gwsignal.models.pyseobnr_model import SEOBNRv5PHM as wf_gen

import bilby_tgr


class TestpSEOBNR(unittest.TestCase):
    def setUp(self):

        self.waveform_kwargs = {
            "waveform_approximant": "SEOBNRv5PHM",
            "mass_1": 60.0,
            "mass_2": 15.0,
            "luminosity_distance": 1429.88270378,
            "a_1": 0.6,
            "a_2": 0.4,
            "tilt_1": 1.4,
            "tilt_2": 0.1,
            "phi_12": 0.76,
            "phi_jl": 1.31,
            "theta_jn": 0.45541468,
            "phase": 0.4675536,
            "minimum_frequency": 15.0,
            "reference_frequency": 20.0,
            "catch_waveform_errors": True,
        }

        self.pseob_params = {
            "domega220",
            "dtau220",
            "domega330",
            "dtau330",
            "domega210",
            "dtau210",
            "domega440",
            "dtau440",
            "domega550",
            "dtau550",
            "domega320",
            "dtau320",
            "domega430",
            "dtau430",
            "dA22",
            "dw22",
            "dA33",
            "dw33",
            "dA21",
            "dw21",
            "dA44",
            "dw44",
            "dA55",
            "dw55",
            "dA32",
            "dw32",
            "dA43",
            "dw43",
            "dTpeak",
            "da6",
            "ddSO",
        }

        self.frequency_array = bilby.core.utils.create_frequency_series(2048, 4)
        self.bad_parameters = copy(self.waveform_kwargs)
        self.bad_parameters["mass_1"] = 300000.0  # this generates an input domain error

    def test_pseob_works_runs_valid_parameters(self):
        self.waveform_kwargs.update({_: 0.0 for _ in self.pseob_params})
        # note: currently if the pSEOB specific variables are all 0, they are removed
        # from the parameters passed to gwsignal.
        # This makes it possible to run those tests without the proper changes in gwsignal
        # and pyseobnr. This is not a desired property and one would need to activate
        # the end-to-end pyseobnr tests, simply by making it so that one of the values
        # is non-zero. This would be possible once the other MR are merged:
        # * https://git.ligo.org/lscsoft/lalsuite/-/merge_requests/2319
        # * https://git.ligo.org/waveforms/software/pyseobnr/-/merge_requests/57

        # example for activating the end-to-end tests from bilby_tgr:
        # self.waveform_kwargs["dTpeak"] = 0.1
        self.assertIsInstance(
            bilby_tgr.pseob.source.gwsignal_binary_black_hole(
                self.frequency_array, **self.waveform_kwargs
            ),
            dict,
        )

    def test_waveform_error_catching(self):
        erroneous_parameters = copy(self.bad_parameters)
        assert erroneous_parameters["catch_waveform_errors"]
        erroneous_parameters.update({_: 0.0 for _ in self.pseob_params})
        self.assertIsNone(
            bilby_tgr.pseob.source.gwsignal_binary_black_hole(
                self.frequency_array, **erroneous_parameters
            )
        )

    def test_waveform_error_raising(self):
        raise_error_parameters = copy(self.waveform_kwargs)
        raise_error_parameters.update(self.bad_parameters)
        raise_error_parameters.update({_: 0.0 for _ in self.pseob_params})
        raise_error_parameters["catch_waveform_errors"] = False
        with self.assertRaises(Exception):
            bilby_tgr.pseob.source.gwsignal_binary_black_hole(
                self.frequency_array, **raise_error_parameters
            )

    def test_report_errors_that_are_not_input_domain_errors(self):
        raise_error_parameters = copy(self.waveform_kwargs)
        raise_error_parameters["mass_1"] = -10
        raise_error_parameters.update({_: 0.0 for _ in self.pseob_params})
        # this should be reported even if we are catching the errors
        raise_error_parameters["catch_waveform_errors"] = True
        with self.assertRaises(Exception):
            bilby_tgr.pseob.source.gwsignal_binary_black_hole(
                self.frequency_array, **raise_error_parameters
            )

    @staticmethod
    def fake_fd_signal(hp, hc, deltaT=1 / 2048):
        # generates fake fd signal
        epoch = lal.LIGOTimeGPS(-300)

        hp_lal = lal.CreateREAL8TimeSeries(
            "hplus", epoch, 0, deltaT, lal.DimensionlessUnit, len(hp)
        )
        hc_lal = lal.CreateREAL8TimeSeries(
            "hcross",
            epoch,
            0,
            1 / 2048,
            lal.DimensionlessUnit,
            len(hp),
        )

        hp_lal.data.data = hp
        hc_lal.data.data = hc

        chirplen = hp_lal.data.length
        chirplen_exp = np.frexp(chirplen)
        chirplen = int(np.ldexp(1, chirplen_exp[1]))
        deltaF = 1.0 / (chirplen * deltaT)

        lal.ResizeREAL8TimeSeries(hp_lal, hp_lal.data.length - chirplen, chirplen)
        lal.ResizeREAL8TimeSeries(hc_lal, hc_lal.data.length - chirplen, chirplen)

        # FFT - Using LAL routines
        hptilde = lal.CreateCOMPLEX16FrequencySeries(
            "FD H_PLUS",
            hp_lal.epoch,
            0.0,
            deltaF,
            lal.DimensionlessUnit,
            int(chirplen / 2.0 + 1),
        )
        hctilde = lal.CreateCOMPLEX16FrequencySeries(
            "FD H_CROSS",
            hc_lal.epoch,
            0.0,
            deltaF,
            lal.DimensionlessUnit,
            int(chirplen / 2.0 + 1),
        )

        plan = lal.CreateForwardREAL8FFTPlan(chirplen, 0)
        lal.REAL8TimeFreqFFT(hctilde, hc_lal, plan)
        lal.REAL8TimeFreqFFT(hptilde, hp_lal, plan)

        frequencies = hptilde.deltaF * np.arange(hptilde.data.length)
        epoch = hptilde.epoch.gpsSeconds + hptilde.epoch.gpsNanoSeconds / 1e9
        return (
            FrequencySeries(
                hptilde.data.data,
                frequencies=frequencies,
                epoch=epoch,
                name="hp",
                unit=u.Unit("strain") * u.s,
            ),
            FrequencySeries(
                hctilde.data.data,
                frequencies=frequencies,
                epoch=epoch,
                name="hc",
                unit=u.Unit("strain") * u.s,
            ),
        )

    def test_pseob_specific_parameters(self):
        # here we test the behaviour of the function "gwsignal_binary_black_hole"
        # until the execution of the "generate_fd_waveform" in gwsignal. In particular
        # the actual generate_fd_waveform is not called and only the parameters passed to
        # the function are checked.
        # The test does not require gwsignal to support pSEOB parameters.

        hp = np.zeros(1000)
        hc = np.zeros(1000)

        return_val_fd_hptilde, return_val_fd_hctilde = self.fake_fd_signal(
            hp, hc, 1 / 2048
        )

        with patch.object(wf_gen, "generate_fd_waveform", autospec=True) as wgen_gen_fd:

            for current_param in self.pseob_params:

                wgen_gen_fd.return_value = return_val_fd_hptilde, return_val_fd_hctilde

                wgen_gen_fd.reset_mock()
                parameters = copy(self.waveform_kwargs)
                parameters.update({_: 0.0 for _ in self.pseob_params})
                parameters.update({current_param: 0.1})

                bilby_tgr.pseob.source.gwsignal_binary_black_hole(
                    self.frequency_array, **parameters
                )

                # checks generated_fd_waveform is called
                wgen_gen_fd.assert_called_once()

                # checks parameters are passed: this is done inside the
                self.assertIsInstance(wgen_gen_fd.call_args_list[0].args[0], wf_gen)

                # check if the currently modified parameter is dA22, dtau32 etc
                for parameter_type in "dA", "domega", "dw", "dtau":
                    if current_param.startswith(parameter_type):
                        param_to_check = f"{parameter_type}_dict"

                        # the trailing part dtau330 -> 330 -> 3,3
                        element_to_check = ",".join(
                            current_param[len(parameter_type) : len(parameter_type) + 2]
                        )

                        self.assertIn(
                            param_to_check, wgen_gen_fd.call_args_list[0].kwargs
                        )
                        self.assertIn(
                            element_to_check,
                            wgen_gen_fd.call_args_list[0].kwargs[param_to_check],
                        )
                        self.assertEqual(
                            wgen_gen_fd.call_args_list[0].kwargs[param_to_check][
                                element_to_check
                            ],
                            0.1,
                        )
                        break

                else:
                    if current_param in ("dTpeak", "da6", "ddSO"):
                        self.assertIn(
                            current_param, wgen_gen_fd.call_args_list[0].kwargs
                        )
                        self.assertEqual(
                            wgen_gen_fd.call_args_list[0].kwargs[current_param], 0.1
                        )
                    else:
                        self.fail(f"Parameter {current_param} not found!")
