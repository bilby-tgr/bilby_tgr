import unittest

import numpy as np
from copy import copy

import bilby
import bilby_tgr


class TestSsbBBH(unittest.TestCase):
    def setUp(self):
        self.parameters = dict(
            mass_1=30.0,
            mass_2=30.0,
            luminosity_distance=400.0,
            a_1=0.4,
            tilt_1=0.2,
            phi_12=1.0,
            a_2=0.8,
            tilt_2=2.7,
            phi_jl=2.9,
            theta_jn=0.3,
            phase=0.0,
            ssb_coeff=8000,
        )
        self.waveform_kwargs = dict(
            waveform_approximant="IMRPhenomPv2",
            reference_frequency=50.0,
            minimum_frequency=20.0,
            catch_waveform_errors=True,
        )
        self.frequency_array = bilby.core.utils.create_frequency_series(2048, 4)
        self.bad_parameters = copy(self.parameters)
        self.bad_parameters["mass_1"] = -30.0

    def tearDown(self):
        del self.parameters
        del self.waveform_kwargs
        del self.frequency_array
        del self.bad_parameters

    def test_ssb_bbh_works_runs_valid_parameters(self):
        self.parameters.update(self.waveform_kwargs)
        self.assertIsInstance(
            bilby_tgr.ssb.source.lal_binary_black_hole(
                self.frequency_array, **self.parameters
            ),
            dict,
        )

    def test_waveform_error_catching(self):
        self.bad_parameters.update(self.waveform_kwargs)
        self.assertIsNone(
            bilby_tgr.ssb.source.lal_binary_black_hole(
                self.frequency_array, **self.bad_parameters
            )
        )

    def test_waveform_error_raising(self):
        raise_error_parameters = copy(self.bad_parameters)
        raise_error_parameters.update(self.waveform_kwargs)
        raise_error_parameters["catch_waveform_errors"] = False
        with self.assertRaises(Exception):
            bilby_tgr.ssb.source.lal_binary_black_hole(
                self.frequency_array, **raise_error_parameters
            )

    def test_ssb_bbh_works_without_waveform_parameters(self):
        self.assertIsInstance(
            bilby_tgr.ssb.source.lal_binary_black_hole(
                self.frequency_array, **self.parameters
            ),
            dict,
        )

    def test_ssb_bbh_xpprecession_version(self):
        self.parameters.update(self.waveform_kwargs)
        self.parameters["waveform_approximant"] = "IMRPhenomXP"

        # Test that we can modify the XP precession version
        out_v223 = bilby_tgr.ssb.source.lal_binary_black_hole(
            self.frequency_array, PhenomXPrecVersion=223, **self.parameters
        )
        out_v102 = bilby_tgr.ssb.source.lal_binary_black_hole(
            self.frequency_array, PhenomXPrecVersion=102, **self.parameters
        )
        self.assertFalse(np.all(out_v223["plus"] == out_v102["plus"]))


if __name__ == "__main__":
    unittest.main()
