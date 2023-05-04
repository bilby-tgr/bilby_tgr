import unittest

import numpy as np
from copy import copy

import bilby
import bilby_tgr


class TestMdrBBH(unittest.TestCase):
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
            alpha=3.0,
            A_eff=1e-20,
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

    def test_mdr_bbh_works_runs_valid_parameters(self):
        self.parameters.update(self.waveform_kwargs)
        self.assertIsInstance(
            bilby_tgr.mdr.source.lal_binary_black_hole(
                self.frequency_array, **self.parameters
            ),
            dict,
        )

    def test_waveform_error_catching(self):
        self.bad_parameters.update(self.waveform_kwargs)
        self.assertIsNone(
            bilby_tgr.mdr.source.lal_binary_black_hole(
                self.frequency_array, **self.bad_parameters
            )
        )

    def test_waveform_error_raising(self):
        raise_error_parameters = copy(self.bad_parameters)
        raise_error_parameters.update(self.waveform_kwargs)
        raise_error_parameters["catch_waveform_errors"] = False
        with self.assertRaises(Exception):
            bilby_tgr.mdr.source.lal_binary_black_hole(
                self.frequency_array, **raise_error_parameters
            )

    def test_mdr_bbh_works_without_waveform_parameters(self):
        self.assertIsInstance(
            bilby_tgr.mdr.source.lal_binary_black_hole(
                self.frequency_array, **self.parameters
            ),
            dict,
        )

    def test_mdr_bbh_xpprecession_version(self):
        self.parameters.update(self.waveform_kwargs)
        self.parameters["waveform_approximant"] = "IMRPhenomXP"

        # Test that we can modify the XP precession version
        out_v223 = bilby_tgr.mdr.source.lal_binary_black_hole(
            self.frequency_array, PhenomXPrecVersion=223, **self.parameters
        )
        out_v102 = bilby_tgr.mdr.source.lal_binary_black_hole(
            self.frequency_array, PhenomXPrecVersion=102, **self.parameters
        )
        self.assertFalse(np.all(out_v223["plus"] == out_v102["plus"]))


class TestMdrBNS(unittest.TestCase):
    def setUp(self):
        self.parameters = dict(
            mass_1=1.4,
            mass_2=1.4,
            luminosity_distance=400.0,
            a_1=0.4,
            a_2=0.3,
            tilt_1=0.2,
            tilt_2=1.7,
            phi_jl=0.2,
            phi_12=0.9,
            theta_jn=1.7,
            phase=0.0,
            lambda_1=100.0,
            lambda_2=100.0,
            alpha=3.0,
            A_eff=1e-20,
        )
        self.waveform_kwargs = dict(
            waveform_approximant="IMRPhenomPv2_NRTidal",
            reference_frequency=50.0,
            minimum_frequency=20.0,
        )
        self.frequency_array = bilby.core.utils.create_frequency_series(2048, 4)

    def tearDown(self):
        del self.parameters
        del self.waveform_kwargs
        del self.frequency_array

    def test_mdr_bns_runs_with_valid_parameters(self):
        self.parameters.update(self.waveform_kwargs)
        self.assertIsInstance(
            bilby_tgr.mdr.source.lal_binary_neutron_star(
                self.frequency_array, **self.parameters
            ),
            dict,
        )

    def test_mdr_bns_works_without_waveform_parameters(self):
        self.assertIsInstance(
            bilby_tgr.mdr.source.lal_binary_neutron_star(
                self.frequency_array, **self.parameters
            ),
            dict,
        )

    def test_fails_without_tidal_parameters(self):
        self.parameters.pop("lambda_1")
        self.parameters.pop("lambda_2")
        self.parameters.update(self.waveform_kwargs)
        with self.assertRaises(TypeError):
            bilby_tgr.mdr.source.lal_binary_neutron_star(
                self.frequency_array, **self.parameters
            )


if __name__ == "__main__":
    unittest.main()
