import unittest

from copy import copy

import bilby
import bilby_tgr


class TestLalBBH(unittest.TestCase):
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
            dchi_Minus2=0.1,
            dchi_0=0.1,
            dchi_1=0.1,
            dchi_2=0.1,
            dchi_3=0.1,
            dchi_3S=0.1,
            dchi_3NS=0.1,
            dchi_4=0.1,
            dchi_4S=0.1,
            dchi_4NS=0.1,
            dchi_5l=0.1,
            dchi_5lS=0.1,
            dchi_5lNS=0.1,
            dchi_6=0.1,
            dchi_6S=0.1,
            dchi_6NS=0.1,
            dchi_6l=0.1,
            dchi_7=0.1,
            dchi_7S=0.1,
            dchi_7NS=0.1,
            dchi_kappaS=0.1,
            dchi_kappaA=0.1,
            f_window_div_f_peak=1.0,
            NumCycles=1,
        )
        self.waveform_kwargs = dict(
            waveform_approximant="IMRPhenomXAS",
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

    def test_lal_bbh_works_runs_valid_parameters(self):
        self.parameters.update(self.waveform_kwargs)
        self.assertIsInstance(
            bilby_tgr.fti.source.lal_binary_black_hole(
                self.frequency_array, **self.parameters
            ),
            dict,
        )

    def test_waveform_error_catching(self):
        self.bad_parameters.update(self.waveform_kwargs)
        self.assertIsNone(
            bilby_tgr.fti.source.lal_binary_black_hole(
                self.frequency_array, **self.bad_parameters
            )
        )

    def test_waveform_error_raising(self):
        raise_error_parameters = copy(self.bad_parameters)
        raise_error_parameters.update(self.waveform_kwargs)
        raise_error_parameters["catch_waveform_errors"] = False
        with self.assertRaises(Exception):
            bilby_tgr.fti.source.lal_binary_black_hole(
                self.frequency_array, **raise_error_parameters
            )

    def test_lal_bbh_works_without_waveform_parameters(self):
        self.assertIsInstance(
            bilby_tgr.fti.source.lal_binary_black_hole(
                self.frequency_array, **self.parameters
            ),
            dict,
        )


class TestLalBNS(unittest.TestCase):
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
            dchi_Minus2=0.1,
            dchi_0=0.1,
            dchi_1=0.1,
            dchi_2=0.1,
            dchi_3=0.1,
            dchi_3S=0.1,
            dchi_3NS=0.1,
            dchi_4=0.1,
            dchi_4S=0.1,
            dchi_4NS=0.1,
            dchi_5l=0.1,
            dchi_5lS=0.1,
            dchi_5lNS=0.1,
            dchi_6=0.1,
            dchi_6S=0.1,
            dchi_6NS=0.1,
            dchi_6l=0.1,
            dchi_7=0.1,
            dchi_7S=0.1,
            dchi_7NS=0.1,
            dchi_kappaS=0.1,
            dchi_kappaA=0.1,
            f_window_div_f_peak=1.0,
            NumCycles=1,
        )
        self.waveform_kwargs = dict(
            waveform_approximant="IMRPhenomXAS_NRTidalv2",
            reference_frequency=50.0,
            minimum_frequency=20.0,
        )
        self.frequency_array = bilby.core.utils.create_frequency_series(2048, 4)

    def tearDown(self):
        del self.parameters
        del self.waveform_kwargs
        del self.frequency_array

    def test_lal_bns_runs_with_valid_parameters(self):
        self.parameters.update(self.waveform_kwargs)
        self.assertIsInstance(
            bilby_tgr.fti.source.lal_binary_neutron_star(
                self.frequency_array, **self.parameters
            ),
            dict,
        )

    def test_lal_bns_works_without_waveform_parameters(self):
        self.assertIsInstance(
            bilby_tgr.fti.source.lal_binary_neutron_star(
                self.frequency_array, **self.parameters
            ),
            dict,
        )

    def test_fails_without_tidal_parameters(self):
        self.parameters.pop("lambda_1")
        self.parameters.pop("lambda_2")
        self.parameters.update(self.waveform_kwargs)
        with self.assertRaises(TypeError):
            bilby_tgr.fti.source.lal_binary_neutron_star(
                self.frequency_array, **self.parameters
            )


if __name__ == "__main__":
    unittest.main()
