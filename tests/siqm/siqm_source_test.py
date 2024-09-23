import unittest
from copy import copy
import bilby
import bilby_tgr


class TestLalSIQM(unittest.TestCase):
    def setUp(self):
        self.parameters = dict(
            mass_1=20.0,
            mass_2=14.0,
            luminosity_distance=400.0,
            a_1=0.4,
            a_2=0.3,
            tilt_1=0.2,
            tilt_2=1.7,
            phi_jl=0.2,
            phi_12=0.9,
            theta_jn=1.7,
            phase=0.0,
            dQuadMon1=100.0,
            dQuadMon2=100.0,
        )
        self.waveform_kwargs = dict(
            waveform_approximant="IMRPhenomPv2",
            reference_frequency=50.0,
            catch_waveform_errors=True,
            minimum_frequency=20.0,
        )
        self.frequency_array = bilby.core.utils.create_frequency_series(2048, 4)
        self.bad_parameters = copy(self.parameters)
        self.bad_parameters["mass_1"] = -30.0

    def tearDown(self):
        del self.parameters
        del self.waveform_kwargs
        del self.frequency_array

    def test_waveform_error_catching_siqm(self):
        self.bad_parameters.update(self.waveform_kwargs)
        self.assertIsNone(
            bilby_tgr.siqm.source.lal_binary_black_hole(
                self.frequency_array, **self.bad_parameters
            )
        )

    def test_waveform_error_raising_siqm(self):
        raise_error_parameters = copy(self.bad_parameters)
        raise_error_parameters.update(self.waveform_kwargs)
        raise_error_parameters["catch_waveform_errors"] = False
        with self.assertRaises(Exception):
            bilby_tgr.siqm.source.lal_binary_black_hole(
                self.frequency_array, **raise_error_parameters
            )

    def test_lal_binary_black_hole_works_without_siqm_parameters(self):
        self.parameters.pop("dQuadMon1")
        self.parameters.pop("dQuadMon2")
        self.parameters.update(self.waveform_kwargs)
        with self.assertRaises(TypeError):
            bilby_tgr.siqm.source.lal_binary_black_hole(
                self.frequency_array, **self.parameters
            )

    def test_lal_binary_black_hole_runs_with_valid_parameters(self):
        self.parameters.update(self.waveform_kwargs)
        self.assertIsInstance(
            bilby_tgr.siqm.source.lal_binary_black_hole(
                self.frequency_array, **self.parameters
            ),
            dict,
        )


if __name__ == "__main__":
    unittest.main()
