import unittest
from copy import copy
import bilby
import bilby_tgr


class Test_LOSA_BBH(unittest.TestCase):
    def setUp(self):
        self.parameters = dict(
            mass_1=10,
            mass_2=8,
            luminosity_distance=200.0,
            a_1=0.1,
            tilt_1=0.1,
            phi_12=1.0,
            a_2=0.05,
            tilt_2=2.7,
            phi_jl=2.9,
            theta_jn=0.3,
            phase=0.0,
            a_by_c=1e-3,
        )
        self.waveform_kwargs = dict(
            waveform_approximant="IMRPhenomXP",
            reference_frequency=20.0,
            minimum_frequency=20.0,
            catch_waveform_errors=True,
        )
        self.frequency_array = bilby.core.utils.create_frequency_series(2048, 4)
        self.bad_parameters = copy(self.parameters)
        self.bad_parameters["mass_1"] = -3.0

    def tearDown(self):
        del self.parameters
        del self.waveform_kwargs
        del self.frequency_array
        del self.bad_parameters

    def test_bbh_losa_works_on_valid_parameters(self):
        self.parameters.update(self.waveform_kwargs)
        self.assertIsInstance(
            bilby_tgr.losa.source.lal_binary_black_hole_losa(
                self.frequency_array, **self.parameters
            ),
            dict,
        )

    def test_waveform_error_catching(self):
        self.bad_parameters.update(self.waveform_kwargs)
        with self.assertRaises(TypeError):
            bilby_tgr.losa.source.lal_binary_black_hole_losa(
                self.frequency_array, **self.bad_parameters
            )

    def test_waveform_error_raising(self):
        raise_error_parameters = copy(self.bad_parameters)
        raise_error_parameters.update(self.waveform_kwargs)
        raise_error_parameters["catch_waveform_errors"] = False
        with self.assertRaises(Exception):
            bilby_tgr.losa.source.lal_binary_black_hole_losa(
                self.frequency_array, **raise_error_parameters
            )

    def test_bbh_losa_works_without_waveform_parameters(self):
        self.assertIsInstance(
            bilby_tgr.losa.source.lal_binary_black_hole_losa(
                self.frequency_array, **self.parameters
            ),
            dict,
        )


class Test_LOSA_BNS(unittest.TestCase):
    def setUp(self):
        self.parameters = dict(
            mass_1=1.4,
            mass_2=1.25,
            luminosity_distance=40.0,
            a_1=0.1,
            tilt_1=0.1,
            phi_12=1.0,
            a_2=0.05,
            tilt_2=2.7,
            phi_jl=2.9,
            theta_jn=0.3,
            phase=0.0,
            lambda_1=150,
            lambda_2=200,
            a_by_c=1e-5,
        )
        self.waveform_kwargs = dict(
            waveform_approximant="IMRPhenomXP_NRTidalv2",
            reference_frequency=20.0,
            minimum_frequency=20.0,
            catch_waveform_errors=True,
        )
        self.frequency_array = bilby.core.utils.create_frequency_series(2048, 4)
        self.bad_parameters = copy(self.parameters)
        self.bad_parameters["mass_1"] = -3.0

    def tearDown(self):
        del self.parameters
        del self.waveform_kwargs
        del self.frequency_array
        del self.bad_parameters

    def test_bns_losa_works_on_valid_parameters(self):
        self.parameters.update(self.waveform_kwargs)
        self.assertIsInstance(
            bilby_tgr.losa.source.lal_binary_neutron_star_losa(
                self.frequency_array, **self.parameters
            ),
            dict,
        )

    def test_waveform_error_catching(self):
        self.bad_parameters.update(self.waveform_kwargs)
        with self.assertRaises(TypeError):
            bilby_tgr.losa.source.lal_binary_neutron_star_losa(
                self.frequency_array, **self.bad_parameters
            )

    def test_waveform_error_raising(self):
        raise_error_parameters = copy(self.bad_parameters)
        raise_error_parameters.update(self.waveform_kwargs)
        raise_error_parameters["catch_waveform_errors"] = False
        with self.assertRaises(Exception):
            bilby_tgr.losa.source.lal_binary_neutron_star_losa(
                self.frequency_array, **raise_error_parameters
            )

    def test_bns_losa_works_without_waveform_parameters(self):
        self.assertIsInstance(
            bilby_tgr.losa.source.lal_binary_neutron_star_losa(
                self.frequency_array, **self.parameters
            ),
            dict,
        )


class Test_LOSA_BBH_Relative_Binning(unittest.TestCase):
    def setUp(self):
        self.parameters = dict(
            mass_1=10,
            mass_2=8,
            luminosity_distance=200.0,
            a_1=0.1,
            tilt_1=0.1,
            phi_12=1.0,
            a_2=0.05,
            tilt_2=2.7,
            phi_jl=2.9,
            theta_jn=0.3,
            phase=0.0,
            a_by_c=1e-3,
            fiducial=1
        )
        self.waveform_kwargs = dict(
            waveform_approximant="IMRPhenomXP",
            reference_frequency=20.0,
            minimum_frequency=20.0,
            catch_waveform_errors=True,
        )
        self.frequency_array = bilby.core.utils.create_frequency_series(2048, 4)
        self.bad_parameters = copy(self.parameters)
        self.bad_parameters["mass_1"] = -3.0

    def tearDown(self):
        del self.parameters
        del self.waveform_kwargs
        del self.frequency_array
        del self.bad_parameters

    def test_bbh_losa_works_on_valid_parameters(self):
        self.parameters.update(self.waveform_kwargs)
        self.assertIsInstance(
            bilby_tgr.losa.source.lal_binary_black_hole_relative_binning_losa(
                self.frequency_array, **self.parameters
            ),
            dict,
        )

    def test_waveform_error_catching(self):
        self.bad_parameters.update(self.waveform_kwargs)
        with self.assertRaises(TypeError):
            bilby_tgr.losa.source.lal_binary_black_hole_relative_binning_losa(
                self.frequency_array, **self.bad_parameters
            )

    def test_waveform_error_raising(self):
        raise_error_parameters = copy(self.bad_parameters)
        raise_error_parameters.update(self.waveform_kwargs)
        raise_error_parameters["catch_waveform_errors"] = False
        with self.assertRaises(Exception):
            bilby_tgr.losa.source.lal_binary_black_hole_relative_binning_losa(
                self.frequency_array, **raise_error_parameters
            )

    def test_bbh_losa_works_without_waveform_parameters(self):
        self.assertIsInstance(
            bilby_tgr.losa.source.lal_binary_black_hole_relative_binning_losa(
                self.frequency_array, **self.parameters
            ),
            dict,
        )


class Test_LOSA_BNS_Relative_Binning(unittest.TestCase):
    def setUp(self):
        self.parameters = dict(
            mass_1=1.4,
            mass_2=1.25,
            luminosity_distance=40.0,
            a_1=0.1,
            tilt_1=0.1,
            phi_12=1.0,
            a_2=0.05,
            tilt_2=2.7,
            phi_jl=2.9,
            theta_jn=0.3,
            phase=0.0,
            lambda_1=150,
            lambda_2=200,
            a_by_c=1e-5,
            fiducial=1
        )
        self.waveform_kwargs = dict(
            waveform_approximant="IMRPhenomXP_NRTidalv2",
            reference_frequency=20.0,
            minimum_frequency=20.0,
            catch_waveform_errors=True,
        )
        self.frequency_array = bilby.core.utils.create_frequency_series(2048, 4)
        self.bad_parameters = copy(self.parameters)
        self.bad_parameters["mass_1"] = -3.0

    def tearDown(self):
        del self.parameters
        del self.waveform_kwargs
        del self.frequency_array
        del self.bad_parameters

    def test_bns_losa_works_on_valid_parameters(self):
        self.parameters.update(self.waveform_kwargs)
        self.assertIsInstance(
            bilby_tgr.losa.source.lal_binary_neutron_star_relative_binning_losa(
                self.frequency_array, **self.parameters
            ),
            dict,
        )

    def test_waveform_error_catching(self):
        self.bad_parameters.update(self.waveform_kwargs)
        with self.assertRaises(TypeError):
            bilby_tgr.losa.source.lal_binary_neutron_star_relative_binning_losa(
                self.frequency_array, **self.bad_parameters
            )

    def test_waveform_error_raising(self):
        raise_error_parameters = copy(self.bad_parameters)
        raise_error_parameters.update(self.waveform_kwargs)
        raise_error_parameters["catch_waveform_errors"] = False
        with self.assertRaises(Exception):
            bilby_tgr.losa.source.lal_binary_neutron_star_relative_binning_losa(
                self.frequency_array, **raise_error_parameters
            )

    def test_bns_losa_works_without_waveform_parameters(self):
        self.assertIsInstance(
            bilby_tgr.losa.source.lal_binary_neutron_star_relative_binning_losa(
                self.frequency_array, **self.parameters
            ),
            dict,
        )


if __name__ == "__main__":
    unittest.main()
