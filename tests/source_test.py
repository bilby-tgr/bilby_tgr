import unittest

import bilby_tgr
import numpy as np


class TestPriorInstantiationWithoutOptionalPriors(unittest.TestCase):
    def setUp(self):
        """ This is called before each test: use to reduce repeated code """
        self.frequency_array = np.linspace(20, 2048, 1000)

    def tearDown(self):
        """ This is called after each test: use to clean up"""
        del self.frequency_array

    def test_lal_binary_black_hole_TIGER(self):
        """ A simple test of the TIGER_binary_black_hole """

        # Call the function with a fixed set of parameters
        polarization_dict = bilby_tgr.source.lal_binary_black_hole_TIGER(
            frequency_array=self.frequency_array,
            mass_1=31,
            mass_2=30,
            luminosity_distance=1000,
            a_1=0.1,
            tilt_1=0.2,
            phi_12=0.2,
            a_2=0.2,
            tilt_2=0.2,
            phi_jl=0.2,
            theta_jn=1.2,
            phase=0.4,
            dchi_0=0.1,
            dchi_1=0.1,
            dchi_2=0.1,
            dchi_3=0.1,
            dchi_4=0.1,
            dchi_5l=0.1,
            dchi_6=0.1,
            dchi_6l=0.1,
            dchi_7=0.1,
            dbeta_2=0.1,
            dbeta_3=0.1,
            dalpha_2=0.1,
            dalpha_3=0.1,
            dalpha_4=0.1,
            dalpha_5=0.1,
            waveform_approximant='IMRPhenomPv2')

        # Run tests to check the output is as expected: a dictionary of arrays
        self.assertTrue(isinstance(polarization_dict, dict))
        self.assertEqual(list(polarization_dict.keys()), ["plus", "cross"])
        self.assertTrue(isinstance(polarization_dict["plus"], np.ndarray))
        self.assertTrue(isinstance(polarization_dict["cross"], np.ndarray))

    def test_lal_binary_neutron_star_TIGER(self):
        """ A simple test of the TIGER_binary_black_hole """

        # Call the function with a fixed set of parameters
        polarization_dict = bilby_tgr.source.lal_binary_neutron_star_TIGER(
            frequency_array=self.frequency_array,
            mass_1=1.4,
            mass_2=1.4,
            luminosity_distance=100,
            a_1=0.1,
            tilt_1=0.2,
            phi_12=0.2,
            a_2=0.05,
            tilt_2=0.2,
            phi_jl=0.2,
            theta_jn=1.2,
            phase=0.4,
            lambda_1=500,
            lambda_2=500,
            dchi_0=0.1,
            dchi_1=0.1,
            dchi_2=0.1,
            dchi_3=0.1,
            dchi_4=0.1,
            dchi_5l=0.1,
            dchi_6=0.1,
            dchi_6l=0.1,
            dchi_7=0.1,
            waveform_approximant='IMRPhenomPv2_NRTidal')

        # Run tests to check the output is as expected: a dictionary of arrays
        self.assertTrue(isinstance(polarization_dict, dict))
        self.assertEqual(list(polarization_dict.keys()), ["plus", "cross"])
        self.assertTrue(isinstance(polarization_dict["plus"], np.ndarray))
        self.assertTrue(isinstance(polarization_dict["cross"], np.ndarray))


if __name__ == "__main__":
    unittest.main()
