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

    def test_non_gr_d_alpha_2_binary_black_hole(self):
        """ A simple test of the non_gr_d_alpha_2_binary_black_hole """

        # Call the function with a fixed set of parameters
        polarization_dict = bilby_tgr.source.non_gr_d_alpha_2_binary_black_hole(
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
            d_alpha_2=0.5)

        # Run tests to check the output is as expected: a dictionary of arrays
        self.assertTrue(isinstance(polarization_dict, dict))
        self.assertEqual(list(polarization_dict.keys()), ["plus", "cross"])
        self.assertTrue(isinstance(polarization_dict["plus"], np.ndarray))
        self.assertTrue(isinstance(polarization_dict["cross"], np.ndarray))


if __name__ == "__main__":
    unittest.main()
