import unittest
import glob

import bilby_tgr


class TestPackagedPriors(unittest.TestCase):
    """ Test that the prepackaged priors load """

    def test_aligned(self):
        filename = "aligned_spins_bbh.prior"
        prior_dict = bilby_tgr.pseob.prior.BBHPriorDict(filename=filename)
        self.assertTrue("chi_1" in prior_dict)
        self.assertTrue("chi_2" in prior_dict)
        self.assertTrue("dA22" in prior_dict)

    def test_aligned_v4(self):
        filename = "aligned_spins_bbh_v4.prior"
        prior_dict = bilby_tgr.pseob.prior.BBHPriorDict(filename=filename)
        self.assertTrue("chi_1" in prior_dict)
        self.assertTrue("chi_2" in prior_dict)
        self.assertFalse("dA22" in prior_dict)

    def test_binary_black_holes(self):
        filename = "precessing_spins_bbh.prior"
        prior_dict = bilby_tgr.pseob.prior.BBHPriorDict(filename=filename)
        self.assertTrue("a_1" in prior_dict)
        self.assertTrue("dA22" in prior_dict)

    def test_all(self):
        prior_files = glob.glob(bilby_tgr.pseob.prior.DEFAULT_PRIOR_DIR + "/*prior")
        for ff in prior_files:
            print("Checking prior file {}".format(ff))
            prior_dict = bilby_tgr.pseob.prior.BBHPriorDict(filename=ff)
            self.assertTrue("chirp_mass" in prior_dict)
            self.assertTrue("mass_ratio" in prior_dict)
            if "precessing" in ff:
                self.assertTrue("a_1" in prior_dict)
            elif "aligned" in ff:
                self.assertTrue("chi_1" in prior_dict)


if __name__ == "__main__":
    unittest.main()
