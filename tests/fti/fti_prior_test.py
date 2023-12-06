from collections import OrderedDict
import unittest
import glob
import os
import pickle

import bilby
import bilby_tgr


class TestBBHPriorDict(unittest.TestCase):
    def setUp(self):
        self.prior_dict = dict()
        self.base_directory = "/".join(
            os.path.dirname(os.path.realpath(__file__)).split("/")[:-2]
        )
        self.filename = os.path.join(
            self.base_directory,
            "bilby_tgr/fti/prior_files/aligned_spins_bbh.prior",
        )
        self.bbh_prior_dict = bilby_tgr.fti.prior.BBHPriorDict(filename=self.filename)
        for key, value in self.bbh_prior_dict.items():
            self.prior_dict[key] = value

    def tearDown(self):
        del self.prior_dict
        del self.filename
        del self.bbh_prior_dict

    def test_create_default_prior(self):
        default = bilby_tgr.fti.prior.BBHPriorDict()
        minima = all(
            [
                self.bbh_prior_dict[key].minimum == default[key].minimum
                for key in default.keys()
            ]
        )
        maxima = all(
            [
                self.bbh_prior_dict[key].maximum == default[key].maximum
                for key in default.keys()
            ]
        )
        names = all(
            [
                self.bbh_prior_dict[key].name == default[key].name
                for key in default.keys()
            ]
        )
        boundaries = all(
            [
                self.bbh_prior_dict[key].boundary == default[key].boundary
                for key in default.keys()
            ]
        )

        self.assertTrue(all([minima, maxima, names, boundaries]))

    def test_create_from_dict(self):
        new_dict = bilby_tgr.fti.prior.BBHPriorDict(dictionary=self.prior_dict)
        for key in self.bbh_prior_dict:
            self.assertEqual(self.bbh_prior_dict[key], new_dict[key])

    def test_redundant_priors_not_in_dict_before(self):
        for prior in [
            "chirp_mass",
            "total_mass",
            "mass_ratio",
            "symmetric_mass_ratio",
            "cos_theta_jn",
            "comoving_distance",
            "redshift",
            "dchi_kappa1",
            "dchi_kappa2",
        ]:
            self.assertTrue(self.bbh_prior_dict.test_redundancy(prior))

    def test_redundant_priors_already_in_dict(self):
        for prior in [
            "mass_1",
            "mass_2",
            "theta_jn",
            "luminosity_distance",
            "dchi_kappaS",
            "dchi_kappaA",
        ]:
            self.assertTrue(self.bbh_prior_dict.test_redundancy(prior))

    def test_correct_not_redundant_priors_masses(self):
        del self.bbh_prior_dict["chirp_mass"]
        for prior in ["chirp_mass", "total_mass", "symmetric_mass_ratio"]:
            self.assertFalse(self.bbh_prior_dict.test_redundancy(prior))

    def test_correct_not_redundant_priors_spin_magnitudes(self):
        del self.bbh_prior_dict["chi_2"]
        self.assertFalse(self.bbh_prior_dict.test_redundancy("chi_2"))

    def test_correct_not_redundant_priors_inclination(self):
        del self.bbh_prior_dict["theta_jn"]
        for prior in ["theta_jn", "cos_theta_jn"]:
            self.assertFalse(self.bbh_prior_dict.test_redundancy(prior))

    def test_correct_not_redundant_priors_distance(self):
        del self.bbh_prior_dict["luminosity_distance"]
        for prior in ["luminosity_distance", "comoving_distance", "redshift"]:
            self.assertFalse(self.bbh_prior_dict.test_redundancy(prior))

    def test_correct_not_redundant_priors_dchi_kappas(self):
        del self.bbh_prior_dict["dchi_kappaA"]
        self.assertFalse(self.bbh_prior_dict.test_redundancy("dchi_kappaA"))

    def test_add_unrelated_prior(self):
        self.assertFalse(self.bbh_prior_dict.test_redundancy("abc"))

    def test_test_has_redundant_priors(self):
        self.assertFalse(self.bbh_prior_dict.test_has_redundant_keys())
        for prior in [
            "mass_1",
            "mass_2",
            "total_mass",
            "symmetric_mass_ratio",
            "cos_theta_jn",
            "comoving_distance",
            "redshift",
            "dchi_kappa1",
            "dchi_kappa2",
        ]:
            self.bbh_prior_dict[prior] = 0
            self.assertTrue(self.bbh_prior_dict.test_has_redundant_keys())
            del self.bbh_prior_dict[prior]

    def test_add_constraint_prior_not_redundant(self):
        self.bbh_prior_dict["dchi_kappaA"] = bilby.prior.Constraint(
            minimum=20, maximum=40, name="dchi_kappaA"
        )
        self.assertFalse(self.bbh_prior_dict.test_has_redundant_keys())

    def test_pickle_prior(self):
        priors = dict(
            chirp_mass=bilby.core.prior.Uniform(10, 20),
            mass_ratio=bilby.core.prior.Uniform(0.125, 1),
        )
        priors = bilby_tgr.fti.prior.BBHPriorDict(priors)
        with open("test.pickle", "wb") as file:
            pickle.dump(priors, file)
        with open("test.pickle", "rb") as file:
            priors_loaded = pickle.load(file)
        self.assertEqual(priors, priors_loaded)
        os.remove("test.pickle")


class TestPackagedPriors(unittest.TestCase):
    """ Test that the prepackaged priors load """

    def test_aligned(self):
        filename = "aligned_spins_bbh.prior"
        prior_dict = bilby_tgr.fti.prior.BBHPriorDict(filename=filename)
        self.assertTrue("chi_1" in prior_dict)
        self.assertTrue("chi_2" in prior_dict)

    def test_binary_black_holes(self):
        filename = "precessing_spins_bbh.prior"
        prior_dict = bilby_tgr.fti.prior.BBHPriorDict(filename=filename)
        self.assertTrue("a_1" in prior_dict)

    def test_all(self):
        prior_files = glob.glob(bilby_tgr.fti.prior.DEFAULT_PRIOR_DIR + "/*prior")
        for ff in prior_files:
            print("Checking prior file {}".format(ff))
            prior_dict = bilby_tgr.fti.prior.BBHPriorDict(filename=ff)
            self.assertTrue("chirp_mass" in prior_dict)
            self.assertTrue("mass_ratio" in prior_dict)
            if "precessing" in ff:
                self.assertTrue("a_1" in prior_dict)
            elif "aligned" in ff:
                self.assertTrue("chi_1" in prior_dict)


class TestBNSPriorDict(unittest.TestCase):
    def setUp(self):
        self.prior_dict = OrderedDict()
        self.base_directory = "/".join(
            os.path.dirname(os.path.realpath(__file__)).split("/")[:-2]
        )
        self.filename = os.path.join(
            self.base_directory,
            "bilby_tgr/fti/prior_files/aligned_spins_bns_tides_on.prior",
        )
        self.bns_prior_dict = bilby_tgr.fti.prior.BNSPriorDict(filename=self.filename)
        for key, value in self.bns_prior_dict.items():
            self.prior_dict[key] = value

    def tearDown(self):
        del self.prior_dict
        del self.filename
        del self.bns_prior_dict
        del self.base_directory

    def test_create_default_prior(self):
        default = bilby_tgr.fti.prior.BNSPriorDict()
        minima = all(
            [
                self.bns_prior_dict[key].minimum == default[key].minimum
                for key in default.keys()
            ]
        )
        maxima = all(
            [
                self.bns_prior_dict[key].maximum == default[key].maximum
                for key in default.keys()
            ]
        )
        names = all(
            [
                self.bns_prior_dict[key].name == default[key].name
                for key in default.keys()
            ]
        )
        boundaries = all(
            [
                self.bns_prior_dict[key].boundary == default[key].boundary
                for key in default.keys()
            ]
        )

        self.assertTrue(all([minima, maxima, names, boundaries]))

    def test_create_from_dict(self):
        new_dict = bilby_tgr.fti.prior.BNSPriorDict(dictionary=self.prior_dict)
        self.assertDictEqual(self.bns_prior_dict, new_dict)

    def test_redundant_priors_not_in_dict_before(self):
        for prior in [
            "chirp_mass",
            "total_mass",
            "mass_ratio",
            "symmetric_mass_ratio",
            "cos_theta_jn",
            "comoving_distance",
            "redshift",
            "lambda_tilde",
            "delta_lambda_tilde",
            "dchi_kappa1",
            "dchi_kappa2",
        ]:
            self.assertTrue(self.bns_prior_dict.test_redundancy(prior))

    def test_redundant_priors_already_in_dict(self):
        for prior in [
            "mass_1",
            "mass_2",
            "chi_1",
            "chi_2",
            "theta_jn",
            "luminosity_distance",
            "lambda_1",
            "lambda_2",
            "dchi_kappaS",
            "dchi_kappaA",
        ]:
            self.assertTrue(self.bns_prior_dict.test_redundancy(prior))

    def test_correct_not_redundant_priors_masses(self):
        del self.bns_prior_dict["chirp_mass"]
        for prior in ["chirp_mass", "total_mass", "symmetric_mass_ratio"]:
            self.assertFalse(self.bns_prior_dict.test_redundancy(prior))

    def test_correct_not_redundant_priors_spin_magnitudes(self):
        del self.bns_prior_dict["chi_2"]
        self.assertFalse(self.bns_prior_dict.test_redundancy("chi_2"))

    def test_correct_not_redundant_priors_inclination(self):
        del self.bns_prior_dict["theta_jn"]
        for prior in ["theta_jn", "cos_theta_jn"]:
            self.assertFalse(self.bns_prior_dict.test_redundancy(prior))

    def test_correct_not_redundant_priors_distance(self):
        del self.bns_prior_dict["luminosity_distance"]
        for prior in ["luminosity_distance", "comoving_distance", "redshift"]:
            self.assertFalse(self.bns_prior_dict.test_redundancy(prior))

    def test_correct_not_redundant_priors_tidal(self):
        del self.bns_prior_dict["lambda_1"]
        for prior in ["lambda_1", "lambda_tilde", "delta_lambda_tilde"]:
            self.assertFalse(self.bns_prior_dict.test_redundancy(prior))

    def test_correct_not_redundant_priors_dchi_kappas(self):
        del self.bns_prior_dict["dchi_kappaA"]
        self.assertFalse(self.bns_prior_dict.test_redundancy("dchi_kappaA"))

    def test_add_unrelated_prior(self):
        self.assertFalse(self.bns_prior_dict.test_redundancy("abc"))

    def test_test_has_redundant_priors(self):
        self.assertFalse(self.bns_prior_dict.test_has_redundant_keys())
        for prior in [
            "mass_1",
            "mass_2",
            "total_mass",
            "symmetric_mass_ratio",
            "cos_theta_jn",
            "comoving_distance",
            "redshift",
            "dchi_kappa1",
            "dchi_kappa2",
        ]:
            self.bns_prior_dict[prior] = 0
            self.assertTrue(self.bns_prior_dict.test_has_redundant_keys())
            del self.bns_prior_dict[prior]

    def test_add_constraint_prior_not_redundant(self):
        self.bns_prior_dict["dchi_kappaA"] = bilby.prior.Constraint(
            minimum=20, maximum=40, name="dchi_kappaA"
        )
        self.assertFalse(self.bns_prior_dict.test_has_redundant_keys())


if __name__ == "__main__":
    unittest.main()
