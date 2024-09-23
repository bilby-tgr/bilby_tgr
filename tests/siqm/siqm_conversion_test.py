import unittest
import bilby_tgr


class TestBasicConversions(unittest.TestCase):
    def setUp(self):
        self.dQuadMon1 = 10.0
        self.dQuadMon2 = 10.0
        self.dQuadMonS = 0.5 * (self.dQuadMon1 + self.dQuadMon2)
        self.dQuadMonA = 0.5 * (self.dQuadMon1 - self.dQuadMon2)

    def test_dQuadMon1_and_dQuadMon2_to_dQuadMonS_and_dQuadMonA(self):
        (
            dQuadMonS,
            dQuadMonA,
        ) = bilby_tgr.siqm.conversion.dQuadMon1_and_dQuadMon2_to_dQuadMonS_and_dQuadMonA(
            self.dQuadMon1, self.dQuadMon2
        )
        self.assertAlmostEqual(self.dQuadMonS, dQuadMonS)
        self.assertAlmostEqual(self.dQuadMonA, dQuadMonA)

    def test_dQuadMonS_and_dQuadMonA_to_dQuadMon1_and_dQuadMon2(self):
        (
            dQuadMon1,
            dQuadMon2,
        ) = bilby_tgr.siqm.conversion.dQuadMonS_and_dQuadMonA_to_dQuadMon1_and_dQuadMon2(
            self.dQuadMonS, self.dQuadMonA
        )
        self.assertAlmostEqual(self.dQuadMon1, dQuadMon1)
        self.assertAlmostEqual(self.dQuadMon2, dQuadMon2)


class TestConvertToLALParams(unittest.TestCase):
    def setUp(self):
        self.search_keys = []
        self.parameters = dict()
        self.component_siqm_parameters = dict(dQuadMon1=10.0, dQuadMon2=10.0)
        self.all_component_pars = self.component_siqm_parameters.copy()
        self.siqm_parameters = self.component_siqm_parameters.copy()
        self.siqm_parameters[
            "dQuadMonS", "dQuadMonA"
        ] = bilby_tgr.siqm.conversion.dQuadMon1_and_dQuadMon2_to_dQuadMonS_and_dQuadMonA(
            **self.all_component_pars
        )

    def tearDown(self):
        del self.search_keys
        del self.parameters

    def siqm_convert(self):
        (
            self.parameters,
            self.added_keys,
        ) = bilby_tgr.siqm.conversion.convert_to_lal_siqm_parameters(self.parameters)

    def _conversion_to_component_siqm(self, keys):
        for key in keys:
            self.parameters[key] = self.siqm_parameters[key]
        self.siqm_convert()
        component_dict = {
            key: self.parameters[key] for key in ["dQuadMon1", "dQuadMon2"]
        }
        self.assertDictEqual(component_dict, self.component_siqm_parameters)

    def test_dQuadMon1_and_dQuadMon2_to_dQuadMonS_and_dQuadMonA(self):
        self._conversion_to_component_siqm(["dQuadMon1", "dQuadMon2"])


class TestGenerateAllParameters(unittest.TestCase):
    def setUp(self):
        self.parameters = dict(
            mass_1=36.0,
            mass_2=29.0,
            a_1=0.4,
            a_2=0.3,
            tilt_1=0.5,
            tilt_2=1.0,
            phi_12=1.7,
            phi_jl=0.3,
            luminosity_distance=2000.0,
            theta_jn=0.4,
            psi=2.659,
            phase=1.3,
            geocent_time=1126259642.413,
            ra=1.375,
            dec=-1.2108,
            dQuadMonS=10,
            dQuadMonA=0,
        )
        self.expected_bbh_keys = [
            "mass_1",
            "mass_2",
            "a_1",
            "a_2",
            "tilt_1",
            "tilt_2",
            "phi_12",
            "phi_jl",
            "luminosity_distance",
            "theta_jn",
            "psi",
            "phase",
            "geocent_time",
            "ra",
            "dec",
            "reference_frequency",
            "waveform_approximant",
            "minimum_frequency",
            "chirp_mass",
            "total_mass",
            "symmetric_mass_ratio",
            "mass_ratio",
            "iota",
            "spin_1x",
            "spin_1y",
            "spin_1z",
            "spin_2x",
            "spin_2y",
            "spin_2z",
            "phi_1",
            "phi_2",
            "chi_eff",
            "chi_1_in_plane",
            "chi_2_in_plane",
            "chi_p",
            "cos_tilt_1",
            "cos_tilt_2",
            "redshift",
            "comoving_distance",
            "mass_1_source",
            "mass_2_source",
            "chirp_mass_source",
            "total_mass_source",
        ]

        self.expected_siqm_keys = [
            "dQuadMon1",
            "dQuadMon2",
            "dQuadMonS",
            "dQuadMonA",
        ]

    def test_generate_all_siqm_parameters(self):
        new_parameters = bilby_tgr.siqm.conversion.generate_all_siqm_parameters(
            self.parameters
        )
        for key in self.expected_bbh_keys + self.expected_siqm_keys:
            self.assertIn(key, new_parameters)


if __name__ == "__main__":
    unittest.main()
