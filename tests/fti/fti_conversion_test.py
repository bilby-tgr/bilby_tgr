import unittest

import pandas as pd

import bilby
import bilby_tgr


class TestBasicConversions(unittest.TestCase):
    def setUp(self):
        self.dchi_kappaS = 150.
        self.dchi_kappaA = -50.
        self.dchi_kappa1 = 100.
        self.dchi_kappa2 = 200.

    def tearDown(self):
        del self.dchi_kappaS
        del self.dchi_kappaA
        del self.dchi_kappa1
        del self.dchi_kappa2

    def test_dkappa1_and_dkappa2_to_dkappaS_and_dkappaA(self):
        dchi_kappaS, dchi_kappaA = bilby_tgr.fti.conversion.dkappa1_and_dkappa2_to_dkappaS_and_dkappaA(
            self.dchi_kappa1, self.dchi_kappa2
        )
        self.assertAlmostEqual(self.dchi_kappaS, dchi_kappaS)
        self.assertAlmostEqual(self.dchi_kappaA, dchi_kappaA)

    def test_dkappaS_and_dkappaA_to_dkappa1_and_dkappa2(self):
        dchi_kappa1, dchi_kappa2 = bilby_tgr.fti.conversion.dkappaS_and_dkappaA_to_dkappa1_and_dkappa2(
            self.dchi_kappaS, self.dchi_kappaA
        )
        self.assertAlmostEqual(self.dchi_kappa1, dchi_kappa1)
        self.assertAlmostEqual(self.dchi_kappa2, dchi_kappa2)


class TestConvertToLALParams(unittest.TestCase):
    def setUp(self):
        self.parameters = dict()

    def tearDown(self):
        del self.parameters

    def bbh_convert(self):
        (
            self.parameters,
            self.added_keys,
        ) = bilby_tgr.fti.conversion.convert_to_lal_binary_black_hole_parameters(self.parameters)

    def test_dkappa1_and_dkappa2_to_dkappaS_and_dkappaA(self):
        self.parameters["dchi_kappa1"] = 100
        self.parameters["dchi_kappa2"] = 200
        dchi_kappaS, dchi_kappaA = bilby_tgr.fti.conversion.dkappa1_and_dkappa2_to_dkappaS_and_dkappaA(
            self.parameters["dchi_kappa1"], self.parameters["dchi_kappa2"]
        )
        self.bbh_convert()
        self.assertAlmostEqual(self.parameters["dchi_kappaS"], dchi_kappaS)
        self.assertAlmostEqual(self.parameters["dchi_kappaA"], dchi_kappaA)


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
            lambda_tilde=1000,
            delta_lambda_tilde=0,
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
            "dchi_Minus2",
            "dchi_0",
            "dchi_1",
            "dchi_2",
            "dchi_3",
            "dchi_3S",
            "dchi_3NS",
            "dchi_4",
            "dchi_4S",
            "dchi_4NS",
            "dchi_5l",
            "dchi_5lS",
            "dchi_5lNS",
            "dchi_6",
            "dchi_6S",
            "dchi_6NS",
            "dchi_6l",
            "dchi_7",
            "dchi_7S",
            "dchi_7NS",
            "dchi_kappaS",
            "dchi_kappaA",
            "dchi_kappa1",
            "dchi_kappa2",
            "f_window_div_f_peak",
            "NumCycles",
        ]
        self.expected_tidal_keys = [
            "lambda_1",
            "lambda_2",
            "lambda_tilde",
            "delta_lambda_tilde",
        ]
        self.data_frame = pd.DataFrame({
            key: [value] * 100 for key, value in self.parameters.items()
        })

    def test_generate_all_bbh_parameters(self):
        self._generate(
            bilby_tgr.fti.conversion.generate_all_bbh_parameters,
            self.expected_bbh_keys,
        )

    def test_generate_all_bns_parameters(self):
        self._generate(
            bilby_tgr.fti.conversion.generate_all_bns_parameters,
            self.expected_bbh_keys + self.expected_tidal_keys,
        )

    def _generate(self, func, expected):
        for values in [self.parameters, self.data_frame]:
            new_parameters = func(values)
            for key in expected:
                self.assertIn(key, new_parameters)

    def test_generate_bbh_parameters_with_likelihood(self):
        priors = bilby_tgr.fti.prior.BBHPriorDict()
        priors["geocent_time"] = bilby.core.prior.Uniform(0.4, 0.6)
        priors["dchi_kappaS"] = bilby.core.prior.Uniform(-500, 500)
        ifos = bilby.gw.detector.InterferometerList(["H1"])
        ifos.set_strain_data_from_power_spectral_densities(duration=1, sampling_frequency=256)
        wfg = bilby.gw.waveform_generator.WaveformGenerator(
            frequency_domain_source_model=bilby_tgr.fti.source.lal_binary_black_hole
        )
        likelihood = bilby.gw.likelihood.GravitationalWaveTransient(
            interferometers=ifos,
            waveform_generator=wfg,
            priors=priors,
            phase_marginalization=True,
            time_marginalization=True,
            reference_frame="H1L1",
        )
        self.parameters["zenith"] = 0.0
        self.parameters["azimuth"] = 0.0
        self.parameters["time_jitter"] = 0.0
        del self.parameters["ra"], self.parameters["dec"]
        self.parameters = pd.DataFrame(self.parameters, index=range(1))
        converted = bilby_tgr.fti.conversion.generate_all_bbh_parameters(
            sample=self.parameters, likelihood=likelihood, priors=priors
        )
        extra_expected = [
            "geocent_time",
            "phase",
            "H1_optimal_snr",
            "H1_matched_filter_snr",
            "ra",
            "dec",
            "dchi_kappa1",
            "dchi_kappa2",
        ]
        for key in extra_expected:
            self.assertIn(key, converted)


if __name__ == "__main__":
    unittest.main()
