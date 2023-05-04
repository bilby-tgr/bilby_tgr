import unittest

import numpy as np
import pandas as pd


import bilby
import bilby_tgr
from bilby_tgr.mdr import conversion


class TestBasicConversionsConsistency(unittest.TestCase):
    # test in back and forth conversion is consistent
    def setUp(self):
        self.alpha = np.array([0., 0.5, 1.5, 2.5, 3.0, 3.5, 4.0])
        self.A_alpha = 1e-20
        self.distance = 1000

    def tearDown(self):
        del self.alpha
        del self.A_alpha
        del self.distance

    def test_A_lambda_conversion(self):
        lam = conversion.A_to_lambda(self.A_alpha, self.alpha)
        A = conversion.lambda_to_A(lam, self.alpha)
        self.assertTrue(
            all(abs(self.A_alpha - A) < 1e-5 * self.A_alpha)
        )

    def test_lambda_lambda_eff_conversion(self):
        lam = conversion.A_to_lambda(self.A_alpha, self.alpha)
        lam_eff = conversion.lambda_to_lambda_eff(lam, self.distance, self.alpha)
        lam_2 = conversion.lambda_eff_to_lambda(lam_eff, self.distance, self.alpha)
        self.assertTrue(
            all(abs(lam - lam_2) < 1e-5 * lam)
        )

    def test_A_lambda_eff_conversion(self):
        lam_eff = conversion.A_to_lambda_eff(self.A_alpha, self.distance, self.alpha)
        A = conversion.lambda_eff_to_A(lam_eff, self.distance, self.alpha)
        self.assertTrue(
            all(abs(self.A_alpha - A) < 1e-5 * self.A_alpha)
        )

    def test_A_A_eff_conversion(self):
        A_eff = conversion.A_to_A_eff(self.A_alpha, self.distance, self.alpha)
        A = conversion.A_eff_to_A(A_eff, self.distance, self.alpha)
        self.assertTrue(
            all(abs(self.A_alpha - A) < 1e-5 * self.A_alpha)
        )


class TestBasicConversions(unittest.TestCase):
    def setUp(self):
        bilby.gw.cosmology.set_cosmology("Planck15")
        self.alpha = 0.0
        self.A_alpha = 1e-20
        self.distance = 2000.
        self.z = 0.3633349489594506
        self.liv_distance = 1481.6623139542735
        self.comoving_distance = 1466.990920476916
        self.lam = 1.239841984055037e+16
        self.lambda_eff = 1.681930332907409e+16
        self.log_lambda_eff = 16.22580800295982
        self.A_eff = 5.43396292703101e-21
        self.sign_A = +1.

    def tearDown(self):
        del self.alpha
        del self.A_alpha
        del self.distance
        del self.z
        del self.liv_distance
        del self.comoving_distance
        del self.lam
        del self.lambda_eff
        del self.log_lambda_eff
        del self.A_eff
        del self.sign_A

    def test_luminosity_distance_to_redshift_scalar(self):
        z = conversion.luminosity_distance_to_redshift_scalar(self.distance)
        self.assertTrue(
            abs(self.z - z) < 1e-5 * self.z
        )

    def test_redshift_to_LIV_distance_scalar(self):
        liv_distance = conversion.redshift_to_LIV_distance_scalar(self.z, self.alpha)
        self.assertTrue(
            abs(self.liv_distance - liv_distance) < 1e-5 * self.liv_distance
        )

    def test_redshift_to_LIV_distance(self):
        liv_distance = conversion.redshift_to_LIV_distance(self.z, self.alpha)
        self.assertTrue(
            abs(self.liv_distance - liv_distance) < 1e-5 * self.liv_distance
        )

    def test_luminosity_distance_to_LIV_distance_scalar(self):
        liv_distance = conversion.luminosity_distance_to_LIV_distance_scalar(self.distance, self.alpha)
        self.assertTrue(
            abs(self.liv_distance - liv_distance) < 1e-5 * self.liv_distance
        )

    def test_luminosity_distance_to_LIV_distance(self):
        liv_distance = conversion.luminosity_distance_to_LIV_distance(self.distance, self.alpha)
        self.assertTrue(
            abs(self.liv_distance - liv_distance) < 1e-5 * self.liv_distance
        )

    def test_comoving_distance_to_LIV_distance(self):
        liv_distance = conversion.comoving_distance_to_LIV_distance(self.comoving_distance, self.alpha)
        self.assertTrue(
            abs(self.liv_distance - liv_distance) < 1e-5 * self.liv_distance
        )

    def test_A_to_lambda(self):
        lam = conversion.A_to_lambda(self.A_alpha, self.alpha)
        self.assertTrue(
            abs(self.lam - lam) < 1e-5 * self.lam
        )

    def test_lambda_to_A(self):
        A_alpha = conversion.lambda_to_A(self.lam, self.alpha)
        self.assertTrue(
            abs(self.A_alpha - A_alpha) < 1e-5 * self.A_alpha
        )

    def test_lambda_to_lambda_eff(self):
        lambda_eff = conversion.lambda_to_lambda_eff(self.lam, self.distance, self.alpha)
        self.assertTrue(
            abs(self.lambda_eff - lambda_eff) < 1e-5 * self.lambda_eff
        )

    def test_lambda_eff_to_lambda(self):
        lam = conversion.lambda_eff_to_lambda(self.lambda_eff, self.distance, self.alpha)
        self.assertTrue(
            abs(self.lam - lam) < 1e-5 * self.lam
        )

    def test_lambda_eff_to_lambda_with_redshift(self):
        lam = conversion.lambda_eff_to_lambda(self.lambda_eff, self.distance, self.alpha, self.z)
        self.assertTrue(
            abs(self.lam - lam) < 1e-5 * self.lam
        )

    def test_A_to_lambda_eff(self):
        lambda_eff = conversion.A_to_lambda_eff(self.A_alpha, self.distance, self.alpha)
        self.assertTrue(
            abs(self.lambda_eff - lambda_eff) < 1e-5 * self.lambda_eff
        )

    def test_lambda_eff_to_A(self):
        A_alpha = conversion.lambda_eff_to_A(self.lambda_eff, self.distance, self.alpha)
        self.assertTrue(
            abs(self.A_alpha - A_alpha) < 1e-5 * self.A_alpha
        )

    def test_lambda_eff_to_A_with_redshift(self):
        A_alpha = conversion.lambda_eff_to_A(self.lambda_eff, self.distance, self.alpha, self.z)
        self.assertTrue(
            abs(self.A_alpha - A_alpha) < 1e-5 * self.A_alpha
        )

    def test_log_lambda_eff_to_A(self):
        A_alpha = conversion.log_lambda_eff_to_A(self.log_lambda_eff, self.distance, self.alpha, self.sign_A)
        self.assertTrue(
            abs(self.A_alpha - A_alpha) < 1e-5 * self.A_alpha
        )

    def test_log_lambda_eff_to_A_with_redshift(self):
        A_alpha = conversion.log_lambda_eff_to_A(self.log_lambda_eff, self.distance, self.alpha, self.sign_A, self.z)
        self.assertTrue(
            abs(self.A_alpha - A_alpha) < 1e-5 * self.A_alpha
        )

    def test_A_to_A_eff(self):
        A_eff = conversion.A_to_A_eff(self.A_alpha, self.distance, self.alpha)
        self.assertTrue(
            abs(self.A_eff - A_eff) < 1e-5 * self.A_eff
        )

    def test_A_eff_to_A(self):
        A_alpha = conversion.A_eff_to_A(self.A_eff, self.distance, self.alpha)
        self.assertTrue(
            abs(self.A_alpha - A_alpha) < 1e-5 * self.A_alpha
        )

    def test_A_eff_to_A_with_redshift(self):
        A_alpha = conversion.A_eff_to_A(self.A_eff, self.distance, self.alpha, self.z)
        self.assertTrue(
            abs(self.A_alpha - A_alpha) < 1e-5 * self.A_alpha
        )


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
            alpha=0.0,
            A_eff=1e-20,
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
            "alpha",
            "A_eff",
            "A_alpha",
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
            conversion.generate_all_bbh_parameters,
            self.expected_bbh_keys,
        )

    def test_generate_all_bns_parameters(self):
        self._generate(
            conversion.generate_all_bns_parameters,
            self.expected_bbh_keys + self.expected_tidal_keys,
        )

    def _generate(self, func, expected):
        for values in [self.parameters, self.data_frame]:
            new_parameters = func(values)
            for key in expected:
                self.assertIn(key, new_parameters)

    def test_generate_bbh_parameters_with_likelihood(self):
        priors = bilby.gw.prior.BBHPriorDict()
        priors["geocent_time"] = bilby.core.prior.Uniform(0.4, 0.6)
        priors["alpha"] = 0.0
        priors["A_eff"] = bilby.core.prior.Uniform(-1e-18, 1e-18)
        ifos = bilby.gw.detector.InterferometerList(["H1"])
        ifos.set_strain_data_from_power_spectral_densities(duration=1, sampling_frequency=256)
        wfg = bilby.gw.waveform_generator.WaveformGenerator(
            frequency_domain_source_model=bilby_tgr.mdr.source.lal_binary_black_hole
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
        converted = conversion.generate_all_bbh_parameters(
            sample=self.parameters, likelihood=likelihood, priors=priors
        )
        extra_expected = [
            "geocent_time",
            "phase",
            "H1_optimal_snr",
            "H1_matched_filter_snr",
            "ra",
            "dec",
        ]
        for key in extra_expected:
            self.assertIn(key, converted)


if __name__ == "__main__":
    unittest.main()
