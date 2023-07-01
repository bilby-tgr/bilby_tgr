import unittest

import numpy as np
import pandas as pd


import bilby
import bilby_tgr
from bilby_tgr.ssb import conversion


class TestBasicConversionsConsistency(unittest.TestCase):
    # test in back and forth conversion is consistent
    def setUp(self):
        self.ssb_coeff = np.array([0., 0.5, 1.5, 2.5, 3.0, 3.5, 4.0]) * 1000
        self.distance = 1000

    def tearDown(self):
        del self.ssb_coeff
        del self.distance

    def test_kv5_ssb_coeff_conversion(self):
        redshift = conversion.luminosity_distance_to_redshift(self.distance)
        kv5 = conversion.kv5_from_ssb_coeff(self.ssb_coeff, redshift)
        ssb_coeff = conversion.ssb_coeff_from_kv5_scalar(kv5, redshift)
        self.assertTrue(
            all(abs(self.ssb_coeff - ssb_coeff) <= 1e-5 * self.ssb_coeff)
        )


class TestBasicConversions(unittest.TestCase):
    def setUp(self):
        bilby.gw.cosmology.set_cosmology("Planck15")
        self.ssb_coeff = 5510.6804959978945
        self.kv5 = 3.103992107409165e-14
        self.distance = 2000.
        self.z = 0.3633349489594506

    def tearDown(self):
        del self.ssb_coeff
        del self.kv5
        del self.distance
        del self.z

    def test_luminosity_distance_to_redshift_scalar(self):
        z = conversion.luminosity_distance_to_redshift_scalar(self.distance)
        self.assertTrue(
            abs(self.z - z) <= 1e-5 * self.z
        )

    def test_kv5_from_ssb_coeff_scalar(self):
        kv5 = conversion.kv5_from_ssb_coeff_scalar(self.ssb_coeff, self.z)
        self.assertTrue(
            abs(self.kv5 - kv5) <= 1e-5 * self.kv5
        )

    def test_ssb_coeff_from_kv5_scalar(self):
        ssb_coeff = conversion.ssb_coeff_from_kv5_scalar(self.kv5, self.z)
        self.assertTrue(
            abs(self.ssb_coeff - ssb_coeff) <= 1e-5 * self.ssb_coeff
        )

    def test_luminosity_distance_to_redshift(self):
        z = conversion.luminosity_distance_to_redshift(self.distance)
        self.assertTrue(
            abs(self.z - z) <= 1e-5 * self.z
        )

    def test_kv5_from_ssb_coeff(self):
        kv5 = conversion.kv5_from_ssb_coeff(self.ssb_coeff, self.z)
        self.assertTrue(
            abs(self.kv5 - kv5) <= 1e-5 * self.kv5
        )

    def test_ssb_coeff_from_kv5(self):
        ssb_coeff = conversion.ssb_coeff_from_kv5(self.kv5, self.z)
        self.assertTrue(
            abs(self.ssb_coeff - ssb_coeff) <= 1e-5 * self.ssb_coeff
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
            ssb_coeff=8000,
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
            "ssb_coeff",
            "kv5",
        ]

        self.data_frame = pd.DataFrame({
            key: [value] * 100 for key, value in self.parameters.items()
        })

    def test_generate_all_bbh_parameters(self):
        self._generate(
            conversion.generate_all_bbh_parameters,
            self.expected_bbh_keys,
        )

    def _generate(self, func, expected):
        for values in [self.parameters, self.data_frame]:
            new_parameters = func(values)
            for key in expected:
                self.assertIn(key, new_parameters)

    def test_generate_bbh_parameters_with_likelihood(self):
        priors = bilby.gw.prior.BBHPriorDict()
        priors["geocent_time"] = bilby.core.prior.Uniform(0.4, 0.6)
        priors["ssb_dim"] = 5
        priors["ssb_coeff"] = bilby.core.prior.Uniform(1e-17, 10000)
        ifos = bilby.gw.detector.InterferometerList(["H1"])
        ifos.set_strain_data_from_power_spectral_densities(duration=1, sampling_frequency=256)
        wfg = bilby.gw.waveform_generator.WaveformGenerator(
            frequency_domain_source_model=bilby_tgr.ssb.source.lal_binary_black_hole
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
