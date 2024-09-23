from bilby.gw.conversion import convert_to_lal_binary_black_hole_parameters
from bilby.gw.conversion import _generate_all_cbc_parameters


def dQuadMonS_and_dQuadMonA_to_dQuadMon1_and_dQuadMon2(dQuadMonS, dQuadMonA):
    """
    Convert symmetric and anti-symmetric combinations of spin-induced quadrupole
    moment parameters to that of binary components.

    Parameters
    ----------
    dQuadMonS: float
            symmetric combination of individual spin-induced quadrupole moment params, dQuadMon1 and dQuadMon2
    dQuadMonA: float
            anti-symmetric combination of individual spin-induced quadrupole moment params, dQuadMon1 and dQuadMon2
    Return
    ------
    dQuadMon1: float
        Spin-induced quadrupole parameter of the heavier object
    dQuadMon2: float
        Spin-induced quadrupole parameter of the lighter object

    """

    dQuadMon1 = dQuadMonS + dQuadMonA
    dQuadMon2 = dQuadMonS - dQuadMonA
    return dQuadMon1, dQuadMon2


def dQuadMon1_and_dQuadMon2_to_dQuadMonS_and_dQuadMonA(dQuadMon1, dQuadMon2):
    """
    Convert spin-induced quadrupole moment params to their symmetric and
    anti symmetric combinations.

    Parameters
    ----------
    dQuadMon1: float
        Spin-induced quadrupole parameter of the heavier object
    dQuadMon2: float
        Spin-induced quadrupole parameter of the lighter object

    Return
    ------
    dQuadMonS: float
        symmetric combination of individual spin-induced quadrupole moment params, dQuadMon1 and dQuadMon2
    dQuadMonA: float
        anti-symmetric combination of individual spin-induced quadrupole moment params, dQuadMon1 and dQuadMon2
    """

    dQuadMonS = 0.5 * (dQuadMon1 + dQuadMon2)
    dQuadMonA = 0.5 * (dQuadMon1 - dQuadMon2)
    return dQuadMonS, dQuadMonA


def convert_to_lal_siqm_parameters(parameters):
    """
    Convert parameters we have into parameters we need.

    This is defined by the parameters of bilby.source.lal_binary_black_hole()


    Mass: mass_1, mass_2
    Spin: a_1, a_2, tilt_1, tilt_2, phi_12, phi_jl
    Extrinsic: luminosity_distance, theta_jn, phase, ra, dec, geocent_time, psi

    This involves popping a lot of things from parameters.
    The keys in added_keys should be popped after evaluating the waveform.

    Parameters
    ----------
    parameters: dict
        dictionary of parameter values to convert into the required parameters

    Return
    ------
    converted_parameters: dict
        dict of the required parameters
    added_keys: list
        keys which are added to parameters during function call
    """

    converted_parameters, added_keys = convert_to_lal_binary_black_hole_parameters(
        parameters
    )

    if not any(
        [
            key in converted_parameters
            for key in ["dQuadMon1", "dQuadMon2", "dQuadMonS", "dQuadMonA"]
        ]
    ):
        converted_parameters["dQuadMon1"] = 0
        converted_parameters["dQuadMon2"] = 0
        added_keys = added_keys + ["dQuadMon1", "dQuadMon2"]
        return converted_parameters, added_keys

    if (
        "dQuadMonS" in converted_parameters.keys()
        and "dQuadMonA" in converted_parameters.keys()
    ):
        (
            converted_parameters["dQuadMon1"],
            converted_parameters["dQuadMon2"],
        ) = dQuadMonS_and_dQuadMonA_to_dQuadMon1_and_dQuadMon2(
            converted_parameters["dQuadMonS"], parameters["dQuadMonA"]
        )

    return converted_parameters, added_keys


def generate_all_siqm_parameters(sample, likelihood=None, priors=None, npool=1):
    """
    From either a single sample or a set of samples fill in all missing
    BBH parameters, in place.

    Parameters
    ==========
    sample: dict or pandas.DataFrame
        Samples to fill in with extra parameters, this may be either an
        injection or posterior samples.
    likelihood: bilby.gw.likelihood.GravitationalWaveTransient, optional
        GravitationalWaveTransient used for sampling, used for waveform and
        likelihood.interferometers.
    priors: dict, optional
        Dictionary of prior objects, used to fill in non-sampled parameters.
    """
    waveform_defaults = {
        "reference_frequency": 50.0,
        "waveform_approximant": "IMRPhenomPv2",
        "minimum_frequency": 20.0,
    }
    output_sample = _generate_all_cbc_parameters(
        sample,
        defaults=waveform_defaults,
        base_conversion=convert_to_lal_siqm_parameters,
        likelihood=likelihood,
        priors=priors,
        npool=npool,
    )
    return output_sample
