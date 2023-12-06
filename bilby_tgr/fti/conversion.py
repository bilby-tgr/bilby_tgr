import bilby.gw.conversion as conversion
from bilby.core.utils import logger


def dkappa1_and_dkappa2_to_dkappaS_and_dkappaA(dchi_kappa1, dchi_kappa2):
    """
    Convert the spin-induced quadrupole moment of the individual objects to
    their symmetric and anti-symmetric combinations.

    Parameters
    ----------
    dchi_kappa1: float
        Spin-induced quadrupole moment of the heavier object
    dchi_kappa2: float
        Spin-induced quadrupole moment of the lighter object

    Return
    ------
    dchi_kappaS: float
        symmetric combination of individual spin-induced quadrupole moments
    dchi_kappaA: float
        anti-symmetric combination of individual spin-induced quadrupole moments
    """

    dchi_kappaS = (dchi_kappa1 + dchi_kappa2) / 2.
    dchi_kappaA = (dchi_kappa1 - dchi_kappa2) / 2.
    return dchi_kappaS, dchi_kappaA


def dkappaS_and_dkappaA_to_dkappa1_and_dkappa2(dchi_kappaS, dchi_kappaA):
    """
    Convert the symmetric and anti-symmetric combinations of spin-induced quadrupole moments to
    the spin-induced quadrupole moments of the individual objects.

    Parameters
    ----------
    dchi_kappaS: float
        symmetric combination of individual spin-induced quadrupole moments
    dchi_kappaA: float
        anti-symmetric combination of individual spin-induced quadrupole moments

    Return
    ------
    dchi_kappa1: float
        Spin-induced quadrupole moment of the heavier object
    dchi_kappa2: float
        Spin-induced quadrupole moment of the lighter object

    """

    dchi_kappa1 = dchi_kappaS + dchi_kappaA
    dchi_kappa2 = dchi_kappaS - dchi_kappaA
    return dchi_kappa1, dchi_kappa2


def generate_fti_parameters(sample, added_keys=[]):
    """
    Generate all SIQM parameters for FTI SIQM test:
        dchi_kappaS, dchi_kappaA, dchi_kappa1, dchi_kappa2,

    Parameters
    ==========
    sample: dict, pandas.DataFrame

    Returns
    =======
    output_sample: dict, pandas.DataFrame
        Updated sample
    """

    output_sample = sample.copy()

    if 'f_window_div_f_peak' not in output_sample.keys():
        output_sample['f_window_div_f_peak'] = 1.0
        added_keys = added_keys + ['f_window_div_f_peak']

    if 'NumCycles' not in output_sample.keys():
        output_sample['NumCycles'] = 1
        added_keys = added_keys + ['NumCycles']

    if not any([key in output_sample for key in
                ['dchi_kappaS', 'dchi_kappaA', 'dchi_kappa1', 'dchi_kappa2']]):
        return output_sample, added_keys

    if (('dchi_kappaS' in output_sample.keys() or 'dchi_kappaA' in output_sample.keys()) and
            ('dchi_kappa1' in output_sample.keys() or 'dchi_kappa2' in output_sample.keys())):
        logger.disabled = False
        logger.warning("Incompatible dkappas provided. This may lead to unexpected behaviour. \
                        Provide (dchi_kappaS and dchi_kappaA) xor (dchi_kappa1 and dchi_kappa2) instead.")

    if 'dchi_kappaS' not in output_sample.keys() and 'dchi_kappaA' not in output_sample.keys():
        if 'dchi_kappa1' not in output_sample.keys():
            output_sample['dchi_kappa1'] = 0
            added_keys = added_keys + ['dchi_kappa1']
        if 'dchi_kappa2' not in output_sample.keys():
            output_sample['dchi_kappa2'] = 0
            added_keys = added_keys + ['dchi_kappa2']
        output_sample['dchi_kappaS'], output_sample['dchi_kappaA'] =\
            dkappa1_and_dkappa2_to_dkappaS_and_dkappaA(
                output_sample['dchi_kappa1'], output_sample['dchi_kappa2'])
        added_keys = added_keys + ['dchi_kappaS', 'dchi_kappaA']

    if 'dchi_kappa1' not in output_sample.keys() and 'dchi_kappa2' not in output_sample.keys():
        if 'dchi_kappaS' not in output_sample.keys():
            output_sample['dchi_kappaS'] = 0
            added_keys = added_keys + ['dchi_kappaS']
        if 'dchi_kappaA' not in output_sample.keys():
            output_sample['dchi_kappaA'] = 0
            added_keys = added_keys + ['dchi_kappaA']
        output_sample['dchi_kappa1'], output_sample['dchi_kappa2'] =\
            dkappaS_and_dkappaA_to_dkappa1_and_dkappa2(
                output_sample['dchi_kappaS'], output_sample['dchi_kappaA'])
        added_keys = added_keys + ['dchi_kappa1', 'dchi_kappa2']

    return output_sample, added_keys


def convert_to_lal_binary_black_hole_parameters(parameters):
    """
    (Adapted from bilby)
    Convert parameters we have into parameters we need.

    This is defined by the parameters of bilby_tgr.fti.source.lal_binary_black_hole()

    Mass: mass_1, mass_2
    Spin: a_1, a_2, tilt_1, tilt_2, phi_12, phi_jl
    Extrinsic: luminosity_distance, theta_jn, phase, ra, dec, geocent_time, psi
    FTI: dchi_Minus2, dchi_0, dchi_1, dchi_2, dchi_3, dchi_3S, dchi_3NS, dchi_4, dchi_4S, dchi_4NS,
        dchi_5l, dchi_5lS, dchi_5lNS, dchi_6, dchi_6S, dchi_6NS, dchi_6l, dchi_7, dchi_7S, dchi_7NS,
        dchi_kappaS, dchi_kappaA, f_window_div_f_peak, NumCycles

    This involves popping a lot of things from parameters.
    The keys in added_keys should be popped after evaluating the waveform.

    Parameters
    ==========
    parameters: dict
        dictionary of parameter values to convert into the required parameters

    Returns
    =======
    converted_parameters: dict
        dict of the required parameters
    added_keys: list
        keys which are added to parameters during function call
    """

    converted_parameters, added_keys = conversion.convert_to_lal_binary_black_hole_parameters(parameters)
    converted_parameters, added_keys = generate_fti_parameters(converted_parameters, added_keys)

    FTI_params = ['dchi_Minus2', 'dchi_0', 'dchi_1', 'dchi_2', 'dchi_3', 'dchi_3S', 'dchi_3NS',
                  'dchi_4', 'dchi_4S', 'dchi_4NS', 'dchi_5l', 'dchi_5lS', 'dchi_5lNS', 'dchi_6',
                  'dchi_6S', 'dchi_6NS', 'dchi_6l', 'dchi_7', 'dchi_7S', 'dchi_7NS',
                  'dchi_kappaS', 'dchi_kappaA', 'dchi_kappa1', 'dchi_kappa2']
    for key in FTI_params:
        if key not in converted_parameters.keys():
            converted_parameters[key] = 0
            added_keys = added_keys + [key]

    return converted_parameters, added_keys


def generate_all_bbh_parameters(sample, likelihood=None, priors=None, npool=1):
    """
    (Adapted from bilby)
    From either a single sample or a set of samples fill in all missing BBH parameters

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

    output_sample = conversion.generate_all_bbh_parameters(sample, likelihood, priors, npool)
    output_sample, _ = generate_fti_parameters(output_sample)

    return output_sample


def convert_to_lal_binary_neutron_star_parameters(parameters):
    """
    (Adapted from bilby)
    Convert parameters we have into parameters we need.

    This is defined by the parameters of bilby_tgr.fti.source.lal_binary_neutron_star()

    Mass: mass_1, mass_2
    Spin: a_1, a_2, tilt_1, tilt_2, phi_12, phi_jl
    Extrinsic: luminosity_distance, theta_jn, phase, ra, dec, geocent_time, psi
    Tidal: lambda_1, lamda_2
    FTI: dchi_Minus2, dchi_0, dchi_1, dchi_2, dchi_3, dchi_3S, dchi_3NS, dchi_4, dchi_4S, dchi_4NS,
        dchi_5l, dchi_5lS, dchi_5lNS, dchi_6, dchi_6S, dchi_6NS, dchi_6l, dchi_7, dchi_7S, dchi_7NS,
        dchi_kappaS, dchi_kappaA, f_window_div_f_peak, NumCycles

    This involves popping a lot of things from parameters.
    The keys in added_keys should be popped after evaluating the waveform.

    Parameters
    ==========
    parameters: dict
        dictionary of parameter values to convert into the required parameters

    Returns
    =======
    converted_parameters: dict
        dict of the required parameters
    added_keys: list
        keys which are added to parameters during function call
    """

    converted_parameters, added_keys = conversion.convert_to_lal_binary_neutron_star_parameters(parameters)
    converted_parameters, added_keys = generate_fti_parameters(converted_parameters, added_keys)

    FTI_params = ['dchi_Minus2', 'dchi_0', 'dchi_1', 'dchi_2', 'dchi_3', 'dchi_3S', 'dchi_3NS',
                  'dchi_4', 'dchi_4S', 'dchi_4NS', 'dchi_5l', 'dchi_5lS', 'dchi_5lNS', 'dchi_6',
                  'dchi_6S', 'dchi_6NS', 'dchi_6l', 'dchi_7', 'dchi_7S', 'dchi_7NS',
                  'dchi_kappaS', 'dchi_kappaA', 'dchi_kappa1', 'dchi_kappa2']
    for key in FTI_params:
        if key not in converted_parameters.keys():
            converted_parameters[key] = 0
            added_keys = added_keys + [key]

    return converted_parameters, added_keys


def generate_all_bns_parameters(sample, likelihood=None, priors=None, npool=1):
    """
    (Adapted from bilby)
    From either a single sample or a set of samples fill in all missing BBH parameters

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

    output_sample = conversion.generate_all_bns_parameters(sample, likelihood, priors, npool)
    output_sample, _ = generate_fti_parameters(output_sample)

    return output_sample
