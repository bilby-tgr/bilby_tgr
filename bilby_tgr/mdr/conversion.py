import bilby
import numpy as np
from bilby.gw.cosmology import get_cosmology
from scipy import constants
from bilby.core.utils.constants import speed_of_light as c, parsec as pc
from scipy.constants import physical_constants

h = physical_constants["Planck constant in eV/Hz"][0] * 1e12

# helper functions from astropy


def isiterable(obj):
    """
    Returns True if the given object is iterable.
    Copied from astropy, as location varies between versions
    """

    try:
        iter(obj)
        return True
    except TypeError:
        return False


def vectorize_if_needed(func, *x):
    """
    Helper function to vectorize functions on array inputs
    Copied here from astropy, because vectorization differs between astropy versions
    """
    if any(map(isiterable, x)):
        return np.vectorize(func)(*x)
    else:
        return func(*x)

# distance conversions - converts between distances and LIV distance measure
# https://arxiv.org/pdf/1110.2720.pdf (eq. 15) written to be of the same form as bilby.gw.conversion functions


def luminosity_distance_to_redshift_scalar(distance, cosmology=None):
    """
    copied from bilby, but non-vectorized
    """
    from astropy import units
    from astropy.cosmology import z_at_value
    cosmology = get_cosmology(cosmology)
    return z_at_value(cosmology.luminosity_distance, distance * units.Mpc).value


def redshift_to_LIV_distance_scalar(z, alpha, cosmology=None):
    """
    adapting from astropy.cosmology integral_comoving_distance
    returns LIV distance measure https://arxiv.org/pdf/1110.2720.pdf (eq. 15)

    cosmology: astropy.cosmology.FLRW, str, dict
        Description of cosmology, one of:
            None - Use DEFAULT_COSMOLOGY from bilby
            Instance of astropy.cosmology.FLRW subclass
            String with name of known Astropy cosmology, e.g., Planck13
            Dictionary with arguments required to instantiate the cosmology
            class.
    """
    from scipy.integrate import quad
    cosmology = get_cosmology(cosmology)

    def integrand(x, *args):
        return cosmology._inv_efunc_scalar(x, *args) * (1 + x)**(alpha - 2.0)

    return (quad(integrand, 0.0, z, args=cosmology._inv_efunc_scalar_args)[0] *
            cosmology.hubble_distance.value * (1 + z)**(1 - alpha))


def redshift_to_LIV_distance(z, alpha, cosmology=None):
    """
    vector version of conversion to LIV distance

    cosmology: astropy.cosmology.FLRW, str, dict
        Description of cosmology, one of:
            None - Use DEFAULT_COSMOLOGY from bilby
            Instance of astropy.cosmology.FLRW subclass
            String with name of known Astropy cosmology, e.g., Planck13
            Dictionary with arguments required to instantiate the cosmology class.
    """
    if isiterable(z):
        z = np.asarray(z)
    if isiterable(alpha):
        alpha = np.asarray(alpha)
    return vectorize_if_needed(redshift_to_LIV_distance_scalar, z, alpha, cosmology)


def luminosity_distance_to_LIV_distance_scalar(distance, alpha, cosmology=None):
    z = luminosity_distance_to_redshift_scalar(distance, cosmology)
    return redshift_to_LIV_distance_scalar(z, alpha, cosmology)


def luminosity_distance_to_LIV_distance(distance, alpha, cosmology=None):
    z = bilby.gw.conversion.luminosity_distance_to_redshift(distance, cosmology)
    return redshift_to_LIV_distance(z, alpha, cosmology)


def comoving_distance_to_LIV_distance(distance, alpha, cosmology=None):
    z = bilby.gw.conversion.comoving_distance_to_redshift(distance, cosmology)
    return redshift_to_LIV_distance(z, alpha, cosmology)


# conversions between LIV parametrizations


def A_to_lambda(A_alpha, alpha, **kwargs):
    """
    Converts lorentz violation parameter expressed in terms of lambda_eff to one expressed in terms of A_alpha

    Parameters:

    A_alpha: np.array or float
        A_alpha in units of peV^(2-alpha)
    alpha: np.array or float
        exponent of lorentz violating momentum term

    cosmology: astropy.cosmology.FLRW, str, dict
        Description of cosmology, one of:
            None - Use DEFAULT_COSMOLOGY from bilby
            Instance of astropy.cosmology.FLRW subclass
            String with name of known Astropy cosmology, e.g., Planck13
            Dictionary with arguments required to instantiate the cosmology
            class.

    Output

    lambda_true: np.array or float
        lambda in units of m
    """

    if isiterable(A_alpha):
        A_alpha = np.asarray(A_alpha)
    if isiterable(alpha):
        alpha = np.asarray(alpha)

    if ((isiterable(alpha) and (2. in alpha)) or (not isiterable(alpha) and alpha == 2.)):
        raise Exception("Error, exponent alpha has value 2, which is degenerate case")
    else:
        # hc in peV*m
        hc = constants.physical_constants["Planck constant in eV/Hz"][0] * constants.c * 1e12
        lambda_true = hc * A_alpha**(1. / (alpha - 2))
        return lambda_true


def lambda_to_A(lambda_true, alpha, **kwargs):
    """
    Parameters:

    lambda_true: np.array or float
        lambda in units of m
    alpha: np.array or float
        exponent of lorentz violating momentum term

    cosmology: astropy.cosmology.FLRW, str, dict
        Description of cosmology, one of:
            None - Use DEFAULT_COSMOLOGY from bilby
            Instance of astropy.cosmology.FLRW subclass
            String with name of known Astropy cosmology, e.g., Planck13
            Dictionary with arguments required to instantiate the cosmology
            class.

    Output

    A_alpha: np.array or float
        A_alpha in units of peV^(2-alpha)
    """

    if isiterable(lambda_true):
        lambda_true = np.asarray(lambda_true)
    if isiterable(alpha):
        alpha = np.asarray(alpha)

    if ((isiterable(alpha) and (2. in alpha)) or (not isiterable(alpha) and alpha == 2.)):
        raise Exception("Error, exponent alpha has value 2, which is degenerate case")
    else:
        # hc in peV*m
        hc = constants.physical_constants["Planck constant in eV/Hz"][0] * constants.c * 1e12
        A_alpha = (lambda_true / hc)**(alpha - 2)
        return A_alpha


def lambda_to_lambda_eff(lambda_true, luminosity_distance, alpha, cosmology=None, **kwargs):
    """
    Parameters:

    lambda_true: np.array or float
        lambda in units of m
    luminosity_distance: np.array or float
        luminosity distance in Mpc
    alpha: np.array or float
        exponent of lorentz violating momentum term

    cosmology: astropy.cosmology.FLRW, str, dict
        Description of cosmology, one of:
            None - Use DEFAULT bilby cosmology
            Instance of astropy.cosmology.FLRW subclass
            String with name of known Astropy cosmology, e.g., Planck13
            Dictionary with arguments required to instantiate the cosmology
            class.

    Output

    lambda_eff: np.array or float
        lambda_eff in units of m
    """

    if isiterable(lambda_true):
        lambda_true = np.asarray(lambda_true)
    if isiterable(luminosity_distance):
        luminosity_distance = np.asarray(luminosity_distance)
    if isiterable(alpha):
        alpha = np.asarray(alpha)

    if ((isiterable(alpha) and (2. in alpha)) or (not isiterable(alpha) and alpha == 2.)):
        raise Exception("Error, exponent alpha has value 2, which is degenerate case")
    else:
        z = bilby.gw.conversion.luminosity_distance_to_redshift(luminosity_distance, cosmology)
        D_alpha = redshift_to_LIV_distance(z, alpha, cosmology)
        lambda_eff = ((1 + z)**(1 - alpha) * luminosity_distance / D_alpha)**(1. / (2 - alpha)) * lambda_true
        return lambda_eff


def lambda_eff_to_lambda(lambda_eff, luminosity_distance, alpha, redshift=None, cosmology=None, **kwargs):
    """
    Parameters:

    lambda_eff: np.array or float
        lambda_eff in units of m
    luminosity_distance: np.array or float
        luminosity distance in Mpc
    alpha: np.array or float
        exponent of lorentz violating momentum term

    cosmology: astropy.cosmology.FLRW, str, dict
        Description of cosmology, one of:
            None - Use default bilby cosmology
            Instance of astropy.cosmology.FLRW subclass
            String with name of known Astropy cosmology, e.g., Planck13
            Dictionary with arguments required to instantiate the cosmology
            class.

    Output

    lambda_true: np.array or float
        lambda in units of m
    """

    if isiterable(lambda_eff):
        lambda_eff = np.asarray(lambda_eff)
    if isiterable(luminosity_distance):
        luminosity_distance = np.asarray(luminosity_distance)
    if isiterable(alpha):
        alpha = np.asarray(alpha)

    if ((isiterable(alpha) and (2. in alpha)) or (not isiterable(alpha) and alpha == 2.)):
        raise Exception("Error, exponent alpha has value 2, which is degenerate case")
    else:
        if redshift is not None:
            z = redshift
        else:
            z = bilby.gw.conversion.luminosity_distance_to_redshift(luminosity_distance, cosmology)
        D_alpha = redshift_to_LIV_distance(z, alpha, cosmology)
        lambda_true = ((1 + z)**(1 - alpha) * luminosity_distance / D_alpha)**(1. / (alpha - 2)) * lambda_eff
        return lambda_true


def A_to_lambda_eff(A_alpha, luminosity_distance, alpha, cosmology=None, **kwargs):
    """
    Parameters:

    A_alpha: np.array or float
        A_alpha in peV^(2-alpha)
    luminosity_distance: np.array or float
        luminosity distance in Mpc
    alpha: np.array or float
        exponent of lorentz violating momentum term

    cosmology: astropy.cosmology.FLRW, str, dict
        Description of cosmology, one of:
            None - Use default bilby cosmology
            Instance of astropy.cosmology.FLRW subclass
            String with name of known Astropy cosmology, e.g., Planck13
            Dictionary with arguments required to instantiate the cosmology
            class.

    Output

    lambda_eff: np.array or float
        lambda_eff in units of m
    """

    lambda_eff = lambda_to_lambda_eff(A_to_lambda(A_alpha, alpha), luminosity_distance, alpha, cosmology)
    return lambda_eff


def lambda_eff_to_A(lambda_eff, luminosity_distance, alpha, redshift=None, cosmology=None, **kwargs):
    """
    Parameters:

    lambda_eff: np.array or float
        lambda_eff in units of m
    luminosity_distance: np.array or float
        luminosity distance in Mpc
    alpha: np.array or float
        exponent of lorentz violating momentum term

    cosmology: astropy.cosmology.FLRW, str, dict
        Description of cosmology, one of:
            None - Use bilby cosmology
            Instance of astropy.cosmology.FLRW subclass
            String with name of known Astropy cosmology, e.g., Planck13
            Dictionary with arguments required to instantiate the cosmology
            class.

    Output

    A_alpha: np.array or float
        A_alpha in peV^(2-alpha)
    """

    A_alpha = lambda_to_A(
        lambda_eff_to_lambda(lambda_eff, luminosity_distance, alpha, redshift=redshift, cosmology=cosmology), alpha)
    return A_alpha


def log_lambda_eff_to_A(log_lambda_eff, luminosity_distance, alpha, sign_A, redshift=None, cosmology=None, **kwargs):
    return sign_A * lambda_eff_to_A(
        10.**log_lambda_eff, luminosity_distance, alpha, redshift=redshift, cosmology=cosmology, **kwargs)


# effective A


def A_to_A_eff(A_alpha, luminosity_distance, alpha, cosmology=None, **kwargs):
    """
    Parameters:

    A_alpha: np.array or float
        A_alpha in peV^(2-alpha)
    luminosity_distance: np.array or float
        luminosity distance in Mpc
    alpha: np.array or float
        exponent of lorentz violating momentum term

    cosmology: astropy.cosmology.FLRW, str, dict
        Description of cosmology, one of:
            None - Use default bilby cosmology
            Instance of astropy.cosmology.FLRW subclass
            String with name of known Astropy cosmology, e.g., Planck13
            Dictionary with arguments required to instantiate the cosmology
            class.

    Output

    A_eff: np.array or float
        A_eff in units of peV^(2-alpha)
    """
    z = bilby.gw.conversion.luminosity_distance_to_redshift(luminosity_distance, cosmology)
    LIV_distance = redshift_to_LIV_distance(z, alpha, cosmology=None)
    return LIV_distance * (1. + z)**(alpha - 1.) / luminosity_distance * A_alpha


def A_eff_to_A(A_eff, luminosity_distance, alpha, redshift=None, cosmology=None, **kwargs):
    """
    Parameters:

    A_eff: np.array or float
        A_eff in peV^(2-alpha)
    luminosity_distance: np.array or float
        luminosity distance in Mpc
    alpha: np.array or float
        exponent of lorentz violating momentum term

    redshift: np.array
        if given, skips redshift calculation step

    cosmology: astropy.cosmology.FLRW, str, dict
        Description of cosmology, one of:
            None - Use default bilby cosmology
            Instance of astropy.cosmology.FLRW subclass
            String with name of known Astropy cosmology, e.g., Planck13
            Dictionary with arguments required to instantiate the cosmology
            class.

    Output

    A_alpha: np.array or float
        A_alpha in units of peV^(2-alpha)
    """
    if redshift is not None:
        z = redshift
    else:
        z = bilby.gw.conversion.luminosity_distance_to_redshift(luminosity_distance, cosmology)
    LIV_distance = redshift_to_LIV_distance(z, alpha, cosmology=None)
    return luminosity_distance * (1. + z)**(1. - alpha) / LIV_distance * A_eff


def mass_graviton_to_A(mass_graviton, **kwargs):
    """
    Parameters:

    mass_graviton: np.array or float
        mass of the graviton in peV

    Output

    A: np.array or float
        A in units of peV^(2-alpha)
    """
    return mass_graviton**2


def A_to_mass_graviton(A_alpha, **kwargs):
    """
    Parameters:

    A: np.array or float
        A in units of peV^(2-alpha)

    Output

    mass_graviton: np.array or float
        mass of the graviton in peV
    """
    return A_alpha**0.5


def mdr_phase_to_A_eff(mdr_phase, luminosity_distance, **kwargs):
    """
    Case alpha = 1 corresponds to constant phase. Conversion for this parametrization
    Parameters:

    mdr_phase: np.array or float
        Constant phase offset introduced by dispersion
    luminosity_distance: np.array or float
        luminosity distance in Mpc

    Output

    A_eff: np.array or float
        A_eff in units of peV^(2-alpha)
    """
    distance = luminosity_distance * pc * 1e6
    A_eff = - mdr_phase * c * h / distance / np.pi
    return A_eff


def A_eff_to_mdr_phase(A_eff, luminosity_distance, **kwargs):
    """
    Case alpha = 1 corresponds to constant phase. Conversion for this parametrization
    Parameters:

    A_eff: np.array or float
        A_eff in units of peV^(2-alpha)
    luminosity_distance: np.array or float
        luminosity distance in Mpc

    Output:

    mdr_phase: np.array or float
        Constant phase offset introduced by dispersion
    """
    distance = luminosity_distance * pc * 1e6
    mdr_phase = -(np.pi * distance * A_eff / h / c)
    return mdr_phase


def generate_all_bbh_parameters(sample, likelihood=None, priors=None, npool=1):
    """
    adapted from bilby
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
    output_sample = bilby.gw.conversion.generate_all_bbh_parameters(sample, likelihood, priors, npool)
    output_sample = generate_liv_parameters(output_sample)
    return output_sample


def generate_all_bns_parameters(sample, likelihood=None, priors=None, npool=1):
    """
    copied from bilby
    From either a single sample or a set of samples fill in all missing
    BNS parameters, in place.

    Since we assume BNS waveforms are aligned, component spins won't be
    calculated.

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
    npool: int, (default=1)
        If given, perform generation (where possible) using a multiprocessing pool
    """
    output_sample = bilby.gw.conversion.generate_all_bns_parameters(sample, likelihood, priors, npool)
    output_sample = generate_liv_parameters(output_sample)
    return output_sample


def generate_liv_parameters(sample, cosmology=None):
    output_sample = sample.copy()
    redshift = None
    if 'redshift' in output_sample:
        redshift = output_sample['redshift']
    if 'alpha' in output_sample:
        if 'mass_graviton_eff' in output_sample:
            output_sample['A_eff'] = mass_graviton_to_A(output_sample["mass_graviton_eff"])
        elif 'mdr_phase' in output_sample:
            output_sample['A_eff'] = mdr_phase_to_A_eff(
                output_sample["mdr_phase"], output_sample['luminosity_distance'])
        if 'A_eff' in output_sample:
            output_sample['A_alpha'] = A_eff_to_A(
                A_eff=output_sample['A_eff'], luminosity_distance=output_sample['luminosity_distance'],
                alpha=output_sample['alpha'], redshift=redshift, cosmology=cosmology)
        elif 'log_lambda_eff' in output_sample and 'sign_A' in output_sample:
            output_sample['A_alpha'] = log_lambda_eff_to_A(
                log_lambda_eff=output_sample['log_lambda_eff'],
                luminosity_distance=output_sample['luminosity_distance'],
                alpha=output_sample['alpha'], sign_A=output_sample['sign_A'],
                redshift=redshift, cosmology=cosmology)

        if np.all(output_sample['alpha'] == 0.0):  # checks if alpha = 0, so if we consider graviton case
            if 'mass_graviton' not in output_sample:
                output_sample['mass_graviton'] = A_to_mass_graviton(output_sample["A_alpha"])
        if np.all(output_sample['alpha'] == 1.0):  # checks if alpha = 1, so constant phase case
            if 'mdr_phase' not in output_sample and 'A_eff' in output_sample:
                output_sample['mdr_phase'] = A_eff_to_mdr_phase(
                    output_sample["A_eff"], output_sample['luminosity_distance'])
    if 'A_alpha' not in output_sample:
        print('Could not generate A_alpha parameter')
    return output_sample


def convert_to_lal_binary_black_hole_parameters(parameters):
    """
    Convert parameters we have into parameters we need.

    This is defined by the parameters of bilby.source.lal_binary_black_hole()


    Mass: mass_1, mass_2
    Spin: a_1, a_2, tilt_1, tilt_2, phi_12, phi_jl
    Extrinsic: luminosity_distance, theta_jn, phase, ra, dec, geocent_time, psi
    MDR: A_eff, alpha

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

    converted_parameters, added_keys = bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters(parameters)
    if 'A_eff' not in converted_parameters:
        if 'mass_graviton_eff' in converted_parameters:
            converted_parameters['A_eff'] = mass_graviton_to_A(converted_parameters["mass_graviton_eff"])
            added_keys.append('A_eff')
        elif 'log_lambda_eff' in converted_parameters and 'sign_A' in converted_parameters:
            converted_parameters['A_eff'] = converted_parameters['sign_A'] * lambda_to_A(
                10**converted_parameters['log_lambda_eff'], converted_parameters['alpha'])
            added_keys.append('A_eff')
        elif 'mdr_phase' in converted_parameters:
            converted_parameters['A_eff'] = mdr_phase_to_A_eff(converted_parameters["mdr_phase"],
                                                               converted_parameters["luminosity_distance"])
            added_keys.append('A_eff')
    return converted_parameters, added_keys


def convert_to_lal_binary_neutron_star_parameters(parameters):
    """
    Convert parameters we have into parameters we need.

    This is defined by the parameters of bilby.source.lal_binary_black_hole()


    Mass: mass_1, mass_2
    Spin: a_1, a_2, tilt_1, tilt_2, phi_12, phi_jl
    Extrinsic: luminosity_distance, theta_jn, phase, ra, dec, geocent_time, psi
    Tidal: lambda_1, lamda_2, lambda_tilde, delta_lambda_tilde
    MDR: A_eff, alpha

    This involves popping a lot of things from parameters.
    The keys in added_keys should be popped after evaluating the waveform.
    For details on tidal parameters see https://arxiv.org/pdf/1402.5156.pdf.

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

    converted_parameters, added_keys = bilby.gw.conversion.convert_to_lal_binary_neutron_star_parameters(parameters)
    if 'A_eff' not in converted_parameters:
        if 'mass_graviton_eff' in converted_parameters:
            converted_parameters['A_eff'] = mass_graviton_to_A(converted_parameters["mass_graviton_eff"])
            added_keys.append('A_eff')
        elif 'log_lambda_eff' in converted_parameters and 'sign_A' in converted_parameters:
            converted_parameters['A_eff'] = converted_parameters['sign_A'] * lambda_to_A(
                10**converted_parameters['log_lambda_eff'], converted_parameters['alpha'])
            added_keys.append('A_eff')
        elif 'mdr_phase' in converted_parameters:
            converted_parameters['A_eff'] = mdr_phase_to_A_eff(converted_parameters["mdr_phase"],
                                                               converted_parameters["luminosity_distance"])
            added_keys.append('A_eff')
    return converted_parameters, added_keys
