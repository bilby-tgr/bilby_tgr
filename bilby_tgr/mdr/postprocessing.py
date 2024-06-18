from .conversion import generate_liv_parameters, redshift_to_LIV_distance
from bilby.gw.conversion import luminosity_distance_to_redshift
import numpy as np


def get_reweighting_weights(posterior, prior_initial, prior_final='default', cosmology=None):
    """
    From Jacobians of transformations, calculates weight for reweighting

    Parameters
    ==========
    posterior: pandas.DataFrame
        Bilby LIV posterior to transform
    prior_initial: string, one of A_eff, log_A_eff, log_lambda_eff
        String specifying prior in which LIV parameter was sampled
    prior_final: string, one of uniform, mass_uniform, default
        String specifying desired prior on A_alpha for reweighting
    cosmology: astropy.cosmology.FLRW, str, dict
        Description of cosmology, one of:
            None - Use default bilby cosmology
            Instance of astropy.cosmology.FLRW subclass
            String with name of known Astropy cosmology, e.g., Planck13
            Dictionary with arguments required to instantiate the cosmology
            class.


    Output

    weighted: pandas.DataFrame
        relative weights of the samples

    """

    assert ((prior_initial in ['A_eff', 'log_A_eff', 'log_lambda_eff', 'mass_graviton_eff', 'mdr_phase']) and
            (prior_final in ['uniform', 'mass_uniform', 'default'])), "Unknown prior"
    posterior_reweighted = posterior.copy()
    if prior_final == 'default':
        if prior_initial == 'mass_graviton_eff':
            prior_final = 'mass_uniform'
        else:
            prior_final = 'uniform'
    if 'A_alpha' not in posterior_reweighted:
        posterior_reweighted = generate_liv_parameters(posterior_reweighted)
    if 'A_alpha' not in posterior_reweighted:
        print('ERROR: Could not find LIV parameter in the posterior')

    # get distance measures for reweighting
    if 'redshift' in posterior_reweighted:
        z = posterior_reweighted['redshift']
    else:
        z = luminosity_distance_to_redshift(
            posterior_reweighted['luminosity_distance'], cosmology=cosmology)
    LIV_distance = redshift_to_LIV_distance(z, posterior_reweighted['alpha'], cosmology=cosmology)

    # distance ratio often useful in reweighting formulas
    distance_ratio = (posterior_reweighted['luminosity_distance'] / LIV_distance *
                      (1 + z)**(1 - posterior_reweighted['alpha']))

    # weights for reweighting from jacobian of the transformation (up to multiplicative constant)
    if prior_final == 'uniform':
        if prior_initial == 'A_eff':
            weights = distance_ratio
        elif prior_initial == 'mass_graviton_eff':
            weights = np.abs(distance_ratio * posterior_reweighted['A_alpha'])**0.5
            # discard negative mass samples (shouldn't be there anyway)
            weights[posterior_reweighted['A_alpha'] < 0] = 0
        elif prior_initial == 'mdr_phase':
            weights = distance_ratio / posterior_reweighted['luminosity_distance']
        else:
            weights = np.abs(posterior_reweighted['A_alpha'])
    elif prior_final == 'mass_uniform':
        assert np.all(posterior_reweighted['alpha'] == 0.0), 'mass_uniform prior is compatible only with alpha=0'
        if prior_initial == 'mass_graviton_eff':
            weights = distance_ratio**0.5
        elif prior_initial == 'A_eff':
            weights = distance_ratio**0.5 / np.abs(posterior_reweighted['A_eff'])**0.5
            # discard negative mass samples (shouldn't be there anyway)
            weights[posterior_reweighted['A_alpha'] < 0] = 0
        else:
            weights = np.abs(posterior_reweighted['A_alpha'])**0.5
    return weights


def reweighted_pdf_1D(posterior, prior_initial, prior_final='uniform',
                      parameter='A_alpha', n_bins=10000, cosmology=None):
    """
    Reweight 1D LIV posterior between prior on effective LIV parameter and prior on the physical
    LIV parameter, and smooth it out with kde.
    Use only with parameters that do NOT have sharp posterior cutoff / periodic bounds

    Parameters
    ==========
    posterior: pandas.DataFrame
        Bilby LIV posterior to transform
    prior_initial: string, one of A_eff, log_A_eff, log_lambda_eff
        String specifying prior in which LIV parameter was sampled
    prior_final: string, one of uniform, mass_uniform,
        String specifying desired prior on A_alpha for reweighting
    parameter: string,
        String specifying the name of the parameter to transform
    n_bins: int,
        Number of interpolation points

    cosmology: astropy.cosmology.FLRW, str, dict
        Description of cosmology, one of:
            None - Use default bilby cosmology
            Instance of astropy.cosmology.FLRW subclass
            String with name of known Astropy cosmology, e.g., Planck13
            Dictionary with arguments required to instantiate the cosmology
            class.


    Output

    interpolator: interpolation function for computing pdf of the given parameter
    p_range: tuple giving the minimum and maximum value of the parameter among the posterior samples

    """
    from scipy.stats import gaussian_kde
    from scipy.interpolate import interp1d

    # Check if the posterior is compatible with LIV analysis
    assert ((prior_initial in ['A_eff', 'log_A_eff', 'log_lambda_eff', 'mass_graviton_eff']) and
            (prior_final in ['uniform', 'mass_uniform', 'default'])), "Unknown prior"
    weights = get_reweighting_weights(posterior, prior_initial, prior_final)

    # make pdf using kde
    samples = posterior[parameter].values.copy()
    pdf = gaussian_kde(samples, weights=weights)
    n_bins = 10000
    p_min = min(samples)
    p_max = max(samples)
    x = np.linspace(2 * min(p_min, 0), 2 * p_max, n_bins)

    return interp1d(x, pdf(x), fill_value=0., bounds_error=False), (p_min, p_max)
