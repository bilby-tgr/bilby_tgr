import bilby
import numpy as np
from bilby.core.utils.constants import speed_of_light as c


def ssb_factor(
        frequency_array, ssb_coeff, **kwargs):
    """
    Compute the modified dispersion of the waveform due to spacetime symmetry breaking (SSB)
    as described in Eq.(12) of arxiv:2210.04481 (only mass-dimension 5 terms in the Lagrangian density are considered)

    Parameters:

    frequency_array: array of floats
        frequencies on which to evaluate the correction [unit: s^{-1}]
    ssb_coeff: float
        modified dispersion coefficient [unit: m*s]

    Output

    correction: np.array
        SSB correction factor beta to apply to cross and plus polarizations
    """

    beta = 2 * np.pi**(1.5) * ssb_coeff / c
    return beta * frequency_array**2


def lal_binary_black_hole(
        frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
        phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, ssb_coeff, **kwargs):
    """
    Binary black hole waveform with modified dispersion due to spacetime symmetry breaking (SSB)
    as described in arxiv:2210.04481. Modified dispersion corresponding to only mass-dimension 5
    terms in the Lagrangian density are considered.

    SSB extra parameters:

    ssb_coeff: float
        modified dispersion coefficient [unit: m*s]
    """

    gr_waveform = bilby.gw.source.lal_binary_black_hole(
        frequency_array=frequency_array, mass_1=mass_1, mass_2=mass_2, luminosity_distance=luminosity_distance,
        a_1=a_1, tilt_1=tilt_1, phi_12=phi_12, a_2=a_2, tilt_2=tilt_2, phi_jl=phi_jl,
        theta_jn=theta_jn, phase=phase, **kwargs)

    ssb_correction = ssb_factor(frequency_array, ssb_coeff, **kwargs)

    if gr_waveform is None:
        return None
    else:
        return dict(plus=gr_waveform['plus'] * np.cos(ssb_correction) - gr_waveform['cross'] * np.sin(ssb_correction),
                    cross=gr_waveform['cross'] * np.cos(ssb_correction) + gr_waveform['plus'] * np.sin(ssb_correction))
