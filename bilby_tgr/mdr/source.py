import bilby
from scipy.constants import physical_constants
import numpy as np
from bilby.core.utils.constants import speed_of_light as c, parsec as pc, solar_mass, gravitational_constant as G


h = physical_constants["Planck constant in eV/Hz"][0] * 1e12  # h in peV*s


def liv_correction_group_velocity(
        frequency_array, luminosity_distance, alpha, A_eff, **kwargs):
    """
    calculates LIV correction using group velocity following https://arxiv.org/pdf/2203.13252.pdf

    Parameters:


    frequency_array: frequencies on which to evaluate the correction, in Hz
    luminosity_distance: float
        luminosity distance in Mpc
    alpha: float
        exponent of lorentz violating momentum term
    A_eff: float
        A_eff in units of peV^(2-alpha)

    Output

    correction: np.array
        LIV phase correction factor, so that the overall modification is exp(1j*correction)"""

    distance = luminosity_distance * pc * 1e6
    zeta = -(np.pi * distance * A_eff * h**(alpha - 2.) / c)
    return zeta * frequency_array**(alpha - 1.)


def liv_correction_phase_velocity(
        frequency_array, mass_1, mass_2, luminosity_distance, alpha, A_eff, **kwargs):
    """
    calculates LIV correction using phase velocity following https://arxiv.org/pdf/1110.2720.pdf,
    but parametrized in A_eff

    Parameters:


    frequency_array: frequencies on which to evaluate the correction, in Hz
    mass_1, mass_2: masses of the binary, in solar mass
    luminosity_distance: float
        luminosity distance in Mpc
    alpha: float
        exponent of lorentz violating momentum term
    A_eff: float
        A_eff in units of peV^(2-alpha)

    Output

    correction: np.array
        LIV phase correction factor, so that the overall modification is exp(1j*correction)"""

    if (alpha == 1.):
        distance = luminosity_distance * pc * 1e6
        zeta = np.pi * distance * A_eff / h / c
        mchirp = bilby.gw.conversion.component_masses_to_chirp_mass(mass_1, mass_2) * solar_mass
        freq_factor = np.pi * mchirp * G / c**3
        return zeta * np.log(freq_factor * frequency_array)
    else:
        distance = luminosity_distance * pc * 1e6
        zeta = (np.pi * distance * A_eff / (alpha - 1.) * h**(alpha - 2.) / c)
        return zeta * frequency_array**(alpha - 1.)


def lal_binary_black_hole(
        frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
        phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, alpha, A_eff, **kwargs):
    """
    Binary black hole model with group velocity LIV correction

    alpha: float
        exponent of lorentz violating momentum term
    A_eff: float
        A_eff in units of peV^(2-alpha)
    """

    gr_waveform = bilby.gw.source.lal_binary_black_hole(
        frequency_array=frequency_array, mass_1=mass_1, mass_2=mass_2, luminosity_distance=luminosity_distance,
        a_1=a_1, tilt_1=tilt_1, phi_12=phi_12, a_2=a_2, tilt_2=tilt_2, phi_jl=phi_jl,
        theta_jn=theta_jn, phase=phase, **kwargs)

    correction = np.zeros_like(frequency_array, dtype=complex)
    if correction[0] == 0:  # avoid raising 0 to negative power
        correction[1:] = np.exp(liv_correction_group_velocity(
            frequency_array[1:], luminosity_distance, alpha, A_eff, **kwargs) * 1j)
    else:
        correction = np.exp(liv_correction_group_velocity(
            frequency_array, luminosity_distance, alpha, A_eff, **kwargs) * 1j)

    if gr_waveform is None:
        return None
    else:
        return dict(plus=gr_waveform['plus'] * correction, cross=gr_waveform['cross'] * correction)


def lal_binary_neutron_star(
        frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
        phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, lambda_1, lambda_2, alpha, A_eff, **kwargs):
    """
    Binary neutron star model with group velocity LIV correction

    alpha: float
        exponent of lorentz violating momentum term
    A_eff: float
        A_eff in units of peV^(2-alpha)
    """

    gr_waveform = bilby.gw.source.lal_binary_neutron_star(
        frequency_array=frequency_array, mass_1=mass_1, mass_2=mass_2, luminosity_distance=luminosity_distance,
        a_1=a_1, tilt_1=tilt_1, phi_12=phi_12, a_2=a_2, tilt_2=tilt_2, phi_jl=phi_jl,
        theta_jn=theta_jn, phase=phase, lambda_1=lambda_1, lambda_2=lambda_2, **kwargs)

    correction = np.zeros_like(frequency_array, dtype=complex)
    if correction[0] == 0:  # avoid raising 0 to negative power
        correction[1:] = np.exp(liv_correction_group_velocity(
            frequency_array[1:], luminosity_distance, alpha, A_eff, **kwargs) * 1j)
    else:
        correction = np.exp(liv_correction_group_velocity(
            frequency_array, luminosity_distance, alpha, A_eff, **kwargs) * 1j)

    if gr_waveform is None:
        return None
    else:
        return dict(plus=gr_waveform['plus'] * correction, cross=gr_waveform['cross'] * correction)


def lal_binary_black_hole_phase_velocity(
        frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
        phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, alpha, A_eff, **kwargs):
    """
    Binary black hole model with phase velocity LIV correction

    alpha: float
        exponent of lorentz violating momentum term
    A_eff: float
        A_eff in units of peV^(2-alpha)
    """

    gr_waveform = bilby.gw.source.lal_binary_black_hole(
        frequency_array=frequency_array, mass_1=mass_1, mass_2=mass_2, luminosity_distance=luminosity_distance,
        a_1=a_1, tilt_1=tilt_1, phi_12=phi_12, a_2=a_2, tilt_2=tilt_2, phi_jl=phi_jl,
        theta_jn=theta_jn, phase=phase, **kwargs)

    correction = np.zeros_like(frequency_array, dtype=complex)
    if correction[0] == 0:  # avoid raising 0 to negative power
        correction[1:] = np.exp(liv_correction_phase_velocity(
            frequency_array[1:], mass_1, mass_2, luminosity_distance, alpha, A_eff, **kwargs) * 1j)
    else:
        correction = np.exp(liv_correction_phase_velocity(
            frequency_array, mass_1, mass_2, luminosity_distance, alpha, A_eff, **kwargs) * 1j)

    if isinstance(gr_waveform, type(None)):
        return None
    else:
        return dict(plus=gr_waveform['plus'] * correction, cross=gr_waveform['cross'] * correction)


def lal_binary_neutron_star_phase_velocity(
        frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
        phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, lambda_1, lambda_2, alpha, A_eff, **kwargs):
    """
    Binary neutron star model with phase velocity LIV correction

    alpha: float
        exponent of lorentz violating momentum term
    A_eff: float
        A_eff in units of peV^(2-alpha)
    """

    gr_waveform = bilby.gw.source.lal_binary_neutron_star(
        frequency_array=frequency_array, mass_1=mass_1, mass_2=mass_2, luminosity_distance=luminosity_distance,
        a_1=a_1, tilt_1=tilt_1, phi_12=phi_12, a_2=a_2, tilt_2=tilt_2, phi_jl=phi_jl,
        theta_jn=theta_jn, phase=phase, lambda_1=lambda_1, lambda_2=lambda_2, **kwargs)

    correction = np.zeros_like(frequency_array, dtype=complex)
    if correction[0] == 0:  # avoid raising 0 to negative power
        correction[1:] = np.exp(liv_correction_phase_velocity(
            frequency_array[1:], mass_1, mass_2, luminosity_distance, alpha, A_eff, **kwargs) * 1j)
    else:
        correction = np.exp(liv_correction_phase_velocity(
            frequency_array, mass_1, mass_2, luminosity_distance, alpha, A_eff, **kwargs) * 1j)

    if isinstance(gr_waveform, type(None)):
        return None
    else:
        return dict(plus=gr_waveform['plus'] * correction, cross=gr_waveform['cross'] * correction)
