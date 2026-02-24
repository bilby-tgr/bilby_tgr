from bilby.gw.source import (lal_binary_black_hole, lal_binary_neutron_star,
                             lal_binary_black_hole_relative_binning, lal_binary_neutron_star_relative_binning)
from .Phase_Amplitude_correction_losa import phase_correction_losa, amplitude_correction_losa
import numpy as np


def lal_binary_black_hole_losa(frequency_array,
                               mass_1,
                               mass_2,
                               luminosity_distance,
                               a_1,
                               tilt_1,
                               phi_12,
                               a_2,
                               tilt_2,
                               phi_jl,
                               theta_jn,
                               phase,
                               a_by_c,
                               **kwargs
                               ):
    '''
    A binary black hole waveform using lal_binary_black_hole with additional -4 PN phase correction due to losa as
    given in arXiv:2302.09651 and 2506.22272.

    Parameters
    ----------
    frequency_array: array_like
        The frequencies at which we want to calculate the strain
    mass_1: float
        The mass of the heavier object in solar masses
    mass_2: float
        The mass of the lighter object in solar masses
    luminosity_distance: float
        The luminosity distance in megaparsec
    a_1: float
        Dimensionless primary spin magnitude
    tilt_1: float
        Primary tilt angle
    phi_12: float
        Azimuthal angle between the component spins
    a_2: float
        Dimensionless secondary spin magnitude
    tilt_2: float
        Secondary tilt angle
    phi_jl: float
        Azimuthal angle between the total and orbital angular momenta
    theta_jn: float
        Orbital inclination
    phase: float
        The phase at coalescence
    a_by_c: float
        The line-of-sight acceleration of the center of the binary in
        s^{-1} (divided by the speed of light)
    kwargs: dict
        Optional keyword arguments.
        e.g.:
        pn_losa_order: int, optional
        The Next to leading order post-newtonian correction to phase due to LOSA.
        Default is -1 (all terms up to 3.5PN order in point-particle, aligned-spin,
                       and 7.5PN in tidal order are considered).

        Note: pn_losa_order should be an integer. Set it to a value 2 times the required PN order.
              For example: pn_losa_order = 3 will represent 1.5PN order.

    Returns
    -------
    dict: A dictionary with the plus and cross polarisation strain modes

    '''

    waveform_kwargs = dict(
        waveform_approximant="IMRPhenomXP",
        reference_frequency=50.0,
        minimum_frequency=20.0,
        maximum_frequency=frequency_array[-1],
        catch_waveform_errors=False,
        pn_spin_order=-1,
        pn_tidal_order=-1,
        pn_phase_order=-1,
        pn_amplitude_order=0,
        pn_losa_order=-1
    )

    waveform_kwargs.update(kwargs)

    waveform_to_return = lal_binary_black_hole(
        frequency_array=frequency_array,
        mass_1=mass_1,
        mass_2=mass_2,
        luminosity_distance=luminosity_distance,
        theta_jn=theta_jn,
        phase=phase,
        a_1=a_1,
        a_2=a_2,
        tilt_1=tilt_1,
        tilt_2=tilt_2,
        phi_12=phi_12,
        phi_jl=phi_jl,
        **waveform_kwargs
    )

    phase_correction = phase_correction_losa(frequency_array[1:], mass_1, mass_2, a_by_c,
                                             a_1=a_1, a_2=a_2, tilt_1=tilt_1, tilt_2=tilt_2,
                                             pn_losa_order=waveform_kwargs["pn_losa_order"])

    amplitude_correction = amplitude_correction_losa(
        frequency_array[1:], mass_1, mass_2, a_by_c)

    for key in ["plus", "cross"]:
        waveform_to_return[key][1:] *= np.exp(-1j *
                                              phase_correction) * (1 + amplitude_correction)

    return waveform_to_return


def lal_binary_neutron_star_losa(frequency_array,
                                 mass_1,
                                 mass_2,
                                 luminosity_distance,
                                 a_1,
                                 tilt_1,
                                 phi_12,
                                 a_2,
                                 tilt_2,
                                 phi_jl,
                                 theta_jn,
                                 phase,
                                 lambda_1,
                                 lambda_2,
                                 a_by_c,
                                 **kwargs
                                 ):
    '''
    Parameters
    ----------
    Parameters of lal_binary_black_hole_losa together with the following:

    lambda_1: float
        Dimensionless tidal deformability of mass_1
    lambda_2: float
        Dimensionless tidal deformability of mass_2
    kwargs: dict
        Optional keyword arguments.
        e.g.:
        pn_losa_order: int, optional
        The Next to leading order post-newtonian correction to phase due to LOSA.
        Default is -1 (all terms up to 3.5PN order in point-particle, aligned-spin,
                       and 7.5PN in tidal order are considered).

        Note: pn_losa_order should be an integer. Set it to a value 2 times the required PN order.
              For example: pn_losa_order = 3 will represent 1.5PN order.

    Returns
    -------
    dict: A dictionary with the plus and cross polarisation strain modes

    '''

    waveform_kwargs = dict(
        waveform_approximant="IMRPhenomXP_NRTidalv2",
        reference_frequency=50.0,
        minimum_frequency=20.0,
        maximum_frequency=frequency_array[-1],
        catch_waveform_errors=False,
        pn_spin_order=-1,
        pn_tidal_order=-1,
        pn_phase_order=-1,
        pn_amplitude_order=0,
        pn_losa_order=-1
    )

    waveform_kwargs.update(kwargs)

    waveform_to_return = lal_binary_neutron_star(
        frequency_array=frequency_array,
        mass_1=mass_1,
        mass_2=mass_2,
        luminosity_distance=luminosity_distance,
        theta_jn=theta_jn,
        phase=phase,
        a_1=a_1,
        a_2=a_2,
        tilt_1=tilt_1,
        tilt_2=tilt_2,
        phi_12=phi_12,
        phi_jl=phi_jl,
        lambda_1=lambda_1,
        lambda_2=lambda_2,
        **waveform_kwargs
    )

    phase_correction = phase_correction_losa(frequency_array[1:], mass_1, mass_2, a_by_c,
                                             a_1=a_1, a_2=a_2, tilt_1=tilt_1, tilt_2=tilt_2,
                                             lambda_1=lambda_1, lambda_2=lambda_2,
                                             pn_losa_order=waveform_kwargs["pn_losa_order"])

    amplitude_correction = amplitude_correction_losa(
        frequency_array[1:], mass_1, mass_2, a_by_c)

    for key in ["plus", "cross"]:
        waveform_to_return[key][1:] *= np.exp(-1j *
                                              phase_correction) * (1 + amplitude_correction)

    return waveform_to_return


def lal_binary_black_hole_relative_binning_losa(
        frequency_array,
        mass_1,
        mass_2,
        luminosity_distance,
        a_1,
        tilt_1,
        phi_12,
        a_2,
        tilt_2,
        phi_jl,
        theta_jn,
        phase,
        a_by_c,
        fiducial,
        **kwargs):
    '''
    Source model to go with RelativeBinningGravitationalWaveTransient likelihood.

    All parameters are the same as in the usual source model: lal_binary_black_hole_losa, except `fiducial`

    fiducial: float
        If fiducial=1, waveform evaluated on the full frequency grid is returned.
        If fiducial=0, waveform evaluated at waveform_kwargs[“frequency_bin_edges”]
        is returned.

    Returns
    -------
    dict: A dictionary with the plus and cross polarisation strain modes

    '''

    waveform_kwargs = dict(
        waveform_approximant='IMRPhenomXP',
        reference_frequency=50.0,
        minimum_frequency=20.0,
        maximum_frequency=frequency_array[-1],
        catch_waveform_errors=False,
        pn_spin_order=-1,
        pn_tidal_order=-1,
        pn_phase_order=-1,
        pn_amplitude_order=0,
        pn_losa_order=-1
    )

    waveform_kwargs.update(kwargs)

    waveform_to_return = lal_binary_black_hole_relative_binning(
        frequency_array=frequency_array, mass_1=mass_1, mass_2=mass_2,
        luminosity_distance=luminosity_distance, theta_jn=theta_jn, phase=phase,
        a_1=a_1, a_2=a_2, tilt_1=tilt_1, tilt_2=tilt_2, phi_jl=phi_jl,
        phi_12=phi_12, fiducial=fiducial, **waveform_kwargs)

    if fiducial == 1:
        phase_correction = phase_correction_losa(frequency_array[1:], mass_1, mass_2, a_by_c,
                                                 a_1=a_1, a_2=a_2, tilt_1=tilt_1, tilt_2=tilt_2,
                                                 pn_losa_order=waveform_kwargs["pn_losa_order"])
        amplitude_correction = amplitude_correction_losa(
            frequency_array[1:], mass_1, mass_2, a_by_c)

        for key in ["plus", "cross"]:
            waveform_to_return[key][1:] *= np.exp(-1j *
                                                  phase_correction) * (1 + amplitude_correction)

    else:
        waveform_kwargs["frequencies"] = waveform_kwargs.pop(
            "frequency_bin_edges")
        phase_correction = phase_correction_losa(waveform_kwargs["frequencies"], mass_1, mass_2, a_by_c,
                                                 a_1=a_1, a_2=a_2, tilt_1=tilt_1, tilt_2=tilt_2,
                                                 pn_losa_order=waveform_kwargs["pn_losa_order"])
        amplitude_correction = amplitude_correction_losa(
            waveform_kwargs["frequencies"], mass_1, mass_2, a_by_c)

        for key in ["plus", "cross"]:
            waveform_to_return[key] *= np.exp(-1j *
                                              phase_correction) * (1 + amplitude_correction)

    return waveform_to_return


def lal_binary_neutron_star_relative_binning_losa(
        frequency_array,
        mass_1,
        mass_2,
        luminosity_distance,
        a_1,
        tilt_1,
        phi_12,
        a_2,
        tilt_2,
        phi_jl,
        lambda_1,
        lambda_2,
        theta_jn,
        phase,
        a_by_c,
        fiducial,
        **kwargs):
    '''
    Source model to go with RelativeBinningGravitationalWaveTransient likelihood.

    All parameters are the same as in the usual source model: lal_binary_neutron_star_losa, except `fiducial`

    fiducial: float
        If fiducial=1, waveform evaluated on the full frequency grid is returned.
        If fiducial=0, waveform evaluated at waveform_kwargs[“frequency_bin_edges”]
        is returned.

    Returns
    -------
    dict: A dictionary with the plus and cross polarisation strain modes

    '''

    waveform_kwargs = dict(
        waveform_approximant='IMRPhenomXP_NRTidalv2',
        reference_frequency=50.0,
        minimum_frequency=20.0,
        maximum_frequency=frequency_array[-1],
        catch_waveform_errors=False,
        pn_spin_order=-1,
        pn_tidal_order=-1,
        pn_phase_order=-1,
        pn_amplitude_order=0,
        pn_losa_order=-1
    )

    waveform_kwargs.update(kwargs)

    waveform_to_return = lal_binary_neutron_star_relative_binning(
        frequency_array=frequency_array, mass_1=mass_1, mass_2=mass_2,
        luminosity_distance=luminosity_distance, theta_jn=theta_jn, phase=phase,
        a_1=a_1, a_2=a_2, tilt_1=tilt_1, tilt_2=tilt_2, phi_12=phi_12,
        phi_jl=phi_jl, lambda_1=lambda_1, lambda_2=lambda_2, fiducial=fiducial, **waveform_kwargs)

    if fiducial == 1:
        phase_correction = phase_correction_losa(frequency_array[1:], mass_1, mass_2, a_by_c,
                                                 a_1=a_1, a_2=a_2, tilt_1=tilt_1, tilt_2=tilt_2,
                                                 lambda_1=lambda_1, lambda_2=lambda_2,
                                                 pn_losa_order=waveform_kwargs["pn_losa_order"])
        amplitude_correction = amplitude_correction_losa(
            frequency_array[1:], mass_1, mass_2, a_by_c)

        for key in ["plus", "cross"]:
            waveform_to_return[key][1:] *= np.exp(-1j *
                                                  phase_correction) * (1 + amplitude_correction)

    else:
        waveform_kwargs["frequencies"] = waveform_kwargs.pop(
            "frequency_bin_edges")
        phase_correction = phase_correction_losa(waveform_kwargs["frequencies"], mass_1, mass_2, a_by_c,
                                                 a_1=a_1, a_2=a_2, tilt_1=tilt_1, tilt_2=tilt_2,
                                                 lambda_1=lambda_1, lambda_2=lambda_2,
                                                 pn_losa_order=waveform_kwargs["pn_losa_order"])
        amplitude_correction = amplitude_correction_losa(
            waveform_kwargs["frequencies"], mass_1, mass_2, a_by_c)

        for key in ["plus", "cross"]:
            waveform_to_return[key] *= np.exp(-1j *
                                              phase_correction) * (1 + amplitude_correction)

    return waveform_to_return
