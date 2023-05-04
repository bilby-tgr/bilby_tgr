import lalsimulation
from bilby.gw.source import _base_lal_cbc_fd_waveform
from lal import CreateDict


def lal_binary_black_hole(
    frequency_array, mass_1, mass_2, luminosity_distance,
    a_1, tilt_1, phi_12, a_2, tilt_2, phi_jl,
    theta_jn, phase,
    domega220, dtau220, domega330, dtau330, domega210, dtau210, domega440, dtau440, domega550, dtau550,
    **kwargs
):
    """
    A binary black hole waveform using lalsimulation with
    parameterized corrections in the TIGER framework.

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
    domega220: float
        The non-GR fractional deviation to the 220 frequency
    dtau220: float
        The non-GR fractional deviation to the 220 damping time
    domega330: float
        The non-GR fractional deviation to the 330 frequency
    dtau330: float
        The non-GR fractional deviation to the 330 damping time
    domega210: float
        The non-GR fractional deviation to the 210 frequency
    dtau210: float
        The non-GR fractional deviation to the 210 damping time
    domega440: float
        The non-GR fractional deviation to the 440 frequency
    dtau440: float
        The non-GR fractional deviation to the 440 damping time
    domega550: float
        The non-GR fractional deviation to the 550 frequency
    dtau550: float
        The non-GR fractional deviation to the 550 damping time
    kwargs: dict
        Optional keyword arguments

    Returns
    -------
    dict: A dictionary with the plus and cross polarisation strain modes
    """
    waveform_kwargs = dict(
        waveform_approximant='pSEOBNRv4HM_PA', reference_frequency=50.0,
        minimum_frequency=20.0, maximum_frequency=frequency_array[-1],
        catch_waveform_errors=False, pn_spin_order=-1, pn_tidal_order=-1,
        pn_phase_order=-1, pn_amplitude_order=0)
    waveform_kwargs.update(kwargs)

    # If no waveform dictionary is passed, create one
    wf_dict = waveform_kwargs.get("lal_waveform_dictionary", CreateDict())

    # Fractional deviation in the 220 frequency and damping time
    if domega220 != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertDOmega220(wf_dict, float(domega220))
    if dtau220 != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertDTau220(wf_dict, float(dtau220))

    # Fractional deviation in the 210 frequency and damping time
    if domega210 != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertDOmega210(wf_dict, float(domega210))
    if dtau210 != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertDTau210(wf_dict, float(dtau210))

    # Fractional deviation in the 330 frequency and damping time
    if domega330 != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertDOmega330(wf_dict, float(domega330))
    if dtau330 != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertDTau330(wf_dict, float(dtau330))

    # Fractional deviation in the 440 frequency and damping time
    if domega440 != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertDOmega440(wf_dict, float(domega440))
    if dtau440 != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertDTau440(wf_dict, float(dtau440))

    # Fractional deviation in the 550 frequency and damping time
    if domega550 != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertDOmega550(wf_dict, float(domega550))
    if dtau550 != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertDTau550(wf_dict, float(dtau550))

    # Update waveform dictionary with TIGER coefficients
    waveform_kwargs["lal_waveform_dictionary"] = wf_dict

    return _base_lal_cbc_fd_waveform(
        frequency_array=frequency_array, mass_1=mass_1, mass_2=mass_2,
        luminosity_distance=luminosity_distance, theta_jn=theta_jn, phase=phase,
        a_1=a_1, a_2=a_2, tilt_1=tilt_1, tilt_2=tilt_2, phi_12=phi_12,
        phi_jl=phi_jl, **waveform_kwargs)
