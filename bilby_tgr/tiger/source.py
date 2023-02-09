import lalsimulation
from bilby.gw.source import _base_lal_cbc_fd_waveform
from lal import CreateDict


def lal_binary_black_hole(
    frequency_array, mass_1, mass_2, luminosity_distance,
    a_1, tilt_1, phi_12, a_2, tilt_2, phi_jl,
    theta_jn, phase,
    dchi_0, dchi_1, dchi_2, dchi_3, dchi_4, dchi_5l, dchi_6, dchi_6l, dchi_7,
    dbeta_2, dbeta_3,
    dalpha_2, dalpha_3, dalpha_4, dalpha_5,
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
    dchi_0: float
        The non-GR inspiral parameter at 0PN
    dchi_1: float
        A non-GR inspiral parameter at 0.5PN
    dchi_2: float
        A non-GR inspiral parameter at 1PN
    dchi_3: float
        A non-GR inspiral parameter at 1.5PN
    dchi_4: float
        A non-GR inspiral parameter at 2PN
    dchi_5l: float
        A non-GR inspiral parameter at 2.5PN, logarithmic correction
    dchi_6: float
        A non-GR inspiral parameter at 3PN
    dchi_6l: float
        A non-GR inspiral parameter at 3PN, logarithmic correction
    dchi_7: float
        A non-GR inspiral parameter at 3.5PN
    dbeta_2: float
        A non-GR intermediate parameter
    dbeta_3: float
        A non-GR intermediate parameter
    dalpha_2: float
        A non-GR merger-ringdown parameter
    dalpha_3: float
        A non-GR merger-ringdown parameter
    dalpha_4: float
        A non-GR merger-ringdown parameter
    dalpha_5: float
        A non-GR merger-ringdown parameter
    kwargs: dict
        Optional keyword arguments

    Returns
    -------
    dict: A dictionary with the plus and cross polarisation strain modes
    """
    waveform_kwargs = dict(
        waveform_approximant='IMRPhenomPv2', reference_frequency=50.0,
        minimum_frequency=20.0, maximum_frequency=frequency_array[-1],
        catch_waveform_errors=False, pn_spin_order=-1, pn_tidal_order=-1,
        pn_phase_order=-1, pn_amplitude_order=0)
    waveform_kwargs.update(kwargs)

    # If no waveform dictionary is passed, create one
    wf_dict = waveform_kwargs.get("lal_waveform_dictionary", CreateDict())

    # Relative shifts for inspiral phase PN coefficients (absolute value for dchi_1)
    if dchi_0 != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertNonGRDChi0(wf_dict, float(dchi_0))
    if dchi_1 != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertNonGRDChi1(wf_dict, float(dchi_1))
    if dchi_2 != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertNonGRDChi2(wf_dict, float(dchi_2))
    if dchi_3 != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertNonGRDChi3(wf_dict, float(dchi_3))
    if dchi_4 != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertNonGRDChi4(wf_dict, float(dchi_4))
    if dchi_5l != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertNonGRDChi5L(wf_dict, float(dchi_5l))
    if dchi_6 != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertNonGRDChi6(wf_dict, float(dchi_6))
    if dchi_6l != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertNonGRDChi6L(wf_dict, float(dchi_6l))
    if dchi_7 != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertNonGRDChi7(wf_dict, float(dchi_7))

    # Relative shifts for intermediate phase coefficients (PhenomD/Pv2)
    if dbeta_2 != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertNonGRDBeta2(wf_dict, float(dbeta_2))
    if dbeta_3 != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertNonGRDBeta3(wf_dict, float(dbeta_3))

    # Relative shifts for merger-ringdown phase coefficients  (PhenomD/Pv2)
    if dalpha_2 != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertNonGRDAlpha2(wf_dict, float(dalpha_2))
    if dalpha_3 != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertNonGRDAlpha3(wf_dict, float(dalpha_3))
    if dalpha_4 != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertNonGRDAlpha4(wf_dict, float(dalpha_4))
    if dalpha_5 != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertNonGRDAlpha5(wf_dict, float(dalpha_5))

    # Update waveform dictionary with TIGER coefficients
    waveform_kwargs["lal_waveform_dictionary"] = wf_dict

    return _base_lal_cbc_fd_waveform(
        frequency_array=frequency_array, mass_1=mass_1, mass_2=mass_2,
        luminosity_distance=luminosity_distance, theta_jn=theta_jn, phase=phase,
        a_1=a_1, a_2=a_2, tilt_1=tilt_1, tilt_2=tilt_2, phi_12=phi_12,
        phi_jl=phi_jl, **waveform_kwargs)


def lal_binary_neutron_star(
        frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
        phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, lambda_1, lambda_2,
        dchi_0, dchi_1, dchi_2, dchi_3, dchi_4, dchi_5l, dchi_6, dchi_6l, dchi_7,
        **kwargs):
    """ A Binary Neutron Star waveform model using lalsimulation with
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
        Azimuthal angle between the two component spins
    a_2: float
        Dimensionless secondary spin magnitude
    tilt_2: float
        Secondary tilt angle
    phi_jl: float
        Azimuthal angle between the total binary angular momentum and the
        orbital angular momentum
    theta_jn: float
        Orbital inclination
    phase: float
        The phase at coalescence
    lambda_1: float
        Dimensionless tidal deformability of mass_1
    lambda_2: float
        Dimensionless tidal deformability of mass_2
    dchi_0: float
        The non-GR inspiral parameter at 0PN
    dchi_1: float
        A non-GR inspiral parameter at 0.5PN
    dchi_2: float
        A non-GR inspiral parameter at 1PN
    dchi_3: float
        A non-GR inspiral parameter at 1.5PN
    dchi_4: float
        A non-GR inspiral parameter at 2PN
    dchi_5l: float
        A non-GR inspiral parameter at 2.5PN, logarithmic correction
    dchi_6: float
        A non-GR inspiral parameter at 3PN
    dchi_6l: float
        A non-GR inspiral parameter at 3PN, logarithmic correction
    dchi_7: float
        A non-GR inspiral parameter at 3.5PN
    kwargs: dict
        Optional keyword arguments
        Supported arguments:
            waveform_approximant
            reference_frequency
            minimum_frequency
            maximum_frequency
            catch_waveform_errors
            pn_spin_order
            pn_tidal_order
            pn_phase_order
            pn_amplitude_order
            mode_array:
                Activate a specific mode array and evaluate the model using those
                modes only.  e.g. waveform_arguments =
                dict(waveform_approximant='IMRPhenomHM', mode_array=[[2,2],[2,-2])
                returns the 22 and 2-2 modes only of IMRPhenomHM.  You can only
                specify modes that are included in that particular model.  e.g.
                waveform_arguments = dict(waveform_approximant='IMRPhenomHM',
                mode_array=[[2,2],[2,-2],[5,5],[5,-5]]) is not allowed because the
                55 modes are not included in this model.  Be aware that some models
                only take positive modes and return the positive and the negative
                mode together, while others need to call both.  e.g.
                waveform_arguments = dict(waveform_approximant='IMRPhenomHM',
                mode_array=[[2,2],[4,-4]]) returns the 22 and 2-2 of IMRPhenomHM.
                However, waveform_arguments =
                dict(waveform_approximant='IMRPhenomXHM', mode_array=[[2,2],[4,-4]])
                returns the 22 and 4-4 of IMRPhenomXHM.

    Returns
    -------
    dict: A dictionary with the plus and cross polarisation strain modes
    """
    waveform_kwargs = dict(
        waveform_approximant='IMRPhenomPv2_NRTidal', reference_frequency=50.0,
        minimum_frequency=20.0, maximum_frequency=frequency_array[-1],
        catch_waveform_errors=False, pn_spin_order=-1, pn_tidal_order=-1,
        pn_phase_order=-1, pn_amplitude_order=0)

    waveform_kwargs.update(kwargs)

    # If no waveform dictionary is passed, create one
    wf_dict = waveform_kwargs.get("lal_waveform_dictionary", CreateDict())

    # Relative shifts for inspiral phase PN coefficients (absolute value for dchi_1)
    if dchi_0 != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertNonGRDChi0(wf_dict, float(dchi_0))
    if dchi_1 != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertNonGRDChi1(wf_dict, float(dchi_1))
    if dchi_2 != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertNonGRDChi2(wf_dict, float(dchi_2))
    if dchi_3 != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertNonGRDChi3(wf_dict, float(dchi_3))
    if dchi_4 != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertNonGRDChi4(wf_dict, float(dchi_4))
    if dchi_5l != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertNonGRDChi5L(wf_dict, float(dchi_5l))
    if dchi_6 != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertNonGRDChi6(wf_dict, float(dchi_6))
    if dchi_6l != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertNonGRDChi6L(wf_dict, float(dchi_6l))
    if dchi_7 != 0.:
        lalsimulation.SimInspiralWaveformParamsInsertNonGRDChi7(wf_dict, float(dchi_7))

    # Update waveform dictionary with TIGER coefficients
    waveform_kwargs["lal_waveform_dictionary"] = wf_dict

    return _base_lal_cbc_fd_waveform(
        frequency_array=frequency_array, mass_1=mass_1, mass_2=mass_2,
        luminosity_distance=luminosity_distance, theta_jn=theta_jn, phase=phase,
        a_1=a_1, a_2=a_2, tilt_1=tilt_1, tilt_2=tilt_2, phi_12=phi_12,
        phi_jl=phi_jl, lambda_1=lambda_1, lambda_2=lambda_2, **waveform_kwargs)
