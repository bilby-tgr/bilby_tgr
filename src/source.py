import lalsimulation
from bilby.gw.source import lal_binary_black_hole
from lal import CreateDict


def non_gr_d_alpha_1(
        frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
        phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, d_alpha_1, **kwargs):
    """ Generate a cbc waveform model with delta alpha 1 using lalsimulation

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
    d_alpha_1: float
        The non-GR parameter
    kwargs: dict
        Optional keyword arguments

    Returns
    -------
    dict: A dictionary with the plus and cross polarisation strain modes
    """
    wf_dict = kwargs.get("lal_waveform_dictionary", CreateDict())
    lalsimulation.SimInspiralWaveformParamsInsertNonGRDAlpha1(wf_dict, float(d_alpha_1))
    kwargs["lal_waveform_dictionary"] = wf_dict
    return lal_binary_black_hole(
        frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
        phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, **kwargs)
