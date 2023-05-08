from lal import CreateDict
import lalsimulation as ls
import lal
import numpy as np
from lal import MSUN_SI

# This gives polarizations from a given hlm mode by adding both the positive and negative m.
# see eqs. (5a) and (5b) of https://arxiv.org/pdf/2203.13937.pdf
# Directly coded from LALSimulation. see here
# https://git.ligo.org/lscsoft/lalsuite/-/blob/master/lalsimulation/lib/LALSimSphHarmMode.c#L175


def AddmodesFD(hlmf, L, M, iota, phiRef):
    hlmf = hlmf.data.data
    Y = lal.SpinWeightedSphericalHarmonic(iota, np.pi / 2 - phiRef, -2, L, M)
    Ymstar = np.conj(lal.SpinWeightedSphericalHarmonic(iota, np.pi / 2 - phiRef, -2, L, -M))
    factorp = 0.5 * (Y + (-1)**(L) * Ymstar)
    factorc = 1j * 0.5 * (Y - (-1)**(L) * Ymstar)
    hpf = factorp * hlmf
    hcf = factorc * hlmf
    return hpf, hcf


def lal_binary_black_hole(
        frequency_array, mass_1, mass_2, luminosity_distance,
        a_1, tilt_1, phi_12, a_2, tilt_2, phi_jl, theta_jn, phase,
        dchi_Minus2, dchi_0, dchi_1, dchi_2, dchi_3, dchi_3S, dchi_3NS, dchi_4, dchi_4S, dchi_4NS,
        dchi_5l, dchi_5lS, dchi_5lNS, dchi_6, dchi_6S, dchi_6NS, dchi_6l, dchi_7, dchi_7S, dchi_7NS,
        dchi_kappaS, dchi_kappaA, f_window_div_f_peak, NumCycles,
        **kwargs):

    """
    A binary black hole waveform using lalsimulation with
    parameterized corrections in the FTI framework.

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
    dchi_Minus2: float
        The non-GR inspiral parameter at -1PN
    dchi_0: float
        The non-GR inspiral parameter at 0PN
    dchi_1: float
        The non-GR inspiral parameter at 0.5PN
    dchi_2: float
        The non-GR inspiral parameter at 1PN
    dchi_3: float
        The non-GR inspiral parameter at 1.5PN
    dchi_3S: float
        The non-GR inspiral parameter at 1.5PN, spin part
    dchi_3NS: float
        The non-GR inspiral parameter at 1.5PN, non-spinning part
    dchi_4: float
        The non-GR inspiral parameter at 2PN
    dchi_4S: float
        The non-GR inspiral parameter at 2PN, spin part
    dchi_4NS: float
        The non-GR inspiral parameter at 2PN, non-spinning part
    dchi_5l: float
        The non-GR inspiral parameter at 2.5PN, logarithmic correction
    dchi_5lS: float
        The non-GR inspiral parameter at 2.5PN, logarithmic correction, spin part
    dchi_5lNS: float
        The non-GR inspiral parameter at 2.5PN, logarithmic correction, non-spinning part
    dchi_6: float
        The non-GR inspiral parameter at 3PN
    dchi_6S: float
        The non-GR inspiral parameter at 3PN, spin part
    dchi_6NS: float
        The non-GR inspiral parameter at 3PN, non-spinning part
    dchi_6l: float
        The non-GR inspiral parameter at 3PN, logarithmic correction
    dchi_7: float
        The non-GR inspiral parameter at 3.5PN
    dchi_7S: float
        The non-GR inspiral parameter at 3.5PN, spin part
    dchi_7NS: float
        The non-GR inspiral parameter at 3.5PN, non-spinning part
    dchi_kappaS: float
        The non-GR inspiral parameter for the symmetric
        combination of the spin-induced quadrupole moments
    dchi_kappaA: float
        The non-GR inspiral parameter for the anti-symmetric
        combination of the spin-induced quadrupole moments
    f_window_div_f_peak: float
        The tapering frequency as fraction of the peak frequency of the 22-mode
    NumCycles: float
        The number of gravitational wave cycles over which the tapering happens
    kwargs: dict
        Optional keyword arguments

    Returns
    -------
    dict: A dictionary with the plus and cross polarisation strain modes
    """

    waveform_kwargs = dict(
        waveform_approximant='SEOBNRv4HM_ROM', reference_frequency=50.0,
        minimum_frequency=20.0, maximum_frequency=frequency_array[-1],
        catch_waveform_errors=False, pn_spin_order=-1, pn_tidal_order=-1,
        pn_phase_order=-1, pn_amplitude_order=0, n_modes=4)
    waveform_kwargs.update(kwargs)

    # If no waveform dictionary is passed, create one
    wf_dict = waveform_kwargs.get("lal_waveform_dictionary", CreateDict())
    nongrParams = lal.CreateDict()

    # Relative shifts for inspiral phase PN coefficients
    if dchi_Minus2 != 0.:
        ls.SimInspiralWaveformParamsInsertNonGRDChiMinus2(wf_dict, float(dchi_Minus2))
        lal.DictInsertREAL8Value(nongrParams, "dchiMinus2", dchi_Minus2)

    if dchi_0 != 0.:
        ls.SimInspiralWaveformParamsInsertNonGRDChi0(wf_dict, float(dchi_0))
        lal.DictInsertREAL8Value(nongrParams, "dchi0", dchi_0)

    if dchi_1 != 0.:
        ls.SimInspiralWaveformParamsInsertNonGRDChi1(wf_dict, float(dchi_1))
        lal.DictInsertREAL8Value(nongrParams, "dchi1", dchi_1)

    if dchi_2 != 0.:
        ls.SimInspiralWaveformParamsInsertNonGRDChi2(wf_dict, float(dchi_2))
        lal.DictInsertREAL8Value(nongrParams, "dchi2", dchi_2)

    if dchi_3 != 0.:
        ls.SimInspiralWaveformParamsInsertNonGRDChi3(wf_dict, float(dchi_3))
        lal.DictInsertREAL8Value(nongrParams, "dchi3", dchi_3)

    if dchi_3S != 0.:
        ls.SimInspiralWaveformParamsInsertNonGRDChi3S(wf_dict, float(dchi_3S))
        lal.DictInsertREAL8Value(nongrParams, "dchi3S", dchi_3S)

    if dchi_3NS != 0.:
        ls.SimInspiralWaveformParamsInsertNonGRDChi3NS(wf_dict, float(dchi_3NS))
        lal.DictInsertREAL8Value(nongrParams, "dchi3NS", dchi_3NS)

    if dchi_4 != 0.:
        ls.SimInspiralWaveformParamsInsertNonGRDChi4(wf_dict, float(dchi_4))
        lal.DictInsertREAL8Value(nongrParams, "dchi4", dchi_4)

    if dchi_4S != 0.:
        ls.SimInspiralWaveformParamsInsertNonGRDChi4S(wf_dict, float(dchi_4S))
        lal.DictInsertREAL8Value(nongrParams, "dchi4S", dchi_4S)

    if dchi_4NS != 0.:
        ls.SimInspiralWaveformParamsInsertNonGRDChi4NS(wf_dict, float(dchi_4NS))
        lal.DictInsertREAL8Value(nongrParams, "dchi4NS", dchi_4NS)

    if dchi_5l != 0.:
        ls.SimInspiralWaveformParamsInsertNonGRDChi5L(wf_dict, float(dchi_5l))
        lal.DictInsertREAL8Value(nongrParams, "dchi5l", dchi_5l)

    if dchi_5lS != 0.:
        ls.SimInspiralWaveformParamsInsertNonGRDChi5LS(wf_dict, float(dchi_5lS))
        lal.DictInsertREAL8Value(nongrParams, "dchi5lS", dchi_5lS)

    if dchi_5lNS != 0.:
        ls.SimInspiralWaveformParamsInsertNonGRDChi5LNS(wf_dict, float(dchi_5lNS))
        lal.DictInsertREAL8Value(nongrParams, "dchi5lNS", dchi_5lNS)

    if dchi_6 != 0.:
        ls.SimInspiralWaveformParamsInsertNonGRDChi6(wf_dict, float(dchi_6))
        lal.DictInsertREAL8Value(nongrParams, "dchi6", dchi_6)

    if dchi_6S != 0.:
        ls.SimInspiralWaveformParamsInsertNonGRDChi6S(wf_dict, float(dchi_6S))
        lal.DictInsertREAL8Value(nongrParams, "dchi6S", dchi_6S)

    if dchi_6NS != 0.:
        ls.SimInspiralWaveformParamsInsertNonGRDChi6NS(wf_dict, float(dchi_6NS))
        lal.DictInsertREAL8Value(nongrParams, "dchi6NS", dchi_6NS)

    if dchi_6l != 0.:
        ls.SimInspiralWaveformParamsInsertNonGRDChi6L(wf_dict, float(dchi_6l))
        lal.DictInsertREAL8Value(nongrParams, "dchi6l", dchi_6l)

    if dchi_7 != 0.:
        ls.SimInspiralWaveformParamsInsertNonGRDChi7(wf_dict, float(dchi_7))
        lal.DictInsertREAL8Value(nongrParams, "dchi7", dchi_7)

    if dchi_7S != 0.:
        ls.SimInspiralWaveformParamsInsertNonGRDChi7S(wf_dict, float(dchi_7S))
        lal.DictInsertREAL8Value(nongrParams, "dchi7S", dchi_7S)

    if dchi_7NS != 0.:
        ls.SimInspiralWaveformParamsInsertNonGRDChi7NS(wf_dict, float(dchi_7NS))
        lal.DictInsertREAL8Value(nongrParams, "dchi7NS", dchi_7NS)

    if dchi_kappaS != 0.:
        ls.SimInspiralWaveformParamsInsertNonGRDChikappaS(wf_dict, float(dchi_kappaS))
        lal.DictInsertREAL8Value(nongrParams, "dchikappaS", dchi_kappaS)

    if dchi_kappaA != 0.:
        ls.SimInspiralWaveformParamsInsertNonGRDChikappaA(wf_dict, float(dchi_kappaA))
        lal.DictInsertREAL8Value(nongrParams, "dchikappaA", dchi_kappaA)

    # Update waveform dictionary with FTI coefficients
    waveform_kwargs["lal_waveform_dictionary"] = wf_dict

    fRef = waveform_kwargs['reference_frequency']
    fMin = waveform_kwargs['minimum_frequency']
    fMax = frequency_array[-1]
    deltaF = frequency_array[1] - frequency_array[0]
    approximant = waveform_kwargs['waveform_approximant']
    nmodes = waveform_kwargs['n_modes']

    m1_SI = mass_1 * MSUN_SI
    m2_SI = mass_2 * MSUN_SI

    iota, S1x, S1y, S1z, S2x, S2y, S2z = ls.SimInspiralTransformPrecessingNewInitialConditions(
        theta_jn, phi_jl, tilt_1, tilt_2, phi_12, a_1, a_2, m1_SI, m2_SI, fRef, phase)

    dL = luminosity_distance * 1.0e6 * ls.lal.PC_SI

    # Treat waveform models without higher modes separately
    if approximant in ['SEOBNRv4_ROM', 'SEOBNRv5_ROM', 'IMRPhenomXAS']:

        approx = ls.GetApproximantFromString(approximant)
        hpf, hcf = ls.SimInspiralChooseFDWaveform(
            m1_SI, m2_SI, 0, 0, S1z, 0, 0, S2z, dL, iota, phase, 0, 0, 0, deltaF, fMin, fMax, fRef, None, approx)

        ls.SimInspiralTestingGRCorrections(
            hpf, 2, 2, m1_SI, m2_SI, S1z, S2z, fMin, fRef, f_window_div_f_peak, NumCycles, nongrParams)
        ls.SimInspiralTestingGRCorrections(
            hcf, 2, 2, m1_SI, m2_SI, S1z, S2z, fMin, fRef, f_window_div_f_peak, NumCycles, nongrParams)

        h_plus = hpf.data.data
        h_cross = hcf.data.data

    elif approximant == 'SEOBNRv4HM_ROM':

        if nmodes == 1:
            approx = ls.GetApproximantFromString('SEOBNRv4_ROM')
            hpf, hcf = ls.SimInspiralChooseFDWaveform(
                m1_SI, m2_SI, 0, 0, S1z, 0, 0, S2z, dL, iota, phase, 0, 0, 0, deltaF, fMin, fMax, fRef, None, approx)

            ls.SimInspiralTestingGRCorrections(
                hpf, 2, 2, m1_SI, m2_SI, S1z, S2z, fMin, fRef, f_window_div_f_peak, NumCycles, nongrParams)
            ls.SimInspiralTestingGRCorrections(
                hcf, 2, 2, m1_SI, m2_SI, S1z, S2z, fMin, fRef, f_window_div_f_peak, NumCycles, nongrParams)

            h_plus = hpf.data.data
            h_cross = hcf.data.data

        elif nmodes in [2, 3, 4, 5]:
            hlm_LAL = ls.SimIMRSEOBNRv4HMROM_Modes(
                phase, deltaF, fMin, fMax, fRef, dL, m1_SI, m2_SI, S1z, S2z, -1, nmodes, use_hybridization=True)

            h_plus = None
            h_cross = None

            modes = [[2, 2], [3, 3], [2, 1], [4, 4], [5, 5]][:nmodes]
            for L, M in modes:
                h = ls.SphHarmFrequencySeriesGetMode(hlm_LAL, L, -M)

                ls.SimInspiralTestingGRCorrections(
                    h, L, M, m1_SI, m2_SI, S1z, S2z, fMin, fRef, f_window_div_f_peak, NumCycles, nongrParams)
                hp_corr_data, hc_corr_data = AddmodesFD(h, L, -M, iota, phase)

                if h_plus is None:
                    h_plus = hp_corr_data
                    h_cross = hc_corr_data
                else:
                    h_plus += hp_corr_data
                    h_cross += hc_corr_data

        else:
            print("ERROR: The number of modes used with SEOBNRv4HM_ROM cannot be n_modes = {}.".format(nmodes))

    else:
        print("ERROR: Waveform approximant {} is not supported. \
               Use SEOBNRv4_ROM, SEOBNRv5_ROM, IMRPhenomXAS, or SEOBNRv4HM_ROM instead.".format(approximant))

    return dict(plus=h_plus, cross=h_cross)
