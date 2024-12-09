import lalsimulation
import numpy as np
from bilby.core import utils
from bilby.core.utils import logger
from bilby.gw.conversion import bilby_to_lalsimulation_spins
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


def gwsignal_binary_black_hole(frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
                               phi_12, a_2, tilt_2, phi_jl, theta_jn, phase,
                               domega220, dtau220, domega330, dtau330, domega210, dtau210, domega440,
                               dtau440, domega550, dtau550, domega320, dtau320, domega430, dtau430,
                               dA22, dw22, dA33, dw33, dA21, dw21, dA44, dw44, dA55, dw55, dA32, dw32, dA43, dw43,
                               dTpeak, da6, ddSO, **kwargs):
    """
    A parametrized SEOB binary black hole waveform model using GWsignal

    Parameters
    ==========
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
        Angle between the total binary angular momentum and the line of sight
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
    domega320: float
        The non-GR fractional deviation to the 320 frequency
    dtau320: float
        The non-GR fractional deviation to the 320 damping time
    domega430: float
        The non-GR fractional deviation to the 430 frequency
    dtau430: float
        The non-GR fractional deviation to the 430 damping time
    dA22: float
        The non-GR fractional deviation to the 22 merger amplitude
    dw22: float
        The non-GR fractional deviation to the 22 merger frequency
    dA33: float
        The non-GR fractional deviation to the 33 merger amplitude
    dw33: float
        The non-GR fractional deviation to the 33 merger frequency
    dA21: float
        The non-GR fractional deviation to the 21 merger amplitude
    dw21: float
        The non-GR fractional deviation to the 21 merger frequency
    dA44: float
        The non-GR fractional deviation to the 44 merger amplitude
    dw44: float
        The non-GR fractional deviation to the 44 merger frequency
    dA55: float
        The non-GR fractional deviation to the 55 merger amplitude
    dw55: float
        The non-GR fractional deviation to the 55 merger frequency
    dA32: float
        The non-GR fractional deviation to the 32 merger amplitude
    dw32: float
        The non-GR fractional deviation to the 32 merger frequency
    dA43: float
        The non-GR fractional deviation to the 43 merger amplitude
    dw43: float
        The non-GR fractional deviation to the 43 merger frequency
    dTpeak: float
        The non-GR additive deviation to the amplitude's peak time
    da6: float
        The non-GR additive deviation to the a6 calibration parameter
    ddSO: float
        The non-GR additive deviation to the dSO calibration parameter
    kwargs: dict
        Optional keyword arguments
        Supported arguments:

        - waveform_approximant
        - reference_frequency
        - minimum_frequency
        - maximum_frequency
        - catch_waveform_errors
        - pn_amplitude_order
        - mode_array:
          Activate a specific mode array and evaluate the model using those
          modes only.  e.g. waveform_arguments =
          dict(waveform_approximant='IMRPhenomHM', mode_array=[[2,2],[2,-2]])
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
    =======
    dict: A dictionary with the plus and cross polarisation strain modes

    Notes
    =====
    This function is a temporary wrapper to the interface that will
    likely be significantly changed or removed in a future release.
    This version is only intended to be used with `SEOBNRv5HM` and `SEOBNRv5PHM` and
    does not have full functionality for other waveform models.
    """

    from lalsimulation.gwsignal import GenerateFDWaveform
    from lalsimulation.gwsignal.models import gwsignal_get_waveform_generator
    import astropy.units as u

    waveform_kwargs = dict(
        waveform_approximant="SEOBNRv5PHM",
        reference_frequency=50.0,
        minimum_frequency=20.0,
        maximum_frequency=frequency_array[-1],
        catch_waveform_errors=False,
        mode_array=None,
        pn_amplitude_order=0,
    )
    waveform_kwargs.update(kwargs)

    waveform_approximant = waveform_kwargs['waveform_approximant']
    if waveform_approximant not in ["SEOBNRv5HM", "SEOBNRv5PHM"]:
        raise ValueError("The only supported approximants are SEOBNRv5PHM and SEOBNRv5HM.")

    reference_frequency = waveform_kwargs['reference_frequency']
    minimum_frequency = waveform_kwargs['minimum_frequency']
    maximum_frequency = waveform_kwargs['maximum_frequency']
    catch_waveform_errors = waveform_kwargs['catch_waveform_errors']
    mode_array = waveform_kwargs['mode_array']
    pn_amplitude_order = waveform_kwargs['pn_amplitude_order']

    if pn_amplitude_order != 0:
        # This is to mimic the behaviour in
        # https://git.ligo.org/lscsoft/lalsuite/-/blob/master/lalsimulation/lib/LALSimInspiral.c#L5542
        if pn_amplitude_order == -1:
            if waveform_approximant in ["SpinTaylorT4", "SpinTaylorT5"]:
                pn_amplitude_order = 3  # Equivalent to MAX_PRECESSING_AMP_PN_ORDER in LALSimulation
            else:
                pn_amplitude_order = 6  # Equivalent to MAX_NONPRECESSING_AMP_PN_ORDER in LALSimulation
        start_frequency = minimum_frequency * 2. / (pn_amplitude_order + 2)
    else:
        start_frequency = minimum_frequency

    # Call GWsignal generator
    wf_gen = gwsignal_get_waveform_generator(waveform_approximant)

    delta_frequency = frequency_array[1] - frequency_array[0]

    frequency_bounds = ((frequency_array >= minimum_frequency) *
                        (frequency_array <= maximum_frequency))

    iota, spin_1x, spin_1y, spin_1z, spin_2x, spin_2y, spin_2z = bilby_to_lalsimulation_spins(
        theta_jn=theta_jn, phi_jl=phi_jl, tilt_1=tilt_1, tilt_2=tilt_2,
        phi_12=phi_12, a_1=a_1, a_2=a_2, mass_1=mass_1 * utils.solar_mass, mass_2=mass_2 * utils.solar_mass,
        reference_frequency=reference_frequency, phase=phase)

    eccentricity = 0.0
    longitude_ascending_nodes = 0.0
    mean_per_ano = 0.0

    # Check if conditioning is needed
    condition = 0
    if wf_gen.metadata["implemented_domain"] == 'time':
        condition = 1

    # Create dict for gwsignal generator
    gwsignal_dict = {'mass1' : mass_1 * u.solMass,
                     'mass2' : mass_2 * u.solMass,
                     'spin1x' : spin_1x * u.dimensionless_unscaled,
                     'spin1y' : spin_1y * u.dimensionless_unscaled,
                     'spin1z' : spin_1z * u.dimensionless_unscaled,
                     'spin2x' : spin_2x * u.dimensionless_unscaled,
                     'spin2y' : spin_2y * u.dimensionless_unscaled,
                     'spin2z' : spin_2z * u.dimensionless_unscaled,
                     'deltaF' : delta_frequency * u.Hz,
                     'f22_start' : start_frequency * u.Hz,
                     'f_max': maximum_frequency * u.Hz,
                     'f22_ref': reference_frequency * u.Hz,
                     'phi_ref' : phase * u.rad,
                     'distance' : luminosity_distance * u.Mpc,
                     'inclination' : iota * u.rad,
                     'eccentricity' : eccentricity * u.dimensionless_unscaled,
                     'longAscNodes' : longitude_ascending_nodes * u.rad,
                     'meanPerAno' : mean_per_ano * u.rad,
                     # 'ModeArray': mode_array,
                     'condition': condition
                     }

    if mode_array is not None:
        gwsignal_dict.update(ModeArray=mode_array)

    if (
        domega220 != 0.0
        or domega330 != 0.0
        or domega210 != 0.0
        or domega440 != 0.0
        or domega550 != 0.0
        or domega320 != 0.0
        or domega430 != 0.0
    ):
        domega_dict = {'2,2': domega220,
                       '2,1': domega210,
                       '3,3': domega330,
                       '3,2': domega320,
                       '4,4': domega440,
                       '4,3': domega430,
                       '5,5': domega550,
                       }
        gwsignal_dict.update(domega_dict=domega_dict)

    if (
        dtau220 != 0.0
        or dtau330 != 0.0
        or dtau210 != 0.0
        or dtau440 != 0.0
        or dtau550 != 0.0
        or dtau320 != 0.0
        or dtau430 != 0.0
    ):
        dtau_dict = {'2,2': dtau220,
                     '2,1': dtau210,
                     '3,3': dtau330,
                     '3,2': dtau320,
                     '4,4': dtau440,
                     '4,3': dtau430,
                     '5,5': dtau550,
                     }
        gwsignal_dict.update(dtau_dict=dtau_dict)
    if (
        dA22 != 0.0
        or dA33 != 0.0
        or dA21 != 0.0
        or dA44 != 0.0
        or dA55 != 0.0
        or dA32 != 0.0
        or dA43 != 0.0
    ):
        dA_dict = {'2,2': dA22,
                   '2,1': dA21,
                   '3,3': dA33,
                   '3,2': dA32,
                   '4,4': dA44,
                   '4,3': dA43,
                   '5,5': dA55,
                   }
        gwsignal_dict.update(dA_dict=dA_dict)

    if (
        dw22 != 0.0
        or dw33 != 0.0
        or dw21 != 0.0
        or dw44 != 0.0
        or dw55 != 0.0
        or dw32 != 0.0
        or dw43 != 0.0
    ):
        dw_dict = {'2,2': dw22,
                   '2,1': dw21,
                   '3,3': dw33,
                   '3,2': dw32,
                   '4,4': dw44,
                   '4,3': dw43,
                   '5,5': dw55,
                   }
        gwsignal_dict.update(dw_dict=dw_dict)

    if dTpeak != 0.0:
        gwsignal_dict.update(dTpeak=dTpeak)

    if da6 != 0.0:
        gwsignal_dict.update(da6=da6)

    if ddSO != 0.0:
        gwsignal_dict.update(ddSO=ddSO)

    # Pass extra waveform arguments to gwsignal
    extra_args = waveform_kwargs.copy()

    for key in [
            "waveform_approximant",
            "reference_frequency",
            "minimum_frequency",
            "maximum_frequency",
            "catch_waveform_errors",
            "mode_array",
            "pn_spin_order",
            "pn_amplitude_order",
            "pn_tidal_order",
            "pn_phase_order",
            "numerical_relativity_file",
    ]:
        if key in extra_args.keys():
            del extra_args[key]

    gwsignal_dict.update(extra_args)

    try:
        hpc = GenerateFDWaveform(gwsignal_dict, wf_gen)
    except Exception as e:
        if not catch_waveform_errors:
            raise
        else:
            EDOM = (
                "Internal function call failed: Input domain error" in e.args[0]
            ) or "Input domain error" in e.args[
                0
            ]
            if EDOM:
                failed_parameters = dict(mass_1=mass_1, mass_2=mass_2,
                                         spin_1=(spin_1x, spin_1y, spin_1z),
                                         spin_2=(spin_2x, spin_2y, spin_2z),
                                         luminosity_distance=luminosity_distance,
                                         iota=iota, phase=phase,
                                         eccentricity=eccentricity,
                                         start_frequency=minimum_frequency)
                logger.warning("Evaluating the waveform failed with error: {}\n".format(e) +
                               "The parameters were {}\n".format(failed_parameters) +
                               "Likelihood will be set to -inf.")
                return None
            else:
                raise

    hplus = hpc.hp
    hcross = hpc.hc

    h_plus = np.zeros_like(frequency_array, dtype=complex)
    h_cross = np.zeros_like(frequency_array, dtype=complex)

    if len(hplus) > len(frequency_array):
        logger.debug("GWsignal waveform longer than bilby's `frequency_array`" +
                     "({} vs {}), ".format(len(hplus), len(frequency_array)) +
                     "probably because padded with zeros up to the next power of two length." +
                     " Truncating GWsignal array.")
        h_plus = hplus[:len(h_plus)]
        h_cross = hcross[:len(h_cross)]
    else:
        h_plus[:len(hplus)] = hplus
        h_cross[:len(hcross)] = hcross

    h_plus *= frequency_bounds
    h_cross *= frequency_bounds

    if condition:
        dt = 1 / hplus.df.value + hplus.epoch.value
        time_shift = np.exp(-1j * 2 * np.pi * dt * frequency_array[frequency_bounds])
        h_plus[frequency_bounds] *= time_shift
        h_cross[frequency_bounds] *= time_shift

    return dict(plus=h_plus, cross=h_cross)
