from __future__ import division, print_function

import bilby
import numpy as np

import bilby_tgr

outdir = "outdir"
label = "siqm_example"
bilby.core.utils.setup_logger(outdir=outdir, label=label)
np.random.seed(88170235)
injection_parameters = dict(
    mass_1=20.0,
    mass_2=10.0,
    a_1=0.7,
    a_2=0.3,
    tilt_1=0.5,
    tilt_2=1.0,
    phi_12=1.7,
    phi_jl=0.3,
    luminosity_distance=200.0,
    theta_jn=0.4,
    psi=2.659,
    phase=1.3,
    geocent_time=1126259642.413,
    ra=1.375,
    dec=-1.2108,
    dQuadMon1=20,
    dQuadMon2=0,
)
injection_parameters
duration = 4
sampling_frequency = 2 * 1024
start_time = injection_parameters["geocent_time"] + 2 - duration
waveform_arguments = dict(
    waveform_approximant="IMRPhenomPv2",
    reference_frequency=50.0,
    minimum_frequency=40.0,
)


waveform_generator = bilby.gw.WaveformGenerator(
    duration=duration,
    sampling_frequency=sampling_frequency,
    frequency_domain_source_model=bilby_tgr.siqm.source.lal_binary_black_hole,
    parameter_conversion=bilby_tgr.siqm.conversion.convert_to_lal_siqm_parameters,
    waveform_arguments=waveform_arguments,
)

# (LIGO-Hanford (H1), LIGO-Livingston (L1), and Virgo (V1)).
# These default to their design sensitivity and start at 40 Hz.

ifos = bilby.gw.detector.InterferometerList(["H1", "L1"])
ifos.set_strain_data_from_power_spectral_densities(
    sampling_frequency=sampling_frequency,
    duration=duration,
    start_time=injection_parameters["geocent_time"] - 3,
)
ifos.inject_signal(
    waveform_generator=waveform_generator, parameters=injection_parameters
)

priors = bilby.gw.prior.BBHPriorDict()
priors["geocent_time"] = bilby.core.prior.Uniform(
    minimum=injection_parameters["geocent_time"] - 1,
    maximum=injection_parameters["geocent_time"] + 1,
    name="geocent_time",
    latex_label="$t_c$",
    unit="$s$",
)
for key in [
    "a_1",
    "a_2",
    "tilt_1",
    "tilt_2",
    "phi_12",
    "phi_jl",
    "psi",
    "ra",
    "dec",
    "geocent_time",
    "phase",
]:
    priors[key] = injection_parameters[key]
priors["dQuadMonS"] = bilby.core.prior.Uniform(0, 500, name="dQuadMonS")
priors["dQuadMonA"] = 0
likelihood = bilby.gw.GravitationalWaveTransient(
    interferometers=ifos, waveform_generator=waveform_generator
)

result = bilby.run_sampler(
    likelihood=likelihood,
    priors=priors,
    sampler="dynesty",
    npoints=1000,
    injection_parameters=injection_parameters,
    outdir=outdir,
    label=label,
)

# Make a corner plot.
result.plot_corner()
