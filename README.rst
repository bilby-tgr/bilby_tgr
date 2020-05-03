|pipeline status| |coverage report|

=========
Bilby_TGR
=========

This is the Bilby Testing General Relativity (TGR) package. It is used to
develop and share analysis done by the LVC TGR group building on the
`Bilby <https://git.ligo.org/lscsoft/bilby>`_ stochastic sampling package.


Installation
------------

To get started developing and using :code:`bilby_tgr`, first clone the repository

.. code-block:: console

   $ git clone git@git.ligo.org:gregory.ashton/bilby_tgr.git
   
To install, enter the cloned directory and run

.. code-block:: console

   $ pip install .
   
To check that you have installed it correctly, open an python prompt and run

.. code-block:: python

   >>> import bilby_tgr
   >>> bilby_tgr.source.non_gr_d_alpha_2_binary_black_hole

If this returns the function, then you have it installed! You can now add new
functions to the sources and access them in the same way.

Running using bilby_pipe
------------------------

Once you have installed :code:`bilby_tgr`, you can use the :code:`bilby_pipe`
package to run stoachastic sampling. For help getting installed and setup with
:code:`bilby_pipe` itself, see `the documentation <https://git.ligo.org/lscsoft/bilby_pipe>`_.
Here, we give an example ini file. Notice that the :code:`frequency-domain-source-model`
is pointing to the :code:`bilby_tgr.source.non_gr_d_alpha_2_binary_black_hole`
function. You can replace this with any new function you care to define in the
:code:`bilby_tgr` package.

.. code-block:: console

   trigger-time = 1186741861.5
   detectors = [H1, L1, V1]
   channel-dict = {H1:DCH-CLEAN_STRAIN_C02, L1:DCH-CLEAN_STRAIN_C02, V1:Hrec_hoft_V1O2Repro2A_16384Hz}
   prior_file = 4s.prior
   time-marginalization=False
   distance-marginalization=True
   phase-marginalization=True
   create-plots=True
   local-generation = True
   psd-dict = {H1:BayesWave_median_PSD_H1.dat, L1:BayesWave_median_PSD_L1.dat, V1:BayesWave_median_PSD_V1.dat}

   label = GW170814
   outdir = d_alpha_2
   accounting = ligo.dev.o3.cbc.pe.lalinference

   duration = 4
   coherence-test = False

   maximum-frequency=1024
   minimum-frequency=20
   sampling-frequency=2048
   reference-frequency = 20
   waveform-approximant = IMRPhenomPV2
   frequency-domain-source-model = bilby_tgr.source.non_gr_d_alpha_2_binary_black_hole

   calibration-model=CubicSpline
   spline-calibration-envelope-dict = {H1:GWTC1_GW170814_H_CalEnv.txt, L1:GWTC1_GW170814_L_CalEnv.txt, V1:GWTC1_GW170814_V_CalEnv.txt}
   spline_calibration-nodes = 10

   deltaT = 0.2
   sampler = dynesty
   sampler-kwargs = {nlive: 1000, sample: rwalk, walks=50, nact=5}
   n-parallel = 4

   transfer-files = False
   
   

.. |pipeline status| image:: https://git.ligo.org/gregory.ashton/bilby_tgr/badges/master/pipeline.svg
   :target: https://git.ligo.org/gregory.ashton/bilby_tgr/commits/master
.. |coverage report| image:: https://gregory.ashton.docs.ligo.org/bilby_tgr/coverage_badge.svg
   :target: https://gregory.ashton.docs.ligo.org/bilby_tgr/htmlcov/
