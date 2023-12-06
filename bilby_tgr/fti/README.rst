========
FTI test
========

The FTI (Flexible Theory-Independent) test of GR adds generic deviations to the PN coefficients of the frequency-domain 
GW phase during inspiral, including deviations of the spin-induced quadrupole moment (SIQM).

For more details regarding the FTI test, see `arXiv:2203.13937 <https://arxiv.org/abs/2203.13937>`_.


Implementation
--------------

Generic deviations to the inspiral are added as a fractional correction to the SPA phase at different PN orders. 
The following corrections coefficients have been implemented in the FTI code: 
-1PN: :code:`dchi_Minus2`; 0PN: :code:`dchi_0`; 0.5PN: :code:`dchi_1`; 1PN: :code:`dchi_2`; 
1.5PN: :code:`dchi_3` / :code:`dchi_3NS` / :code:`dchi_3S`; 2PN: :code:`dchi_4` / :code:`dchi_4NS` / :code:`dchi_4S`; 
2.5PN: :code:`dchi_5l` / :code:`dchi_5lNS` / :code:`dchi_5lS`; 3PN: :code:`dchi_6` / :code:`dchi_6NS` / :code:`dchi_6S`, :code:`dchi_6l`; 
3.5PN: :code:`dchi_7` / :code:`dchi_7NS` / :code:`dchi_7S`,
which, from 1.5PN onwards, can be normalized to include effects in the full/nonspinning/spin-only sectors respectively.

FTI can also be used for spin-induced quadrupole moment tests via the symmetric and anti-symmetric combinations 
of the kappa-coefficients (:code:`dchi_kappaS`, :code:`dchi_kappaA`). The SIQM corrections are included in the 
inspiral phase at 2PN and 3PN (see `arXiv:1701.06318 <https://arxiv.org/abs/1701.06318>`_ for the coefficients).

The corrections to the inspiral phase are smoothly tapered of towards merger. The tapering function is controlled by 
two parameters: :code:`f_window_div_f_peak`, the frequency at which is the tapering happens as fraction of the peak 
frequency of the 22-mode, and :code:`NumCycles`, the number of inspiral cycles over which the tapering happens.

The FTI test is implemented for non-precessing frequency domain waveforms. More precisely, the following approximants 
are supported for BBH: :code:`SEOBNRv4HM_ROM`, :code:`SEOBNRv4_ROM`, :code:`SEOBNRv5_ROM`, and :code:`IMRPhenomXAS`. 
For :code:`SEOBNRv4HM_ROM`, it is possible to specify how many modes are used in the waveform using the waveform argument :code:`n_modes` (default is 4). 
This takes the first :code:`n_modes` in the list of modes :code:`[[2, 2], [3, 3], [2, 1], [4, 4], [5, 5]]`. 
For BNS and NSBH the following approximants are supported: 
:code:`SEOBNRv4_ROM_NRTidalv2`, :code:`SEOBNRv4T_surrogate`, :code:`IMRPhenomXAS_NRTidalv2`, and :code:`SEOBNRv4_ROM_NRTidalv2_NSBH`.

Note that for the ROM waveforms additional ROM files are needed and the environment variable :code:`LAL_DATA_PATH` 
needs to be set to point to their location. The ROM files can be downloaded from the :code:`lalsuite-extra` package 
`here <https://git.ligo.org/lscsoft/lalsuite-extra/-/tree/master/data/lalsimulation>`_. They are also available on OSG clusters through cvmfs.
To use this, one has to set the :code:`LAL_DATA_PATH` environment variable before submitting a job:

.. code-block:: console
   
   $ export LAL_DATA_PATH=/cvmfs/singularity.opensciencegrid.org/.images/4f/8299f8f688d0175846129a95b353b07c452dd708a4c5fb8427a5f867610196/opt/lal_data:/cvmfs/singularity.opensciencegrid.org/.images/4c/6d00a32836cc022a244b3683316670784eecb59993be1c6a3f6a74f6f86a3f/opt/lal_data


Usage
-----

After installing :code:`bilby_tgr`, you can use the :code:`bilby_pipe` package to run an analysis 
(see the :code:`bilby_pipe` `documentation <https://git.ligo.org/lscsoft/bilby_pipe>`_). 
To run using the FTI test, a few lines need to be changed or added to the :code:`bilby_pipe` configuration file. 
The things that need to be in a :code:`config.ini` file to run FTI are:

* :code:`frequency-domain-source-model=bilby_tgr.fti.source.lal_binary_black_hole`
* :code:`default-prior=bilby_tgr.fti.prior.BBHPriorDict`
* :code:`conversion-function=bilby_tgr.fti.conversion.convert_to_lal_binary_black_hole_parameters`
* :code:`generation-function=bilby_tgr.fti.conversion.generate_all_bbh_parameters`
* :code:`waveform-approximant` needs to be one of :code:`SEOBNRv4HM_ROM`, :code:`SEOBNRv4_ROM`, :code:`SEOBNRv5_ROM`, :code:`IMRPhenomXAS`
* The number of modes used with the :code:`SEOBNRv4HM_ROM` waveform can be specified using the line :code:`waveform-arguments-dict={'n_modes': 4}`

The other thing that is needed is a prior file. The prior file is a standard aligned spin CBC prior 
plus the FTI parameters one wants to sample over.

An example with :code:`.ini` and :code:`.prior` files can be found `here <https://git.ligo.org/lscsoft/bilby_tgr/-/blob/master/examples/fti>`_. 
To run the example, enter the directory with the config and prior files and run

.. code-block:: console

   $ bilby_pipe GW170814_dchi0.ini --submit
