pSEOBNRv5 model
===============

``pSEOBNRv5`` is an extension to the ``SEOBNRv5HM`` and ``SEOBNRv5PHM`` models adding parametrized and generic GR deviations to the respective underlying models. 
For more details regarding the ``pSEOBNR`` model, see `arXiv:2205.05132 <https://arxiv.org/pdf/2205.05132>`_.

Implementation
--------------

The ``SEOBNRv5PHM`` waveform model describes binary-black-holes with spin precession, includes several higher harmonics, notably the ``(2,1), (3,3), (3,2), (4,4), (4,3), (5,5)``, in addition to the dominant ``(2,2)`` mode, and is calibrated to numerical-relativity simulations of binary-black-holes with aligned-spins.
Because this analysis directly modifies parameters in an inspiral-merger-ringdown waveform model, it takes full advantage of the signal length and its SNR.

The parameters of the deviation can independently be applied on:

* the quasi-normal-mode frequency for each mode as a fractional deviation
* the quasi-normal-mode damping time for each mode as a fractional deviation
* the merger amplitude for each mode as a fractional deviation
* the merger frequency for each mode as a fractional deviation
* the amplitude's peak time as an additive deviation
* the model's calibration parameter (``a6`` and ``dSO``) as an additive deviation
* A setting ``omega_prec_deviation`` (set to True by default) includes the fractional deviation in the J-frame QNM frequencies also into the precession rate computation (see Eq. 13 `arXiv:2301.06558 <https://arxiv.org/abs/2301.06558>`_)

The pipeline has been reviewed when including deviation parameters to the quasi-normal-mode frequency and dampling time of the fundamental ``(2,2,0)`` mode. The review wiki is available `here <https://git.ligo.org/cbc-testinggr/reviews/pseobnrv5phm>`_.
For this analysis, in addition to the general TGR selection criteria, an SNR greater than 8 is required in both the pre- and post-inspiral parts of the signal, since a reasonable inspiral SNR is needed to constrain the remnant properties expected in the GR prediction and break the degeneracy between the fundamental ringdown frequency deviation parameter and the remnant mass.

Usage
-----

After installing ``bilby_tgr``, you can use the ``bilby_pipe`` package to run an analysis 
(see the ``bilby_pipe`` `documentation <https://git.ligo.org/lscsoft/bilby_pipe>`_). 

To run using the ``pSEOBNRv5`` model test, a few lines need to be changed or added to the ``bilby_pipe`` configuration file. 
The things that need to be in a ``config.ini`` file to run pSEOB are:

* ``frequency-domain-source-model=bilby_tgr.pseob.source.gwsignal_binary_black_hole``
* ``waveform-approximant`` needs to be one of ``SEOBNRv5PHM`` or ``SEOBNRv5HM``
* Any additional waveform arguments, for example ``waveform-arguments-dict={'lmax_nyquist':2, 'omega_prec_deviation': True}``

A prior file is also needed: it is based on e.g. a standard precessing or aligned spin CBC prior 
plus the ``pSEOBNR`` deviations specific parameters one wants to sample over.

An example with ``.ini`` and ``.prior`` files can be found `here <https://git.ligo.org/lscsoft/bilby_tgr/-/blob/master/examples/pseob>`_. 
To run the example, enter the directory with the config and prior files and run

.. code-block:: console

   $ bilby_pipe pSEOB_config.ini --submit

