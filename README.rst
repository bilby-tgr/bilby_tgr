|pipeline status| |coverage report|

=========
Bilby_TGR
=========

This is the Bilby Testing General Relativity (TGR) package. It is used for
analyses done by the LVK TGR group building on the
`Bilby <https://github.com/bilby-dev/bilby>`_ stochastic sampling package.

The package includes the following tests of GR:

* FTI
* LOSA
* MDR
* pSEOB
* SIQM
* SSB
* TIGER

Installation
------------

To install the latest released version of :code:`bilby_tgr` from `pypi <https://pypi.org/project/bilby-tgr>`_, simply run

.. code-block:: console

   $ pip install bilby_tgr

To check that you have installed it correctly, open a python prompt and run

.. code-block:: python

   >>> import bilby_tgr
   >>> bilby_tgr.tiger.source.lal_binary_black_hole_TIGER_PhenomP

If this returns the function, then you have it installed!

To install :code:`bilby_tgr` from source, first clone the repository:

.. code-block:: console

   $ git clone https://github.com/bilby-tgr/bilby_tgr.git

Then enter the cloned directory and run

.. code-block:: console

    $ cd bilby_tgr/
    $ pip install .


Running using bilby_pipe
------------------------

Once you have installed :code:`bilby_tgr`, you can use the :code:`bilby_pipe`
package to run stochastic sampling. For help getting installed and setup with
:code:`bilby_pipe` itself, see `the documentation <https://git.ligo.org/lscsoft/bilby_pipe>`_.

In the `examples` directory, we give example ini files for the different testing GR pipelines.
Notice that the :code:`frequency-domain-source-model` is pointing to a source function in this package,
e.g. :code:`bilby_tgr.tiger.source.lal_binary_black_hole_TIGER_PhenomP`.
You can replace this with any new function you care to define in the :code:`bilby_tgr` package.

The other thing you need is a prior file (4s.prior in the above ini).
This will be a standard CBC prior, plus any new parameters.


Acknowledging/Citing Bilby TGR
------------------------------

If you have used Bilby TGR in your scientific work, please acknowledge us in your papers/proposals.
Please cite the appropriate version from https://doi.org/10.5281/zenodo.10940209.

Additionally, please cite :code:`bilby` appropriately as described at https://lscsoft.docs.ligo.org/bilby/citing-bilby.html. 
If you used the automation module :code:`bilby_tgr.asimov`, please additionally cite asimov and cbcflow.


.. |pipeline status| image:: https://git.ligo.org/lscsoft/bilby_tgr/badges/master/pipeline.svg
   :target: https://git.ligo.org/lscsoft/bilby_tgr/commits/master
.. |coverage report| image:: https://lscsoft.docs.ligo.org/bilby_tgr/coverage_badge.svg
   :target: https://lscsoft.docs.ligo.org/bilby_tgr/htmlcov/
