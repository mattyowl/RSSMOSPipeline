RSSMOSPipeline uses the `Robot Framework <http://robotframework.org/>`_ for tests.

.. note::  The files needed to run the test suite are found in the
           `tests <https://github.com/mattyowl/RSSMOSPipeline/tree/master/tests>`_
           directory of the **RSSMOSPipeline** source distribution.

To run a set of tests that should take about X minutes:

.. code-block::

   robot tests.robot

To run a single test in a test suite do, e.g., 

.. code-block::

   robot -t "MOS reduction using slit file" tests.robot

Check the ``plots/`` directory for output from some tests.

