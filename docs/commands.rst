.. _Usage:

=======================
RSSMOSPipeline Commands
=======================

The **RSSMOSPipeline** package includes a number of command-line programs, each of which is described below.


.. _reducerCommand:
    
rss_mos_reducer
---------------

.. argparse::
   :filename: ../bin/rss_mos_reducer
   :func: makeParser
   :prog: rss_mos_reducer
   
   :program:`rss_mos_reducer` reduces SALT RSS data taken in either longslit or multi-object spectroscopy (MOS) modes
   to extracted, sky subtracted spectra.


.. _createArcModelCommand:

rss_mos_create_arc_model
------------------------

.. argparse::
   :filename: ../bin/rss_mos_create_arc_model
   :func: makeParser
   :prog: rss_mos_create_arc_model
   
   :program:`rss_mos_create_arc_model` is a tool for creating reference models for wavelength calibration from arc spectra.


.. _inspectArcModelCommand:

rss_mos_inspect_arc_model
-------------------------

.. argparse::
   :filename: ../bin/rss_mos_inspect_arc_model
   :func: makeParser
   :prog: rss_mos_inspect_arc_model

   :program:`rss_mos_inspect_arc_model` is a tool for inspecting reference models used for wavelength calibration from arc spectra.

