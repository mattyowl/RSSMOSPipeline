
=====
Usage
=====

Running the pipeline
====================

1. Download and unpack the archive containing your SALT data. This pipeline will operate on the
   contents of the ``product`` directory.


2. Run the code, e.g.::

    rss_mos_reducer product reducedDir maskName

   ``maskName`` is the combination of ``OBJECT_MASKID`` from the corresponding FITS header keywords.
   So, for example if your data for an object named J0058.1+0031 were taken with a mask called
   P001788N07 and live in a directory called ``product`` then you would run the pipeline with::

    rss_mos_reducer product reduced J0058.1+0031\_P001788N07

   If you want to just list the masks found under a directory called ``product`` (without doing anything
   with the data) you can use::

    rss_mos_reducer product reduced list

   And if you're super confident that nothing can go wrong, you can ask the pipeline to process
   data for all the masks it can find (in this case, under the ``product`` directory)::

    rss_mos_reducer product reduced all

   The processed data from each mask will be found in its own subdirectory under ``reduced/`` in the above
   example.

   You can run the code with the ``-h`` switch to see a help message with more options. The most
   important are ``-t``, for setting the threshold used for identifying MOS slitlets in the flat
   field frames (only needed if you find slits are missing \- use the DS9 ``.reg`` files to check),
   and ``-i``, which switches on the iterative sky subtraction method for the spectral extraction.
   This currently seems to work better than the default extraction method, but is slower.


3. After the pipeline has finished, you will find data at various stages of reduction under the
   ``reducedData`` dir (``reduced`` for the above example). At the moment, this is not cleaned up.

   Two methods are used for extracting spectra:

   1. Individual 1D spectra are extracted from each science frame and then stacked. This might be
      desirable if object traces shift from frame-to-frame. These are placed under::

        `reducedDir/OBJECT_MASKID/1DSpec_extractAndStack/`

   2. The 2D spectra are stacked and then a 1d spectrum is extracted. These are placed under::

        `reducedDir/OBJECT_MASKID/1DSpec_2DSpec_stackAndExtract/`.

      The stacked 2D spectra are also placed in this directory.

   Currently method 2 seems to work best.

   If the iterative sky subtraction method is enabled (using the ``-i`` switch), then
   ``_iterative`` is appended to the name of the output directory (e.g., ``1DSpec_2DSpec_stackAndExtract_iterative``).

   File names include the slit number (e.g., ``*_SLIT10.fits``). You can check which slit
   (identified from the flat frame) corresponds with which extracted spectrum by examining the
   ``masterFlat_?.fits`` file(s) and the corresponding DS9 region file(s) (ending ``.reg``) in
   the ``reduced/OBJECT_MASKID`` directory.

   The 1D spectra are in FITS table format (stored in the extension ``1D_SPECTRUM``) with the following
   columns:

   * ``SPEC`` - the object spectrum
   * ``SKYSPEC`` - the spectrum of the sky
   * ``LAMBDA`` - the wavelength scale corresponding to both of the above (in Angstroms)

   The 2D spectra are FITS images, which can be examined with DS9. They include the wavelength
   solution as the WCS stored in the header.


4. The ``rss_mos_visual_inspector`` tool can be used to look at the 1D spectra, compare to template
   spectra from SDSS, and measure redshifts. This needs the ``tkinter`` module to be installed, and
   for the ``matplotlib` backend to be set to ``TkAgg``.

   Run this using, e.g.::

    rss_mos_visual_inspector 1DSpec_2DSpec_stackAndExtract/1D_*.fits results

   In this example, each spectrum will be displayed in turn, and a text file containing the results
   of redshift measurements (and any comments by the user) will be written in the ``results/``
   directory, together with a plot of each spectrum. These are recorded from the state when the user
   clicks the ``Done - show next`` button in the top right of the graphical user interface. Hopefully
   it should be fairly obvious how to use most of the plotting controls. Note that you must click
   ``Redraw plot`` after making changes to either the template redshift or spectral smoothing.

   .. note:: The ``XC Galaxies``, ``XC LRGs`` and ``XC QSOs`` buttons will not work unless you have IRAF
             and PyRAF installed (and even then the canned settings may not work well on all setups).


Wavelength calibration
======================

For wavelength calibration, a reference model needs to be made. So far this has been done for
the following grating/lamp/detector binning combinations:

* pg0900 Ar (2x2 binning)
* pg0900 Ne (1x2 binning)
* pg0900 Ne (2x2 binning)
* pg0900 Xe (2x2 binning)
* pg1300 Ar (2x2 binning)
* pg1300 Ar (1x2 binning)
* pg1300 Ne (1x2 binning)
* pg1300 Xe (1x2 binning)
* pg1300 CuAr (2x2 binning)
* pg1800 Ne (2x2 binning)

More can be added relatively easily (see the :ref:`createArcModelCommand` script, which will
print instructions when you run it).

The quickest way to check the wavelength calibration is to inspect the ``arcTransformTest*.png`` and
``skyCheck*SLIT*.png`` files, found under ``reducedDir/OBJECT_MASKID/diagnostics/``. The former
show the transformed arc spectrum compared to the reference model, while the latter plot the sky
signal extracted from each science frame for each slit, with the reference wavelengths of prominent
sky lines shown as vertical dashed lines for comparison. If these don't line up, then an additional
(or revised) reference model is needed.

The size in Angstroms of any offset from prominent sky lines in each final stacked spectrum is logged in
``reducedDir/OBJECT_MASKID/diagnostics/skyWavelengthCalibCheck.csv``. As well as the median offset
of identified sky lines with respect to their reference wavelengths, the results averaged across all
slits are listed at the bottom of this file (look for "all slits median" and "all slits RMS"). This
test shows that wavelength calibration using the pg0900 Ne (2x2 binning) reference model is good to
< 1 Angstrom. Other reference models have not been tested extensively yet.


Longslit mode
=============

The pipeline can also run on longslit data. The code detects object traces in each science frame,
and creates "pseudo-slitlets" around each detected object. The processing steps are otherwise identical to those
for MOS data. If necessary, you can use the ``--longslit-threshold`` option to change the detection threshold used.


Manually specifying slit (or object spectrum) locations
=======================================================

If you know the position of the slits in the FITS file, and/or the automated object/slit finding is not
working well on your data for some reason, you can specify them in a plain-text file.
To do this, you can do this using the ``-F`` argument, e.g.::

    rss_mos_reducer -F slit_loc.txt product reduced all

This bypasses the automated slit/object finding routine and uses the input locations given in the text file.
Be sure to structure the slit location files with three columns: ``slitno`` (i.e., slit number),
``ystart`` (i.e., y-coord coinciding with beginning of slit), and ``yend`` (i.e., y-coord coinciding with end of slit).

Below is a very simple example of a manual slit location file, for just one object, found in rows 988--1026
of some SALT RSS image::

    slitno   ystart   yend
    10       988      1026


No flats? No problem!
=====================

If you don't have flat fields for your observations, you can skip flat fielding altogether with the ``-n`` argument.
If no flat fields are found for your observations, the pipeline will attempt to drop into this mode automatically.

Note that if you do not have flats for MOS data, you *must* also specify the location of your slits using
the ``-F`` argument (see above).


Things which can/should be improved/added
===========================================

* Extraction (not optimal at the moment)
* Spectrophotometric calibration using a standard (not implemented)
