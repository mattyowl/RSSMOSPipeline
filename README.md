# RSSMOSPipeline
Pipeline for reducing both longslit and multi-object spectroscopy from the Robert Stobie Spectrograph
on [SALT](https://www.salt.ac.za/).

* **Documentation**: Link to be added
* **License**: [GPL v3](https://github.com/simonsobs/n/blob/main/LICENSE)
* **Authors**: Matt Hilton, Melissa Moris
* **Installation**: `pip install RSSMOSPipeline`
* **Support**: Please use the [GitHub issues page](https://github.com/mattyowl/RSSMOSPipeline/issues), or contact [Matt Hilton](mailto:matt.hilton@wits.ac.za)

Please note this software is under development at the moment, and the instructions in this 
README file may not always be up to date.

## Important notes
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

More can be added relatively easily (see the `rss_mos_create_arc_model` script, which will 
print instructions when you run it).

The quickest way to check the wavelength calibration is to inspect the `arcTransformTest*.png` and 
`skyCheck_SLIT*.png` files, found under `reducedDir/OBJECT_MASKID/diagnostics/`. The former 
show the transformed arc spectrum compared to the reference model, while the latter plot the sky 
signal extracted from each science frame for each slit, with the reference wavelengths of prominent
sky lines shown as vertical dashed lines for comparison. If these don't line up, then an additional
(or revised) reference model is needed. 

The size in Angstroms of any offset from prominent sky lines in each final stacked spectrum is logged in 
`reducedDir/OBJECT_MASKID/diagnostics/skyWavelengthCalibCheck.csv`. As well as the median offset
of identified sky lines with respect to their reference wavelengths, the results averaged across all
slits are listed at the bottom of this file (look for "all slits median" and "all slits RMS"). This 
test shows that wavelength calibration using the pg0900 Ne (2x2 binning) reference model is good to 
< 1 Angstrom. Other reference models have not been tested extensively yet.

The pipeline has not yet been optimised for speed. On an Intel Core i5-3320M (2.60 GHz), it currently 
(Oct 2016) takes 30 minutes to run a mask with 28 slits using the iterative sky subtraction method,
or 12 minutes with the non-iterative sky subtraction method.

## Software needed
The pipeline is written in pure python (only 3.x is supported now). It needs the following modules to 
be installed:

* numpy
* scipy
* astropy
* matplotlib
* IPython

IPython is used for debugging, but isn't really needed to run the pipeline. The install script (see 
below) should install the needed modules automatically if they are not already on your system.

## Installation
From PyPI:

```
pip install RSSMOSPipeline
```

(note: if using, e.g., Ubuntu, then you would need to add `$HOME/.local/bin` to `$PATH`, as
this is where scripts like `rss_mos_reducer` are installed).

Or, as root, from the source code archive:
    
```
sudo python setup.py install
```

Or, in your home directory:
    
```
python setup.py install --prefix=$HOME/.local
```

Then add `$HOME/.local/bin` to `$PATH`, and e.g., `$HOME/.local/lib/python3.6/site-packages` to $PYTHONPATH.

```
export PATH=$HOME/.local/bin:$PATH    
export PYTHONPATH=$HOME/.local/lib/python3.6/site-packages:$PYTHONPATH
```

Or, possibly:

```
python setup.py install --user
```

(this works on Ubuntu).

## How to run

1. Download and unpack the archive containing your SALT data. This pipeline will operate on the
   contents of the `product` directory.

2.  Run the code, e.g.,

    ```
    rss_mos_reducer product reducedDir maskName
    ```

    `maskName` is the combination of `OBJECT_MASKID` from the corresponding FITS header keywords.
    
    So, for example if your data for an object named J0058.1+0031 were taken with a mask called 
    P001788N07 and live in a directory called `product` then you would run the pipeline with:

    ```
    rss_mos_reducer product reduced J0058.1+0031_P001788N07
    ```
    
    If you want to just list the masks found under a directory called `product` (without doing anything
    with the data) you can use:
    
    ```
    rss_mos_reducer product reduced list
    ```
    
    And if you're super confident that nothing can go wrong, you can ask the pipeline to process
    data for all the masks it can find (in this case, under the `product` directory)
    
    ```
    rss_mos_reducer product reduced all
    ```
    
    The processed data from each mask will be found in its own subdirectory under `reduced/` in the above
    example.
    
    You can run the code with the `-h` switch to see a help message with more options. The most 
    important are `-t`, for setting the threshold used for identifying MOS slitlets in the flat
    field frames (only needed if you find slits are missing - use the DS9 `.reg` files to check), 
    and `-i`, which switches on the iterative sky subtraction method for the spectral extraction.
    This currently seems to work better than the default extraction method, but is slower.

3.  After the pipeline has finished, you will find data at various stages of reduction under the `reducedData` 
    dir (`reduced` for the above example). At the moment, this is not cleaned up.

    Two methods are used for extracting spectra:
    
    1.  Individual 1d spectra are extracted from each science frame and then stacked. This might be 
        desirable if object traces shift from frame-to-frame. These are placed under 
        `reducedDir/OBJECT_MASKID/1DSpec_extractAndStack/`
        
    2.  The 2d spectra are stacked and then a 1d spectrum is extracted. These are placed under 
        `reducedDir/OBJECT_MASKID/1DSpec_2DSpec_stackAndExtract/`. The stacked 2d spectra are also 
        placed in this directory.
        
    Currently method (2) works best (Oct 2016).
    
    If the iterative sky subtraction method is enabled (using the `-i` switch), then 
    `_iterative` is appended to the name of the output directory 
    (e.g., `1DSpec_2DSpec_stackAndExtract_iterative`).
    
    File names include the slit number (e.g., `*_SLIT10.fits`). You can check which slit (identified 
    from the flat frame) corresponds with which extracted spectrum by examining the 
    `masterFlat_?.fits` file(s) and the corresponding DS9 region file(s) (ending `.reg`) in the `reduced/OBJECT_MASKID` directory.

    The 1d spectra are in FITS table format (stored in the extension `1D_SPECTRUM`) with the following 
    columns:

    * `SPEC` 		- the object spectrum
    * `SKYSPEC` 	- the spectrum of the sky
    * `LAMBDA`	- the wavelength scale corresponding to both of the above (in Angstroms)
    
    The 2d spectra are `.fits` images, which can be examined with DS9. They include the wavelength
    solution as the WCS stored in the header.

4.  The `rss_mos_visual_inspector` tool can be used to look at the 1D spectra, compare to template
    spectra from SDSS, and measure redshifts. This needs the `tkinter` module to be installed, and
    for the `matplotlib` backend to be set to `TkAgg`.

    Run this using, e.g.,

    ```
    rss_mos_visual_inspector 1DSpec_2DSpec_stackAndExtract/1D_*.fits results
    ```

    In this example, each spectrum will be displayed in turn, and a text file containing the results
    of redshift measurements (and any comments by the user) will be written in the `results/`
    directory, together with a plot of each spectrum. These are recorded from the state when the user
    clicks the `Done - show next` button in the top right of the graphical user interface. Hopefully
    it should be fairly obvious how to use most of the plotting controls. Note that you must click
    `Redraw plot` after making changes to either the template redshift or spectral smoothing.
    
    Note that the `XC Galaxies`, `XC LRGs` and `XC QSOs` buttons will not work unless you have IRAF
    and PyRAF installed (and even then the canned settings may not work well on all setups).

## Longslit mode
The pipeline can also run on longslit data. The code detects object traces in each science frame,
and creates "pseudo-slitlets" around each detected object. The processing steps are otherwise identical to those
for MOS data. If necessary, you can use the --longslit-threshold option to change the detection threshold used.

## Additional flags

### Manual slit locations
If you know the position of the slits in the fits file, you can specify them in an ascii file that will be read 
in using astropy Tables and used to extract slits from the data. To do this, you can do this using the `-F`
argument, e.g.,

```
rss_mos_reducer -F slit_loc.txt product reduced all
```

This ignores the automated slit finding routine and uses the input location. Be sure to structure the slit
location files with three columns: `slitno` (i.e. slit number), `ystart` (i.e. y-coord coinciding with beginning
of slit), `yend` (i.e. y-coord coinciding with end of slit).

### No flats? No problem!
If you don't have flat fields for your observations, you can skip flat fielding altogether with the `-n` argument. To do
this, you *must* also specify the location of your slits using the `-F` argument (see above).

## Things which can/should be improved/added
* Extraction (not optimal at the moment)
* Spectrophotometric calibration using a standard (not implemented)

## Comments, bug reports, help, suggestions etc..
Please contact [Matt Hilton](mailto:matt.hilton@wits.ac.za), or [open an issue on GitHub](https://github.com/mattyowl/RSSMOSPipeline/issues).
