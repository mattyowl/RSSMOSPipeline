# RSSMOSPipeline
Pipeline for reducing multi-object spectroscopy from the Robert Stobie Spectrograph on SALT.

## Important notes
For wavelength calibration, a reference model needs to be made. So far this has been done for
the following grating/lamp/detector binning combinations:

* pg0900 Ar (2x2 binning)
* pg0900 Ne (2x2 binning)
* pg0900 Xe (2x2 binning)
* pg1300 Ar (2x2 binning)
* pg1800 Ne (2x2 binning)

More can be added relatively easily (see the modelArcSpectra dir), but no documentation on 
this yet.

## Software needed

The pipeline is written in pure python (2.7.x). It needs the following modules to be installed:

* numpy
* scipy
* atpy
* pyfits
* matplotlib

IPython is used for debugging, but isn't needed to run the pipeline. So you can comment out
```python
import IPython
``` 
from the `rss_mos_reducer.py` script if you don't want to install IPython.

Instead of pyfits and atpy, astropy can potentially be used, but this isn't implemented yet.

## How to run

1. Unpack the archive containing the pipeline code (`rss_mos_reducer.py`) and the `modelArcSpectra`
folder to some directory.

2. Copy your raw data to the same directory, placing it in a folder (e.g., `rawData/`). This 
should contain the contents of the `product/` dir, as provided by SALT.

3.  Run the code, e.g.,

    ```
    python rss_mos_reducer.py rawData reducedData maskName
    ```

    You can find maskName from the FITS headers for your data - look for the `MASKID` keyword.
    
    So, for example if your data were taken with a mask called P000807N01 and live in a directory
    called `ACTTest` then you would run the pipeline with:

    ```
    python rss_mos_reducer.py ACTTest reducedACTTest P000807N01
    ```
    
    If you want to just list the masks found under a directory called `product` (without doing anything
    with the data) you can use:
    
    ```
    python rss_mos_reducer.py product reduced list
    ```
    
    And if you're super confident that nothing can go wrong, you can ask the pipeline to process
    data for all the masks it can find (in this case, under the `product` directory)
    
    ```
    python rss_mos_reducer.py product reduced all
    ```
    
    The processed data from each mask will be found in its own subdirectory under `reduced/` in the above
    example.

4.  After the pipeline has finished, you will find data at various stages of reduction under the
`   reducedData` dir (`reducedACTTest` for the above example). At the moment, this is not cleaned up.

    The extracted 1d spectra can be found under (for the example above) `reducedACTTest/P000807N01/1DSpec/`.
    They include the slit number in the file name (e.g., `*_SLIT10.fits`). You can check which slit
    (identified from the flat frame) corresponds with which extracted spectrum by examining the
    `masterFlat_?.fits` file(s) and the corresponding DS9 region file(s) (ending .reg) in the `reducedData/`
    directory.

    The 1d spectra are in FITS table format (stored in the extension `1D_SPECTRUM`) with the following 
    columns:

    * `SPEC` 		- the object spectrum
    * `SKYSPEC` 	- the spectrum of the sky
    * `LAMBDA`	- the wavelength scale corresponding to both of the above (in Angstroms)

## Things which can/should be improved/added
* flat fielding (polynomial fit probably not optimal)
* extraction (not optimal at the moment)
* spectrophotometric calibration using a standard (not implemented)

## Comments, bug reports, help, suggestions etc..
Please contact Matt Hilton <hiltonm@ukzn.ac.za>. I am happy to tweak this to work with as many
gratings combinations as needed, but will need data to work with (and time).

