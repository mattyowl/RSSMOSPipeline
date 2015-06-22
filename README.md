# RSSMOSPipeline
Pipeline for reducing multi-object spectroscopy from the Robert Stobie Spectrograph on SALT.

## Important notes
For wavelength calibration, a reference model needs to be made. So far this has been done for
the following grating/lamp/detector binning combinations:

* pg0900 Ar (2x2 binning)
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
'import IPython' from the rss\_mos\_reducer.py script if you don't want to install IPython.

Instead of pyfits and atpy, astropy can potentially be used, but this isn't implemented yet.

## How to run

1. Unpack the archive containing the pipeline code (rss\_mos\_reducer.py) and the modelArcSpectra
folder to some directory.

2. Copy your raw data to the same directory, placing it in a folder (e.g., rawData/). This 
should contain the contents of the product/ dir, as provided by SALT.

3. Run the code, e.g.,

'''
python rss\_mos\_reducer.py rawData reducedData maskName
'''

You can find maskName from the FITS headers for your data - look for the MASKID keyword.

So, for example if your data were taken with a mask called P000807N01 and live in a directory
called ACTTest then you would run the pipeline with:

'''
python rss\_mos\_reducer.py ACTTest reducedACTTest P000807N01
'''

4. After the pipeline has finished, you will find data at various stages of reduction under the
reducedData dir (reducedACTTest for the above example). At the moment, this is not cleaned up.

The extracted 1d spectra can be found under (for the example above) reducedACTTest/P000807N01/1DSpec/.
They include the slit number in the file name (e.g., \_SLIT10.fits). You can check which slit
(identified from the flat frame) corresponds with which extracted spectrum by examining the
masterFlat\_0.fits file and the corresponding DS9 region file (ending .reg) in the reducedData/
directory.

The 1d spectra are in FITS table format (stored in the extension 1D\_SPECTRUM with the following 
columns:

SPEC 	- the object spectrum
SKYSPEC - the spectrum of the sky
LAMBDA	- the wavelength scale corresponding to both of the above (in Angstroms)
 






