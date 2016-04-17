"""Create reference spectrum model.

Takes as input 2d arc spectrum image.

If we don't find a plain text file of coords matching a row to wavelength in angstroms, 
we display a plot so the user can make one.

Then we fit and save a model.

"""

import os
import sys
import pyfits
import atpy
from astLib import *
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage
from scipy import interpolate
from scipy import optimize
import optparse
import pickle
import IPython
#plt.matplotlib.interactive(True)

#-------------------------------------------------------------------------------------------------------------
def parseApproxArcCoordsFile(fileName):
    """Returns a dictionary indexed by x coord of line wavelengths in Angstroms.
    
    """
    
    inFile=open(fileName, "r")
    lines=inFile.readlines()
    inFile.close()
    
    coordsDict={}
    for line in lines:
        if len(line) > 3 and line[0] != "#":
            bits=line.split()
            coordsDict[int(bits[1])]=float(bits[0])
    
    return coordsDict

#-------------------------------------------------------------------------------------------------------------
def detectLines(data, sigmaCut = 3.0, thresholdSigma = 5.0, featureMinPix = 30):
    """Detect lines in a 2d arc spectrum. Uses the central row of the 2d spectrum only.
    
    Returns: featureTable, segmentationMap
    
    """
    
    # Detect arc lines
    mean=0
    sigma=1e6
    for i in range(20):
        nonZeroMask=np.not_equal(data, 0)
        mask=np.less(abs(data-mean), sigmaCut*sigma)
        mask=np.logical_and(nonZeroMask, mask)
        mean=np.mean(data[mask])
        sigma=np.std(data[mask])
    detectionThreshold=thresholdSigma*sigma
    mask=np.greater(data-mean, detectionThreshold)

    # Get feature positions, number of pixels etc.
    # Find features in 2d, match to wavelength coord in centre row
    segmentationMap, numObjects=ndimage.label(mask)
    sigPixMask=np.equal(mask, 1)
    objIDs=np.unique(segmentationMap)
    objNumPix=ndimage.sum(sigPixMask, labels = segmentationMap, index = objIDs)
    objPositions_centreRow=ndimage.center_of_mass(data[data.shape[0]/2], labels = segmentationMap, index = objIDs)
    objPositions_centreRow=np.array(objPositions_centreRow).flatten()
    minPixMask=np.greater(objNumPix, featureMinPix)
    featureTable=atpy.Table()
    featureTable.add_column('id', objIDs[minPixMask])
    featureTable.add_column('x_centreRow', objPositions_centreRow[minPixMask])
    featureTable.add_column('y_centreRow', [data.shape[0]/2]*len(featureTable))
    featureTable.add_column('amplitude', data[data.shape[0]/2, np.array(np.round(featureTable['x_centreRow']), dtype = int)])

    return featureTable, segmentationMap

#-------------------------------------------------------------------------------------------------------------
def tagWavelengthFeatures(featureTable, approxCoordsDict):
    """This adds a wavelength column to featureTable, tagging features which are nearest in x-coord to the 
    contents of features in approxCoordsDict. This will only work on the reference arc spectrum.
    
    Removes features from featureTable which are not tagged with wavelengths.
    
    Returns featureTable
    
    """

    featureTable.add_column('wavelength', np.zeros(len(featureTable)))
    for x in approxCoordsDict.keys():
        wavelength=approxCoordsDict[x]
        rowNumber=np.argmin(abs(x-featureTable['x_centreRow']))
        featureTable['wavelength'][rowNumber]=wavelength

    featureTable=featureTable.where(featureTable['wavelength'] != 0)
    
    return featureTable

#-------------------------------------------------------------------------------------------------------------
def makeModelArcSpectrum(data, approxCoordsDict, outFileName, yRow, sigmaCut = 3.0, thresholdSigma = 5.0, 
                         featureMinPix = 30):
    """Make reference model arc spectrum. This has wavelengths of features identified in a table. We also
    save the middle row of the spectrum.
    
    """
    
    # Detect and tag features with known wavelengths in reference spectrum
    featureTable, segMap=detectLines(data)
    featureTable=tagWavelengthFeatures(featureTable, approxCoordsDict)
    data_centreRow=data[yRow]
    
    # Save reference model as a pickled dictionary
    refModelDict={'featureTable': featureTable, 'arc_centreRow': data_centreRow}
    pickleFile=file(outFileName, "wb")
    pickler=pickle.Pickler(pickleFile)
    pickler.dump(refModelDict)
    pickleFile.close()
        
#-------------------------------------------------------------------------------------------------------------
def findWavelengthCalibration(arcFileName, modelFileName, sigmaCut = 3.0, thresholdSigma = 5.0, 
                              featureMinPix = 50, order = 6):
    """Find wavelength calibration for .fits image arcFileName containing a 2d arc spectrum.
    
    modelFileName is the path to a model made by makeModelArcSpectrum.
    
    Returns an array of polynomial fit coefficients that can be fed into wavelengthCalibrateAndRectify
    
    NOTE: CHECK ORDER, ADD STEPPIX FOR LONGSLIT
    
    """
    
    # Load reference model
    pickleFile=file(modelFileName, "rb")
    unpickler=pickle.Unpickler(pickleFile)
    refModelDict=unpickler.load()
    pickleFile.close()
    
    # First need to find arc features
    arcImg=pyfits.open(arcFileName)
    arcData=arcImg['SCI'].data
    arcFeatureTable, arcSegMap=detectLines(arcData)
    
    # Use cross correlation to get initial guess at shift between arc and model
    arc_centreRow=arcData[arcData.shape[0]/2]
    shift=np.argmax(np.correlate(refModelDict['arc_centreRow'], arc_centreRow, mode = 'full'))-refModelDict['arc_centreRow'].shape[0]
    
    # Find wavelength dependent scale change by brute force - takes ~2 sec.
    # We tried using optimize.minimize but it wasn't robust and didn't work as well as this
    # NOTE: We may need to set min/max scales somehow
    data_x=np.arange(0, arc_centreRow.shape[0])
    scales=np.linspace(-0.05, 0.05, 2000)
    overlaps=[]
    for s in scales:
        tck=interpolate.splrep((data_x+shift)+s*data_x, arc_centreRow)
        arc_centreRow_shifted=interpolate.splev(data_x, tck, ext = 1)
        overlap=np.trapz(abs(refModelDict['arc_centreRow']-arc_centreRow_shifted))
        overlaps.append(overlap)
    bestFitScale=scales[np.argmin(overlaps)]
    #plt.plot(scales, overlaps, 'r.')
    
    # Sanity check: coord transformation arc -> ref model coords
    x=np.arange(0, len(arc_centreRow))
    arc_x_shifted=(x+shift)+x*bestFitScale
    plt.plot(x, refModelDict['arc_centreRow'], 'b-')
    plt.plot(arc_x_shifted, arc_centreRow, 'r-')
    plt.close()
    
    # Tag wavelength features in arc
    arcFeatureTable.add_column('wavelength', np.zeros(len(arcFeatureTable)))
    maxDistancePix=20   # maximum error allowed in our coord transformation
    for row in refModelDict['featureTable']:
        dist=abs(row['x_centreRow']-(arcFeatureTable['x_centreRow']+shift+arcFeatureTable['x_centreRow']*bestFitScale))
        if dist.min() < maxDistancePix:
            index=np.argmin(dist)
            arcFeatureTable['wavelength'][index]=row['wavelength']
    arcFeatureTable=arcFeatureTable.where(arcFeatureTable['wavelength'] != 0)
    
    # Sanity check: tagged features
    plt.plot(arcData[arcData.shape[0]/2], 'b-')
    for row in arcFeatureTable:
        plt.text(row['x_centreRow'], row['amplitude'], row['wavelength'])
    plt.close()
    #plt.show()
      
    # Find 2d wavelength solution which we can use for rectification/wavelength calibration
    # Fit functions for how feature x positions change with y
    # For longslit, we have a border with no data - here is a rough estimate of size in %-age terms
    border=int(round(arcData.shape[0]*0.05))
    stepPix=10
    ys=np.arange(border, arcData.shape[0]-border, stepPix)
    order=6
    for o in range(order+1):
        arcFeatureTable.add_column('poly%d' % (o), np.zeros(len(arcFeatureTable)))
    count=0
    #plt.matplotlib.interactive(True)
    for row in arcFeatureTable:
        count=count+1
        print "... %d/%d ..." % (count, len(arcFeatureTable))
        xs=np.zeros(ys.shape)
        for i in range(len(ys)):
            objPositions=ndimage.center_of_mass(arcData[ys[i]], labels = arcSegMap, index = arcFeatureTable['id'])
            xs[i]=objPositions[np.where(arcFeatureTable['id'] == row['id'])[0]][0]
        # For long slit, this is not linear it is also, obviously, double valued!
        validMask=np.equal(np.isnan(xs), False)
        coeffs=np.polyfit(ys[validMask], xs[validMask], 6)
        poly=np.poly1d(coeffs)
        xsFit=poly(ys)
        for o in range(order+1):
            row['poly%d' % (o)]=coeffs[o]
        #plt.plot(xs[validMask], ys[validMask], 'r.')
        #plt.plot(xsFit, ys, 'k--')
        #IPython.embed()
        #sys.exit()

    # Wavelength calibration and model with arbitrary order polynomials - get coeffs for each row
    # We could potentially fit these coeffs as fn. of y - they are all well behaved
    # This array should be all we need to wavelength calibrate + rectify
    # KEEP TRACK OF ORDER HERE! NEED TO GET ALL COEFFS, BUT SUBSQUENT FIT CAN BE 2ND ORDER
    wavelengths=arcFeatureTable['wavelength']
    fitCoeffsArr=[]
    for y in range(arcData.shape[0]):
        xs=[]
        for row in arcFeatureTable:
            coeffs=[]
            for o in range(order+1):
                coeffs.append(row['poly%d' % (o)])
            poly=np.poly1d(coeffs)
            xs.append(poly(y))
        xs=np.array(xs)
        fitCoeffsArr.append(np.polyfit(xs, wavelengths, 2))
    fitCoeffsArr=np.array(fitCoeffsArr)
    
    # Sanity check: wavelength calibration model with tagged features
    wavelengthCalibPoly=np.poly1d(fitCoeffsArr[arcData.shape[0]/2])
    wavelengths=wavelengthCalibPoly(np.arange(arcData.shape[1]))
    plt.plot(wavelengths, arcData[arcData.shape[0]/2], 'r-')
    plt.plot(arcFeatureTable['wavelength'], arcFeatureTable['amplitude'], 'bo')
    for row in arcFeatureTable:
        plt.text(row['wavelength'], row['amplitude'], row['wavelength'])
    plt.close()
    
    return fitCoeffsArr

#-------------------------------------------------------------------------------------------------------------
def wavelengthCalibrateAndRectify(inFileName, outFileName, fitCoeffsArr):
    """Applies the wavelength calibration, and rectification, to 2d spectrum in inFileName, writing 
    output to outFileName. The wavelength calibration is provided in fitCoeffsArr (see 
    findWavelengthCalibration)
    
    """
    
    img=pyfits.open(inFileName)
    data=img[0].data
    header=img[0].header
    
    # Using above, make an array containing wavelengths
    wavelengthsMap=np.zeros(data.shape)
    for y in range(data.shape[0]):
        wavelengthCalibPoly=np.poly1d(fitCoeffsArr[y])
        wavelengthsMap[y]=wavelengthCalibPoly(np.arange(data.shape[1]))
    #astImages.saveFITS("wavelengthsMap.fits", wavelengthsMap, None)
    
    # How we would want our wavelength map to look after applying some transformation
    # To make things easier later, make a linear wavelength scale
    wavelengths_centreRow=wavelengthsMap[wavelengthsMap.shape[0]/2]
    maxWavelength=wavelengths_centreRow.max()
    minWavelength=wavelengths_centreRow.min()
    linearWavelengthRange=np.linspace(minWavelength, maxWavelength, data.shape[1])
    FITSWavelengthScale=linearWavelengthRange[1]-linearWavelengthRange[0]
    FITSRefLambda=linearWavelengthRange[0]
    FITSRefPixel=1                              # Remember index from 1 is FITS convention
    rectWavelengthsMap=np.array([linearWavelengthRange]*data.shape[0])
    #astImages.saveFITS("rectWavelengthsMap.fits", rectWavelengthsMap, None)

    # Remap the data to our preferred linear wavelength scale
    # Assume we can treat each row independently
    # Save linear spectral WCS in header
    rectifiedData=np.zeros(data.shape)
    for y in range(data.shape[0]):
        tck=interpolate.splrep(wavelengthsMap[y], data[y])
        rectifiedData[y]=interpolate.splev(rectWavelengthsMap[y], tck)
    header['CTYPE1']='LINEAR'
    header['DISPAXIS']=1
    header['CRVAL1']=FITSRefLambda
    header['CRPIX1']=FITSRefPixel
    header['CD1_1']=FITSWavelengthScale
    header['CDELT1']=FITSWavelengthScale
    header['CUNIT1']='Angstroms'
    newImg=pyfits.HDUList()
    hdu=pyfits.PrimaryHDU(None, header)   
    hdu.data=rectifiedData
    newImg.append(hdu)
    newImg.writeto(outFileName, clobber = True)
    
    # Sanity check plot: linear wavelength scale
    plt.plot(rectWavelengthsMap[data.shape[0]/2], rectifiedData[data.shape[0]/2], 'r-')
    plt.close()
    
#-------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    usage="""
    %prog arcFrame2d.fits

    """

    parser = optparse.OptionParser(usage)
    parser.add_option("-y", "--y", dest="yRow", help="row from which to extract arc spectrum")
    options, args = parser.parse_args()

    if len(args) != 1:
        parser.error("What the...?!?! You fool! That's the wrong number of arguments!")    

    arcFileName=args[0]
    
    #if os.path.exists(modelsDir) == False:
        #os.makedirs(modelsDir)
    
    img=pyfits.open(arcFileName)
    data=img['SCI'].data
    header=img[0].header

    if options.yRow != None:
        row=int(options.yRow)
    else:
        row=data.shape[0]/2
            
    coordsFileName=header['GRATING']+"_"+header['LAMPID']+"_"+header['CCDSUM'].replace(" ", "x")+".txt"
    modelFileName="RefModel_"+coordsFileName.replace(".txt", ".pickle")
    
    if os.path.exists(coordsFileName) == False:
        print "Coordinates file not found - you need to make one..."
        print "1. Use the plot window to record approximate pixel coords (x-axis) for spectral lines identified in a calibrated arc plot (see http://pysalt.salt.ac.za/lineatlas/lineatlas.html and look under the heading 'Arc Lamp Plots')."
        print "2. Save a text file (name it %s) in the current directory, which has columns: Wavelengh (Angstroms, float), pixel coord (integer) on x-axis." % (coordsFileName)
        print "3. Close the plot window and re-run this script to generate the model."
        plt.plot(data[row])
        plt.title(coordsFileName)
        plt.show()
        sys.exit()
    
    approxCoordsDict=parseApproxArcCoordsFile(coordsFileName)

    # Makes the reference model
    makeModelArcSpectrum(data, approxCoordsDict, modelFileName, row)
    sys.exit()
    
    # Test wavelength calibration and rectification - this will make it into pipeline code once ready
    #testFileName="testArc1.fits"
    testFileName=arcFileName
    #testFileName="testArc3.fits"
    wavelengthCalibCoeffs=findWavelengthCalibration(testFileName, modelFileName)
    wavelengthCalibrateAndRectify(testFileName, "wavelengthCalib_"+testFileName, wavelengthCalibCoeffs)
