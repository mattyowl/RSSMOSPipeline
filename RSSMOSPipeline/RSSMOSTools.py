"""

    Copyright 2014-2022 Matt Hilton
    Copyright 2023 Matt Hilton, Melissa Morris

    This file is part of RSSMOSPipeline.

    RSSMOSPipeline is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    RSSMOSPipeline is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with RSSMOSPipeline.  If not, see <http://www.gnu.org/licenses/>.
    
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import glob
import time
import datetime
import astropy.table as atpy
import argparse
from scipy import interpolate
from scipy import ndimage
from scipy import optimize
from scipy import stats
import pickle
import logging
import RSSMOSPipeline
# import IPython
from astropy.table import Table
#plt.matplotlib.interactive(True)

#-------------------------------------------------------------------------------------------------------------
if sys.version_info.major == 3:
    PICKLE_OPTIONS={'encoding': 'latin1'}
else:
    PICKLE_OPTIONS={}

LOGFILE=None
REF_MODEL_DIR=RSSMOSPipeline.__path__[0]+os.path.sep+"data"+os.path.sep+"modelArcSpectra"

# For checking wavelength calibration accuracy
# Useful plot: http://www.astro.keele.ac.uk/jkt/GrSpInstructions/skylines1.jpg
checkSkyLines=[
               #4980.7,
               #5461.0,
               5577.0, 
               5893.0, 
               6300.304,
               6363.708,
               6863.955,
               #6923.220,
               #7276.405,
               #7316.282,
               #7340.885,
               #7750.640,
               #7794.112,
               #7821.503,
               #7913.708,
               #7993.332,
               #8344.602,
               #8399.170,
               #8430.174,
               #8761.314,
               #8767.912,
               #8778.333
               ]

logger=logging.getLogger('RSSMOSPipeline')
    
#-------------------------------------------------------------------------------------------------------------
def listToString(inList, prefix = "", delimiter=","):
    """Converts a list into a comma delimited string, can add prefixes if needed.
    
    """
    outString=""
    for i in inList:
        if outString != "":
            outString=outString+delimiter
        outString=outString+prefix+str(i)
    return outString

#-------------------------------------------------------------------------------------------------------------
def listToFile(inList, fileName, extension = None):
    """Converts a list of file names into a file that can be used in iraf with the @ syntax.
    
    Returns the file name of the list file, with @ prefix added already.
    
    """
    
    listFile=open(fileName, "w")
    for inFile in inList:
        if extension == None:
            listFile.write(inFile+"\n")
        else:
            listFile.write(inFile+extension+"\n")
    listFile.close()

    return "@"+fileName

#-------------------------------------------------------------------------------------------------------------
def makeOutputFileNameList(inFileNameList, prefix, outDir):
    """Generates output file names in the form outDir/<prefix>filename.fits for all file names in the 
    inFileNameList
    
    """
    
    outList=[]
    for inFileName in inFileNameList:
        outList.append(makeOutputFileName(inFileName, prefix, outDir))
    return outList
    
#-------------------------------------------------------------------------------------------------------------
def makeOutputFileName(inputFileName, prefix, outDir):
    """Given a raw file name in the form outDir/filename.fits, generate output filename in the form
    reduced/<prefix>filename.fits
    
    """
    
    return outDir+os.path.sep+prefix+os.path.split(inputFileName)[-1]

#-------------------------------------------------------------------------------------------------------------
def splitMEF(MEFFileName, rootOutFileName):
    """This splits a MEF.
    
    """
    
    img=pyfits.open(MEFFileName)
    for i in range(len(img)):
        if img[i].name == 'SCI':
            newImg=pyfits.HDUList()
            hdu=pyfits.PrimaryHDU(None, img[i].header)   
            hdu.data=img[i].data
            newImg.append(hdu)
            outFileName=rootOutFileName.replace(".fits", "_%d.fits" % (i))
            if os.path.exists(outFileName) == True:
                os.remove(outFileName)
            newImg.writeto(outFileName, overwrite = True)    

#-------------------------------------------------------------------------------------------------------------
def getImageInfo(rawDir, modelArcsDir = None):
    """Sorts through all .fits files, making lists of biases, flats, science frames, arc frames, for each 
    night of observations, grating, position combo. Returns a dictionary of lists.
    
    NOTE: Explicitly avoiding standard stars for now 
    
    """
 
    pickleFileName=rawDir+os.path.sep+"imageInfo.pkl"
    logger.info("Reading image headers (cache file location: %s)" % (pickleFileName))    
    previousFiles=[]
    infoDict={}
    if os.path.exists(pickleFileName) == True and modelArcsDir is None:
        with open(pickleFileName, "rb") as pickleFile:
            unpickler=pickle.Unpickler(pickleFile)
            infoDict=unpickler.load()
            previousFiles=unpickler.load()
        
    # Now organised by objectName and meta info -> flats, arcs, object frames
    # This is so we can use one pipeline for MOS and longslit
    files=glob.glob(rawDir+os.path.sep+"mbxgp*.fits")
    newFiles=[]
    for f in files:
        if f not in previousFiles:
            newFiles.append(f)
    
    # First, get object names, and assemble global list of arcs and flats
    arcsList=[]
    flatsList=[]
    for f in newFiles:
        img=pyfits.open(f)
        header=img[0].header
        if header['OBSMODE'] == 'SPECTROSCOPY':
            obsType=header['CCDTYPE']
            maskID=header['MASKID']
            if obsType == 'OBJECT':
                maskName=header['OBJECT'].replace('"', "")+"_"+maskID
                infoDict[maskName]={}
                infoDict[maskName][maskID]={}
                infoDict[maskName]['maskID']=maskID             # Clunky, but convenient
                infoDict[maskName]['maskType']=header['MASKTYP']  
                infoDict[maskName]['objName']=header['OBJECT'].replace("'", "").replace('"', "")  
            elif obsType == 'FLAT':
                flatsList.append(f)
            elif obsType == 'ARC':
                arcsList.append(f)
    
    # Now add flats etc.
    for maskName in infoDict.keys():
        
        for f in newFiles:
            img=pyfits.open(f)
            header=img[0].header
            if header['OBSMODE'] == 'SPECTROSCOPY':# and header['MASKTYP'] == 'MOS':
    
                dateObs=header['DATE-OBS']
                timeObs=header['TIME-OBS']
                maskID=header['MASKID']
                obsType=header['CCDTYPE']
                objName=header['OBJECT'].replace("'", "").replace('"', "")  

                if maskID == infoDict[maskName]['maskID']:
                                    
                    if obsType not in infoDict[maskName][maskID].keys():
                        infoDict[maskName][maskID][obsType]=[]
                    if obsType == 'OBJECT' and objName == infoDict[maskName]['objName']:
                        infoDict[maskName][maskID][obsType].append(f)
                        # Just so we can track this later in output 1d spectra
                        infoDict[maskName][maskID]['RA']=header['RA']
                        infoDict[maskName][maskID]['DEC']=header['DEC']
                        # Add matching arcs and flats
                        for o in ['FLAT', 'ARC', 'modelFileNames', 'modelExists']:
                            if o not in infoDict[maskName][maskID].keys():
                                infoDict[maskName][maskID][o]=[]
                        possFlats=findMatchingFilesByTime(f, flatsList)
                        for p in possFlats:
                            matchesSettings=False
                            with pyfits.open(p) as pImg:
                                if pImg[0].header['GRATING'] == header['GRATING'] and pImg[0].header['CCDSUM'] == header['CCDSUM']:
                                    matchesSettings=True
                            if p not in infoDict[maskName][maskID]['FLAT'] and matchesSettings == True:
                                infoDict[maskName][maskID]['FLAT'].append(p)
                        possArcs=findMatchingFilesByTime(f, arcsList)
                        for p in possArcs:
                            if p not in infoDict[maskName][maskID]['ARC']:
                                matchesSettings=False
                                with pyfits.open(p) as pImg:
                                    if pImg[0].header['GRATING'] == header['GRATING'] and pImg[0].header['CCDSUM'] == header['CCDSUM']:
                                        matchesSettings=True
                                    binning=pImg[0].header['CCDSUM'].replace(" ", "x")
                                    grating=pImg[0].header['GRATING']
                                    lampid=pImg[0].header['LAMPID']
                                    modelFileNamesGlob=REF_MODEL_DIR+os.path.sep+"RefModel_"+grating+"_"+lampid+"_*.pickle"
                                modelFileNames=glob.glob(modelFileNamesGlob)
                                if modelArcsDir is not None:
                                    modelFileNames=modelFileNames+glob.glob(modelArcsDir+os.path.sep+"RefModel_"+grating+"_"+lampid+"_*.pickle")
                                for modelFileName in modelFileNames:
                                    if modelFileName not in infoDict[maskName][maskID]['modelFileNames'] and matchesSettings == True and lampid != "NONE":
                                        infoDict[maskName][maskID]['modelFileNames'].append(modelFileName)
                                        if os.path.exists(modelFileName) == True:
                                            infoDict[maskName][maskID]['modelExists'].append(True)
                                        else:
                                            infoDict[maskName][maskID]['modelExists'].append(False)                                    
                                        infoDict[maskName][maskID]['ARC'].append(p)

    # For subsequent runs
    with open(pickleFileName, "wb") as pickleFile:
        previousFiles=files
        pickler=pickle.Pickler(pickleFile)
        pickler.dump(infoDict)
        pickler.dump(previousFiles)
            
    return infoDict

#-------------------------------------------------------------------------------------------------------------
def makeMasterFlats(maskDict, outDir, deltaHours = 0.5):
    """Make master flats from files in 'FLAT' key of maskDict. Automatically group flats taken within
    deltaHours. Adds paths to 'masterFlat' key.
    
    """
    
    flatLists=groupFilesListByTime(maskDict['FLAT'])
        
    maskDict['masterFlats']=[]

    logger.info("Making masterFlats (it is a good idea to check alignment with object spectra at the cutting stage, and remove any flats which aren't aligned)")
    
    for i in range(len(flatLists)):
        flatFiles=flatLists[i]
        masterFlatPath=outDir+os.path.sep+"masterFlat_%d.fits" % (i)
        logger.info("making %s (%s)" % (masterFlatPath, flatFiles))
        if os.path.exists(masterFlatPath) == False:
            flatCube=[]
            for f in flatFiles:
                img=pyfits.open(f)
                flatCube.append(img['SCI'].data)
            flatCube=np.array(flatCube)
            flatData=np.median(flatCube, axis = 0)
            img['SCI'].data=flatData
            if os.path.exists(masterFlatPath) == True:
                os.remove(masterFlatPath)
            img.writeto(masterFlatPath, overwrite = True)
        maskDict['masterFlats'].append(masterFlatPath)
        
#-------------------------------------------------------------------------------------------------------------
def getCTimeFromHeader(fileName):
    """Get unix ctime from header DATE-OBS and TIME-OBS keywords.
    
    Returns ctime (accurate to second level only)
    
    """
    
    img=pyfits.open(fileName)
    header=img[0].header
    dateObs=header['DATE-OBS']
    timeObs=header['TIME-OBS']
    y, m, d=dateObs.split("-")
    y=int(y)
    m=int(m)
    d=int(d)
    h, mn, s=timeObs.split(":")
    h=int(h)
    mn=int(mn)
    sec=int(s.split(".")[0])
    usec=int(s.split(".")[-1])*1000
    dateTimeObs=datetime.datetime(y, m, d, h, mn, sec, usec)
    ctime=int(dateTimeObs.strftime("%s"))
    
    return ctime

#-------------------------------------------------------------------------------------------------------------
def groupFilesListByTime(filesList, deltaHours = 0.5):
    """Given a list of files, splits it such that they are in groups covered by deltaHours.
    
    """
    
    ctimes=[]
    for f in filesList:
        ctime=getCTimeFromHeader(f)
        ctimes.append(ctime)  # 1 second accuracy is enough for us
        
    outFileLists=[]
    outCTimeLists=[]
    for fileName, ctime in zip(filesList, ctimes):
        foundList=False
        for cl, fl in zip(outCTimeLists, outFileLists):
            for c, f in zip(cl, fl):
                if abs(c-ctime) < 3600.0*deltaHours:
                    if ctime not in cl and fileName not in fl:
                        cl.append(ctime)
                        fl.append(fileName)
                    foundList=True
        if foundList == False:
            outFileLists.append([])
            outCTimeLists.append([])
            outFileLists[-1].append(fileName)
            outCTimeLists[-1].append(ctime)
            
    return outFileLists
#-------------------------------------------------------------------------------------------------------------
def writeDS9SlitRegions(regFileName, slitsDict, imageFileName):
    """Write a DS9 .reg file showing location of the slits. The image corresponding to imageFileName is only
    used to set the horizontal coordinates of the slit locations (i.e., so they are centred).
    
    """
    
    # Write out a .reg file so we can match slits to objects
    img=pyfits.open(imageFileName)
    centreColumn=int(img['SCI'].header['NAXIS1']/2)
    img.close()
    outFile=open(regFileName, "w")
    outFile.write("# DS9 region file\n")
    outFile.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    outFile.write("image\n")
    for key in slitsDict.keys():
        outFile.write("point(%.3f,%.3f) # point=boxcircle text={SLIT%d}\n" % (centreColumn, (slitsDict[key]['yMax']+slitsDict[key]['yMin'])/2.0, key))            
        outFile.write("box(%.1f,%.1f,%.1f,%.1f)\n" % (centreColumn, (slitsDict[key]['yMax']+slitsDict[key]['yMin'])/2.0+1, centreColumn*2.0, (slitsDict[key]['yMax']-slitsDict[key]['yMin'])))
    outFile.close()        
    
#-------------------------------------------------------------------------------------------------------------
def cutIntoSlitLets(maskDict, outDir, threshold = 0.1, slitFileName = None, noFlat = False):
    """Cuts files into slitlets, making MEF files. 
            
    threshold is the parameter used by findSlits
    
    """

    if len(maskDict['masterFlats']) == 0:
        raise Exception("No flat field through the slit mask is provided, so slit locations cannot be found - you will need to re-run using the -n and -F switches (use rss_mos_reducer -h to get help on the different options).")

    maskDict['slitsDicts']={}

    # If specified, creates slits based on input locations. Also skips flat fielding.
    if slitFileName is not None and noFlat:
        masterFlatPath = 'noflat'
        slitsDict=slitsFromFile(slitFileName)
        maskDict['slitsDicts'][masterFlatPath]=slitsDict

    # If specified, creates slits based on input locations. Does not skip flat fielding.
    elif slitFileName is not None and noFlat != True:
        for i in range(len(maskDict['masterFlats'])):
            masterFlatPath=maskDict['masterFlats'][i]
            cutMasterFlatPath=masterFlatPath.replace("masterFlat", "cmasterFlat")
            slitsDict=slitsFromFile(slitFileName)
            maskDict['slitsDicts'][masterFlatPath]=slitsDict

    # We find slits using master flats, and store slit locations in a dictionary indexed by masterFlatPath
    # We have a routine to pull out a corresponding master flat (and hence slitsDict) for every image.
    else:
        for i in range(len(maskDict['masterFlats'])):
            masterFlatPath=maskDict['masterFlats'][i]
            cutMasterFlatPath=masterFlatPath.replace("masterFlat", "cmasterFlat")
            slitsDict=findSlits(masterFlatPath, threshold = threshold)
            maskDict['slitsDicts'][masterFlatPath]=slitsDict

        # To avoid problems with occasional missing slits (if we treat each flat separately), use
        # the slits found from the first flat as a reference, and find the y-shifts between them
        refDict=maskDict['slitsDicts'][maskDict['masterFlats'][0]]
        img=pyfits.open(maskDict['masterFlats'][0])
        height=img[1].data.shape[0]
        ref=np.zeros(height)
        for key in refDict.keys():
            ref[int(round(refDict[key]['yCentre']))]=1
        shiftsDict={}
        for key in maskDict['slitsDicts'].keys():
            if key != maskDict['masterFlats'][0]:
                slitsDict=maskDict['slitsDicts'][key]
                g=np.zeros(height)
                for skey in slitsDict.keys():
                    g[int(round(slitsDict[skey]['yCentre']))]=1
                corr, corrMax, shift=fftCorrelate(ref, g)
                shiftsDict[key]=shift
            else:
                shiftsDict[key]=0.

        # Remake all slits dictionaries, based on reference minus shift
        newSlitsDicts={}
        for i in range(len(maskDict['masterFlats'])):
            objPath=maskDict['masterFlats'][i]
            newSlitsDicts[objPath]={}
            for key in refDict:
                newSlitsDicts[objPath][key]={}
                newSlitsDicts[objPath][key]['yMin']=refDict[key]['yMin']-int(round(shiftsDict[objPath]))
                newSlitsDicts[objPath][key]['yMax']=refDict[key]['yMax']-int(round(shiftsDict[objPath]))
        maskDict['slitsDicts']=newSlitsDicts
        
    # ^^^ Tidy all the above up later
        
    # Cut arcs, matched here with appropriate object frames
    toCutList=maskDict['OBJECT']#+maskDict['ARC']
    outCutList=makeOutputFileNameList(toCutList, "c", outDir)
    logger.info("It is a good idea to check that for the below the corresponding DS9 .reg file aligns with the slits and that object spectra are actually centred in the slits")
    maskDict['cutFlatDict']={}
    maskDict['cutArcDict']={}
    for f, outFileName in zip(toCutList, outCutList):
        if noFlat == False:
            flatFileName=findMatchingFilesByTime(f, maskDict['masterFlats'], timeInterval = None)[0]
            slitsDict=maskDict['slitsDicts'][flatFileName]
            logger.info("cutting %s (and arcs, flats) using %s for slit definition ..." % (f, flatFileName))
            label=os.path.split(flatFileName)[-1].replace(".fits", "")
            # Flat
            cutFlatFileName=makeOutputFileName(flatFileName, "c"+label, outDir)
            cutSlits(flatFileName, cutFlatFileName, slitsDict)
            maskDict['cutFlatDict'][f]=cutFlatFileName
        else:
            label='masterFlat_0'
            flatFileName = 'masterFlat_0'
        # Object
        cutSlits(f, outFileName, slitsDict)
        # Arc
        arcFileName=findMatchingFilesByTime(f, maskDict['ARC'], timeInterval = None)[0]
        cutArcFileName=makeOutputFileName(arcFileName, "c"+label, outDir)
        cutSlits(arcFileName, cutArcFileName, slitsDict)
        maskDict['cutArcDict'][f]=cutArcFileName
        # DS9 regions
        regFileName=flatFileName.replace(".fits", "_slitLocations.reg")
        writeDS9SlitRegions(regFileName, slitsDict, f)
        
#-------------------------------------------------------------------------------------------------------------
def cutIntoPseudoSlitLets(maskDict, outDir, thresholdSigma = 3.0, noFlat = False):
    """For longslit data. Finds objects, and then cuts into pseudo-slitlets: we take some region +/- Y pixels
    around the object trace and pretend that is a MOS slitlet. Outputs MEF files.
            
    """
        
    # Find object traces in OBJECT frames and cut +/- some distance in Y around them
    maskDict['slitsDicts']={}
    for i in range(len(maskDict['OBJECT'])):
        objPath=maskDict['OBJECT'][i]
        slitsDict=findPseudoSlits(objPath, thresholdSigma = thresholdSigma)
        maskDict['slitsDicts'][objPath]=slitsDict
    
    # There can be significant offsets between object traces in longslit frames
    # So... compare all the pseudo-slits we assigned and measure a y-offset from the i=0 dict
    # This will be used by cutSlits, if present
    refDict=maskDict['slitsDicts'][maskDict['OBJECT'][0]]
    img=pyfits.open(maskDict['OBJECT'][0])
    height=img[1].data.shape[0]
    ref=np.zeros(height)
    for key in refDict.keys():
        ref[int(round(refDict[key]['yCentre']))]=1
    shiftsDict={}
    for key in maskDict['slitsDicts'].keys():
        if key != maskDict['OBJECT'][0]:
            slitsDict=maskDict['slitsDicts'][key]
            g=np.zeros(height)
            for skey in slitsDict.keys():
                g[int(round(slitsDict[skey]['yCentre']))]=1
            corr, corrMax, shift=fftCorrelate(ref, g)
            shiftsDict[key]=shift
        else:
            shiftsDict[key]=0.
    
    # Remake all pseudo-slits dictionaries, based on reference minus shift
    newSlitsDicts={}
    for i in range(len(maskDict['OBJECT'])):
        objPath=maskDict['OBJECT'][i]
        newSlitsDicts[objPath]={}
        for key in refDict:
            newSlitsDicts[objPath][key]={}
            newSlitsDicts[objPath][key]['yMin']=refDict[key]['yMin']-int(round(shiftsDict[objPath]))
            newSlitsDicts[objPath][key]['yMax']=refDict[key]['yMax']-int(round(shiftsDict[objPath]))
    maskDict['slitsDicts']=newSlitsDicts
    
    # ^^^ Tidy all the above up later
    
    # Cut arcs, flats, matched here with appropriate object frames
    toCutList=maskDict['OBJECT']#+maskDict['ARC']
    outCutList=makeOutputFileNameList(toCutList, "c", outDir)
    logger.info("It is a good idea to check from the corresponding DS9 .reg file that object spectra are actually centred in the pseudo-slits")
    maskDict['cutFlatDict']={}
    maskDict['cutArcDict']={}
    for f, outFileName in zip(toCutList, outCutList):
        if noFlat == False:
            logger.info("cutting %s (and arcs, flats) using %s for slit definition ..." % (f, f))
            slitsDict=maskDict['slitsDicts'][f]
            label=os.path.split(f)[-1].replace(".fits", "")
            flatFileName=findMatchingFilesByTime(f, maskDict['masterFlats'], timeInterval = None)[0]
            cutFlatFileName=makeOutputFileName(flatFileName, "c"+label, outDir)
            cutSlits(flatFileName, cutFlatFileName, slitsDict)
            maskDict['cutFlatDict'][f]=cutFlatFileName
        else:
            label='masterFlat_0'
            flatFileName='masterFlat_0'
        # Object
        cutSlits(f, outFileName, slitsDict)
        # Arc
        arcFileName=findMatchingFilesByTime(f, maskDict['ARC'], timeInterval = None)[0]
        cutArcFileName=makeOutputFileName(arcFileName, "c"+label, outDir)
        cutSlits(arcFileName, cutArcFileName, slitsDict)
        maskDict['cutArcDict'][f]=cutArcFileName
    
        # Write out a .reg file so we can match slits to objects
        img=pyfits.open(f)
        centreColumn=int(img['SCI'].header['NAXIS1']/2)
        img.close()
        regFileName=outFileName.replace(".fits", "_slitLocations.reg")
        outFile=open(regFileName, "w")
        outFile.write("# DS9 region file\n")
        outFile.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
        outFile.write("image\n")
        for key in slitsDict.keys():
            outFile.write("point(%.3f,%.3f) # point=boxcircle text={SLIT%d}\n" % (centreColumn, (slitsDict[key]['yMax']+slitsDict[key]['yMin'])/2.0, key))            
            outFile.write("box(%.1f,%.1f,%.1f,%.1f)\n" % (centreColumn, (slitsDict[key]['yMax']+slitsDict[key]['yMin'])/2.0+1, centreColumn*2.0, (slitsDict[key]['yMax']-slitsDict[key]['yMin'])))
        outFile.close() 
        
#-------------------------------------------------------------------------------------------------------------
def cutSlits(inFileName, outFileName, slitsDict):
    """Makes a MEF file containing slitlets.
    
    """
            
    img=pyfits.open(inFileName)
    d=img['SCI'].data
            
    newImg=pyfits.HDUList()
    hdu=pyfits.PrimaryHDU(None, img[0].header)
    newImg.append(hdu)
    for slitKey in slitsDict.keys():
        slitData=d[slitsDict[slitKey]['yMin']:slitsDict[slitKey]['yMax']]
        hdu=pyfits.ImageHDU(data = slitData, header = None, name = 'SLIT%d' % (slitKey))
        newImg.append(hdu)
        
    if os.path.exists(outFileName):
        os.remove(outFileName)
        
    newImg.writeto(outFileName)
    newImg.close()

#-------------------------------------------------------------------------------------------------------------
def slitsFromFile(slitFileName):
    """
    Creates dictionary that defines slits, using definitions from a file.
    File has 3 columns: slitno, ystart, and yend
    """
    dat = Table.read(slitFileName,format='ascii')
    slitsDict = {}
    for row in dat:
        slitsDict[row['slitno']] = {'yMin':row['ystart'],'yMax':row['yend'],'yCentre': (row['ystart']+row['yend'])/2.}
    return slitsDict
    
#-------------------------------------------------------------------------------------------------------------
def findSlits(flatFileName, minSlitHeight = 5, threshold = 0.1):
    """Find the slits, without using any info from the mask design file...
    
    minSlitHeight is used to throw out any problematic weird too-narrow slits (if any)
    
    threshold is used to identify slit edges. Values in the range 0.1 - 0.4 work best.
    
    Returns a dictionary which can be fed into cutSlits
    
    """ 
    
    img=pyfits.open(flatFileName)
    d=img['SCI'].data
    
    # Take out spectrum of flat lamp (approx)
    a=np.median(d, axis = 0)
    a=np.array([a]*d.shape[0])
    zeroMask=np.equal(a, 0)
    nonZeroMask=np.not_equal(a, 0)
    d[nonZeroMask]=d[nonZeroMask]/a[nonZeroMask]
    d[zeroMask]=0.
    #d[np.isnan(d)]=0.0

    # Use grad to find edges
    prof=np.median(d, axis = 1)
    grad=np.gradient(prof)
    plusMask=np.greater(grad, threshold)
    minusMask=np.less(grad, threshold*-1)

    # This looks for alternating +/-, but will merge slits which butt up against each other
    slitsDict={}
    lookingFor=1
    yMin=None
    yMax=None
    slitCount=0
    for i in range(len(plusMask)):
        if lookingFor == 1:
            if plusMask[i] == True:
                yMin=i
                lookingFor=0
        if lookingFor == 0:
            if minusMask[i] == True:
                yMax=i
                lookingFor=1
        if yMin != None and yMax != None and (yMax - yMin) > minSlitHeight:
            slitCount=slitCount+1
            slitsDict[slitCount]={'yMin': yMin, 'yMax': yMax, 'yCentre': (yMax+yMin)/2.}    
            yMin=None
            yMax=None
    
    ## Debugging
    #print("Check slitsDict")
    #IPython.embed()
    #sys.exit()
    
    # Slits can be bendy: measure the bendiness 
    # above routine misses large chunks of red end of bendy slits at top of mask
    # if we can just cut out a larger slitlet, the rectification / wavelength calibration can handle unbending
    # Note: below doesn't work...
    #xSlitsDict={}
    
    #xBinStep=50
    #xBinEdges=[]
    #for i in range(d.shape[1]/xBinStep):
        #xBinEdges.append(i*xBinStep)
    
    #for key in slitsDict.keys():
        #xSlitsDict[key]={}
        #xSlitsDict[key]['x']=[]
        #xSlitsDict[key]['yMin']=[]
        #xSlitsDict[key]['yMax']=[]
        #for i in range(len(xBinEdges)-1):
            
            #xMin=xBinEdges[i]
            #xMax=xBinEdges[i+1]
            #x=(xMin+xMax)/2
            
            ## Use grad to find edges, skip chip gaps
            #prof=np.median(d[:, xMin:xMax], axis = 1)
            #if np.sum(prof) > 0:
                #xSlitsDict[key]['x'].append(x)
                #grad=np.gradient(prof)
                ##y=(slitsDict[key]['yMax']+slitsDict[key]['yMin'])/2
                #diffMin=abs(np.where(grad > 0.1)[0]-slitsDict[key]['yMin'])
                #yMin=np.where(grad > 0.1)[0][np.where(diffMin == diffMin.min())[0][0]]
                #diffMax=abs(np.where(grad < -0.1)[0]-slitsDict[key]['yMax'])
                #yMax=np.where(grad < -0.1)[0][np.where(diffMax == diffMax.min())[0][0]]            
                #xSlitsDict[key]['yMin'].append(yMin)
                #xSlitsDict[key]['yMax'].append(yMax)
    
    ## For now, just cut within max(yMax), min(yMin) and rely on rectify step to unbend
    #for key in xSlitsDict:
        #xSlitsDict[key]['yMin']=min(xSlitsDict[key]['yMin'])
        #xSlitsDict[key]['yMax']=max(xSlitsDict[key]['yMax'])
        #del xSlitsDict[key]['x']
    #slitsDict=xSlitsDict
    
    #---
    # Below works on A3827 - delete once tested with new algorithm
    #img=pyfits.open(flatFileName)
    #d=img['SCI'].data
    #prof=np.median(d, axis = 1)
    #threshold=np.median(prof)*2
        
    #sigPix=np.array(np.greater(prof, threshold), dtype=int)
    #sigPixMask=np.equal(sigPix, 1)
    #segmentationMap, numObjects=ndimage.label(sigPix)
    #maxSlits=segmentationMap.max()

    #slitsDict={}
    #for i in range(1, maxSlits+1):
        #indices=np.where(segmentationMap == i)[0]
        #slitsDict[i]={'yMin': indices.min(), 'yMax': indices.max()}

    ## Sanity check plot
    #plt.plot(prof)
    #plt.plot([np.median(prof)*2]*len(prof), 'r-')
    #IPython.embed()
    #sys.exit()

    return slitsDict

#-------------------------------------------------------------------------------------------------------------
def findPseudoSlits(objFileName, skyRows = 20, minSlitHeight = 10., thresholdSigma = 3., minTraceWidth = 5,\
                    halfBlkSize = 100, filterPix = 50, maskEdgePix = 250):
    """Finds object traces in longslit data, defines regions +/- skyRows around them, so we can treat in 
    the same way as MOS slitlets.
    
    objects are detected as peaks in the SNR profile across the slit. Use minTraceWidth to set the number
    of pixels in the SNR profile that must be above thresholdSigma for an object to be detected.
    
    halfBlkSize is used for measuring local background. This value is given in ubinned pixels and will be
    adjusted accordingly according to binning in vertical direction along the detector CCDs.
    
    filterPix is used for filtering out large scale gradient in the profile in the spatial direction.
    Given in unbinned pixels.

    maskEdgePix is used to throw out the edges (top and bottom in spatial direction) which are useless.
    Given in unbinned pixels.

    Returns a dictionary which can be fed into cutSlits
    
    """ 
    
    with pyfits.open(objFileName) as img:
        d=img['SCI'].data

    halfBlkSize=int(halfBlkSize/img['SCI'].header['CDELT2']) # So input is in unbinned pixels
    filterPix=int(filterPix/img['SCI'].header['CDELT2'])
    maskEdgePix=int(maskEdgePix/img['SCI'].header['CDELT2'])

    # Take out spectrum of flat lamp (approx)
    a=np.median(d, axis = 0)
    a=np.array([a]*d.shape[0])
    zeroMask=np.equal(a, 0)
    nonZeroMask=np.not_equal(a, 0)
    d[nonZeroMask]=d[nonZeroMask]/a[nonZeroMask]
    d[zeroMask]=0.
    
    # Profile may not be flattened anyway due to diff illumination - filter it to flatten gradient
    prof=np.median(d, axis = 1)
    unfiltered=prof
    prof[np.less(prof, 0)]=0.
    gBig=ndimage.gaussian_filter1d(prof, filterPix)
    prof=prof-gBig
    # prof[np.less(prof, 0)]=0.
    prof[:maskEdgePix]=0
    prof[-maskEdgePix:]=0

    # Old ----
    # # Find local background, noise (running clipped mean)
    # prof=np.median(d, axis = 1)
    # prof[np.less(prof, 0)]=0.
    # sigmaCut=3.0
    # bck=np.zeros(prof.shape)
    # sig=np.zeros(prof.shape)
    # for y in range(prof.shape[0]):
    #     yMin=y-halfBlkSize
    #     yMax=y+halfBlkSize
    #     if yMin < 0:
    #         yMin=0
    #     if yMax > prof.shape[0]-1:
    #         yMax=prof.shape[0]-1
    #     mean=0
    #     sigma=1e6
    #     for i in range(20):
    #         nonZeroMask=np.not_equal(prof[yMin:yMax], 0)
    #         mask=np.less(abs(prof[yMin:yMax]-mean), sigmaCut*sigma)
    #         mean=np.mean(prof[yMin:yMax][mask])
    #         sigma=np.std(prof[yMin:yMax][mask]-mean)
    #     bck[y]=mean
    #     sig[y]=sigma
    # Old ^^^^

    # Estimate noise in filtered profile (ignoring that pixels are correlated etc.)
    sigmaCut=3.0
    sigma=1e6
    mean=0
    for i in range(20):
        mask=np.less(abs(prof-mean), sigmaCut*sigma)
        mean=np.mean(prof[mask])
        sigma=np.std(prof[mask])
    bck=0
    sig=sigma
   
    # Detect peaks
    profSNR=(prof-bck)/sig    
    profSNR[:maskEdgePix]=0
    profSNR[-maskEdgePix:]=0
    mask=np.greater(profSNR, thresholdSigma)    
    segmentationMap, numObjects=ndimage.label(mask)
    sigPixMask=np.equal(mask, 1)
    objIDs=np.unique(segmentationMap)
    objNumPix=ndimage.sum(sigPixMask, labels = segmentationMap, index = objIDs)
    objPositions=ndimage.center_of_mass(prof, labels = segmentationMap, index = objIDs)
    objPositions=np.array(objPositions).flatten()
    minPixMask=np.greater(objNumPix, minTraceWidth)
    
    # Define pseudo slits, including sky rows next to each object
    # We need yCentre here for using cross correlation to get the shift between images
    slitsDict={}
    slitCount=0
    for objID, yPos, traceWidth in zip(objIDs[minPixMask], objPositions[minPixMask], objNumPix[minPixMask]):
        yMin=int(round(yPos-(traceWidth/2.+skyRows)))
        yMax=int(round(yPos+(traceWidth/2.+skyRows)))
        if (yMax - yMin) > minSlitHeight:
            slitCount=slitCount+1
            slitsDict[slitCount]={'yMin': yMin, 'yMax': yMax, 'yCentre': yPos}    

    if slitsDict == {}:
        raise Exception("didn't find any object traces - try adjusting --longslit-threshold")
    
    return slitsDict

#-------------------------------------------------------------------------------------------------------------
def makeChipGapMask(data, numGaps = 2, thresholdValue = 1e-3):
    """Given an image array (data), find and mask the chip gaps.
    
    Returns a mask where 1 = chip gap, 0 otherwise.
    
    """
    
    # This has to work with any image - flats or background subtracted object spectra
    lowMaskValue=thresholdValue
    med=np.median(data, axis = 0)
    chipGapMask=np.array(np.less(abs(med), lowMaskValue), dtype = float)  # flags chip gaps as noise
    segmentationMap, numObjects=ndimage.label(chipGapMask)
    sigPixMask=np.equal(chipGapMask, 1)
    objIDs=np.unique(segmentationMap)
    if len(objIDs) > 0:
        objNumPix=ndimage.sum(sigPixMask, labels = segmentationMap, index = objIDs)
        newChipGapMask=np.zeros(med.shape)
        for i in range(numGaps):
            newChipGapMask[segmentationMap == objIDs[np.argmax(objNumPix)]]=1
            objNumPix[np.argmax(objNumPix)]=0
    chipGapMask=np.array([newChipGapMask]*data.shape[0])

    # Old
    # lowMaskValue=2.0
    # minPix=1000
    # chipGapMask=np.array(np.less(data, lowMaskValue), dtype = float)  # flags chip gaps as noise
    # segmentationMap, numObjects=ndimage.label(chipGapMask)
    # sigPixMask=np.equal(chipGapMask, 1)
    # objIDs=np.unique(segmentationMap)
    # if len(objIDs) > 0:
    #     objNumPix=ndimage.sum(sigPixMask, labels = segmentationMap, index = objIDs)
    #     for objID, nPix in zip(objIDs, objNumPix):
    #         if nPix < minPix:
    #             chipGapMask[np.equal(segmentationMap, objID)]=0.0
    
    return chipGapMask

#-------------------------------------------------------------------------------------------------------------
def applyFlatField(maskDict, outDir):
    """Applies the flat field correction. Let's do this in place...
    
    """

    logger.info("Applying flat field")
    
    toFlatList=makeOutputFileNameList(maskDict['OBJECT'], "c", outDir)
    for rawFileName, f in zip(maskDict['OBJECT'], toFlatList):
        cutMasterFlatPath=maskDict['cutFlatDict'][rawFileName]
        img=pyfits.open(f)
        flatImg=pyfits.open(cutMasterFlatPath)
        extensionsList=[]
        for hdu in img:
            if "SLIT" in hdu.name:
                extensionsList.append(hdu.name)
        
        for extension in extensionsList:
            data=img[extension].data
            flatfield=flatImg[extension].data
            med=np.median(flatfield, axis = 0)
                        
            # Find the chip gaps and make a mask
            chipGapMask=makeChipGapMask(flatfield, thresholdValue = 2)

            # Invert it as we want weights for interpolator
            gapsMask=np.equal(np.median(chipGapMask, axis = 0), False)
            
            # Fit and remove spectrum of flat lamp
            # Polynomial fit
            #x=np.arange(len(med))
            #poly=np.poly1d(np.polyfit(x[gapsMask], med[gapsMask], 40)) # was 10
            #mod=np.array([poly(x)]*data.shape[0])
            #flatfield=flatfield/mod
            #plt.matplotlib.interactive(True)
            #plt.plot(x, med)
            #plt.plot(x, poly(x))

            # Bit of extra masking - sudden changes in gradient
            # (could be around chip gaps, or at edges)
            # grad=np.gradient(med*gapsMask)
            # gapsMask=gapsMask*(abs(grad) < 500)

            # Extrapolate (rather than set to zero) on left/right edges
            # (keeps spline well behaved near edges)
            testSum=0
            for i in range(1, len(med)):
                testSum=testSum+med[-i]
                if testSum > 0:
                    med[-i:]=med[-i]
                    break
            
            # Spline fit - fall back to smaller number of knots if some problem
            # If this fails, it should be caught in the rss_mos_reducer script
            # Don't want too many knots - trying to bridge the chip gaps
            knotsToTry=[80, 40, 20]
            for numKnots in knotsToTry:
                x=np.arange(len(med))
                spacing=x.max()/numKnots
                w=gapsMask # weights are normally 1/stdev
                t=np.linspace(spacing, x.max()-spacing, numKnots-2)
                tck=interpolate.splrep(x, med, w = w, t = t, s = 0)
                if np.any(np.isnan(tck[1])) == True:
                    logger.warning('spline fit to flatfield contains NaNs - trying lower number of knots')
                    continue
                else:
                    break
            if np.any(np.isnan(tck[1])) == True:
                raise Exception("Making flat field model failed - NaNs in the spline fit.")
            mod=np.array([interpolate.splev(x, tck)]*data.shape[0])
            flatfield=flatfield/mod
            #plt.matplotlib.interactive(True)
            #plt.plot(x, med)
            #plt.plot(x, interpolate.splev(x, tck), 'r-')
            
            zeroMask=np.equal(flatfield, 0)
            nonZeroMask=np.not_equal(flatfield, 0)
            data[nonZeroMask]=data[nonZeroMask]/flatfield[nonZeroMask]
            data[zeroMask]=0.0
            img[extension].data=data
        
        if os.path.exists(f) == True:
            os.remove(f)
            
        img.writeto(f, overwrite = True)

#-------------------------------------------------------------------------------------------------------------
def findMatchingFilesByTime(inputFileName, possibleFilesList, timeInterval = 3600.0):
    """Returns a list of filenames in possibleFilesList that were obtained within timeInterval (in seconds) of
    the inputFileName. Use to find corresponding arcs, flats.
    
    Set timeInterval = None, to retrieve only the nearest file in time.
    
    Returns list of fileNames
    
    """
        
    ctimes=[]
    for f in possibleFilesList:
        ctime=getCTimeFromHeader(f)
        ctimes.append(ctime)
    ctimes=np.array(ctimes)
        
    fileCTime=getCTimeFromHeader(inputFileName)
    
    tab=atpy.Table()
    tab.add_column(atpy.Column(possibleFilesList, 'fileNames'))
    tab.add_column(atpy.Column(abs(fileCTime-ctimes), 'dt'))
    tab.sort('dt')
    
    if timeInterval is not None:
        matchedList=tab['fileNames'][np.where(tab['dt'] < timeInterval)].tolist()
    else:
        matchedList=tab['fileNames'][np.where(tab['dt'] == tab['dt'].min())].tolist()
    
    return matchedList

#-------------------------------------------------------------------------------------------------------------
def detectLines(data, sigmaCut = 3.0, thresholdSigma = 2.0, featureMinPix = 30, numBins = 1):
    """Detect lines in a 2d (or 1d) arc spectrum. If 2d, uses the central row of the 2d spectrum only.
    
    Returns: featureTable, segmentationMap
    
    """
    
    #---
    # More complicated arc line detection, which allows varying background along dispersion direction
    # This reduces to the old version when numBins = 1
    if data.ndim == 2:
        binEdges=np.round(np.linspace(0, data.shape[1], numBins+1))
    elif data.ndim == 1:
        binEdges=np.round(np.linspace(0, data.shape[0], numBins+1))

    mean=np.zeros(binEdges.shape[0]-1)
    sigma=1e6*np.ones(binEdges.shape[0]-1)
    for i in range(20):
        for k in range(binEdges.shape[0]-1):
            if data.ndim == 2:
                dataSlice=data[:, int(binEdges[k]):int(binEdges[k+1])]
            elif data.ndim == 1:
                dataSlice=data[int(binEdges[k]):int(binEdges[k+1])]
            nonZeroMask=np.not_equal(dataSlice, 0)
            mask=np.less(abs(dataSlice-mean[k]), sigmaCut*sigma[k])
            mask=np.logical_and(nonZeroMask, mask)
            mean[k]=np.mean(dataSlice[mask])
            sigma[k]=np.std(dataSlice[mask])
    detectionThreshold=thresholdSigma*sigma
    
    mask=np.zeros(data.shape)
    for k in range(binEdges.shape[0]-1):
        if data.ndim == 2:
            dataSlice=data[:, int(binEdges[k]):int(binEdges[k+1])]
            mask[:, int(binEdges[k]):int(binEdges[k+1])]=np.greater(dataSlice-mean[k], detectionThreshold[k])
        elif data.ndim == 1:
            dataSlice=data[int(binEdges[k]):int(binEdges[k+1])]
            mask[int(binEdges[k]):int(binEdges[k+1])]=np.greater(dataSlice-mean[k], detectionThreshold[k])
    
    #---
    # Old, global method
    # Detect arc lines
    #mean=0
    #sigma=1e6
    #for i in range(20):
        #nonZeroMask=np.not_equal(data, 0)
        #mask=np.less(abs(data-mean), sigmaCut*sigma)
        #mask=np.logical_and(nonZeroMask, mask)
        #mean=np.mean(data[mask])
        #sigma=np.std(data[mask])
    #detectionThreshold=thresholdSigma*sigma
    #mask=np.greater(data-mean, detectionThreshold)
    #---
    
    # Get feature positions, number of pixels etc.
    # Find features in 2d, match to wavelength coord in centre row
    segmentationMap, numObjects=ndimage.label(mask)
    sigPixMask=np.equal(mask, 1)
    objIDs=np.unique(segmentationMap)
    objNumPix=ndimage.sum(sigPixMask, labels = segmentationMap, index = objIDs)
    if data.ndim == 2:
        #objPositions_centreRow=ndimage.center_of_mass(data[data.shape[0]/2], labels = segmentationMap, index = objIDs)
        objPositions_centreRow=ndimage.maximum_position(data[int(data.shape[0]/2)], labels = segmentationMap[int(data.shape[0]/2)], index = objIDs)
        objAmplitudes_centreRow=ndimage.maximum(data[int(data.shape[0]/2)], labels = segmentationMap[int(data.shape[0]/2)], index = objIDs)
    elif data.ndim == 1:
        # ndmage.centre_of_mass can be led astray... just use local maximum
        #objPositions_centreRow=ndimage.center_of_mass(data, labels = segmentationMap, index = objIDs)
        objPositions_centreRow=ndimage.maximum_position(data, labels = segmentationMap, index = objIDs)
        objAmplitudes_centreRow=ndimage.maximum(data, labels = segmentationMap, index = objIDs)

    objPositions_centreRow=np.array(objPositions_centreRow).flatten()
    objAmplitudes_centreRow=np.array(objAmplitudes_centreRow).flatten()
    minPixMask=np.greater(objNumPix, featureMinPix)
    featureTable=atpy.Table()
    featureTable.add_column(atpy.Column(objIDs[minPixMask], 'id'))
    featureTable.add_column(atpy.Column(objPositions_centreRow[minPixMask], 'x_centreRow'))
    if data.ndim == 2:
        featureTable.add_column(atpy.Column([int(data.shape[0]/2)]*len(featureTable), 'y_centreRow'))
        featureTable.add_column(atpy.Column(data[int(data.shape[0]/2), np.array(np.round(featureTable['x_centreRow']), dtype = int)], 'amplitude')) 
    elif data.ndim == 1:
        featureTable.add_column(atpy.Column(objAmplitudes_centreRow[minPixMask], 'amplitude'))

    return featureTable, segmentationMap

#-------------------------------------------------------------------------------------------------------------
def fftCorrelate(f, g):
    """Does zero-padded fft correlation between arrays f, g.
    
    Returns corr, corrMax, shift
    
    """
    
    # Upsampling doesn't help unless go factor of several
    upSample=10.0
    fUp=ndimage.zoom(f, upSample)
    gUp=ndimage.zoom(g, upSample)
    
    # Zero padding
    numPaddedSamples=len(fUp)*2
    n=None
    for i in range(1, 30):
        if 2**i > numPaddedSamples:
            n=i
            break
    if n == None:
        raise Exception("Wavelength range covered is too big!")
    fPadded=np.zeros(2**n)
    fPadded[int(fUp.shape[0]/2):int(fUp.shape[0]/2+fUp.shape[0])]=fUp[:]
    gPadded=np.zeros(2**n)
    gPadded[int(gUp.shape[0]/2):int(gUp.shape[0]/2+gUp.shape[0])]=gUp[:]
    
    # FFT correlate
    fTemplate=np.fft.fft(fPadded)
    fSpectrum=np.fft.fft(np.fft.fftshift(gPadded))
    fxCorr=fTemplate*fSpectrum.conj()
    corr=np.fft.ifft(fxCorr).real
    
    # Get shift, accounting for zero padding and upsampling
    corrMax=corr.max()
    corrMaxIndex=float(np.argmax(corr))
    shift=(corrMaxIndex-float(fPadded.shape[0])/2.-1)/upSample
    
    return corr, corrMax, shift

#-------------------------------------------------------------------------------------------------------------
def minFunc_findScale(s, shift, arcRow, normRefModel, data_x):
    """For optimize.minimise_scalar.
    
    """
    
    tck=interpolate.splrep((data_x+shift)+s*data_x, arcRow)
    arcRow_scaled=interpolate.splev(data_x, tck, ext = 1)
    arcMean=np.mean(arcRow_scaled)
    arcStd=np.std(arcRow_scaled)
    

    # Old style: minimize what we call overlap below - this behaves well with optimize.minimize
    # (overlap vs scale is a function with a clear minimum)
    overlap=np.trapz(abs(normRefModel[:data_x.shape[0]]-(arcRow_scaled-arcMean)/arcStd))
        
    return overlap

#-------------------------------------------------------------------------------------------------------------
def minFunc_findShift(shift, scale, arcRow, normRefModel, data_x):
    """For optimize.minimise_scalar.
    
    """
    
    tck=interpolate.splrep((data_x+shift)+scale*data_x, arcRow)
    arcRow_scaled=interpolate.splev(data_x, tck, ext = 1)
    arcMean=np.mean(arcRow_scaled)
    arcStd=np.std(arcRow_scaled)
    
    # Old style: minimize what we call overlap below - this behaves well with optimize.minimize
    # (overlap vs scale is a function with a clear minimum)
    overlap=np.trapz(abs(normRefModel[:data_x.shape[0]]-(arcRow_scaled-arcMean)/arcStd))
        
    return overlap

#-------------------------------------------------------------------------------------------------------------
def minFunc_findShiftAndScale(x, arcRow, normRefModel, data_x):
    """For optimize.minimise - return 1/corrMax
    
    """
    
    #t0=time.time()
    shift=x[0]
    s=x[1]
    tck=interpolate.splrep((data_x+shift)+s*data_x, arcRow)
    arcRow_scaled=interpolate.splev(data_x, tck, ext = 1)
    arcMean=np.mean(arcRow_scaled)
    arcStd=np.std(arcRow_scaled)
    
    # Old style: minimize what we call overlap below - this behaves well with optimize.minimize
    # (overlap vs scale is a function with a clear minimum)
    overlap=np.trapz(abs(normRefModel[:data_x.shape[0]]-(arcRow_scaled-arcMean)/arcStd))
    #t1=time.time()
    #print "minFunc_findShiftAndScale took %.3f sec" % (t1-t0)

    return overlap

#-------------------------------------------------------------------------------------------------------------
def findScaleAndShift(arcRow, refModelDict, numScales = 101):
    """Find best fit stretch and scale to transform arcRow to reference model.
    
    Returns maxmimum of cross correlation, and corresponding best fit scale and shift
        
    """
 
    tab, arcRowSegMap=detectLines(arcRow)
    arcRowSegMap[np.greater(arcRowSegMap, 0)]=1.
    tab, refModelSegMap=detectLines(refModelDict['arc_centreRow'])
    refModelSegMap[np.greater(refModelSegMap, 0)]=1.
    
    # fftCorrelate
    data_x=np.arange(0, arcRowSegMap.shape[0])  
    corrMaxList=[]
    scalesArr=np.linspace(-0.2, 0.2, numScales)
    shiftsList=[]
    #overlapsList=[]
    for fftScale in scalesArr:
        tck=interpolate.splrep(data_x+fftScale*data_x, arcRowSegMap)
        arcRow_scaled=interpolate.splev(data_x, tck, ext = 1)
        corr, corrMax, fftShift=fftCorrelate(refModelSegMap, arcRow_scaled) 
        #overlap=minFunc_findShiftAndScale([fftShift, fftScale], arcRowSegMap, refModelSegMap, data_x)
        #overlapsList.append(overlap)
        corrMaxList.append(corrMax)
        shiftsList.append(fftShift)
    index=corrMaxList.index(max(corrMaxList))
    shift=shiftsList[index]
    corrMax=corrMaxList[index]
    s=scalesArr[index]
    #minOverlap=overlapsList[index]
    #corrMax=1./minOverlap
    #print "... using fftCorrelate shift = %.3f, scale = %.3f ..." % (shift, s)
    
    return corrMax, s, shift
  
#-------------------------------------------------------------------------------------------------------------
def selectBestRefModel(modelFileNameList, arcData, thresholdSigma = 2.0):
    """Returns the reference arc model which maximises the cross correlation with the arc data.
    
    """

    bestCorrMaxList=[]
    bestFitShiftList=[]
    bestFitScaleList=[]
    refModelDictsList=[]
    arcFeatureTablesList=[]
    arcSegMapsList=[]
    fitDictList=[]
    for modelFileName in modelFileNameList:
    
        # Load reference model
        pickleFile=open(modelFileName, "rb")
        unpickler=pickle.Unpickler(pickleFile, **PICKLE_OPTIONS)
        refModelDict=unpickler.load()
        refModelDictsList.append(refModelDict)
        pickleFile.close()
        
        # Helps for tracking later
        refModelDict['modelFileName']=os.path.split(modelFileName)[-1]
        
        # First need to find arc features
        arcFeatureTable, arcSegMap=detectLines(arcData, thresholdSigma = thresholdSigma)
        arcFeatureTablesList.append(arcFeatureTable)
        arcSegMapsList.append(arcSegMap)
        
        # Replaced np.correlate with fft based correlation
        # Find shift and wavelength dependent scale change (stretch, then shift)
        arcRow=arcData[int(arcData.shape[0]/2)]
        bestCorrMax, bestFitScale, bestFitShift=findScaleAndShift(arcRow, refModelDict, numScales = 22)
        bestCorrMaxList.append(bestCorrMax)
        bestFitScaleList.append(bestFitScale)
        bestFitShiftList.append(bestFitShift)
 
    # So... which arc model worked best? Use it...
    bestModelIndex=np.argmax(bestCorrMaxList)
    refModelDict=refModelDictsList[bestModelIndex]
    arcFeatureTable=arcFeatureTablesList[bestModelIndex]
    arcSegMap=arcSegMapsList[bestModelIndex]
    bestFitScale=bestFitScaleList[bestModelIndex]
    bestFitShift=bestFitShiftList[bestModelIndex]
    
    logger.info("using refModel = %s" % (modelFileNameList[bestModelIndex]))
    
    return refModelDict, arcFeatureTable, arcSegMap
    
#-------------------------------------------------------------------------------------------------------------
def findWavelengthCalibration(arcData, modelFileName, sigmaCut = 3.0, thresholdSigma = 2.0, 
                              featureMinPix = 50, order = 2, diagnosticsDir = None, diagnosticsLabel = None):
    """Find wavelength calibration for .fits image arcFileName containing a 2d arc spectrum.
    
    modelFileName is the path to a model made by makeModelArcSpectrum.
    
    Returns a dictionary with keys:
        fitCoeffsArr                - polynomial fit coefficients for use with wavelengthCalibrateAndRectify
        refModelFileName            - the reference model that was selected
        numArcFeaturesIdentified    - the number of arc lines found
    
    """

    # Track some useful info here - will get written into the diagnostics dir later
    diagnosticDict={}
    
    # We now allow multiple reference models for each grating/lamp/binning config
    # This is useful if e.g., we have some MOS slit which is way to the red/blue end of the detector
    # First select the model to use based on centre row only (saves much time)
    # Choose best from maximum cross correlation
    # NOTE: This now searches arcs with different spatial binning too
    modelFileNameList=glob.glob(modelFileName.split(".pickle")[0][:-1]+"*.pickle")
    refModelDict, arcFeatureTable, arcSegMap=selectBestRefModel(modelFileNameList, arcData, 
                                                                thresholdSigma = thresholdSigma)
    diagnosticDict['refModelFileName']=refModelDict['modelFileName']
    
    # Check sizes match - sometimes they differ by a pixel, why is a bit of a mystery...
    # It's safe to truncate at the right
    maxLength=min(refModelDict['arc_centreRow'].shape[0], arcData.shape[1])
    arcData=arcData[:, :maxLength]
    arcSegMap=arcSegMap[:, :maxLength]
    refModelDict['arc_centreRow']=refModelDict['arc_centreRow'][:maxLength]
    
    # Continue with previous 2d wavelength calib method
    # This step is important, but not crucial (although when it does go wrong, it is spectacular)
    # This just affects which lines are cross-identified between the arc and the reference model
    yIndex=int(arcData.shape[0]/2)
    arcRow=arcData[yIndex]
    bestCorrMax, bestFitScale, bestFitShift=findScaleAndShift(arcRow, refModelDict, numScales = 401)
    arc_centreRow=arcRow
    
    # Arc transform test plot
    #plt.matplotlib.interactive(True)
    data_x=np.arange(0, arc_centreRow.shape[0])        
    x=np.arange(0, len(arc_centreRow))
    arc_x_shifted=x*(1+bestFitScale)+bestFitShift
    plt.figure(figsize=(12,8))
    plt.plot(x, refModelDict['arc_centreRow'][:data_x.shape[0]]/refModelDict['arc_centreRow'][:data_x.shape[0]].mean(), 'b-', label ='ref model')
    plt.plot(arc_x_shifted, arc_centreRow/arc_centreRow.mean(), 'r-', label = 'transformed arc')
    plt.title("shift = %.3f, scale = %.3f" % (bestFitShift, bestFitScale))
    plt.semilogy()
    plt.legend()
    plt.savefig(diagnosticsDir+os.path.sep+"arcTransformTest_"+diagnosticsLabel+".png")
    plt.close()    
    
    # Tag features by transforming model coords to arc coords and looking for closest match
    # Set maxDistancePix high enough to get enough matches, but not so high bring in rubbish
    # Weak link may be centroiding done in detectLines
    arcFeatureTable.add_column(atpy.Column(np.zeros(len(arcFeatureTable)), 'wavelength'))    
    maxDistancePix=10.0
    for row in refModelDict['featureTable']:
        transformed_model_x=(row['x_centreRow']-bestFitShift)/(1+bestFitScale)
        dist=abs(arcFeatureTable['x_centreRow']-transformed_model_x)
        if dist.min() < maxDistancePix:
            index=np.argmin(dist)
            arcFeatureTable['wavelength'][index]=row['wavelength']
    arcFeatureTable=arcFeatureTable[np.where(arcFeatureTable['wavelength'] != 0)]
    diagnosticDict['numArcFeaturesIdentified']=len(arcFeatureTable)
    if len(arcFeatureTable) < 5:
        logger.info("aborted - less than 5 features identified in the arc")
        diagnosticDict['fitCoeffsArr']=None
        return diagnosticDict
    
    # Tagging check in pixel space
    #plt.matplotlib.interactive(True)
    plt.figure(figsize=(12, 8))
    plt.plot(arcRow, 'k-')
    plt.plot(arcFeatureTable['x_centreRow'], arcFeatureTable['amplitude'], 'bo')
    for row in arcFeatureTable:
        plt.text(row['x_centreRow'], row['amplitude'], row['wavelength'])
    plt.xlabel("x (pixels)")
    plt.ylabel("Relative Flux")
    plt.title("Features tagged")
    plt.semilogy()
    plt.savefig(diagnosticsDir+os.path.sep+"taggedFeaturesPixelCoords_"+diagnosticsLabel+".png")
    plt.close()
        
    # Find 2d wavelength solution which we can use for rectification/wavelength calibration
    # Fit functions for how feature x positions change with y
    ys=np.arange(arcData.shape[0])
    arcFeatureTable.add_column(atpy.Column(np.zeros(len(arcFeatureTable)), 'slope'))
    arcFeatureTable.add_column(atpy.Column(np.zeros(len(arcFeatureTable)), 'intercept'))
    for row in arcFeatureTable:
        xs=np.zeros(arcData.shape[0])
        for i in range(len(ys)):
            objPositions=ndimage.center_of_mass(arcData[ys[i]], labels = arcSegMap, index = arcFeatureTable['id'])
            xs[i]=float(objPositions[np.where(arcFeatureTable['id'] == row['id'])[0][0]][0])
        # Linear fit should allow us to work out shear - here, given y, we want x
        # We probably don't want all this info (ys, xs), but keep for now
        # We could use polynomial instead of linear (see below)
        try:
            slope, intercept=np.polyfit(ys, xs, 1)
        except:
            raise Exception("polyfit failed")
        row['slope']=slope
        row['intercept']=intercept

    # Wavelength calibration and model with arbitrary order polynomials - get coeffs for each row
    # Now iterating, because can blow up if point on edge is way off
    # On each iteration, we mask the highest residual above a threshold of delta wavelength > 100 Angstroms
    wavelengths=arcFeatureTable['wavelength']
    fitCoeffsArr=[]
    for y in range(arcData.shape[0]):
        xs=[]
        for row in arcFeatureTable:
            xs.append(row['slope']*y + row['intercept'])
        xs=np.array(xs)
        mask=np.array([True]*len(xs))
        removedPoint=True
        while removedPoint == True:
            coeffs=np.polyfit(xs[mask], wavelengths[mask], order)
            testPoly=np.poly1d(coeffs)
            testWavelengths=testPoly(xs)
            res=abs(wavelengths-testWavelengths)
            if res[mask].max() > 100.:
                removedPoint=True
                mask[np.equal(res, res[mask].max())]=False
            else:
                removedPoint=False
        fitCoeffsArr.append(coeffs)
    fitCoeffsArr=np.array(fitCoeffsArr)
    
    # Check of wavelength calibration fit
    wavelengthCalibPoly=np.poly1d(fitCoeffsArr[yIndex])
    testWavelengths=wavelengthCalibPoly(np.arange(arcRow.shape[0]))
    plt.figure(figsize=(12,8))
    plt.title("refModelFileName = %s, numArcFeaturesIdentified = %d" % (diagnosticDict['refModelFileName'], diagnosticDict['numArcFeaturesIdentified']))
    plt.plot(np.arange(arcRow.shape[0]), testWavelengths, 'k--')
    plt.plot(arcFeatureTable['x_centreRow'], arcFeatureTable['wavelength'], 'r.')
    plt.xlabel("x (pixels)")
    plt.ylabel("Wavelength (Angstroms)")
    plt.savefig(diagnosticsDir+os.path.sep+"wavelengthCalibModelFit_"+diagnosticsLabel+".png")
    plt.close()

    # Sanity check: tagged features in wavelength calibrated spectrum
    plt.figure(figsize=(12, 8))
    plt.title("Wavelength calibration check")
    plt.plot(testWavelengths, arcRow, 'k-')
    plt.plot(arcFeatureTable['wavelength'], arcFeatureTable['amplitude'], 'bo')
    for row in arcFeatureTable:
        plt.text(row['wavelength'], row['amplitude'], row['wavelength'])
    plt.xlabel("Wavelength (Angstroms)")
    plt.ylabel("Relative Flux")
    plt.semilogy()
    plt.savefig(diagnosticsDir+os.path.sep+"taggedFeaturesWavelengthCoords_"+diagnosticsLabel+".png")
    plt.close()
        
    resultDict={'fitCoeffsArr': fitCoeffsArr}
    for key in diagnosticDict:
        resultDict[key]=diagnosticDict[key]
    
    return resultDict

#-------------------------------------------------------------------------------------------------------------
def wavelengthCalibrateAndRectify(inFileName, outFileName, wavelengthCalibDict, extensionsList = "all",
                                  makeDiagnosticPlots = False):
    """Applies the wavelength calibration, and rectification, to all extensions of inFileName, writing 
    output to outFileName. The wavelength calibration is provided in wavelengthCalibDict, where each key
    corresponds to each extension number (see findWavelengthCalibration)
    
    """
    
    logger.info("Applying wavelength calibration and rectifying (%s)" % (inFileName))
    img=pyfits.open(inFileName)
    if extensionsList == "all":
        extensionsList=[]
        for hdu in img:
            if "SLIT" in hdu.name:
                extensionsList.append(hdu.name)
                
    for extension in extensionsList:
    
        logger.info("%s" % (extension))
        
        data=img[extension].data
        header=img[extension].header
        fitCoeffsArr=wavelengthCalibDict[extension]['fitCoeffsArr']
        
        # Can carry on if wavelength calib fails for a slit... fix later...
        if np.any(fitCoeffsArr) != None:
            
            # Using above, make an array containing wavelengths
            wavelengthsMap=np.zeros(data.shape)
            for y in range(data.shape[0]):
                wavelengthCalibPoly=np.poly1d(fitCoeffsArr[y])
                wavelengthsMap[y]=wavelengthCalibPoly(np.arange(data.shape[1]))
            #astImages.saveFITS("wavelengthsMap.fits", wavelengthsMap, None)
            
            # How we would want our wavelength map to look after applying some transformation
            # To make things easier later, make a linear wavelength scale
            wavelengths_centreRow=wavelengthsMap[int(wavelengthsMap.shape[0]/2)]
            maxWavelength=wavelengths_centreRow.max()
            minWavelength=wavelengths_centreRow.min()
            linearWavelengthRange=np.linspace(minWavelength, maxWavelength, data.shape[1])
            FITSWavelengthScale=linearWavelengthRange[1]-linearWavelengthRange[0]
            FITSRefLambda=linearWavelengthRange[0]
            FITSRefPixel=1                              # Remember index from 1 is FITS convention
            rectWavelengthsMap=np.array([linearWavelengthRange]*data.shape[0])
            #astImages.saveFITS("rectWavelengthsMap.fits", rectWavelengthsMap, None)
            
            # Debugging
            #if os.path.split(inFileName)[-1] == "cmbxgpP201211100057.fits":
                #print "check mapping of spectra to .fits"
                #IPython.embed()
                #sys.exit()
                
            # Remap the data to our preferred linear wavelength scale
            # Assume we can treat each row independently
            # Save linear spectral WCS in header
            rectifiedData=np.zeros(data.shape)
            for y in range(data.shape[0]):
                try:
                    tck=interpolate.splrep(wavelengthsMap[y], data[y])
                    rectifiedData[y]=interpolate.splev(rectWavelengthsMap[y], tck, ext = 1)
                except:
                    logger.warning("wavelength calibration/rectification error: this slit will be blank")
                    break
            img[extension].data=rectifiedData
            header['CTYPE1']='LINEAR'
            header['DISPAXIS']=1
            header['CRVAL1']=FITSRefLambda
            header['CRPIX1']=FITSRefPixel
            header['CD1_1']=FITSWavelengthScale
            header['CDELT1']=FITSWavelengthScale
            header['CUNIT1']='Angstroms'
                
            ## Sanity check plot: linear wavelength scale
            #if makeDiagnosticPlots == True:
                #diagnosticsDir=os.path.split(outFileName)[0]+os.path.sep+"diagnostics"
                #if os.path.exists(diagnosticsDir) == False:
                    #os.makedirs(diagnosticsDir)
                #plt.plot(rectWavelengthsMap[data.shape[0]/2], rectifiedData[data.shape[0]/2], 'k-')
                #plt.xlabel("Wavelength (Angstroms)")
                #plt.ylabel("Relative Intensity")
                #plt.title("%s - %s" % (os.path.split(inFileName)[-1], extension))
                #plt.savefig(diagnosticsDir+os.path.sep+"wavelengthCalibCheck_%s_%s.png" % (os.path.split(outFileName)[-1].replace(".fits", ""), extension))
                #plt.close()

    # Write output
    if os.path.exists(outFileName) == True:
        os.remove(outFileName)
        
    img.writeto(outFileName, overwrite = True)
    
#-------------------------------------------------------------------------------------------------------------
def wavelengthCalibration2d(maskDict, outDir, extensionsList = "all", modelArcsDir = None):
    """Finds 2d wavelength calibration from arc frames, applies to arc frames and object frames, rectifying
    them and also interpolating to a linear wavelength scale to make life easier later.
    
    Should be fully automatic, assuming a suitable reference model is available.
    
    """
       
    diagnosticsDir=outDir+os.path.sep+"diagnostics"
    if os.path.exists(diagnosticsDir) == False:
        os.makedirs(diagnosticsDir)
        
    logger.info("Finding 2d wavelength solution")
    maskDict['wavelengthCalib']={}
    for key in maskDict['cutArcDict'].keys():
        
        cutArcPath=maskDict['cutArcDict'][key]
        
        logger.info("arc = %s" % (cutArcPath))
              
        img=pyfits.open(cutArcPath)
        
        if extensionsList == "all":
            extensionsList=[]
            for hdu in img:
                if "SLIT" in hdu.name:
                    extensionsList.append(hdu.name)

        binning=img[0].header['CCDSUM'].replace(" ", "x")
        grating=img[0].header['GRATING']
        lampid=img[0].header['LAMPID']
        
        # Now we don't care about spatial binning for finding arcs
        binning=binning[:-1]+"?"
        modelFileNames=glob.glob(REF_MODEL_DIR+os.path.sep+"RefModel_"+grating+"_"+lampid+"_"+binning+".pickle")
        if modelArcsDir is not None:
            modelFileNames=modelFileNames+glob.glob(modelArcsDir+os.path.sep+"RefModel_"+grating+"_"+lampid+"_"+binning+".pickle")

        if cutArcPath not in maskDict['wavelengthCalib'].keys():
            maskDict['wavelengthCalib'][cutArcPath]={}
            for extension in extensionsList:
                logger.info("extension = %s" % (extension))
                arcData=img[extension].data
                foundModelFile=False
                for modelFileName in modelFileNames:
                    if os.path.exists(modelFileName) == True:
                        foundModelFile=True
                        break
                if foundModelFile == False:
                    print("No reference model exists for grating %s, lamp %s, with binning %s." % (grating, lampid, binning))
                    print("Use rss_mos_create_arc_model to make a reference model and then re-run.")
                    print("arcFileName: %s" % (cutArcPath))
                    sys.exit()
                diagnosticsLabel=os.path.split(cutArcPath)[-1].replace(".fits", "")+"_"+extension

                maskDict['wavelengthCalib'][cutArcPath][extension]=findWavelengthCalibration(arcData, modelFileName, diagnosticsDir = diagnosticsDir, diagnosticsLabel = diagnosticsLabel)
    
    # Write out some diagnostic info (ref model used, number of arc lines found
    outFile=open(diagnosticsDir+os.path.sep+"wavelengthCalibDiagnostics.csv", "w")
    outFile.write("#cutArcFileName\textension\trefModelFileName\tnumArcFeaturesIdentified\n") 
    for arcKey in maskDict['wavelengthCalib'].keys():
        slitKeys=list(maskDict['wavelengthCalib'][arcKey].keys())
        slitKeys.sort()
        for slitKey in slitKeys:
            refModelFileName=maskDict['wavelengthCalib'][arcKey][slitKey]['refModelFileName']
            numArcFeaturesIdentified=maskDict['wavelengthCalib'][arcKey][slitKey]['numArcFeaturesIdentified']
            outFile.write("%s\t%s\t%s\t%d\n" % (os.path.split(arcKey)[-1], slitKey, 
                                                refModelFileName, numArcFeaturesIdentified))
    outFile.close()
    
    # Apply the calibration to the arc frames (diagnostic purposes)
    # 'rw' prefix => rectified and wavelength calibrated
    for key in maskDict['cutArcDict'].keys():
        cutArcPath=maskDict['cutArcDict'][key]                
        rectArcPath=makeOutputFileName(cutArcPath, "rw", outDir)
        wavelengthCalibrateAndRectify(cutArcPath, rectArcPath, maskDict['wavelengthCalib'][cutArcPath],
                                      extensionsList = extensionsList, makeDiagnosticPlots = True)
    
    # Apply the calibration to object spectra           
    for fileName in maskDict['OBJECT']:
        cutArcPath=maskDict['cutArcDict'][fileName]                
        cutPath=makeOutputFileName(fileName, "c", outDir)
        rectPath=makeOutputFileName(fileName, "rwc", outDir)
        wavelengthCalibrateAndRectify(cutPath, rectPath, maskDict['wavelengthCalib'][cutArcPath], 
                                      extensionsList = extensionsList)        
   
#-------------------------------------------------------------------------------------------------------------
def measureProfile(data, mask, minTraceWidth = 4., halfBlkSize = 50, sigmaCut = 3.):
    """Used in the spectral extraction to fit the object profile in the y direction.
    
    Now using masked arrays, which makes this robust when iterating.
    
    """
    
    ## Find local background, noise (running clipped mean)
    ## NOTE: direct from mos pipeline for finding pseudo slits
    #d=data
    #prof=np.median(d, axis = 1)    
    #prof[np.less(prof, 0)]=0.     
    #bck=np.zeros(prof.shape)
    #sig=np.zeros(prof.shape)
    #for y in range(prof.shape[0]):
        #yMin=y-halfBlkSize
        #yMax=y+halfBlkSize
        #if yMin < 0:
            #yMin=0
        #if yMax > prof.shape[0]-1:
            #yMax=prof.shape[0]-1
        #mean=0
        #sigma=1e6
        #for i in range(20):
            #nonZeroMask=np.not_equal(prof[yMin:yMax], 0)
            #mask=np.less(abs(prof[yMin:yMax]-mean), sigmaCut*sigma)
            #mean=np.mean(prof[yMin:yMax][mask])
            #sigma=np.std(prof[yMin:yMax][mask])            
        #bck[y]=mean
        #if sigma > 0:
            #sig[y]=sigma
        #else:
            #sig[y]=np.std(prof)

    # Non-running version of the above ^^^
    # Tweaked to not break in the case of complete signal domination
    d=np.ma.masked_array(data, mask)
    prof=np.ma.median(d, axis = 1)    
    prof[np.less(prof, 0)]=0.     
    mean=np.mean(prof)
    sigma=np.std(prof)
    for i in range(10):
        nonZeroMask=np.not_equal(prof, 0)
        mask=np.less(abs(prof-mean), sigmaCut*sigma)
        mean=np.mean(prof[mask])
        sigma=np.std(prof[mask])            
    bck=mean
    sig=sigma
            
    # Not sure if this will work all the time...
    profSNR=(prof-bck)/sig    
    prof=profSNR
    prof[np.less(prof, 0)]=0.
    try:
        prof=prof/prof.max()
    except:
        logger.warning("profile measurement failed")
        prof=np.zeros(data.shape[0])
    
    if np.any(np.isnan(prof)) == True:
        raise Exception("nans in object profile")
        
    # Sometimes we get rubbish in slit edges being counted as signal - spot that, and remove it
    segmentationMap, numObjects=ndimage.label(prof)
    if numObjects > 1:
        sigPixMask=np.greater(prof, 0)
        objIDs=np.unique(segmentationMap)
        objNumPix=ndimage.sum(sigPixMask, labels = segmentationMap, index = objIDs)
        signalObjID=objIDs[np.where(objNumPix == objNumPix.max())][0]
        for objID in objIDs:
            if objID != signalObjID:
                mask=np.equal(segmentationMap, objID)
                prof[mask]=0
            
    return prof
    
#-------------------------------------------------------------------------------------------------------------
def iterativeWeightedExtraction(data, maxIterations = 1000, subFrac = 0.8, runningProfile = None, \
                                iterateProfile = False, throwAwayRows = 2):
    """Extract 1d spectrum of object, sky, and find noisy pixels affected by cosmic rays while we're at it.
    This is somewhat similar to the Horne optimal extraction algorithm. We solve:

    ws1*s + k + wn1*n1 = v1
    ws2*s + k + wn2*n2 = v2
    ws3*s + k + wn3*n3 = v3

    |ws1 1 wn1 0   0   | x | s  | = | v1 |
    |ws2 1 0   wn2 0   |   | k  | = | v2 |
    |ws3 1 0   0   wn3 |   | n1 | = | v3 |
                           | n2 |
                           | n3 |

    where ws = signal weight, wk = sky weight == 1, wn = noise weight, s = signal, k = sky, n = noise.

    Sky weight has to be 1, because we have signal + sky everywhere, and the sky level should be the same 
    across all rows. Noise weights are simply 0 or 1 and are used to mask cosmic rays. This seems pretty 
    effective.

    We do all this column by column. The sky estimate is subtracted iteratively.
        
    Returns extracted signal, sky (both 1d) and 2d masked array of data with CRs flagged
    We also fill CR-flagged pixels with median sky (in case want to regrid later for using
    other extraction methods)
        
    """

    logger.info("extracting spectrum")

    # Throw away rows at edges as these often contain noise
    if throwAwayRows > 0:
        data=data[throwAwayRows:-throwAwayRows]
    
    if data.shape[0] == 0:
        logger.warning("all data cut in iterativeWeightedExtraction")
        return np.zeros(data.shape[1]), np.zeros(data.shape[1]), np.zeros(data.shape)
    
    # Find the chip gaps and make a mask
    chipGapMask=makeChipGapMask(data)

    # First measurement of the profile of the object
    if np.any(runningProfile) == None:
        prof=measureProfile(data, np.zeros(data.shape))
    
    # Iterative sky subtraction
    wn2d=np.zeros(data.shape)+chipGapMask               # Treat chip gaps as noise
    skySub=np.zeros(data.shape)+data
    signalArr=np.zeros([maxIterations, data.shape[1]])
    diff=1e9
    tolerance=1e-5
    k=0
    skyTotal=np.zeros(data.shape[1]) # we need to add to this each iteration
    
    # For ignoring junk at the edges
    skyMask2d=identifySky(data)

    warned=False
    while diff > tolerance or k > maxIterations:
        t0=time.time()
        xArr=[]
        if iterateProfile == True and runningProfile == None:
            prof=measureProfile(skySub, wn2d)   # non-running prof
        for i in range(data.shape[1]):
            if np.any(runningProfile) != None:
                prof=runningProfile[:, i]
            # Running prof - this behaves strangely...
            #runWidth=200
            #iMin=i-runWidth
            #iMax=i+runWidth
            #if iMin < 0:
                #iMin=0
            #if iMax > data.shape[1]-1:
                #iMax=data.shape[1]-1
            #prof=measureProfile(skySub[:, iMin:iMax], wn2d[:, iMin:iMax])
            # As you were...
            b=skySub[:, i]
            A=np.zeros([b.shape[0]+2, b.shape[0]])
            A[0]=prof
            # Sky should be same everywhere... if something stands out, it isn't sky, so zap
            # This should take care of any rubbish at slit edges causing oversubtracted sky
            skyMask=maskNoisyData(b)
            skyMask=np.equal(skyMask, False)
            A[1]=np.array(skyMask, dtype = int)#1.#(1.-prof) 
            #A[1]=skyMask2d[:, i]
            # CR masking
            wn=wn2d[:, i]
            for j in range(b.shape[0]):
                A[2+j, j]=wn[j] # noise weights - if 1, zap that pixel (CR or bad);i f we flag a CR, set signal and sky in that pixel to zero weight
                if wn[j] == 1.0:
                    A[0, j]=0.0
                    A[1, j]=0.0
            try:
                x, R=optimize.nnls(A.transpose(), b)
            except:
                if warned == False:
                    logger.warning("nnls problem - setting column to zero")
                warned=True
                x=np.zeros(b.shape[0]+2)
            xArr.append(x)
        
        # Below here as usual
        xArr=np.array(xArr).transpose()
        sky=xArr[1]
        sky2d=np.array([sky]*data.shape[0])
        try:
            skySub=skySub-subFrac*sky2d 
        except:
            raise Exception("skySub failure")
        signal=xArr[0]
        signal[np.less(signal, 0)]=0.
        signalArr[k]=signal
        # CR rejection
        arr=skySub
        thresholdSigma=30.0 # was 30
        sigmaCut=3.0
        mean=0
        sigma=1e6
        lowMaskValue=2.0
        for i in range(20):
            gtrZeroMask=np.greater(arr, lowMaskValue)
            mask=np.less(abs(arr-mean), sigmaCut*sigma)
            mask=np.logical_and(gtrZeroMask, mask)
            mean=np.mean(arr[mask])
            sigma=np.std(arr[mask])
        detectionThreshold=thresholdSigma*sigma
        wn2d=np.array(np.greater(arr-mean, detectionThreshold), dtype = float)
        wn2d[np.less(arr, 0)]=0.0
        wn2d=wn2d+chipGapMask                   # Add in the mask for chip gaps
        wn2d[np.greater(wn2d, 1)]=1.0
        # Grow the CRMask
        wn2d=ndimage.uniform_filter(wn2d, 2)
        wn2d[np.greater(wn2d, 0.4)]=1.0
        
        # Did we converge?
        t1=time.time()
        if k > 0:
            diff=np.sum((signalArr[k]-signalArr[k-1])**2)
            logger.info("iteration %d (diff = %.3e, time taken = %.3f sec)" % (k, diff, t1-t0))
        k=k+1
        # Keep track of sky, we don't want to report just the residual
        skyTotal=skyTotal+subFrac*sky

    # Or... let's just use all of the above for CR-flagging
    mData=np.ma.masked_array(data, wn2d)
    
    return signal, skyTotal, mData

#-------------------------------------------------------------------------------------------------------------
def identifySky(data):
    """Scan the data by column, and flag pixels as being within +/- 3 sigma of the mean in the column.
    
    Returns a 2d array of sky weights, for use by weightedExtraction and iterativeWeightedExtraction.
    
    """

    skyMask2d=np.zeros(data.shape)
    for i in range(data.shape[1]):
        skyMask2d[:, i]=maskNoisyData(data[:, i])
        skyMask2d[:, i]=np.equal(skyMask2d[:, i], 0)
    
    return skyMask2d

#-------------------------------------------------------------------------------------------------------------
def weightedExtraction(data, medColumns = 10, thresholdSigma = 30.0, sigmaCut = 3.0, profSigmaPix = 4.0, \
                       throwAwayRows = 2):
    """Extract 1d spectrum of object, sky, and find noisy pixels affected by cosmic rays while we're at it.
    This was (supposed) to be similar to the Horne optimal extraction. We solve:

    ws1*s + k + wn1*n1 = v1
    ws2*s + k + wn2*n2 = v2
    ws3*s + k + wn3*n3 = v3

    |ws1 1 wn1 0   0   | x | s  | = | v1 |
    |ws2 1 0   wn2 0   |   | k  | = | v2 |
    |ws3 1 0   0   wn3 |   | n1 | = | v3 |
                           | n2 |
                           | n3 |

    where ws = signal weight, wk = sky weight == 1, wn = noise weight, s = signal, k = sky, n = noise.

    Sky weight has to be 1, because we have signal + sky everywhere, and the sky level should be the same 
    across all rows. Noise weights are simply 0 or 1 and are used to mask cosmic rays. This seems pretty 
    effective.

    We do all this column by column. 
    
    It turns out this is a great way to get the cosmic rays, but not great for sky line subtraction.
    
    Returns extracted signal, sky (both 1d) and 2d masked array of data with CRs flagged.
    We also fill CR-flagged pixels with median sky (in case want to regrid later for using
    other extraction methods)

    """

    # Throw away rows at edges as these often contain noise
    if throwAwayRows > 0:
        data=data[throwAwayRows:-throwAwayRows]
    
    # Find the chip gaps and make a mask
    chipGapMask=makeChipGapMask(data)

    # All of the below is really just CR rejection now...
    # Assume one object per slit and a fixed width
    #traceHalfWidth=4
    #prof=np.median(data, axis = 1)
    #peakIndex=np.where(prof == prof.max())[0]
    #x=np.arange(len(prof))
    #xMin=peakIndex-traceHalfWidth
    #try:
        #if xMin < 0:
            #xMin=0
    #except:
        #print "xMin seems to be array: check peakIndex"
        #IPython.embed()
        #sys.exit()
    #xMax=peakIndex+traceHalfWidth
    #if xMax > len(prof)-1:
        #xMax=len(prof)-1
    #prof[:xMin]=0.0
    #prof[xMax:]=0.0  

    # First measurement of the profile of the object
    prof=measureProfile(data, np.zeros(data.shape))
    if prof.sum() == 0:
        logger.warning("profile measurement is zero everywhere - returning empty spectrum")
        return np.zeros(data.shape[1]), np.zeros(data.shape[1]), np.zeros(data.shape)
        
    # Which rows to use as sky?
    skyMask2d=identifySky(data)
    
    # Fit/extract signal, sky (but really this is just a good way to find and mask cosmic rays)
    wn2d=np.zeros(data.shape)+chipGapMask               # Treat chip gaps as noise
    ws2d=np.zeros(data.shape)
    wk2d=np.zeros(data.shape)
    recSky=np.zeros(data.shape)
    for k in range(10):
        sky=np.zeros(data.shape[1])
        signal=np.zeros(data.shape[1])
        for i in range(data.shape[1]):
            # See testWeightedExtraction2.py, testWeightedExtraction4.py for alternatives using SVD, pseudoInverse etc. 
            v=data[:, i]
            w=np.zeros((v.shape[0], 2+v.shape[0]))
            w[:, 0]=prof.reshape(w[:, 0].shape)     # signal weight - varies across rows
            # Sky should be same everywhere... if something stands out, it isn't sky, so zap
            # This should take care of any rubbish at slit edges causing oversubtracted sky
            w[:, 1]=skyMask2d[:, i]#np.array(skyMask, dtype = int)#1.#(1.-prof) 
            #w[:, 1]=1                               # sky weight - needs to be the same across all rows
            wn=wn2d[:, i]
            for j in range(v.shape[0]):
                w[j, 2+j]=wn[j]                     # noise weights - if 1, zap that pixel (CR or bad)
                # Weights must sum to 1
                if wn[j] == 1.0:
                    w[j, 0]=0.0
                    w[j, 1]=0.0
            ws2d[:, i]=w[:, 0]
            wk2d[:, i]=w[:, 1]
            w=np.matrix(w)
            w=np.nan_to_num(w)
            x, R=optimize.nnls(w, v)
            signal[i]=x[0]
            sky[i]=x[1]

        # Detect noisy pixels
        # Construct a residual image and spot columns/pixels affected by cosmic ray hits etc.
        # Where we have CR hits, have -ve pixels in those columns - ignore those
        # To find, mask -ve in res, find 3sigma clipped median, sigma, and then mask those pixels, then iterate
        sky2d=np.array([sky]*data.shape[0])*wk2d
        signal2d=np.array([signal]*data.shape[0])*ws2d
        res2d=data-sky2d-signal2d
        
        arr=res2d
        thresholdSigma=30.0
        sigmaCut=3.0
        mean=0
        sigma=1e6
        lowMaskValue=2.
        for i in range(20):
            lastSigma=sigma
            gtrZeroMask=np.greater(arr, lowMaskValue)
            mask=np.less(abs(arr-mean), sigmaCut*sigma)
            mask=np.logical_and(gtrZeroMask, mask)
            mean=np.mean(arr[mask])
            sigma=np.std(arr[mask])
            if lastSigma == sigma:
                break
        detectionThreshold=thresholdSigma*sigma
        wn2d=np.array(np.greater(arr-mean, detectionThreshold), dtype = float)
        wn2d[np.less(arr, 0)]=0.0
        
        # Insert code (perhaps) to grow mask around CR hits (code is in testWeightedExtraction4)
        # However, experiments showed this sometimes made things worse - fix later...
        
        # Add in the mask for chip gaps
        wn2d=wn2d+chipGapMask
        wn2d[np.greater(wn2d, 1)]=1.0
            
    # We'll return masked array of the data - this has only CRs flagged (with 1s)
    mData=np.ma.masked_array(data, wn2d)

    # Fix masked values by filling with median sky
    # This works better than sklearn.imputer
    # skyFill masked values correspond only to the chip gaps
    skyFill=np.ma.median(mData, axis = 0)
    for i in range(mData.shape[0]):
        mData[i][mData[i].mask]=skyFill[mData[i].mask]
    
    return signal, sky, mData

#-------------------------------------------------------------------------------------------------------------
def checkWavelengthCalibUsingSky(sky, wavelengths, featureMinPix = 5, mismatchLimit = 50.):
    """Checks the wavelength calibration using the sky spectrum, matching against prominent lines.
    
    Returns the median wavelength offset (in Angstroms), and the number of lines used in the measurement.
    
    """
        
    # Find the chip gaps and make a mask - put this in its own routine...
    lowMaskValue=2.0
    minPix=1000
    chipGapMask=np.array(np.less(sky, lowMaskValue), dtype = float)  # flags chip gaps as noise
    segmentationMap, numObjects=ndimage.label(chipGapMask)
    sigPixMask=np.equal(chipGapMask, 1)
    objIDs=np.unique(segmentationMap)
    objNumPix=ndimage.sum(sigPixMask, labels = segmentationMap, index = objIDs)
    gapsList=[]
    for objID, nPix in zip(objIDs, objNumPix):    
        if nPix > 10 and nPix < minPix:
            chipGapMask[np.equal(segmentationMap, objID)]=0.0
            indices=np.where(segmentationMap == objID)[0]
            gapsList.append([wavelengths[min(indices)], wavelengths[max(indices)]])

    # To aid in finding sky lines, subtract off any large scale trend first...
    try:
        smoothSky=ndimage.uniform_filter1d(sky, 30)
    except:
        raise Exception("smoothSky problem")
    smoothSky=ndimage.uniform_filter1d(sky, int(sky.shape[0]/3))
    bckSubSky=sky-smoothSky
    
    # If detectLines misses the lines we intend to check, and picks up others, then we will get crazy answers
    # The simplest way to deal with this is to throw out crazy answers (> 50 Angstroms off say)
    cdelt=wavelengths[1]-wavelengths[0]
    crval1=wavelengths[0]
    featureTable, segMap=detectLines(bckSubSky, featureMinPix = featureMinPix, numBins = 8)
    featureTable.add_column(atpy.Column(featureTable['x_centreRow']*cdelt+crval1, 'wavelength'))
    diffs=[]
    usedSkyLines=[]
    usedArcLines=[]
    for c in checkSkyLines:
        inGap=False
        for g in gapsList:
            if c > g[0] and c < g[1]:
                inGap=True
        if inGap == False and c > wavelengths.min()+20 and c < wavelengths.max()-20:
            diff=featureTable['wavelength']-c
            mask=np.equal(np.abs(diff), np.abs(diff).min())
            if abs(diff[mask]) < mismatchLimit:
                diffs.append(diff[mask])
                usedSkyLines.append(c)
    
    return np.median(diffs), len(diffs)

#-------------------------------------------------------------------------------------------------------------
def extractAndStackSpectra(maskDict, outDir, extensionsList = "all", iterativeMethod = False, subFrac = 0.4):
    """Extracts and stacks spectra from science frames which have already been wavelength calibrated.
        
    Two methods are used for extracting spectra:
    
    1.  Individual 1d spectra are extracted from each science frame and then stacked. This might be 
        desirable if object traces shift from frame-to-frame. These are placed under 
        outDir/1DSpec_extractAndStack/
        
    2.  The 2d spectra are stacked and then a 1d spectrum is extracted. These are placed under 
        outDir/1DSpec_2DSpec_stackAndExtract/. The stacked 2d spectra are also placed in this
        directory.
        
    It seems like method (2) works best (Oct 2016).
    
    If iterativeMethod = True, then an iterative method for spectral extraction is used. In which
    case _iterative is added to the name of the output directory.
    
    subFrac is the fraction of the sky to be subtracted in each iteration of iterativeWeightedExtraction,
    if iterativeMethod = True.
        
    Cosmic rays are flagged in the first method, and masked out (replaced with sky) for the second method.
    In both cases, the median is used to stack the frames (or 1D spectra), and of course all data are
    projected onto a common wavelength scale first.
    
    Output 1d spectra are in .fits table format, with columns SPEC, SKYSPEC, and LAMBDA.
    Output 2d spectra (written as outDir/1DSpec_2DSpec/2D_*.fits) are fits images.
        
    """

    diagnosticsDir=outDir+os.path.sep+"diagnostics"
    if os.path.exists(diagnosticsDir) == False:
        os.makedirs(diagnosticsDir)
     
    extractStackSpecDir=outDir+os.path.sep+"1DSpec_extractAndStack"
    if iterativeMethod == True:
        extractStackSpecDir=extractStackSpecDir+"_iterative"
    if os.path.exists(extractStackSpecDir) == False:
        os.makedirs(extractStackSpecDir)
        
    stackExtractSpecDir=outDir+os.path.sep+"1DSpec_2DSpec_stackAndExtract"
    if iterativeMethod == True:
        stackExtractSpecDir=stackExtractSpecDir+"_iterative"
    if os.path.exists(stackExtractSpecDir) == False:
        os.makedirs(stackExtractSpecDir)    

    # Log checks of wavelength calibration
    skyWavelengthCalibCheckList=[]
        
    # Get list of extensions
    cutArcPath=maskDict['cutArcDict'][maskDict['OBJECT'][0]]                
    img=pyfits.open(cutArcPath)
    if extensionsList == "all":
        extensionsList=[]
        for hdu in img:
            if "SLIT" in hdu.name:
                extensionsList.append(hdu.name)
    
    # The way we stack... identify signal dominated rows and average them to a 1d spectrum, then stack all 1d
    # Do same for sky rows
    # They may all be projected onto slightly different wavelength coordinates system... have to deal with that also
    toStackList=makeOutputFileNameList(maskDict['OBJECT'], 'rwc', outDir)
    logger.info("Extracting and stacking")
    for extension in extensionsList:
        logger.info("%s" % (extension))
        signalList=[]
        skyList=[]
        wavelengthsList=[]
        headersList=[]
        CRMaskedDataCube=[]
        for fileName in toStackList:
            img=pyfits.open(fileName)
            foundExtension=False
            try:
                data=img[extension].data
                testCDELT1=img[extension].header['CDELT1']
                foundExtension=True
            except KeyError:
                logger.warning("skipping %s in %s" % (extension, fileName))
                foundExtension=False
            
            if foundExtension == True:
                
                logger.info("%s" % (fileName))
                
                header=img[extension].header

                # Extract calibrated wavelength scale, assuming left most pixel corresponds to CRVAL1
                # NOTE: if this fails due to missing CDELT1, it probably means the slit we found is a wacky shape - check the .reg file                   
                w=np.arange(data.shape[1])*header['CDELT1']+header['CRVAL1']
                if w[0] != header['CRVAL1']:
                    raise Exception("wavelength of pixel 0 doesn't correspond to CRVAL1 - what happened?")
                
                # Extract signal, sky and CR-flagged 2d spectrum data
                # If blank slit (which it would be if we skipped over something failing earlier), insert blank row
                if np.nonzero(data)[0].shape[0] > 0:
                    
                    if iterativeMethod == True:
                        signal, sky, mData=iterativeWeightedExtraction(data, subFrac = subFrac)
                    else:
                        t0=time.time()
                        signal, sky, mData=weightedExtraction(data)
                        t1=time.time()
                    CRMaskedDataCube.append(mData)
                    signalList.append(signal)
                    skyList.append(sky)
                    wavelengthsList.append(w)
                    headersList.append(header)
                    
                    # Diagnostic plot of smoothed spectrum
                    plt.figure()
                    plt.plot(w, ndimage.uniform_filter1d(signal, 15))
                    plt.xlabel("wavelength")
                    plt.ylabel("relative flux")
                    plotFileName=diagnosticsDir+os.path.sep+os.path.split(fileName)[-1].replace(".fits", "_%s_signal.png" % (extension))
                    plt.savefig(plotFileName)
                    plt.close()
                    
                else:
                    logger.warning("empty slit")
                    #signalList.append(np.zeros(data.shape[1]))
                    #skyList.append(np.zeros(data.shape[1]))
        
        CRMaskedDataCube=np.ma.masked_array(CRMaskedDataCube)
        
        if len(signalList) == 0:
            logger.warning("failed to construct CRMaskedDataCube - signalList empty - skipping %s" % (extension))
            continue
        
        # Wavelength calib diagnostic: does the sky line up from each science frame?
        fluxMax=0
        for sky, wavelengths, header in zip(skyList, wavelengthsList, headersList):
            if sky.shape[0] > 0:
                flux=sky/np.median(sky)
                if flux.max() > fluxMax:
                    fluxMax=flux.max()
                plt.plot(wavelengths, flux)
                medianOffset, numLines=checkWavelengthCalibUsingSky(sky, wavelengths, featureMinPix = 5)
                logger.info("individual frame - sky wavelength calib check: medianOffset = %.3f Angstroms, numLines = %d" %(medianOffset, numLines))
        for l in checkSkyLines:
            plt.plot([l]*10, np.linspace(0, fluxMax, 10), 'k--')
        plt.ylim(0, fluxMax)
        plt.xlabel("Wavelength (Angstroms)")
        plt.ylabel("Relative Flux")
        plt.savefig(diagnosticsDir+os.path.sep+"skyCheck_"+extension+".png")
        plt.close()
                    
        # Make stacked spectrum - interpolate onto common wavelength scale, then take median
        # We could make this fancier (noise weighting etc.)...
        signalArr=np.array(signalList)
        skyArr=np.array(skyList)
        wavelengthsArr=np.array(wavelengthsList)
        regrid_wavelengths=np.median(wavelengthsArr, axis = 0)
        regrid_signalArr=np.zeros(signalArr.shape)
        regrid_skyArr=np.zeros(skyArr.shape)
        for i in range(signalArr.shape[0]):
            tck=interpolate.splrep(wavelengthsArr[i], signalArr[i])
            regrid_signalArr[i]=interpolate.splev(regrid_wavelengths, tck, ext = 1)
            tck=interpolate.splrep(wavelengthsArr[i], skyArr[i])
            regrid_skyArr[i]=interpolate.splev(regrid_wavelengths, tck, ext = 1) 
        signal=np.median(regrid_signalArr, axis = 0)
        sky=np.median(regrid_skyArr, axis = 0)
        try:
            if sky.shape[0] > 0:
                medianOffset, numLines=checkWavelengthCalibUsingSky(sky, regrid_wavelengths, featureMinPix = 5)
                logger.info("extractAndStack - sky wavelength calib check: medianOffset = %.3f Angstroms, numLines = %d" % (medianOffset, numLines))
                outFileName=extractStackSpecDir+os.path.sep+"1D_"+maskDict['objName'].replace(" ", "_")+"_"+maskDict['maskID']+"_"+extension+".fits"
                write1DSpectrum(signal, sky, regrid_wavelengths, outFileName, maskDict['RA'], maskDict['DEC'])
        except:
            print("sky has no shape")
            IPython.embed()
            sys.exit()
        
        # 2d combined/extracted, projecting CR-masked arrays to same wavelength grid
        # NOTE: These all have different wavelength calibrations stored in wavelengthsList
        # NOTE: This assumes that the wavelength calibrations have similar accuracy... if [0] has a problem, then we do!
        projDataCube=[]
        refWavelengths=wavelengthsList[0]
        refHeader=headersList[0]
        for data, wavelengths in zip(CRMaskedDataCube, wavelengthsList):
            if projDataCube == []:
                projDataCube.append(data)
            else:
                projData=np.zeros(projDataCube[0].shape)
                projMask=np.zeros(projDataCube[0].shape)
                for i in range(min([projData.shape[0], data.shape[0]])):
                    tck=interpolate.splrep(wavelengths, data.data[i])
                    projData[i]=interpolate.splev(refWavelengths, tck, ext = 1)
                    if len(CRMaskedDataCube.mask.shape) > 0:    # sometimes data.mask[i] is simply False, and not a mask
                        tck=interpolate.splrep(wavelengths, data.mask[i])
                        projMask[i]=interpolate.splev(refWavelengths, tck, ext = 1)
                        projMask[i][np.less(projMask[i], 0.1)]=0.
                        projMask[i]=np.array(projMask[i], dtype = bool)
                projDataCube.append(np.ma.masked_array(projData, projMask))
        projDataCube=np.ma.masked_array(projDataCube)       
        # Very slow if use the masked array below...
        med=np.array(np.ma.median(projDataCube, axis = 0))   
        if iterativeMethod == True:
            signal, sky, mData=iterativeWeightedExtraction(med, subFrac = subFrac)
        else:
            signal, sky, mData=weightedExtraction(med)       
        # Put chipGapMask into 1d spectra
        chipGapMask=np.median(makeChipGapMask(med), axis = 0)
        outFileName=stackExtractSpecDir+os.path.sep+"1D_"+maskDict['objName'].replace(" ", "_")+"_"+maskDict['maskID']+"_"+extension+".fits"
        write1DSpectrum(signal, sky, refWavelengths, outFileName, maskDict['RA'], maskDict['DEC'],
                        mask = chipGapMask)
        
        # Experimenting with a method that will handle running profile
        #t0=time.time()
        #signal, sky, skySubbed2d=finalExtraction(med, subFrac = subFrac)
        #t1=time.time()
        #print "... final extraction (took %.3f sec) ..." % (t1-t0)
        #outFileName=stackExtractSpecDir+os.path.sep+"1D_"+maskDict['objName'].replace(" ", "_")+"_"+maskDict['maskID']+"_"+extension+"_testFinal.fits"
        #write1DSpectrum(signal, sky, refWavelengths, outFileName, maskDict['RA'], maskDict['DEC'])        
        #print "final extract again"
        #IPython.embed()
        #sys.exit()
        
        # Write 2d combined spectrum
        outFileName=stackExtractSpecDir+os.path.sep+"2D_"+maskDict['objName'].replace(" ", "_")+"_"+maskDict['maskID']+"_"+extension+".fits"
        newImg=pyfits.HDUList()
        hdu=pyfits.PrimaryHDU(None, refHeader)   
        hdu.data=np.array(med)
        newImg.append(hdu)
        if os.path.exists(outFileName) == True:
            os.remove(outFileName)
        newImg.writeto(outFileName, overwrite = True)
        newImg.close()
        
        # Quantify wavelength calibration accuracy using sky
        if sky.sum() > 0:
            medianOffset, numLines=checkWavelengthCalibUsingSky(sky, refWavelengths, featureMinPix = 5)
            logger.info("stackAndExtract - sky wavelength calib check: medianOffset = %.3f Angstroms, numLines = %d" % (medianOffset, numLines))
            skyWavelengthCalibCheckList.append([extension, medianOffset, numLines])
        
        # Write 2d, sky subtracted combined spectrum
        #outFileName=stackExtractSpecDir+os.path.sep+"2D_noSky_"+maskDict['objName'].replace(" ", "_")+"_"+maskDict['maskID']+"_"+extension+".fits"
        #newImg=pyfits.HDUList()
        #hdu=pyfits.PrimaryHDU(None, refHeader)   
        #hdu.data=np.array(skySubbed2d)
        #newImg.append(hdu)
        #newImg.writeto(outFileName, clobber = True)
        #newImg.close()

    # Write table of sky wavelength calib test results
    logFile=open(diagnosticsDir+os.path.sep+"skyWavelengthCalibCheck.csv", "w")
    logFile.write("#extension\tmedianOffsetAngstroms\tnumLines\n")
    offsets=[]
    for row in skyWavelengthCalibCheckList:
        extension=row[0]
        medianOffset=row[1]
        numLines=row[2]
        logFile.write("%s\t%.3f\t%d\n" % (extension, medianOffset, numLines))
        offsets.append(medianOffset)
    offsets=np.array(offsets)
    RMS=np.sqrt(np.mean(offsets**2))
    logFile.write("all slits median\t%.3f\t%d\n" % (np.median(offsets), len(offsets)))
    logFile.write("all slits RMS\t%.3f\t%d\n" % (RMS, len(offsets)))
    logFile.close()
    
#-------------------------------------------------------------------------------------------------------------
def maskNoisyData(data, sigmaCut = 3.0):
    """Return a mask with zeros for data with values < sigma*sigmaCut selected. 
    
    """

    mean=0
    sigma=1e6
    for i in range(20):
        lastSigma=sigma
        nonZeroMask=np.not_equal(data, 0)
        mask=np.less(abs(data-mean), sigmaCut*sigma)
        mask=np.logical_and(nonZeroMask, mask)
        mean=np.mean(data[mask])
        sigma=np.std(data[mask])
        if sigma == lastSigma:
            break            
    mask=np.greater(abs(data-mean), sigmaCut*sigma)

    return mask
    
#-------------------------------------------------------------------------------------------------------------
def fitProfile(data, mask, borderPix = 4):
    """Fit the profile in the y-direction, with a truncated Gaussian. We downweight pixels towards the edge
    of the slit, because sometimes there is junk there (from cutting it). This means we might miss some
    objects where the trace is towards the edge of the slit.
    
    mask is used to mask noise (or at least chip gaps)
    
    Returns the centre of the Gaussian profile, and its sigma. If everything was masked, returns
    [-99, -99], so these can easily be masked out.
    
    """
    
    # Prior - upweight pixels near centre, for finding the centre of the object trace
    prior=np.zeros(data.shape[0])
    prior[borderPix:-borderPix]=1
    prior=ndimage.uniform_filter1d(prior, borderPix)
    
    # First, flag noisy pixels
    noiseMask=maskNoisyData(data)
    mData=np.ma.masked_array(data, np.logical_or(mask, noiseMask))#NoisyData(data, sigmaCut = 2.0)    
    
    prof=np.ma.median(mData, axis = 1)
    
    if prof.mask.sum() == prof.shape[0]:
        # All is masked
        return -99, -99
    
    prof=prof*prior
    prof=prof/prof.max()
    x0=np.where(prof == prof.max())[0][0]
    x=np.arange(prof.shape[0])
    sigmaRange=np.linspace(1, int(prof.shape[0]/2), prof.shape[0]*2)
    resArr=np.zeros(sigmaRange.shape[0])
    for i in range(len(sigmaRange)):
        s=sigmaRange[i]
        gauss=np.exp(-((x-x0)**2)/(2*s**2))
        resArr[i]=((prof-gauss)**2).sum()
        #plt.plot(gauss)
    sigma=sigmaRange[np.where(resArr == resArr.min())]/2
    #fittedProf=np.exp(-((x-x0)**2)/(2*sigma**2))
    
    return x0, sigma
    
#-------------------------------------------------------------------------------------------------------------
def finalExtraction(data, subFrac = 0.8):
    """This fits for the object trace, so it can vary in y-position along the slit. Does a very simple
    extraction, which still appears to be working better for sky subtraction than the fancier attempts.
    
    Returns signal, sky, and 2d sky-subtracted spectrum 
    
    """

    # Throw away rows at edges as these often contain noise
    #throwAwayRows=2
    #data=data[throwAwayRows:-throwAwayRows]
    
    # Find the chip gaps and make a mask
    lowMaskValue=2.0
    minPix=1000
    chipGapMask=np.array(np.less(data, lowMaskValue), dtype = float)  # flags chip gaps as noise
    segmentationMap, numObjects=ndimage.label(chipGapMask)
    sigPixMask=np.equal(chipGapMask, 1)
    objIDs=np.unique(segmentationMap)
    objNumPix=ndimage.sum(sigPixMask, labels = segmentationMap, index = objIDs)
    for objID, nPix in zip(objIDs, objNumPix):    
        if nPix < minPix:
            chipGapMask[np.equal(segmentationMap, objID)]=0.0
    
    # Initial mask
    wn2d=np.zeros(data.shape)+chipGapMask               # Treat chip gaps as noise
        
    # Running profile measurement
    profCentres=np.zeros(data.shape[1])
    profSigmas=np.zeros(data.shape[1])
    for i in range(data.shape[1]):
        runWidth=100
        iMin=i-runWidth
        iMax=i+runWidth
        if iMin < 0:
            iMin=0
        if iMax > data.shape[1]-1:
            iMax=data.shape[1]-1
        profCentres[i], profSigmas[i]=fitProfile(data[:, iMin:iMax], wn2d[:, iMin:iMax])

    # Fit for trace centre, just use median for trace width sigma (doesn't vary by that much)
    x=np.arange(data.shape[1])
    mask=np.greater(profCentres, 0) # fitProfile returns -99 for completely masked data
    result=stats.linregress(x[mask], profCentres[mask])
    traceCentre=x*result.slope+result.intercept
    traceSigma=np.median(profSigmas[mask])
    
    # Make 2d running profile
    runningProf=np.zeros(data.shape)
    x=np.arange(data.shape[0])
    for i in range(data.shape[1]):
        runningProf[:, i]=np.exp(-((x-traceCentre[i])**2)/(2*traceSigma**2))

    # New ---
    # Iterative sky subtraction, with the object trace we found
    signal, sky, mData=iterativeWeightedExtraction(data, subFrac = subFrac, runningProfile = runningProf, iterateProfile = False, throwAwayRows = 0)

    # Make 2d sky-subtracted spectrum
    sky2d=np.array([sky]*data.shape[0])
    skySubbed2d=data-sky2d
    
    # Old ---
    # Fancy weighted extraction, assuming we have no noise now (masked CRs and replaced with sky earlier)
    #signal=np.zeros(data.shape[1])
    #sky=np.zeros(data.shape[1])
    #for i in range(data.shape[1]):
        #A=np.mat([runningProf[:, i], np.ones(data.shape[0])]).transpose()
        #b=data[:, i]
        ##X=np.dot(np.linalg.pinv(A), b)
        ##signal[i]=X[0, 0]
        ##sky[i]=X[0, 1]
        #X, R=optimize.nnls(A, b)
        #signal[i]=X[0]
        #sky[i]=X[1]
    # Or... this simple method seems to work better on sky subtraction
    # But won't work if signal dominated?
    #sky=np.median(data, axis = 0)   # This won't work if signal dominated.
    
    # Final answer...
    #sky2d=np.array([sky]*data.shape[0])
    #skySubbed2d=data-sky2d
    #signal=np.average(skySubbed2d, weights = runningProf, axis = 0)

    ## This just makes things work better in visualTemplateRedshift
    #signal[np.less(signal, 0)]=0. 
       
    #---

    # Zap the chip gaps (again)
    #chipGapMask1d=np.greater(np.median(chipGapMask, axis = 0), 0)
    #signal[chipGapMask1d]=0.
    #sky[chipGapMask1d]=0.
    #skySubbed2d[np.where(chipGapMask == 1)]=0.
    
    return signal, sky, skySubbed2d

#-------------------------------------------------------------------------------------------------------------
def write1DSpectrum(signal, sky, wavelength, outFileName, maskRA, maskDec, mask = None):
    """Writes 1D spectrum to .fits table file.
    
    """
    
    specColumn=pyfits.Column(name='SPEC', format='D', array=signal)
    skyColumn=pyfits.Column(name='SKYSPEC', format='D', array=sky)
    lambdaColumn=pyfits.Column(name='LAMBDA', format='D', array=wavelength)
    cols=[specColumn, skyColumn, lambdaColumn]
    if np.any(mask) != None:
        maskColumn=pyfits.Column(name='MASK', format='D', array=mask)
        cols.append(maskColumn)
    tabHDU=pyfits.BinTableHDU.from_columns(cols)
    tabHDU.name='1D_SPECTRUM'
    HDUList=pyfits.HDUList([pyfits.PrimaryHDU(), tabHDU])
    HDUList[0].header['MASKRA']=maskRA
    HDUList[0].header['MASKDEC']=maskDec
    if os.path.exists(outFileName) == True:
        os.remove(outFileName)
    HDUList.writeto(outFileName, overwrite = True)   
    

    
    
