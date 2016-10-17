#!/usr/bin/env python

"""Pipeline for reducing SALT RSS MOS data, using the stuff that comes in the product/ dir.

"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pyfits
import glob
import time
import datetime
import atpy
import argparse
#from astLib import *
from scipy import interpolate
from scipy import ndimage
from scipy import optimize
from scipy import stats
import pickle
import IPython
#plt.matplotlib.interactive(True)

#-------------------------------------------------------------------------------------------------------------
LOGFILE=None
REF_MODEL_DIR="modelArcSpectra"

#-------------------------------------------------------------------------------------------------------------
def trace(message, verbose = True, logFile = None):
    """Simple trace function. Prints to screen if verbose == True, prints to logFile if logFile is not None.
    Automatically adds newline characters to end of message.
    
    """  
    if verbose == True:
        print message
    if logFile != None:
        logFile.write(message+"\n")
        
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
    
    listFile=file(fileName, "w")
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
            newImg.writeto(rootOutFileName.replace(".fits", "_%d.fits" % (i)), clobber = True)    

#-------------------------------------------------------------------------------------------------------------
def getImageInfo(rawDir):
    """Sorts through all .fits files, making lists of biases, flats, science frames, arc frames, for each 
    night of observations, grating, position combo. Returns a dictionary of lists.
    
    NOTE: Explicitly avoiding standard stars for now 
    
    """
 
    # Now organised by objectName and meta info -> flats, arcs, object frames
    # This is so we can use one pipeline for MOS and longslit
    infoDict={}
    files=glob.glob(rawDir+os.path.sep+"mbxgp*.fits")  # Files either GMOS N or S depending on mask name
    
    # First, get object names
    for f in files:
        img=pyfits.open(f)
        header=img[0].header
        if header['OBSMODE'] == 'SPECTROSCOPY':
            obsType=header['CCDTYPE']
            maskID=header['MASKID']
            if obsType == 'OBJECT':
                maskName=header['OBJECT']+"_"+maskID
                infoDict[maskName]={}
                infoDict[maskName][maskID]={}
                infoDict[maskName]['maskID']=maskID             # Clunky, but convenient
                infoDict[maskName]['maskType']=header['MASKTYP']  
                infoDict[maskName]['objName']=header['OBJECT'].replace("'", "").replace('"', "")  

    # Now add flats etc.
    for maskName in infoDict.keys():
        
        for f in files:
            img=pyfits.open(f)
            header=img[0].header
            if header['OBSMODE'] == 'SPECTROSCOPY':# and header['MASKTYP'] == 'MOS':
    
                dateObs=header['DATE-OBS']
                timeObs=header['TIME-OBS']
                maskID=header['MASKID']
                obsType=header['CCDTYPE']
                objName=header['OBJECT'].replace("'", "").replace('"', "")  

                if maskID == infoDict[maskName]['maskID']:
                                    
                    # NOTE: we could add some checks for e.g. grating, camera station etc. here
                    # i.e., match flats, arcs to object frames here
                    if obsType not in infoDict[maskName][maskID].keys():
                        infoDict[maskName][maskID][obsType]=[]
                    if obsType != 'OBJECT':
                        infoDict[maskName][maskID][obsType].append(f)
                    elif obsType == 'OBJECT' and objName == infoDict[maskName]['objName']:
                        infoDict[maskName][maskID][obsType].append(f)
                        # Just so we can track this later in output 1d spectra
                        infoDict[maskName][maskID]['RA']=header['RA']
                        infoDict[maskName][maskID]['DEC']=header['DEC']
                    
    return infoDict

#-------------------------------------------------------------------------------------------------------------
def makeMasterFlats(maskDict, outDir, deltaHours = 0.5):
    """Make master flats from files in 'FLAT' key of maskDict. Automatically group flats taken within
    deltaHours. Adds paths to 'masterFlat' key.
    
    """
    
    flatLists=groupFilesListByTime(maskDict['FLAT'])
        
    maskDict['masterFlats']=[]

    print ">>> Making masterFlats (it is a good idea to sin bin any that aren't aligned with object spectra at cutting stage)"
    
    for i in range(len(flatLists)):
        flatFiles=flatLists[i]
        masterFlatPath=outDir+os.path.sep+"masterFlat_%d.fits" % (i)
        print "... making %s (%s) ..." % (masterFlatPath, flatFiles)
        if os.path.exists(masterFlatPath) == False:
            flatCube=[]
            for f in flatFiles:
                img=pyfits.open(f)
                flatCube.append(img['SCI'].data)
            flatCube=np.array(flatCube)
            flatData=np.median(flatCube, axis = 0)
            img['SCI'].data=flatData
            img.writeto(masterFlatPath, clobber = True)
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
    centreColumn=img['SCI'].header['NAXIS1']/2
    img.close()
    outFile=file(regFileName, "w")
    outFile.write("# DS9 region file\n")
    outFile.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    outFile.write("image\n")
    for key in slitsDict.keys():
        outFile.write("point(%.3f,%.3f) # point=boxcircle text={SLIT%d}\n" % (centreColumn, (slitsDict[key]['yMax']+slitsDict[key]['yMin'])/2.0, key))            
        outFile.write("box(%.1f,%.1f,%.1f,%.1f)\n" % (centreColumn, (slitsDict[key]['yMax']+slitsDict[key]['yMin'])/2.0+1, centreColumn*2.0, (slitsDict[key]['yMax']-slitsDict[key]['yMin'])))
    outFile.close()        
    
#-------------------------------------------------------------------------------------------------------------
def cutIntoSlitLets(maskDict, outDir, threshold = 0.1):
    """Cuts files into slitlets, making MEF files. 
            
    threshold is the parameter used by findSlits
    
    """

    # We find slits using master flats, and store slit locations in a dictionary indexed by masterFlatPath
    # We have a routine to pull out a corresponding master flat (and hence slitsDict) for every image.
    maskDict['slitsDicts']={}
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
        ref[refDict[key]['yCentre']]=1
    shiftsDict={}
    for key in maskDict['slitsDicts'].keys():
        if key != maskDict['masterFlats'][0]:
            slitsDict=maskDict['slitsDicts'][key]
            g=np.zeros(height)
            for skey in slitsDict.keys():
                g[slitsDict[skey]['yCentre']]=1
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
    print ">>> It is a good idea to check that for the below the corresponding DS9 .reg file aligns with the slits and that object spectra are actually centred in the slits..."
    maskDict['cutFlatDict']={}
    maskDict['cutArcDict']={}
    for f, outFileName in zip(toCutList, outCutList):
        flatFileName=findMatchingFileByTime(f, maskDict['masterFlats'])
        slitsDict=maskDict['slitsDicts'][flatFileName]
        print "... cutting %s (and arcs, flats) using %s for slit definition ..." % (f, flatFileName)
        label=os.path.split(flatFileName)[-1].replace(".fits", "")
        # Object
        cutSlits(f, outFileName, slitsDict)
        # Arc
        arcFileName=findMatchingFileByTime(f, maskDict['ARC'])
        cutArcFileName=makeOutputFileName(arcFileName, "c"+label, outDir)
        cutSlits(arcFileName, cutArcFileName, slitsDict)
        maskDict['cutArcDict'][f]=cutArcFileName
        # Flat
        cutFlatFileName=makeOutputFileName(flatFileName, "c"+label, outDir)
        cutSlits(flatFileName, cutFlatFileName, slitsDict)
        maskDict['cutFlatDict'][f]=cutFlatFileName
        # DS9 regions
        regFileName=flatFileName.replace(".fits", "_slitLocations.reg")
        writeDS9SlitRegions(regFileName, slitsDict, f)
        
#-------------------------------------------------------------------------------------------------------------
def cutIntoPseudoSlitLets(maskDict, outDir):
    """For longslit data. Finds objects, and then cuts into pseudo-slitlets: we take some region +/- Y pixels
    around the object trace and pretend that is a MOS slitlet. Outputs MEF files.
            
    """
        
    # Find object traces in OBJECT frames and cut +/- some distance in Y around them
    maskDict['slitsDicts']={}
    for i in range(len(maskDict['OBJECT'])):
        objPath=maskDict['OBJECT'][i]
        slitsDict=findPseudoSlits(objPath)
        maskDict['slitsDicts'][objPath]=slitsDict
    
    # There can be significant offsets between object traces in longslit frames
    # So... compare all the pseudo-slits we assigned and measure a y-offset from the i=0 dict
    # This will be used by cutSlits, if present
    refDict=maskDict['slitsDicts'][maskDict['OBJECT'][0]]
    img=pyfits.open(maskDict['OBJECT'][0])
    height=img[1].data.shape[0]
    ref=np.zeros(height)
    for key in refDict.keys():
        ref[refDict[key]['yCentre']]=1
    shiftsDict={}
    for key in maskDict['slitsDicts'].keys():
        if key != maskDict['OBJECT'][0]:
            slitsDict=maskDict['slitsDicts'][key]
            g=np.zeros(height)
            for skey in slitsDict.keys():
                g[slitsDict[skey]['yCentre']]=1
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
    print ">>> It is a good idea to check from the corresponding DS9 .reg file that object spectra are actually centred in the pseudo-slits..."
    maskDict['cutFlatDict']={}
    maskDict['cutArcDict']={}
    for f, outFileName in zip(toCutList, outCutList):
        print "... cutting %s (and arcs, flats) using %s for slit definition ..." % (f, f)
        slitsDict=maskDict['slitsDicts'][f]
        label=os.path.split(f)[-1].replace(".fits", "")
        # Object
        cutSlits(f, outFileName, slitsDict)
        # Arc
        arcFileName=findMatchingFileByTime(f, maskDict['ARC'])
        cutArcFileName=makeOutputFileName(arcFileName, "c"+label, outDir)
        cutSlits(arcFileName, cutArcFileName, slitsDict)
        maskDict['cutArcDict'][f]=cutArcFileName
        # Flat
        flatFileName=findMatchingFileByTime(f, maskDict['masterFlats'])
        cutFlatFileName=makeOutputFileName(flatFileName, "c"+label, outDir)
        cutSlits(flatFileName, cutFlatFileName, slitsDict)
        maskDict['cutFlatDict'][f]=cutFlatFileName
    
        # Write out a .reg file so we can match slits to objects
        img=pyfits.open(f)
        centreColumn=img['SCI'].header['NAXIS1']/2
        img.close()
        regFileName=outFileName.replace(".fits", "_slitLocations.reg")
        outFile=file(regFileName, "w")
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
    d=d/a
    d[np.isnan(d)]=0.0

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
    
    # Debugging
    #print "Check slitsDict"
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
def findPseudoSlits(objFileName, skyRows = 20, minSlitHeight = 10., thresholdSigma = 3., minTraceWidth = 5):
    """Finds object traces in longslit data, defines regions +/- skyRows around them, so we can treat in 
    the same way as MOS slitlets.
    
    objects are detected as peaks in the SNR profile across the slit. Use minTraceWidth to set the number
    of pixels in the SNR profile that must be above thresholdSigma for an object to be detected.
    
    Returns a dictionary which can be fed into cutSlits
    
    """ 
    
    img=pyfits.open(objFileName)
    d=img['SCI'].data
    
    # Take out spectrum of flat lamp (approx)
    a=np.median(d, axis = 0)
    d=d/a
    d[np.isnan(d)]=0.0
        
    # Find local background, noise (running clipped mean)
    prof=np.median(d, axis = 1)    
    prof[np.less(prof, 0)]=0.
    halfBlkSize=50
    sigmaCut=3.0        
    bck=np.zeros(prof.shape)
    sig=np.zeros(prof.shape)
    for y in range(prof.shape[0]):
        yMin=y-halfBlkSize
        yMax=y+halfBlkSize
        if yMin < 0:
            yMin=0
        if yMax > prof.shape[0]-1:
            yMax=prof.shape[0]-1
        mean=0
        sigma=1e6
        for i in range(20):
            nonZeroMask=np.not_equal(prof[yMin:yMax], 0)
            mask=np.less(abs(prof[yMin:yMax]-mean), sigmaCut*sigma)
            mean=np.mean(prof[yMin:yMax][mask])
            sigma=np.std(prof[yMin:yMax][mask])            
        bck[y]=mean
        sig[y]=sigma
    
    # Detect peaks
    profSNR=(prof-bck)/sig    
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
            
    return slitsDict

#-------------------------------------------------------------------------------------------------------------
def applyFlatField(maskDict, outDir):
    """Applies the flat field correction. Let's do this in place...
    
    """

    print ">>> Applying flat field..."
    
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
            
            # Find chip gaps
            threshold=200
            grad=np.gradient(med)
            plusMask=np.greater(grad, threshold)
            minusMask=np.less(grad, -1*threshold) 
            gapsDict={}
            lookingFor=0
            xMin=None
            xMax=None
            gapsCount=0
            for i in range(len(plusMask)):
                if lookingFor == 1:
                    if plusMask[i] == True:
                        xMax=i
                        lookingFor=0
                if lookingFor == 0:
                    if minusMask[i] == True:
                        xMin=i
                        lookingFor=1
                if xMin != None and xMax != None:
                    gapsCount=gapsCount+1
                    gapsDict[gapsCount]={'xMin': xMin-1, 'xMax': xMax+1}    
                    xMin=None
                    xMax=None
            gapsMask=np.ones(len(med), dtype = bool)
            for key in gapsDict.keys():
                gapsMask[gapsDict[key]['xMin']:gapsDict[key]['xMax']]=False
            
            # Fit and remove spectrum of flat lamp
            x=np.arange(len(med))
            poly=np.poly1d(np.polyfit(x[gapsMask], med[gapsMask], 10))
            mod=np.array([poly(x)]*data.shape[0])
            flatfield=flatfield/mod
            #plt.plot(x, med)
            #plt.plot(x, poly(x))
            #IPython.embed()
            #sys.exit()
            
            data=data/flatfield
            data[np.isnan(data)]=0.0
            img[extension].data=data
        
        img.writeto(f, clobber = True)

#-------------------------------------------------------------------------------------------------------------
def findMatchingFileByTime(inputFileName, possibleFilesList):
    """Identify the file name in possibleFilesList that is closest to inputFileName in terms of time.
    Use to find corresponding arcs, flats.
    
    Returns fileName
    
    """
        
    ctimes=[]
    for f in possibleFilesList:
        ctime=getCTimeFromHeader(f)
        ctimes.append(ctime)
    ctimes=np.array(ctimes)
    
    fileCTime=getCTimeFromHeader(inputFileName)
    bestMatchIndex=np.where(abs(fileCTime-ctimes) == abs(fileCTime-ctimes).min())[0][0]
    
    return possibleFilesList[bestMatchIndex]

#-------------------------------------------------------------------------------------------------------------
def flattenAndRectify(maskDict, outDir, subtractSky = False):
    """Apply the appropriate flat field and rectification transform to all arc and science frames.
    
    """
        
    # Flatten arcs and science frames, applying rectification as we go
    # We need to work out the appropriate transform to use on a file-by-file basis
    # Do this by matching to nearest arc in time
    toFlatList=maskDict['ARC']+maskDict['OBJECT']
    toFlatList=makeOutputFileNameList(toFlatList, "c", outDir)
    outFlattenedList=makeOutputFileNameList(toFlatList, "f", outDir)
    for f, outFileName in zip(toFlatList, outFlattenedList):

        # Find nearest masterFlat in time
        masterFlatPath=findMatchingFileByTime(f, maskDict['masterFlats'])
        cutMasterFlatPath=makeOutputFileName(masterFlatPath, "c", outDir)
        try:
            masterFlatImg=pyfits.open(cutMasterFlatPath)
        except:
            print "Hmm - can't find cutMasterFlatPath"
            IPython.embed()
            sys.exit()
        
        # Find nearest arc in time
        arcPath=findMatchingFileByTime(f, maskDict['ARC'])        
        
        toFlatImg=pyfits.open(f)
        for i in range(len(toFlatImg)):
            if "SLIT" in toFlatImg[i].name:
                
                toFlatData=toFlatImg[i].data
                flatData=masterFlatImg[i].data
                
                # Use response instead of vanilla flat fielding
                # If we don't remove these files, we end up appending lots of extensions
                if os.path.exists("resp_in.fits") == True:
                    os.remove("resp_in.fits")
                if os.path.exists("resp_out.fits") == True:
                    os.remove("resp_out.fits")
                newImg=pyfits.HDUList()
                hdu=pyfits.PrimaryHDU(None, masterFlatImg[i].header)   
                hdu.header.update("DISPAXIS", 1)
                hdu.data=flatData
                newImg.append(hdu)
                newImg.writeto("resp_in.fits", clobber = True)  
                iraf.flpr() # give time to write resp_in.fits 
                longslit.response(calibration="resp_in.fits", normalization="resp_in.fits", 
                                    response="resp_out.fits", interactive="no")
                respData=pyfits.getdata("resp_out.fits")
                toFlatImg[i].data=(toFlatData/respData)#*(np.mean(flatData))
                toFlatImg[i].data=np.nan_to_num(toFlatImg[i].data)
                
                # Rectify
                xCorrTransform=maskDict['xCorrTransforms'][arcPath][toFlatImg[i].name]
                toFlatImg[i].data=applyXCorrRectification(toFlatImg[i].data, xCorrTransform)    
                
                # Subtract sky here on science frames - should make CR rejection easier
                # Take sky as everything below median + 1 sigma
                if toFlatImg[0].header['CCDTYPE'] == 'OBJECT' and subtractSky == True:
                    # Need to do something a little more sophisticated to handle chip gaps
                    gapMask=np.less(toFlatImg[i].data, 1.0)
                    gapMask=np.greater(ndimage.uniform_filter(np.array(gapMask, dtype = float), 5), 0.1) # dilation
                    prof=np.median(toFlatImg[i].data, axis = 1)
                    med=np.median(prof)
                    std=np.std(prof)
                    yMask=np.less(prof, med+std)
                    sky=np.median(toFlatImg[i].data[yMask], axis = 0)
                    skyImage=np.array([sky]*toFlatImg[i].data.shape[0])
                    toFlatImg[i].data=toFlatImg[i].data-skyImage
                    toFlatImg[i].data[gapMask]=0.0
                                    
        toFlatImg.writeto(outFileName, clobber = True)

#-------------------------------------------------------------------------------------------------------------
def detectLines(data, sigmaCut = 3.0, thresholdSigma = 2.0, featureMinPix = 30):
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
        raise Exception, "Wavelength range covered is too big!"
    fPadded=np.zeros(2**n)
    fPadded[fUp.shape[0]/2:fUp.shape[0]/2+fUp.shape[0]]=fUp[:]
    gPadded=np.zeros(2**n)
    gPadded[gUp.shape[0]/2:gUp.shape[0]/2+gUp.shape[0]]=gUp[:]
    
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
    
    shift=x[0]
    s=x[1]
    tck=interpolate.splrep((data_x+shift)+s*data_x, arcRow)
    arcRow_scaled=interpolate.splev(data_x, tck, ext = 1)
    arcMean=np.mean(arcRow_scaled)
    arcStd=np.std(arcRow_scaled)
    
    # Old style: minimize what we call overlap below - this behaves well with optimize.minimize
    # (overlap vs scale is a function with a clear minimum)
    overlap=np.trapz(abs(normRefModel[:data_x.shape[0]]-(arcRow_scaled-arcMean)/arcStd))
        
    return overlap

#-------------------------------------------------------------------------------------------------------------
def findScaleAndShift(arcRow, refModelDict, quick = False):
    """Find best fit stretch and scale to transform arcRow to reference model
    
    Use quick = True to skip the basinhopping brute force algorithm (needed for accurate wavelength 
    calibration, but not for the model selection step)
    
    """
    
    # Use cross correlation to get initial guess at shift between arc and model
    shift=np.argmax(np.correlate(refModelDict['arc_centreRow'], arcRow, mode = 'full'))-refModelDict['arc_centreRow'].shape[0]
    
    # Sometimes the refModel is shorter by a pixel than arcRow (we can handle the reverse case)
    # NOTE: added when adapting MOS pipeline to work on longslit
    if arcRow.shape[0] > refModelDict['arc_centreRow'].shape[0]:
        arcRow=arcRow[:refModelDict['arc_centreRow'].shape[0]]
    
    # New optimize based method (robust when used with overlap method, but not with xcorr)
    data_x=np.arange(0, arcRow.shape[0])        
    modelStd=np.std(refModelDict['arc_centreRow'])
    modelMean=np.mean(refModelDict['arc_centreRow'])
    normRefModel=(refModelDict['arc_centreRow']-modelMean)/modelStd
    
    result=optimize.minimize_scalar(minFunc_findScale, bounds = (-0.04, 0.04), method = 'Bounded', 
                                    args = (shift, arcRow, normRefModel, data_x))
    s=result['x']
    
    # xcorr is still best to use for selecting between models with corrMax
    tck=interpolate.splrep((data_x+shift)+s*data_x, arcRow)
    arcRow_scaled=interpolate.splev(data_x, tck, ext = 1)
    arcMean=np.mean(arcRow_scaled)
    arcStd=np.std(arcRow_scaled)
    corr, corrMax, extraShift=fftCorrelate(normRefModel, (arcRow_scaled-arcMean)/arcStd)    

    # This brute force approach is robust... optimize.minimize is not
    # And we checked... it is better to do it in this order!
    if quick == False:
        result=optimize.basinhopping(minFunc_findShiftAndScale, [shift, s], minimizer_kwargs = {'args': (arcRow, normRefModel, data_x), 'bounds': [[shift-100, shift+100], [-0.04, 0.04]]})
        shift=result['x'][0]
        s=result['x'][1]
   
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
        pickleFile=file(modelFileName, "rb")
        unpickler=pickle.Unpickler(pickleFile)
        refModelDict=unpickler.load()
        refModelDictsList.append(refModelDict)
        pickleFile.close()
        
        # First need to find arc features
        arcFeatureTable, arcSegMap=detectLines(arcData, thresholdSigma = thresholdSigma)
        arcFeatureTablesList.append(arcFeatureTable)
        arcSegMapsList.append(arcSegMap)
        
        # Replaced np.correlate with fft based correlation
        # Find shift and wavelength dependent scale change (stretch, then shift)
        arcRow=arcData[arcData.shape[0]/2]
        bestCorrMax, bestFitScale, bestFitShift=findScaleAndShift(arcRow, refModelDict, quick = True)
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
    
    print "... using refModel = %s ..." % (modelFileNameList[bestModelIndex])
    
    return refModelDict, arcFeatureTable, arcSegMap
    
#-------------------------------------------------------------------------------------------------------------
def findWavelengthCalibration(arcData, modelFileName, sigmaCut = 3.0, thresholdSigma = 2.0, 
                              featureMinPix = 50, order = 2, diagnosticsDir = None, diagnosticsLabel = None):
    """Find wavelength calibration for .fits image arcFileName containing a 2d arc spectrum.
    
    modelFileName is the path to a model made by makeModelArcSpectrum.
    
    Returns an array of polynomial fit coefficients that can be fed into wavelengthCalibrateAndRectify
    
    """

    # We now allow multiple reference models for each grating/lamp/binning config
    # This is useful if e.g., we have some MOS slit which is way to the red/blue end of the detector
    # First select the model to use based on centre row only (saves much time)
    # Choose best from maximum cross correlation
    modelFileNameList=glob.glob(modelFileName.split(".pickle")[0]+"*.pickle")
    refModelDict, arcFeatureTable, arcSegMap=selectBestRefModel(modelFileNameList, arcData, 
                                                                thresholdSigma = thresholdSigma)

    # Continue with previous 2d wavelength calib method
    yIndex=arcData.shape[0]/2
    arcRow=arcData[yIndex]
    bestCorrMax, bestFitScale, bestFitShift=findScaleAndShift(arcRow, refModelDict, quick = False)
    arc_centreRow=arcRow
    
    # Sanity check plot - this is key, if this is off, all is off
    #plt.matplotlib.interactive(True)
    data_x=np.arange(0, arc_centreRow.shape[0])        
    x=np.arange(0, len(arc_centreRow))
    arc_x_shifted=x*(1+bestFitScale)+bestFitShift
    plt.figure(figsize = (15, 10))
    plt.plot(x, refModelDict['arc_centreRow'][:data_x.shape[0]]/refModelDict['arc_centreRow'][:data_x.shape[0]].mean(), 'b-', label ='ref model')
    plt.plot(arc_x_shifted, arc_centreRow/arc_centreRow.mean(), 'r-', label = 'transformed arc')
    plt.semilogy()
    plt.legend()
    plt.savefig(diagnosticsDir+os.path.sep+"arcTransformTest_"+diagnosticsLabel+".png")
    plt.close()
    
    # New - based on J0018 skyCheck.png, the old method is better...
    ## Calibration from the reference model
    ## Then use bestFitScale, bestFitShift to transform
    ##plt.matplotlib.interactive(True)
    #transformed_model_x=(refModelDict['featureTable']['x_centreRow']-bestFitShift)/(1+bestFitScale)
    #fitCoeffs=np.polyfit(transformed_model_x, refModelDict['featureTable']['wavelength'], order)
    #wavelengthCalibPoly=np.poly1d(fitCoeffs)
    #wavelengths=wavelengthCalibPoly(transformed_model_x)
    #plt.plot(transformed_model_x, refModelDict['featureTable']['wavelength'], 'r.')
    #plt.plot(transformed_model_x, wavelengths, 'b-')
    #plt.xlabel("Transformed x (pix)")
    #plt.ylabel("Wavelength (Angstroms)")
    #plt.savefig(diagnosticsDir+os.path.sep+"transformedFeaturesFit_"+diagnosticsLabel+".png")
    #plt.close()

    ## Heck, let's just get best fit shifts over the slit
    ## Fix the scale to that found for centreRow
    #shifts=[]
    #data_x=np.arange(0, arcRow.shape[0])        
    #modelStd=np.std(refModelDict['arc_centreRow'])
    #modelMean=np.mean(refModelDict['arc_centreRow'])
    #normRefModel=(refModelDict['arc_centreRow']-modelMean)/modelStd
    #fitCoeffsArr=[]
    #for y in range(arcData.shape[0]):
        #t0=time.time()
        #result=optimize.minimize_scalar(minFunc_findShift, bounds = (bestFitShift-50, bestFitShift+50), 
                                        #method = 'Bounded', args = (bestFitScale, arcData[y], 
                                                                    #normRefModel, data_x))
        #shift=result['x']
        #transformed_model_x=(refModelDict['featureTable']['x_centreRow']-shift)/(1+bestFitScale)
        #fitCoeffs=np.polyfit(transformed_model_x, refModelDict['featureTable']['wavelength'], order)
        #fitCoeffsArr.append(fitCoeffs)
    #fitCoeffsArr=np.array(fitCoeffsArr)
        
    # Old ---    
    # Tag features by transforming model coords to arc coords and looking for closest match
    # Looking at above sanity plot, seems like the weak link here is centroiding done for 
    # x_centreRow in arcFeatureTable?
    arcFeatureTable.add_column('wavelength', np.zeros(len(arcFeatureTable)))    
    maxDistancePix=20.0
    for row in refModelDict['featureTable']:
        transformed_model_x=(row['x_centreRow']-bestFitShift)/(1+bestFitScale)
        dist=abs(arcFeatureTable['x_centreRow']-transformed_model_x)
        if dist.min() < maxDistancePix:
            index=np.argmin(dist)
            arcFeatureTable['wavelength'][index]=row['wavelength']
    arcFeatureTable=arcFeatureTable.where(arcFeatureTable['wavelength'] != 0)
    if len(arcFeatureTable) == 0:
        raise Exception, "No features identified in arc"    
        
    # Sanity check: tagged features
    #plt.plot(arcData[arcData.shape[0]/2], 'b-')
    #for row in arcFeatureTable:
        #plt.text(row['x_centreRow'], row['amplitude'], row['wavelength'])
    ##plt.close()
    #print "sanity check tag"
    #IPython.embed()
    #sys.exit()
    
    # Fit wavelength solution on centre row to check order of fit
    # We can't go to higher order if we get something nonsensical (double valued)
    # NOTE: 4th order is worse with spectra we've done before and worked ok...
    #order=4
    #orderOkay=False
    #while orderOkay == False:
        #fitCoeffs=np.polyfit(arcFeatureTable['x_centreRow'], arcFeatureTable['wavelength'], order)
        #poly=np.poly1d(fitCoeffs)
        #try:
            #tck=interpolate.splrep(poly(data_x), data_x)
            #orderOkay=True
        #except:
            #order=order-1
    #diff=arcFeatureTable['wavelength']-poly(arcFeatureTable['x_centreRow'])
    ##plt.plot(arcFeatureTable['x_centreRow'], diff, 'r.')
    #print "... order = %d, median residual = %.3f Angstroms" % (order, np.median(diff))
    #plt.plot(arcFeatureTable['x_centreRow'], arcFeatureTable['wavelength'], 'r.')
    #plt.plot(data_x, poly(data_x), 'k--')
    
    # Find 2d wavelength solution which we can use for rectification/wavelength calibration
    # Fit functions for how feature x positions change with y
    ys=np.arange(arcData.shape[0])
    arcFeatureTable.add_column('slope', np.zeros(len(arcFeatureTable)))
    arcFeatureTable.add_column('intercept', np.zeros(len(arcFeatureTable)))
    for row in arcFeatureTable:
        xs=np.zeros(arcData.shape[0])
        for i in range(len(ys)):
            objPositions=ndimage.center_of_mass(arcData[ys[i]], labels = arcSegMap, index = arcFeatureTable['id'])
            xs[i]=objPositions[np.where(arcFeatureTable['id'] == row['id'])[0]][0]
        # Linear fit should allow us to work out shear - here, given y, we want x
        # We probably don't want all this info (ys, xs), but keep for now
        # We could use polynomial instead of linear (see below)
        try:
            slope, intercept=np.polyfit(ys, xs, 1)
        except:
            print "polyfit failed"
            IPython.embed()
            sys.exit()
        row['slope']=slope
        row['intercept']=intercept

    # Wavelength calibration and model with arbitrary order polynomials - get coeffs for each row
    # We could potentially fit these coeffs as fn. of y - they are all well behaved
    # This array should be all we need to wavelength calibrate + rectify
    wavelengths=arcFeatureTable['wavelength']
    fitCoeffsArr=[]
    for y in range(arcData.shape[0]):
        xs=[]
        for row in arcFeatureTable:
            xs.append(row['slope']*y + row['intercept'])
        xs=np.array(xs)
        try:
            fitCoeffsArr.append(np.polyfit(xs, wavelengths, order))
        except:
            print "WARNING: wavelength calib failed"
            return None
    fitCoeffsArr=np.array(fitCoeffsArr)

    # Sanity check: wavelength calibration model with tagged features
    if diagnosticsLabel != None and diagnosticsDir != None:
        yCheck=arcData.shape[0]/2
        wavelengthCalibPoly=np.poly1d(fitCoeffsArr[yCheck])
        wavelengths=wavelengthCalibPoly(np.arange(arcData.shape[1]))
        plt.plot(wavelengths, arcData[arcData.shape[0]/2], 'r-')
        plt.plot(arcFeatureTable['wavelength'], arcFeatureTable['amplitude'], 'bo')
        for row in arcFeatureTable:
            plt.text(row['wavelength'], row['amplitude'], row['wavelength'])
        plt.semilogy()
        plt.savefig(diagnosticsDir+os.path.sep+"taggedFeatures_%s.png" % (diagnosticsLabel))
        plt.close()

    # Debugging particular slits
    #if diagnosticsLabel == 'SLIT15':
        #IPython.embed()
        #sys.exit()
    #---
    
    return fitCoeffsArr

#-------------------------------------------------------------------------------------------------------------
def wavelengthCalibrateAndRectify(inFileName, outFileName, wavelengthCalibDict, extensionsList = "all",
                                  makeDiagnosticPlots = False):
    """Applies the wavelength calibration, and rectification, to all extensions of inFileName, writing 
    output to outFileName. The wavelength calibration is provided in wavelengthCalibDict, where each key
    corresponds to each extension number (see findWavelengthCalibration)
    
    """
    
    print ">>> Applying wavelength calibration and rectifying (%s) ..." % (inFileName)
    img=pyfits.open(inFileName)
    if extensionsList == "all":
        extensionsList=[]
        for hdu in img:
            if "SLIT" in hdu.name:
                extensionsList.append(hdu.name)
                
    for extension in extensionsList:
    
        print "... %s ..." % (extension)
        
        data=img[extension].data
        header=img[extension].header
        fitCoeffsArr=wavelengthCalibDict[extension]
        
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
            wavelengths_centreRow=wavelengthsMap[wavelengthsMap.shape[0]/2]
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
                    print "WARNING: splrep error, this slit will be blank"
            img[extension].data=rectifiedData
            header['CTYPE1']='LINEAR'
            header['DISPAXIS']=1
            header['CRVAL1']=FITSRefLambda
            header['CRPIX1']=FITSRefPixel
            header['CD1_1']=FITSWavelengthScale
            header['CDELT1']=FITSWavelengthScale
            header['CUNIT1']='Angstroms'
                
            # Sanity check plot: linear wavelength scale
            if makeDiagnosticPlots == True:
                diagnosticsDir=os.path.split(outFileName)[0]+os.path.sep+"diagnostics"
                if os.path.exists(diagnosticsDir) == False:
                    os.makedirs(diagnosticsDir)
                plt.plot(rectWavelengthsMap[data.shape[0]/2], rectifiedData[data.shape[0]/2], 'k-')
                plt.xlabel("Wavelength (Angstroms)")
                plt.ylabel("Relative Intensity")
                plt.title("%s - %s" % (os.path.split(inFileName)[-1], extension))
                plt.savefig(diagnosticsDir+os.path.sep+"wavelengthCalibCheck_%s_%s.png" % (os.path.split(outFileName)[-1].replace(".fits", ""), extension))
                plt.close()

    # Write output
    img.writeto(outFileName, clobber = True)
    
#-------------------------------------------------------------------------------------------------------------
def wavelengthCalibration2d(maskDict, outDir, extensionsList = "all"):
    """Finds 2d wavelength calibration from arc frames, applies to arc frames and object frames, rectifying
    them and also interpolating to a linear wavelength scale to make life easier later.
    
    Should be fully automatic, assuming a suitable reference model is available.
    
    """
       
    diagnosticsDir=outDir+os.path.sep+"diagnostics"
    if os.path.exists(diagnosticsDir) == False:
        os.makedirs(diagnosticsDir)
        
    print ">>> Finding 2d wavelength solution ..."
    maskDict['wavelengthCalib']={}
    for key in maskDict['cutArcDict'].keys():
        
        cutArcPath=maskDict['cutArcDict'][key]
        
        print "--> arc = %s ..." % (cutArcPath)
              
        img=pyfits.open(cutArcPath)
        
        if extensionsList == "all":
            extensionsList=[]
            for hdu in img:
                if "SLIT" in hdu.name:
                    extensionsList.append(hdu.name)

        binning=img[0].header['CCDSUM'].replace(" ", "x")
        grating=img[0].header['GRATING']
        lampid=img[0].header['LAMPID']
        modelFileName=REF_MODEL_DIR+os.path.sep+"RefModel_"+grating+"_"+lampid+"_"+binning+".pickle"

        maskDict['wavelengthCalib'][cutArcPath]={}
        for extension in extensionsList:
            print "... extension = %s ..." % (extension)
            arcData=img[extension].data
            if os.path.exists(modelFileName) == False:
                print "No reference model exists for grating %s, lamp %s, with binning %s" % (grating, lampid, binning)
                print "Use createModelArcSpectrum.py script under modelArcSpectra dir"
                print "(arcFileName: %s)" % (arcFileName)
                sys.exit()
            maskDict['wavelengthCalib'][cutArcPath][extension]=findWavelengthCalibration(arcData, modelFileName, diagnosticsDir = diagnosticsDir, diagnosticsLabel = extension)
    
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
        print "profile measurement fail - no object detected?"
        print "add some sensible default profile for this case"
        IPython.embed()
        sys.exit()
    
    if np.any(np.isnan(prof)) == True:
        print "nans in object profile"
        IPython.embed()
        sys.exit()
        
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
def iterativeWeightedExtraction(data, maxIterations = 1000, subFrac = 0.8, runningProfile = None,
                                iterateProfile = False):
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

    print "... extracting spectrum ..."
    
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
                raise Exception, "nnls problem - full of noise?"
                ## This column must be completely full of noise
                #x=np.zeros(data.shape[0]+2)
            xArr.append(x)
        
        # Below here as usual
        xArr=np.array(xArr).transpose()
        sky=xArr[1]
        sky2d=np.array([sky]*data.shape[0])
        skySub=skySub-subFrac*sky2d 
        signal=xArr[0]
        signal[np.less(signal, 0)]=0.
        signalArr[k]=signal
        # CR rejection
        arr=skySub
        thresholdSigma=30.0 # was 30
        sigmaCut=3.0
        mean=0
        sigma=1e6
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
            print "... iteration %d ( diff = " % (k), diff, "time taken = %.3f sec" % (t1-t0), ")"
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
def weightedExtraction(data, medColumns = 10, thresholdSigma = 30.0, sigmaCut = 3.0, profSigmaPix = 4.0):
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
    print ">>> Extracting and stacking..."
    for extension in extensionsList:
        print "... %s ..." % (extension)
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
                foundExtension=True
            except KeyError:
                print "... WARNING: missing %s in %s ..." % (extension, fileName)
                foundExtension=False
            
            if foundExtension == True:
                
                print "... %s ..." % (fileName)
                
                header=img[extension].header
                
                # Extract calibrated wavelength scale, assuming left most pixel corresponds to CRVAL1
                # NOTE: if this fails due to missing CDELT1, it probably means the slit we found is a wacky shape - check the .reg file
                w=np.arange(data.shape[1])*header['CDELT1']+header['CRVAL1']
                    
                if w[0] != header['CRVAL1']:
                    raise Exception, "wavelength of pixel 0 doesn't correspond to CRVAL1 - what happened?"
                
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
                    print "WARNING: empty slit"
                    signalList.append(np.zeros(data.shape[1]))
                    skyList.append(np.zeros(data.shape[1]))

        # Wavelength calib diagnostic: does the sky line up from each science frame?
        checkSkyLines=[5577, 5893, 6300]  # [OI], NaI, [OI] sky lines, for checking wavelength calib
        fluxMax=0
        for sky, wavelengths in zip(skyList, wavelengthsList):
            flux=sky/np.median(sky)
            if flux.max() > fluxMax:
                fluxMax=flux.max()
            plt.plot(wavelengths, flux)
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
        wavelength=np.median(wavelengthsArr, axis = 0)
        regrid_signalArr=np.zeros(signalArr.shape)
        regrid_skyArr=np.zeros(skyArr.shape)
        for i in range(signalArr.shape[0]):
            tck=interpolate.splrep(wavelengthsArr[i], signalArr[i])
            regrid_signalArr[i]=interpolate.splev(wavelength, tck, ext = 1)
            tck=interpolate.splrep(wavelengthsArr[i], skyArr[i])
            regrid_skyArr[i]=interpolate.splev(wavelength, tck, ext = 1) 
        signal=np.median(regrid_signalArr, axis = 0)
        sky=np.median(regrid_skyArr, axis = 0) 
        outFileName=extractStackSpecDir+os.path.sep+"1D_"+maskDict['objName'].replace(" ", "_")+"_"+maskDict['maskID']+"_"+extension+".fits"
        write1DSpectrum(signal, sky, wavelength, outFileName)

        # 2d combined/extracted, projecting CR-masked arrays to same wavelength grid
        # NOTE: These all have different wavelength calibrations stored in wavelengthsList
        # NOTE: This assumes that the wavelength calibrations have similar accuracy... if [0] has a problem, then we do!
        CRMaskedDataCube=np.ma.masked_array(CRMaskedDataCube)
        projDataCube=[]
        refWavelengths=wavelengthsArr[0]
        refHeader=headersList[0]
        for data, wavelengths, header in zip(CRMaskedDataCube, wavelengthsList, headersList):
            if projDataCube == []:
                projDataCube.append(data)
            else:
                projData=np.zeros(projDataCube[0].shape)
                projMask=np.zeros(projDataCube[0].shape)
                for i in range(min([projData.shape[0], data.shape[0]])):
                    tck=interpolate.splrep(w, data.data[i])
                    projData[i]=interpolate.splev(refWavelengths, tck, ext = 1)
                    tck=interpolate.splrep(w, data.mask[i])
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
        # Experimenting with a method that will handle running profile
        #t0=time.time()
        #signal, sky, skySubbed2d=finalExtraction(med, subFrac = subFrac)
        #t1=time.time()
        #print "... final extraction (took %.3f sec) ..." % (t1-t0)
        
        outFileName=stackExtractSpecDir+os.path.sep+"1D_"+maskDict['objName'].replace(" ", "_")+"_"+maskDict['maskID']+"_"+extension+".fits"
        write1DSpectrum(signal, sky, refWavelengths, outFileName)
        
        # Write 2d combined spectrum
        outFileName=stackExtractSpecDir+os.path.sep+"2D_"+maskDict['objName'].replace(" ", "_")+"_"+maskDict['maskID']+"_"+extension+".fits"
        newImg=pyfits.HDUList()
        hdu=pyfits.PrimaryHDU(None, refHeader)   
        hdu.data=np.array(med)
        newImg.append(hdu)
        newImg.writeto(outFileName, clobber = True)
        newImg.close()

        # Write 2d, sky subtracted combined spectrum
        #outFileName=stackExtractSpecDir+os.path.sep+"2D_noSky_"+maskDict['objName'].replace(" ", "_")+"_"+maskDict['maskID']+"_"+extension+".fits"
        #newImg=pyfits.HDUList()
        #hdu=pyfits.PrimaryHDU(None, refHeader)   
        #hdu.data=np.array(skySubbed2d)
        #newImg.append(hdu)
        #newImg.writeto(outFileName, clobber = True)
        #newImg.close()

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
    sigmaRange=np.linspace(1, prof.shape[0]/2, prof.shape[0]*2)
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
    signal, sky, mData=iterativeWeightedExtraction(data, subFrac = subFrac, runningProfile = runningProf, iterateProfile = False)

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
def write1DSpectrum(signal, sky, wavelength, outFileName):
    """Writes 1D spectrum to .fits table file.
    
    """
    
    specColumn=pyfits.Column(name='SPEC', format='D', array=signal)
    skyColumn=pyfits.Column(name='SKYSPEC', format='D', array=sky)
    lambdaColumn=pyfits.Column(name='LAMBDA', format='D', array=wavelength)
    tabHDU=pyfits.new_table([specColumn, skyColumn, lambdaColumn])
    tabHDU.name='1D_SPECTRUM'
    HDUList=pyfits.HDUList([pyfits.PrimaryHDU(), tabHDU])
    HDUList[0].header['MASKRA']=maskDict['RA']
    HDUList[0].header['MASKDEC']=maskDict['DEC']
    HDUList.writeto(outFileName, clobber=True)   
    
#-------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    parser = argparse.ArgumentParser("rss_mos_reducer.py")
    parser.add_argument("rawDir", help="The directory that holds the partially-processed data. Usually called product.")
    parser.add_argument("reducedDir", help="The directory where the reduced data will be written, along with any diagnostic data. This directory will be created if it doesn't already exist.")
    parser.add_argument("maskName", help="Use 'all' to reduce all data found under rawDir/, or 'list' to list all masks (by object) found under rawDir/. The maskName is made from the keyword combination OBJECT_MASKID found in the .fits headers.") 
    parser.add_argument("-t", "--threshold", dest="threshold", default=0.1, help="Threshold used for the slit finding algorithm - check that all slits are being found using e.g. masterflat_0.fits and the .reg files produced by the pipeline. Values in the range 0.1-0.4 work best (default=0.2; increase this value if some slits are getting divided; decrease if slits missing). This option only applies to MOS masks.")
    parser.add_argument("-i", "--iterative-extraction", dest="iterativeMethod", action='store_true', help = "Enable the iterative spectral extraction method.")
    parser.add_argument("-f", "--skysub-fraction", dest="subFrac", default=0.8, help="Fraction of the sky background to be subtracted in each iteration of the iterative spectral extraction algorithm (default=0.8; increase this value for faster convergence). This only has an effect if the -i switch is also used.")
    parser.add_argument("-e", "--exclude-masks", dest="excludeMasks", default="", help="Names of masks to exclude (if using maskName = 'all'). Separate mask names with , but no spaces (e.g., -e M1,M2). Useful for avoiding inclusion of e.g., standard stars.")
    parser.add_argument("-s", "--slits", dest="extensionsList", default="all", help="Reduce the data corresponding to the given slit names only. Separate slit names with , but no spaces (e.g., -e SLIT9,SLIT15).")

    args = parser.parse_args()

    rawDir=args.rawDir
    baseOutDir=args.reducedDir
    maskName=args.maskName
    
    threshold=float(args.threshold)
    subFrac=float(args.subFrac)
    iterativeMethod=args.iterativeMethod
    excludeMasks=args.excludeMasks.split(",")
    extensionsList=args.extensionsList
    if extensionsList != "all":
        extensionsList=extensionsList.split(",")
    
    if os.path.exists(baseOutDir) == False:
        os.makedirs(baseOutDir)
        
    # Sort out what's what...
    infoDict=getImageInfo(rawDir)
    
    if maskName == 'list':
        print "Masks found: %s" % (str(infoDict.keys()))
        sys.exit()
    elif maskName != 'all':
        shortDict={}
        for key in infoDict.keys():
            if key == maskName:
                shortDict[key]=infoDict[key]
        infoDict=shortDict
        
    if maskName != 'all' and maskName not in infoDict.keys():
        print "ERROR: maskName not found. Try using 'list' to see available maskNames."
        sys.exit()

    masksToReduce=[]
    for key in infoDict.keys():
        if key not in excludeMasks:
            masksToReduce.append(key)
            
    # We're organised by object name, reduce each in turn
    for maskName in masksToReduce:
        
        print ">>> Mask: %s" % (maskName)
        
        outDir=baseOutDir+os.path.sep+maskName
        if os.path.exists(outDir) == False:
            os.makedirs(outDir)
    
        # Log what options we're running with
        logFile=file(outDir+os.path.sep+"rss_mos_reducer.log", "w")
        logFile.write("started: %s\n" % (datetime.datetime.now().isoformat()))
        for key in args.__dict__.keys():
            logFile.write("%s: %s\n" % (key, str(args.__dict__[key])))
        logFile.close()
    
        # Tied ourselves in knots a bit here...
        maskDict=infoDict[maskName][infoDict[maskName]['maskID']]
        maskDict['maskID']=infoDict[maskName]['maskID']
        maskDict['objName']=infoDict[maskName]['objName']
        maskType=infoDict[maskName]['maskType']
        
        makeMasterFlats(maskDict, outDir)

        if maskType == 'MOS':
            cutIntoSlitLets(maskDict, outDir, threshold = threshold)
        elif maskType == 'LONGSLIT':
            cutIntoPseudoSlitLets(maskDict, outDir)
        
        applyFlatField(maskDict, outDir)
        
        wavelengthCalibration2d(maskDict, outDir, extensionsList = extensionsList)

        extractAndStackSpectra(maskDict, outDir, extensionsList = extensionsList, iterativeMethod = iterativeMethod, subFrac = subFrac)
 
    
    
