"""

Library for running RSSMOSPipeline tests using Robot Framework

"""

import os
import sys
import subprocess
import shutil
import numpy as np
from astLib import *
import astropy.io.fits as pyfits
import astropy.table as atpy
from astropy.coordinates import SkyCoord, match_coordinates_sky
import pylab as plt

#------------------------------------------------------------------------------------------------------------
class RSSMOSPipelineTests(object):

    def __init__(self):
        """Basic set-up for running tests. For those that need datasets to be downloaded, a setup_<name>
        routine should be added (see, e.g., setup_longslit).
        
        """

        self._status = ''

        self.cacheDir="testsCache"
        if os.path.exists(self.cacheDir) == False:
            os.makedirs(self.cacheDir)
        
        # Why don't we just run all tests from the cache dir?
        self.runDir=os.path.abspath(self.cacheDir)
            
        self.plotsDir="plots"
        if os.path.exists(self.plotsDir) == False:
            os.makedirs(self.plotsDir)
    
    
    def setup_MOS(self):
        """Set-up for tests that use MOS data - downloads them if not found.
        
        """

        thisDir=os.getcwd()
        self.productDir="product"
        if os.path.exists(self.cacheDir+os.path.sep+self.productDir+os.path.sep+"mbxgpP201708180091.fits") == False:
            print(">>> Downloading MOS data [this will be cached under %s so only needs to be done once] ..." % (self.cacheDir))
            os.chdir(self.cacheDir)
            os.system("wget 'https://www.dropbox.com/scl/fi/1w72qdjrzn7pvjl93zhd2/mos-test-data.tar.gz?dl=0&rlkey=zxh7ddok0mimazf65yw6qx8sw'")
            os.system("tar -zxvf 'mos-test-data.tar.gz?dl=0&rlkey=zxh7ddok0mimazf65yw6qx8sw'")
            os.remove("mos-test-data.tar.gz?dl=0&rlkey=zxh7ddok0mimazf65yw6qx8sw")
            os.chdir(thisDir)


    def setup_longslit(self):
        """Set-up for tests that use longslit data - downloads them if not found.

        """

        thisDir=os.getcwd()
        self.productDir="product_longslit"
        if os.path.exists(self.cacheDir+os.path.sep+self.productDir+os.path.sep+"mbxgpP202102040049.fits") == False:
            os.chdir(self.cacheDir)
            self._get_longslit_data()
            os.chdir(thisDir)


    def setup_no_flats_longslit(self):
        """Set-up for tests that use longslit data - downloads them if not found.

        """

        thisDir=os.getcwd()
        self.productDir="product_longslit_noflats"
        if os.path.exists(self.cacheDir+os.path.sep+self.productDir+os.path.sep+"mbxgpP202102040049.fits") == False:
            os.chdir(self.cacheDir)
            self._get_longslit_data()
            os.chdir(thisDir)


    def _get_longslit_data(self):
            print(">>> Downloading longslit data [this will be cached under %s so only needs to be done once] ..." % (self.cacheDir))
            os.system("wget 'https://www.dropbox.com/scl/fi/w42hj8x0xa46aqpkwdugf/longslit-test-data.tar.gz?rlkey=j52i0zkqnax9xgx94ujck1wkc&dl=0'")
            # os.system("wget 'https://www.dropbox.com/scl/fi/1w72qdjrzn7pvjl93zhd2/mos-test-data.tar.gz?dl=0&rlkey=zxh7ddok0mimazf65yw6qx8sw'")
            os.system("tar -zxvf 'longslit-test-data.tar.gz?rlkey=j52i0zkqnax9xgx94ujck1wkc&dl=0'")
            os.remove("longslit-test-data.tar.gz?rlkey=j52i0zkqnax9xgx94ujck1wkc&dl=0")


    def makeSlitFile(self, slitNum, yStart, yEnd):
        """Makes a file called manualSlitLoc.txt that contains location of a single slit, manually specified.

        """
        with open(self.cacheDir+os.path.sep+"manualSlitLoc.txt", "w") as outFile:
            outFile.write("slitno ystart yend\n")
            outFile.write("%s     %d     %d\n" % (int(slitNum), int(yStart), int(yEnd)))


    def run_reduction(self, reducedDir = None, maskName = None):
        args=['rss_mos_reducer', self.productDir, reducedDir, maskName]
        self._run_command(args)


    def run_reduction_on_selected_slits(self, reducedDir = None, maskName = None, slits = None):
        args=['rss_mos_reducer', self.productDir, reducedDir, maskName, '-s', slits]
        self._run_command(args)


    def run_no_flats_reduction_on_selected_slits(self, reducedDir = None, maskName = None, slits = None):
        args=['rss_mos_reducer', self.productDir, reducedDir, maskName, '-n', '-s', slits]
        self._run_command(args)


    def run_reduction_using_slit_file(self, reducedDir = None, maskName = None, slitFileName = None):
        args=['rss_mos_reducer', self.productDir, reducedDir, maskName, '-F', slitFileName]
        self._run_command(args)


    def cross_match(self, inCatalogFileName, outCatalogFileName, radiusArcmin = 1.0):
        """Cross matches input and output source catalogs.
        
        """
        inTab=atpy.Table().read(inCatalogFileName)
        outTab=atpy.Table().read(outCatalogFileName)
        self.inTab, self.outTab, self.rDeg=catalogs.crossMatch(inTab, outTab, 
                                                               radiusArcmin = radiusArcmin)
    
    
    def getRADecKeys(self, tab):
        """Returns the column names in which RA, dec coords are stored, after trying a couple of variations.
        
        """
        RAKeysToTry=['ra', 'RADeg']
        decKeysToTry=['dec', 'decDeg']
        RAKey, decKey=None, None
        for key in RAKeysToTry:
            if key in tab.keys():
                RAKey=key
                break
        for key in decKeysToTry:
            if key in tab.keys():
                decKey=key
                break
        if RAKey == None or decKey == None:
            raise Exception("couldn't identify RA, dec columns in the supplied table")
        
        return RAKey, decKey
        
        
    def check_recovered_ratio(self, inKey, outKey, toleranceSigma = 1.0, expectedRatio = 1.0, SNRCut = 4,
                              SNRKey = 'fixed_SNR', errInKey = None, errOutKey = None,
                              plotLabel = None, plotsDir = "plots"):
        """Catalogs must have been cross matched before this can be run.
        Calculate the ratio of the columns pointed to by inKey, outKey. If tolerance (defined as
        1 - average ratio) is exceeded, the test is failed. SNRCut is applied to SNRKey in the output catalog.
        
        """

        inTab=self.inTab
        outTab=self.outTab
        mask=outTab[SNRKey] > SNRCut
        x=inTab[inKey]
        y=outTab[outKey]
        meanRatio=np.mean(y[mask])/np.mean(x[mask])
        bsRatios=[]
        for i in range(5000):
            indices=np.random.randint(0, len(x[mask]), len(x[mask]))
            bsX=x[mask][indices]
            bsY=y[mask][indices]
            bsRatios.append(np.mean(bsY)/np.mean(bsX))
        meanRatioErr=np.percentile(abs(bsRatios-meanRatio), 68.3)

        label="<input %s>/<output %s> = %.3f ± %.3f (%s > %.1f)" % (inKey, outKey, meanRatio, meanRatioErr, SNRKey, SNRCut)
        print("%s" % (label))

        if plotLabel is not None:
            plotMin=0.1*min([inTab[inKey].min(), outTab[outKey].min()])
            plotMax=10*max([inTab[inKey].max(), outTab[outKey].max()])
            plotRange=np.linspace(plotMin, plotMax, 100)
            plt.figure(figsize = (10, 8))
            if errInKey is not None and errOutKey is not None:
                plt.errorbar(inTab[inKey], outTab[outKey], yerr = outTab[errOutKey],
                     xerr = inTab[errInKey], elinewidth = 1.5, ecolor = '#AAAAFF',
                     fmt = 'D', ms = 6, label = None)
            else:
                plt.plot(inTab[inKey], outTab[outKey], 'D')
            plt.xlabel("input "+inKey)
            plt.ylabel("output "+outKey)
            plt.plot(plotRange, plotRange, 'k--')
            plt.xlim(plotMin, plotMax)
            plt.ylim(plotMin, plotMax)
            plt.loglog()
            plt.title(label, fontdict = {'size': plotTitleSize})
            plt.savefig(plotsDir+os.path.sep+plotLabel+"_XvY.png")
            plt.close()
        if abs((expectedRatio-meanRatio)/meanRatioErr) > toleranceSigma:
            self._status="FAILED"
        else:
            self._status="SUCCESS"
            
            
    def check_recovered_positions(self, toleranceArcsec = 12.0, SNRKey = 'SNR', SNRMax = 10.0, 
                                  plotLabel = None, plotsDir = "plots"):
        """Blah
        
        """
        inTab=self.inTab
        outTab=self.outTab
        rDeg=self.rDeg
        SNRs=outTab[SNRKey]
        medOffsetArcmin=np.median(rDeg)*60
        SNRMask=np.less(SNRs, SNRMax)
        medOffsetArcmin_SNR10=np.median(rDeg[SNRMask])*60#[SNRMask])*60
        print('... median recovered position offset = %.2f" (full sample)' % (medOffsetArcmin*60))
        label='median recovered position offset = %.2f" (SNR < 10)' % (medOffsetArcmin_SNR10*60)
        print("... %s" % (label))
        if plotLabel is not None:
            plt.figure(figsize = (10, 8))
            plt.plot(SNRs, rDeg*3600., 'r.')
            plt.semilogx()
            plt.xlabel("SNR")
            plt.ylabel('Position offset (")')
            plt.title(label, fontdict = {'size': plotTitleSize})
            plt.savefig(plotsDir+os.path.sep+plotLabel+"_posRec.png")
            plt.close()
        if medOffsetArcmin_SNR10*60 > toleranceArcsec:
            self._status="FAILED"
        else:
            self._status="SUCCESS"
    
    
    def subtract_maps(self, map1FileName, map2FileName, outMapFileName):
        with pyfits.open(map1FileName) as img1:
            d1=img1[0].data
            wcs=astWCS.WCS(img1[0].header, mode = 'pyfits')
        with pyfits.open(map2FileName) as img2:
            d2=img2[0].data
        maps.saveFITS(outMapFileName, d1-d2, wcs)


    def check_map_sigma(self, map1FileName, expectedSigma, tol = 1.0):
        with pyfits.open(map1FileName) as img1:
            d1=img1[0].data
        print("args:", map1FileName, expectedSigma)
        diff=abs(np.std(d1.flatten()) - float(expectedSigma))
        print("diff from expected sigma", diff)
        if diff > tol:
            self._status="FAILED"
        else:
            self._status="SUCCESS"
    
            
    def status_should_be(self, expected_status):
        if expected_status != self._status:
            raise AssertionError("Expected status to be '%s' but was '%s'."
                                 % (expected_status, self._status))


    def _run_command(self, args):
        thisDir=os.getcwd()
        os.chdir(self.runDir)
        print(args)
        process=subprocess.run(args, universal_newlines=True, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
        if process.returncode != 0:
             print(process.stdout)
             raise AssertionError("Return code of '%s' is non-zero." % (str(args)))
        os.chdir(thisDir)
