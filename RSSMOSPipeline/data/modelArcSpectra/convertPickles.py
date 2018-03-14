"""

Convert pickles from atpy to astropy.table

"""

import os
import sys
import glob
import atpy
import astropy.table as atab
import pickle
import IPython

inFileNames=glob.glob("*.pickle")

for modelFileName in inFileNames:
    
    pickleFile=file(modelFileName, "rb")
    unpickler=pickle.Unpickler(pickleFile)
    refModelDict=unpickler.load()
    pickleFile.close()
    
    newTab=atab.Table()
    tab=refModelDict['featureTable']
    for k in tab.keys():
        newTab.add_column(atab.Column(tab[k], k))
    refModelDict['featureTable']=newTab
    
    pickleFile=file("converted_"+modelFileName, "wb")
    pickler=pickle.Pickler(pickleFile)
    pickler.dump(refModelDict)
    pickleFile.close()
