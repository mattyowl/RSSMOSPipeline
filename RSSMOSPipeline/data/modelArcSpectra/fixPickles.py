"""Fix all pickles by reading in and writing out again in binary format (otherwise problems
for e.g., Mac users?)

"""

import os
import sys
import glob
import pickle
import IPython

inFileNames=glob.glob("*.pickle")
for f in inFileNames:
    pickleFile=file(f, "r")
    unpickler=pickle.Unpickler(pickleFile)
    refModelDict=unpickler.load()
    pickleFile.close()
    
    pickleFile=file(f, "wb")
    pickler=pickle.Pickler(pickleFile)
    pickler.dump(refModelDict)
    pickleFile.close()
    
    #IPython.embed()
    #sys.exit()
    