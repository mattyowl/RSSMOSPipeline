# -*- coding: iso-8859-1 -*-
#
# RSSMOSPipeline install script

import os
import glob
from setuptools import setup
from setuptools import Extension
#from distutils.core import setup
#from distutils.extension import Extension
#from Cython.Distutils import build_ext
import numpy
#import popen2

setup(name='RSSMOSPipeline',
      version="1.0.dev0",
      url="https://github.com/mattyowl/RSSMOSPipeline",
      author='Matt Hilton',
      author_email='hiltonm@ukzn.ac.za',
      classifiers=[],
      description='Pipeline for reducing both longslit and multi-object spectroscopic data from the Robert Stobie Spectrograph on SALT.',
      long_description="""Pipeline for reducing both longslit (added May 2016) and multi-object spectroscopy from the Robert Stobie Spectrograph on SALT.""",
      packages=['RSSMOSPipeline'],
      package_data={'RSSMOSPipeline': ['data/*', 'data/modelArcSpectra/*']},
      scripts=['bin/rss_mos_reducer', 'bin/rss_mos_create_arc_model', 'bin/rss_mos_inspect_arc_model'],
      #cmdclass={'build_ext': build_ext},
      #ext_modules=[Extension("nemoCython", ["nemo/nemoCython.pyx"], include_dirs=[numpy.get_include()])]
      install_requires=["astropy >= 3.2",
                        "numpy >= 1.10",
                        "matplotlib >= 2.0",
                        "scipy >= 1.0"]
)
