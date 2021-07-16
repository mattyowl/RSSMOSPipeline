# -*- coding: iso-8859-1 -*-
#
# RSSMOSPipeline install script

import os
import glob
from setuptools import setup
from setuptools import Extension
import versioneer
import numpy

setup(name='RSSMOSPipeline',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      url="https://github.com/mattyowl/RSSMOSPipeline",
      author='Matt Hilton',
      author_email='hiltonm@ukzn.ac.za',
      classifiers=[],
      description='Pipeline for reducing both longslit and multi-object spectroscopic data from the Robert Stobie Spectrograph on SALT.',
      long_description="""Pipeline for reducing both longslit and multi-object spectroscopy from the Robert Stobie Spectrograph on SALT.""",
      packages=['RSSMOSPipeline'],
      package_data={'RSSMOSPipeline': ['data/*', 'data/modelArcSpectra/*', 'data/templateSpectra/*']},
      scripts=['bin/rss_mos_reducer', 'bin/rss_mos_create_arc_model', 'bin/rss_mos_inspect_arc_model', 'bin/rss_mos_visual_inspector'],
      install_requires=["astropy >= 3.2",
                        "numpy >= 1.10",
                        "matplotlib >= 2.0",
                        "scipy >= 1.0",
                        "astLib >= 0.11.6"]
)
