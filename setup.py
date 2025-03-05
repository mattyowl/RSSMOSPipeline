# -*- coding: iso-8859-1 -*-
#
# RSSMOSPipeline install script

import os
import glob
from setuptools import setup
from setuptools import Extension
import numpy

setup(name='RSSMOSPipeline',
      version='0.7.0',
      author='Matt Hilton + RSSMOSPipeline Contributors',
      author_email='matt.hilton@wits.ac.za',
      packages=['RSSMOSPipeline'],
      package_data={'RSSMOSPipeline': ['data/*', 'data/modelArcSpectra/*', 'data/templateSpectra/*', 'data/templateSpectra/TremontiStarburstTemplate/*']},
      scripts=['bin/rss_mos_reducer', 'bin/rss_mos_create_arc_model', 'bin/rss_mos_inspect_arc_model', 'bin/rss_mos_visual_inspector']
)
