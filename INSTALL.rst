**RSSMOSPipeline** is written in pure Python (only 3.x is supported now). It needs the following modules to be installed:

* numpy
* scipy
* astropy
* matplotlib
* IPython

The install script (see below) should install the needed modules automatically if they are not already on your system.
IPython is used for debugging, but isn't really needed to run the pipeline.

The latest tagged version of **RSSMOSPipeline** can be installed using ``pip``:

.. code-block::

   pip install RSSMOSPipeline

You may also install using the standard ``setup.py`` script, e.g., as root:

.. code-block::

   sudo python setup.py install

Alternatively,

.. code-block::

   python setup.py install --user

will install ``RSSMOSPipeline`` under ``$HOME/.local`` (on Ubuntu), and in some other default location on Mac.

You can also use the ``--prefix`` option, e.g.,

.. code-block::

   python setup.py install --prefix=$HOME/local

and then add ``$HOME/local/bin`` to $PATH, and e.g., ``$HOME/local/lib/python3.6/site-packages`` to
$PYTHONPATH (adjust the path according to your Python version number).

.. code-block::

   export PATH=$HOME/local/bin:$PATH
   export PYTHONPATH=$HOME/local/lib/python3.6/site-packages:$PYTHONPATH

If **RSSMOSPipeline** has installed correctly, then you should find its command line tools are
available, for example,

.. code-block::

   rss_mos_reducer -h

should display a helpful message about the command-line options for the main ``rss_mos_reducer`` command.
