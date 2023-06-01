.. _Development:

===================================
Contributing to Further Development
===================================

RSSMOSPipeline is hosted on `Github <https://github.com/mattyowl/RSSMOSPipeline>`_ and is available under a free
software license. Help is appreciated in its development.


Reporting Issues
----------------

Found a bug? Please tell us, either on the `issues page <https://github.com/mattyowl/RSSMOSPipeline/issues>`_
or by `email <matt.hilton@wits.ac.za>`_.


Contributing Code
-----------------

Want to add a feature or fix something? The preferred method of contributing code is to clone
the `repository <https://github.com/mattyowl/RSSMOSPipeline>`_  and work on your new feature
in your own branch::

    git clone https://github.com/mattyowl/RSSMOSPipeline.git
    git checkout -b the-name-of-your-branch

When the time comes to commit your changes, please contact `Matt Hilton <matt.hilton@wits.ac.za>`_  with your
GitHub username, in order to be granted write access to the repository. You only need to do this once.

Alternatively, you can work on changes in a fork of the repository.

When you are ready for your changes to be added to the ``master`` branch, please 
issue a Pull Request for your changes to be reviewed. 


Style
^^^^^

When adding code, please adhere to the style used throughout RSSMOSPipeline where possible.

* As you may notice, RSSMOSPipeline uses `camelCase <https://en.wikipedia.org/wiki/Camel_case>`_ throughout -
  please keep it that way.

* Indent with 4 spaces.

* The maximum line length is 110 characters (sometimes it makes sense to break this).

* Docstrings *should* follow the `Google style <https://www.sphinx-doc.org/en/master/usage/extensions/example_google.html>`_.
  This has only been partially done (so far) in the existing code. At the very least, every function
  should have a docstring of some kind that describes what it does, even if it is re-formatted later.


Testing
^^^^^^^

RSSMOSPipeline uses the Robot framework for tests (see :ref:`TestingPage`). These are integration tests,
rather than unit tests (at least at the moment), and can take a while to run. While more work 
(and more tests) need to be added, you should check that these tests still pass (or at least, 
don't crash) before committing your changes.

