============
Installation
============

Supported platforms & Python versions
=====================================

* Linux or Mac OS X
* Python 3.3 or greater

Basic installation
==================

#. Ensure you have Python 3 installed.

#. Install the ``cython`` and ``numpy`` packages:
   
   One way to do this is:

   .. code-block:: bash
                
      $ pip3 install cython numpy
   
   But the exact commands depend on your platform. To install numpy
   via pip, you must have the Python headers installed. On
   Debian-based systems, like Ubuntu, this is the "python3-dev" package.
   
   For more help installing numpy, see
   http://www.scipy.org/scipylib/download.html

#. Install dcaf from PyPI, or directly from git for the bleeding-edge version:
   
   .. code-block:: bash
       
      $ pip3 install dcaf
   
   OR
   
   .. code-block:: bash
        
      $ pip3 install git+ssh://git@bitbucket.org/wrenlab/dcaf.git
      
   The remaining dependencies should be installed for you.
