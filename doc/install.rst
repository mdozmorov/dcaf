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
   
   The exact commands depend on your platform. To install numpy via
   pip, you must have the Python headers installed. On Debian-based
   systems, like Ubuntu, this is the "python3-dev" package.
   
   Many Linux distributions also have a pre-compiled version of numpy
   available through the system package manager.
   
   For more help installing numpy, see
   http://www.scipy.org/scipylib/download.html

#. Install dcaf from PyPI, or directly from git for the bleeding-edge version:
   
   .. code-block:: bash
       
      $ pip3 install dcaf
   
   OR
   
   .. code-block:: bash
        
      $ pip3 install git+ssh://git@bitbucket.org/wrenlab/dcaf.git
      
   The remaining required dependencies should be installed for you.

Optional dependencies
=====================

There are some (non-Python) packages that may be required for full
functionality. These are generally third-party bioinformatics programs
written in a variety of languages. The ``dcaf`` installer will check
for these dependencies and warn if any are missing, but will continue
with the installation.

The dependencies, categorized by task, are:

High-throughput sequencing
--------------------------

- `FASTQC - A quality control tool for high-throughput sequence data <http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_

Genome analysis
---------------

- `Jim Kent's source utilities (from the UCSC Genome Browser) <http://genomewiki.ucsc.edu/index.php/The_source_tree>`_

Expression analysis
-------------------

- `URSA (Unveiling RNA Sample Annotation) from Troyanskaya Lab <https://bitbucket.org/youngl/ursa_backend>`_
