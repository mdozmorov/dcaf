Welcome to dcaf's documentation!
================================

Overview
--------

``dcaf`` is a Python package containing various bioinformatics code
developed in the Wren Lab at the `Oklahoma Medical Research Foundation
<http://omrf.org/>`_.

This package contains code and scripts for various bioinformatics tasks:

- Wrappers for sequence and variant analysis, principally RNA-seq
- Tools for managing a PostgreSQL database of taxon, gene, expression,
  ontology, and textual (MEDLINE) data
- Readers for various bioinformatics file formats, including OBO_,
  SOFT_, and BBI_ (BigWig/BigBED)
- Text mining, enrichment analysis, and other aspects of expression
  analysis
- Meta-analysis of gene expression data, including gene
  function prediction and imputation of experimental sample parameters
  such as tissue, age, and gender

The principal developers are `Cory Giles <mailto:mail@corygil.es>`_
and `Mikhail Dozmorov <mailto:dozmorovm@omrf.org>`_, and the package
is licensed under the `GNU Affero General Public License version 3
<http://www.gnu.org/licenses/agpl-3.0.html>`_ (or any later version
thereof, at your option).

General Information
-------------------

.. toctree::
   :maxdepth: 3
   
   install
   configuration
   tutorial
   text
   api
   todo
             
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. _OBO: http://www.geneontology.org/GO.format.obo-1_2.shtml
.. _SOFT: http://www.ncbi.nlm.nih.gov/geo/info/soft2.html
.. _BBI: http://bioinformatics.oxfordjournals.org/content/26/17/2204.long
