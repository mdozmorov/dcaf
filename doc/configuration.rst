.. _configuration:

=============
Configuration
=============

Many aspects of ``dcaf`` functionality can be altered by creating a
file called ``dcaf.cfg`` or ``.dcaf.cfg``, and placing it in your
system's ``PATH`` or in your ``HOME`` directory. By using this
configuration file, you can alter, e.g.:

- The connection parameters to the ``dcaf`` PostgreSQL database
- The connection parameters to the UCSC Genome Browser MySQL database (if,
  for example, you have a local mirror)
- The directory and settings for caching downloaded and intermediate files
- Logging settings

Default configuration
---------------------

The following are the default configuration settings, with explanatory
comments. Override these settings as needed in your ``dcaf.cfg``:

.. literalinclude:: ../data/defaults.cfg
