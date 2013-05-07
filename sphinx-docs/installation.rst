Installation
============

Each machine that runs the PFMFind client or FSIndex requires the PFMFind Python module and its prerequisites. A machine running the database needs to be have PostgreSQL installed. Here, we only describe installation of the Python module; for PostgreSQL-related instruction please refer to the appropriate `PostgreSQL manual <http://www.postgresql.org/docs/>`_.

Prerequisites
-------------

PFMFind requires Python2 (ideally version 2.7) with Tkinter. Before installing PFMFind you also need to download and install the following Python modules (the versions used for testing are indicated in square brackets):

* `Numpy <http://www.numpy.org/>`_ [1.6.2]

* `Scipy <http://www.scipy.org/>`_ [0.11]

* `Biopython <http://www.biopython.org/>`_ [1.58]

* `Pmw megawidgets toolkit <http://pmw.sourceforge.net/>`_ [1.3.3]

* `psycopg2 <http://initd.org/psycopg/>`_ [2.4.5]

Please consult the documentation of individual packages for information of how to install them.

Installing on UNIX platforms
----------------------------

PFMFind is installed in the standard `distutils way <http://docs.python.org/inst/inst.html>`_: Extract the source distribution, go to the top distribution directory and type::

  $ python setup.py install


Please read the `Installing Python Modules <http://docs.python.org/inst/inst.html>`_ for customisation options. Note that you need a C compiler (such as gcc) to build the C extensions included in PFMFind.

Other platforms
---------------

Version 0.5 is currently not supported for platforms other than UNIX. You will need to build the entire package yourself (including C extensions) and install it as appropriate for your case.
