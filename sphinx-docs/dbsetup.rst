Database Setup
==============

The first step in setting up a working PFMFind system is to populate a database of sequences and their annotations.


Quick guide
-----------

* Download the NCBI Taxonomy database from `<ftp://ftp.ncbi.nih.gov/pub/taxonomy/>`_ (this can be skipped and the taxonomy downloaded automatically by the ``PFMFsetupdb.py`` script.

* Download the appropriate Uniprot, Uniref and Interpro files:

  * `UniprotKB files <http://www.uniprot.org/downloads>`_ in text (.dat) format;

  * `Uniref files <http://www.uniprot.org/downloads>`_ in XML format;

  * Interpro file ``protein2ipr.dat`` from `<http://www.ebi.ac.uk/interpro/download.html>`_.

*  Prepare the configuration file (see :ref:`sec-dbconfig` below).

* Run the ``PFMFsetupdb.py`` script (for example with configuration file ``dbsetup.xml``::

  $ PFMFsetupdb.py dbsetup.xml

* Wait until the database is populated. If you have several configuration files associated with different steps, repeat the step 4. as needed.


``PFMFsetupdb.py`` script
-------------------------


``PFMFsetupdb.py`` creates and loads the sequence database schema in the following sequence:

* Create a database schema and BioSQL tables;

* Load NCBI taxonomy information;

* Load Uniprot sequences and annotations;

* Load Uniref cluster information;

* Load Interpro domain informations.

Depending on the configuration file, each step can be done separately or all together. The sequence however should be preserved except that steps 4. and 5. may be interchanged (or even totally omitted if only sequence/Uniprot data is required).

.. tip::
   While PFMFind requires certain file formats, it does not require that the data originates from the Uniprot, Uniref or Interpro databases. For example, you may supply your own custom sequence database in Uniprot format, or a set of sequence features in Interpro format.


.. _sec-dbconfig:

Configuration file
------------------

An example of a ``PFMFsetupdb.py`` configuration file in XML format is shown below::

   <?xml version="1.0" encoding="UTF-8"?>
   <!DOCTYPE PFMF_db_setup SYSTEM "PFMFdb.dtd">
   <PFMF_db_setup>
     <Database driver="psycopg2" user="aleksand" db="PFMFind"/>
     <Schema name="SwissProt02082005" create="1"/>
     <Sql_scripts sql_start="biosqldb-pg-nocnstr.sql" sql_end="biosqldb-pg-cnstr.sql"/>
     <Taxonomy copy="1">
        <Taxon_dir>/home/aleksand/data/bio/taxon_tables</Taxon_dir>
     </Taxonomy>
     <Uniprot_file namespace="SwissProt">uniprot_sprot.dat</Uniprot_file>
     <Uniref_file namespace="SwissProt">uniref50.xml</Uniref_file>
     <InterPro_file namespace="SwissProt">protein2interpro.dat</InterPro_file>
   </PFMF_db_setup>

All configuration tags are under ``<PFMF_db_setup>`` tag and the specification is given in the file ``PFMFdb.dtd`` located in the ``examples/setup_config`` directory of the PFMFind distribution. Hence, each configuration file can be checked for validity before being passed to ``PFMFsetupdb.py``. Below is the description of each element and its attributes.


``<Database>``
^^^^^^^^^^^^^^

This is a mandatory element that is empty and that must start the configuration options. It has six possible attributes:

* ``driver`` (mandatory): the PostgreSQL Python driver, usually *psycopg2* but others may work as well;

* ``db`` (mandatory): the name of the PostgreSQL database to connect;

* ``user`` (optional): database user name;

* ``password`` (optional): database password;

* ``host`` (optional): the host the database is running on (if left out, *localhost* is assumed);

* ``port`` (optional): the port on the host the database is listening to.


``<Schema>``
^^^^^^^^^^^^
An optional empty element (the default PostgreSQL schema is used if omitted) describing the schema that will contain the dataset and annotations. Two attributes:

* ``name`` (mandatory): the name of the schema;

* ``create`` (optional): whether to create the schema. Should be set to *1* if yes and to*0* if the schema already exists. The default is *0*.

``<Sql_dir>``
^^^^^^^^^^^^^

An optional element whose value is the path to the SQL scripts specified by the ``<Sql_scripts>`` element described below. The default is the current working directory.

``<Sql_scripts>``
^^^^^^^^^^^^^^^^^^

An optional element indicating the SQL scripts to be run to create the required tables. One script (attribute ``sql_start`` is run at the beginning and the other (attribute ``sql_end`` at the end of the run of ``PFMFsetupdb.py``. If the tag is omitted, the default scripts are run. To skip a default script (for example, if the tables are already created), set its corresponding attribute to empty string.

``<Taxonomy>``
^^^^^^^^^^^^^^

Optional element with one optional attribute (``copy``). If specified, it must include the element ``<Taxon_dir>``, which indicates the storage directory for the taxon data from NCBI. The ``copy`` attribute can be set to *0* (default), *1* or *2*. If it is *2*, only the tab-separated tables will be loaded into database (this can be used to reuse the tables from several PostgreSQL schemas); if it is *1* the tables will in addition be created from the NCBI taxon data; if it is *0* the taxon data will be downloaded from the NCBI ftp site.

``<Uniprot_file>``, ``<Uniref_file>``, ``<InterPro_file>``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These three tags describe the path to files containing the Uniprot, Uniref and InterPro data, respectively. Each can be repeated as many times as necessary (or omitted). The loaders for Uniref and InterPro consider only the annotations for those (Uniprot) sequences that were stored before - hence, the order of tags is important. All three have the same attributes: ``namespace`` (mandatory, dataset identifier), ``sql_start`` and ``sql_end`` (optional SQL scripts to run before and after loading the data).
