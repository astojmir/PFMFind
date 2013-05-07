Index Setup
===========

After all sequences are loaded, it is time to set up search indexes using the ``PFMFsetupix.py``. The configuration parameters are once again supplied through an XML file, such as the one shown below::

      <?xml version="1.0" encoding="UTF-8"?>
      <!DOCTYPE PFMF_index_setup SYSTEM "PFMFix.dtd">
      <PFMF_index_setup>
        <Database driver="psycopg" host="130.195.61.38" port="5432"
                  user="aleksander" password="aleksander" db="PFMFind"/>
        <Index_dir>/home/aleksander/data/uniprot-3.5/FSindex</Index_dir>
        <Dataset name="sprot" schema="PFMFind02" namespace="SwissProt"
                 max_residues="80000000">
          <Index length="10">
            <Partition>STAN#ILVM#KRDEQ#WFYH#GPC</Partition>
          </Index>
          <Index length="12">
            <Partition>STAN#ILVM#KRDEQ#WFYHGPC</Partition>
          </Index>
        </Dataset>
      </PFMF_index_setup>


Configuration file
------------------

The index configuration file is very similar to the :ref:`database configuration file <sec-dbconfig>`. All tags are under the ``<PFMF_index_setup>`` tag and the specification is given in the file ``PFMFix.dtd`` located in the ``examples/setup_config`` directory of the PFMFind distribution. Below is the description of each element and its attributes.

``<Database>``
^^^^^^^^^^^^^^

These are database connection details. Same as for database setup described :ref:`earlier <sec-dbconfig>`.

``<Index_dir>``
^^^^^^^^^^^^^^^

A mandatory element describing the path that will contain the sequence files in FASTA format and their indexes.

``<Dataset>``
^^^^^^^^^^^^^

An element containing details about each dataset to be indexed. There can be multiple ``<Dataset>`` elements, each specifying multiple indexes.

It has four possible attributes:

* ``name`` (mandatory): the dataset identifier used as a prefix for all files related to this dataset;

* ``schema`` (optional): the PostgreSQL schema containing the dataset. If omitted, the setup script attempts to use the schemas from the default PostgreSQL path.

* ``namespace`` (optional): the namespace associated with the dataset. Each schema can contain multiple datasets, which are distinguished using the namespace identifier. If omitted, all sequences from the schema are retrieved.

* ``max_residues`` (optional): the maximum number of residues a single sequence file can contain. If the dataset contains more than that number of residues, multiple files and indexes will be created (Note that each part may contain slightly more residues because the full sequences are not broken up). These parts can be loaded by slave index servers and the whole search distributed by a master server.

Each ``<Dataset>`` element can contain one or more ``<Index>`` elements, whose only (mandatory) attribute is ``length``, giving the fragment length to be used to create the index. In turn, each ``<Index>`` element contains one or more ``<Partition>`` elements, having no attributes and containing a string specifying the amino acid alphabet partitions.

Alphabet partitions are separated by ``#`` (e.g. ``STAN#ILVM#KRDEQ#WFYH#GPC``) and each ``<Partition>`` element specifies an alphabet partition for a single position in the fragment (up to the fragment length). If there are fewer ``<Partition>`` elements than the specified fragment length, the last ``<Partition>`` element is used for all remaining positions. Hence, in order to have the same alphabet partitions at all positions, it is sufficient to specify a single ``<Partition>`` element.
