This directory contains SQL scripts to be used when populating the BioSQL
database.

* biosqldb-pg.sql

The original BioSQL Postges script -- not to be used.

* biosqldb-pg-nocnstr

The BioSQL Postges script with all constraints and indexes removed from tables.
To be used at start.

* biosqldb-pg-cnstr

Constraints and indexes added back.
To be used after all data is loaded.

* biosqldb-pg-fk.sql

Foreign key constraints
