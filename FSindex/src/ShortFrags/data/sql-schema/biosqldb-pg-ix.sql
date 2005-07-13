-- CREATE INDEX db_auth on biodatabase ( authority ); 
-- CREATE INDEX taxparent ON taxon ( parent_taxon_id ); 
-- CREATE INDEX taxnametaxonid ON taxon_name ( taxon_id ); 
-- CREATE INDEX taxnamename ON taxon_name ( name ); 
-- CREATE INDEX term_ont ON term ( ontology_id ); 
-- CREATE INDEX trmdbxref_dbxrefid ON term_dbxref ( dbxref_id ); 

-- relationship between controlled vocabulary / ontology term 
-- we use subject/predicate/object but this could also 
-- be thought of as child/relationship-type/parent. 
-- the subject/predicate/object naming is better as we 
-- can think of the graph as composed of statements. 
-- 
-- we also treat the relationshiptypes / predicates as 
-- controlled terms in themselves; this is quite useful 
-- as a lot of systems (eg GO) will soon require 
-- ontologies of relationship types (eg subtle differences 
-- in the partOf relationship) 
-- 
-- this table probably won't be filled for a while, the core 
-- will just treat ontologies as flat lists of terms 
CREATE SEQUENCE term_relationship_pk_seq;
CREATE TABLE term_relationship ( 
	 term_relationship_id INTEGER DEFAULT nextval ( 'term_relationship_pk_seq' ) NOT NULL , 
	 subject_term_id INTEGER NOT NULL , 
	 predicate_term_id INTEGER NOT NULL , 
	 object_term_id INTEGER NOT NULL , 
	 ontology_id INTEGER NOT NULL ) ; 
-- 	 PRIMARY KEY ( term_relationship_id ) , 
-- 	 UNIQUE ( subject_term_id , predicate_term_id , object_term_id , ontology_id ) ) ; 

-- CREATE INDEX trmrel_predicateid ON term_relationship ( predicate_term_id ); 
-- CREATE INDEX trmrel_objectid ON term_relationship ( object_term_id ); 
-- CREATE INDEX trmrel_ontid ON term_relationship ( ontology_id ); 
-- CONFIG: you may want to add this if you can't get the optimizer to
-- use the composite index for the initial keys
--CREATE INDEX trmrel_subjectid ON term_relationship(subject_term_id); 

-- This lets one associate a single term with a term_relationship 
-- effecively allowing us to treat triples as 1st class terms.
-- 
-- At this point this table is only supported in Biojava. If you want
-- to know more about the rationale and idea behind it, read the
-- following article that Mat Pocock posted to the mailing list:
-- http://www.open-bio.org/pipermail/biosql-l/2003-October/000455.html
CREATE TABLE term_relationship_term (
        term_relationship_id INTEGER NOT NULL,
        term_id              INTEGER NOT NULL ) ;
--         PRIMARY KEY ( term_relationship_id ),
--         UNIQUE ( term_id ) 
-- );

-- the infamous transitive closure table on ontology term relationships 
-- this is a warehouse approach - you will need to update this regularly 
-- 
-- the triple of (subject, predicate, object) is the same as for ontology 
-- relationships, with the exception of predicate being the greatest common 
-- denominator of the relationships types visited in the path (i.e., if 
-- relationship type A is-a relationship type B, the greatest common 
-- denominator for path containing both types A and B is B) 
-- 
-- See the GO database or Chado schema for other (and possibly better 
-- documented) implementations of the transitive closure table approach. 
CREATE SEQUENCE term_path_pk_seq;
CREATE TABLE term_path ( 
         term_path_id INTEGER DEFAULT nextval ( 'term_path_pk_seq' ) NOT NULL ,
	 subject_term_id INTEGER NOT NULL , 
	 predicate_term_id INTEGER NOT NULL , 
	 object_term_id INTEGER NOT NULL , 
	 ontology_id INTEGER NOT NULL , 
	 distance INTEGER ) ; 
-- 	 PRIMARY KEY (term_path_id),
-- 	 UNIQUE ( subject_term_id , predicate_term_id , object_term_id , ontology_id , distance ) ) ; 

-- CREATE INDEX trmpath_predicateid ON term_path ( predicate_term_id ); 
-- CREATE INDEX trmpath_objectid ON term_path ( object_term_id ); 
-- CREATE INDEX trmpath_ontid ON term_path ( ontology_id ); 
-- CONFIG: you may want to add this if you can't get the optimizer to
-- use the composite index for the initial keys
--CREATE INDEX trmpath_subjectid ON term_path(subject_term_id); 

-- we can be a bioentry without a biosequence, but not visa-versa 
-- most things are going to be keyed off bioentry_id 
--
-- accession is the stable id, display_id is a potentially volatile, 
-- human readable name. 
--
-- Version may be unknown, may be undefined, or may not exist for a certain
-- accession or database (namespace). We require it here to avoid RDBMS-
-- dependend enforcement variants (version is in a compound alternative key),
-- and to simplify query construction for UK look-ups. If there is no version
-- the convention is to put 0 (zero) here. Likewise, a record with a version
-- of zero means the version is to be interpreted as NULL.
--
-- not all entries have a taxon, but many do. 
--
-- one bioentry only has one taxon! (weirdo chimerias are not handled. tough) 
--
-- Name maps to display_id in bioperl. We have a different column name 
-- here to avoid confusion with the naming convention for foreign keys. 
CREATE SEQUENCE bioentry_pk_seq;
CREATE TABLE bioentry ( 
	 bioentry_id INTEGER DEFAULT nextval ( 'bioentry_pk_seq' ) NOT NULL , 
	 biodatabase_id INTEGER NOT NULL , 
	 taxon_id INTEGER , 
	 name VARCHAR ( 40 ) NOT NULL , 
	 accession VARCHAR ( 40 ) NOT NULL , 
	 identifier VARCHAR ( 40 ) , 
	 division VARCHAR ( 6 ) , 
	 description TEXT , 
	 version INTEGER NOT NULL ) ; 
-- 	 PRIMARY KEY ( bioentry_id ) , 
-- 	 UNIQUE ( accession , biodatabase_id , version ) , 
-- CONFIG: uncomment one (and only one) of the two lines below. The
-- first puts a uniqueness constraint on the identifier column alone;
-- the other one puts a uniqueness constraint on identifier only
-- within a namespace.
--	 UNIQUE ( identifier ) 
--	 UNIQUE ( identifier , biodatabase_id ) 
-- ) ; 

-- CREATE INDEX bioentry_name ON bioentry ( name ); 
-- CREATE INDEX bioentry_db ON bioentry ( biodatabase_id ); 
-- CREATE INDEX bioentry_tax ON bioentry ( taxon_id ); 

-- 
-- bioentry-bioentry relationships: these are typed 
-- 
CREATE SEQUENCE bioentry_relationship_pk_seq;
CREATE TABLE bioentry_relationship ( 
	 bioentry_relationship_id INTEGER DEFAULT nextval ( 'bioentry_relationship_pk_seq' ) NOT NULL , 
	 object_bioentry_id INTEGER NOT NULL , 
	 subject_bioentry_id INTEGER NOT NULL , 
	 term_id INTEGER NOT NULL , 
	 rank INTEGER ) ; 
-- 	 PRIMARY KEY ( bioentry_relationship_id ) , 
-- 	 UNIQUE ( object_bioentry_id , subject_bioentry_id , term_id ) ) ; 

-- CREATE INDEX bioentryrel_trm ON bioentry_relationship ( term_id ); 
-- CREATE INDEX bioentryrel_child ON bioentry_relationship (subject_bioentry_id);
-- CONFIG: you may want to add this if you can't get the optimizer to
-- use the composite index for the initial keys 
--CREATE INDEX bioentryrel_parent ON bioentry_relationship(object_bioentry_id);

-- for deep (depth > 1) bioentry relationship trees we need a transitive 
-- closure table too 
CREATE TABLE bioentry_path ( 
	 object_bioentry_id INTEGER NOT NULL , 
	 subject_bioentry_id INTEGER NOT NULL , 
	 term_id INTEGER NOT NULL , 
	 distance INTEGER ) ;
-- 	 UNIQUE ( object_bioentry_id , subject_bioentry_id , term_id , distance ) ) ; 

-- CREATE INDEX bioentrypath_trm ON bioentry_path ( term_id ); 
-- CREATE INDEX bioentrypath_child ON bioentry_path ( subject_bioentry_id ); 
-- you may want to add this for mysql because MySQL often is broken with 
-- respect to using the composite index for the initial keys 
--CREATE INDEX bioentrypath_parent ON bioentry_path(object_bioentry_id); 

-- some bioentries will have a sequence 
-- biosequence because sequence is sometimes a reserved word 
CREATE TABLE biosequence ( 
	 bioentry_id INTEGER NOT NULL , 
	 version INTEGER , 
	 length INTEGER , 
	 alphabet VARCHAR ( 10 ) , 
	 seq TEXT ) ;
-- 	 PRIMARY KEY ( bioentry_id ) ) ; 

-- CONFIG: add these only if you want them: 
-- ALTER TABLE biosequence ADD COLUMN ( isoelec_pt NUMERIC(4,2) ); 
-- ALTER TABLE biosequence ADD COLUMN (	mol_wgt DOUBLE PRECISION ); 
-- ALTER TABLE biosequence ADD COLUMN ( perc_gc DOUBLE PRECISION ); 

-- database cross-references (e.g., GenBank:AC123456.1) 
--
-- Version may be unknown, may be undefined, or may not exist for a certain
-- accession or database (namespace). We require it here to avoid RDBMS-
-- dependend enforcement variants (version is in a compound alternative key),
-- and to simplify query construction for UK look-ups. If there is no version
-- the convention is to put 0 (zero) here. Likewise, a record with a version
-- of zero means the version is to be interpreted as NULL.
--
CREATE SEQUENCE dbxref_pk_seq;
CREATE TABLE dbxref ( 
	 dbxref_id INTEGER DEFAULT nextval ( 'dbxref_pk_seq' ) NOT NULL , 
	 dbname VARCHAR ( 40 ) NOT NULL , 
	 accession VARCHAR ( 40 ) NOT NULL , 
	 version INTEGER NOT NULL ) ; 
--	 PRIMARY KEY ( dbxref_id ) , 
--	 UNIQUE ( accession , dbname , version ) ) ; 

-- CREATE INDEX dbxref_db ON dbxref ( dbname ); 

-- for roundtripping embl/genbank, we need to have the "optional ID" 
-- for the dbxref. 
-- 
-- another use of this table could be for storing 
-- descriptive text for a dbxref. for example, we may want to 
-- know stuff about the interpro accessions we store (without 
-- importing all of interpro), so we can attach the text 
-- description as a synonym 
CREATE TABLE dbxref_qualifier_value ( 
	 dbxref_id INTEGER NOT NULL , 
	 term_id INTEGER NOT NULL , 
	 rank INTEGER NOT NULL DEFAULT 0 , 
	 value TEXT ) ; 
-- 	 PRIMARY KEY ( dbxref_id , term_id , rank ) ) ; 

-- CREATE INDEX dbxrefqual_dbx ON dbxref_qualifier_value ( dbxref_id ); 
-- CREATE INDEX dbxrefqual_trm ON dbxref_qualifier_value ( term_id ); 

-- Direct dblinks. It is tempting to do this 
-- from bioentry_id to bioentry_id. But that wont work 
-- during updates of one database - we will have to edit 
-- this table each time. Better to do the join through accession 
-- and db each time. Should be almost as cheap 
CREATE TABLE bioentry_dbxref ( 
	 bioentry_id INTEGER NOT NULL , 
	 dbxref_id INTEGER NOT NULL , 
	 rank INTEGER ) ; 
-- 	 PRIMARY KEY ( bioentry_id , dbxref_id ) ) ; 

-- CREATE INDEX dblink_dbx ON bioentry_dbxref ( dbxref_id ); 

-- We can have multiple references per bioentry, but one reference 
-- can also be used for the same bioentry. 
-- 
-- No two references can reference the same reference database entry 
-- (dbxref_id). This is where the MEDLINE id goes: PUBMED:123456. 
CREATE SEQUENCE reference_pk_seq;
CREATE TABLE reference ( 
	 reference_id INTEGER DEFAULT nextval ( 'reference_pk_seq' ) NOT NULL , 
	 dbxref_id INTEGER , 
	 location TEXT NOT NULL , 
	 title TEXT , 
	 authors TEXT , 
	 crc VARCHAR ( 32 ) ) ; 
-- 	 PRIMARY KEY ( reference_id ) , 
-- 	 UNIQUE ( dbxref_id ) , 
-- 	 UNIQUE ( crc ) ) ; 

-- bioentry to reference associations 
CREATE TABLE bioentry_reference ( 
	 bioentry_id INTEGER NOT NULL , 
	 reference_id INTEGER NOT NULL , 
	 start_pos INTEGER , 
	 end_pos INTEGER , 
	 rank INTEGER NOT NULL DEFAULT 0 ) ; 
-- 	 PRIMARY KEY ( bioentry_id , reference_id , rank ) ) ; 

-- CREATE INDEX bioentryref_ref ON bioentry_reference ( reference_id ); 

-- We can have multiple comments per seqentry, and 
-- comments can have embedded '\n' characters 
CREATE SEQUENCE comment_pk_seq;
CREATE TABLE comment ( 
	 comment_id INTEGER DEFAULT nextval ( 'comment_pk_seq' ) NOT NULL , 
	 bioentry_id INTEGER NOT NULL , 
	 comment_text TEXT NOT NULL , 
	 rank INTEGER NOT NULL DEFAULT 0 ) ; 
-- 	 PRIMARY KEY ( comment_id ) , 
-- 	 UNIQUE ( bioentry_id , rank ) ) ; 

-- tag/value and ontology term annotation for bioentries goes here
CREATE TABLE bioentry_qualifier_value ( 
	 bioentry_id INTEGER NOT NULL , 
	 term_id INTEGER NOT NULL , 
	 value TEXT , 
	 rank INTEGER NOT NULL DEFAULT 0 ) ; 
--	 UNIQUE ( bioentry_id , term_id , rank ) ) ; 

-- CREATE INDEX bioentryqual_trm ON bioentry_qualifier_value ( term_id ); 

-- feature table. We cleanly handle 
--   - simple locations 
--   - split locations 
--   - split locations on remote sequences 
CREATE SEQUENCE seqfeature_pk_seq;
CREATE TABLE seqfeature ( 
	 seqfeature_id INTEGER DEFAULT nextval ( 'seqfeature_pk_seq' ) NOT NULL , 
	 bioentry_id INTEGER NOT NULL , 
	 type_term_id INTEGER NOT NULL , 
	 source_term_id INTEGER NOT NULL , 
	 display_name VARCHAR ( 64 ) , 
	 rank INTEGER NOT NULL DEFAULT 0 ) ;
-- 	 PRIMARY KEY ( seqfeature_id ) , 
-- 	 UNIQUE ( bioentry_id , type_term_id , source_term_id , rank ) ) ; 

-- CREATE INDEX seqfeature_trm ON seqfeature ( type_term_id ); 
-- CREATE INDEX seqfeature_fsrc ON seqfeature ( source_term_id ); 
-- CONFIG: you may want to add this if you can't get the optimizer to
-- use the composite index for the initial keys 
--CREATE INDEX seqfeature_bioentryid ON seqfeature(bioentry_id); 

-- seqfeatures can be arranged in containment hierarchies. 
-- one can imagine storing other relationships between features, 
-- in this case the term_id can be used to type the relationship 
CREATE SEQUENCE seqfeature_relationship_pk_seq;
CREATE TABLE seqfeature_relationship ( 
	 seqfeature_relationship_id INTEGER DEFAULT nextval ( 'seqfeature_relationship_pk_seq' ) NOT NULL , 
	 object_seqfeature_id INTEGER NOT NULL , 
	 subject_seqfeature_id INTEGER NOT NULL , 
	 term_id INTEGER NOT NULL , 
	 rank INTEGER ) ; 
-- 	 PRIMARY KEY ( seqfeature_relationship_id ) , 
-- 	 UNIQUE ( object_seqfeature_id , subject_seqfeature_id , term_id ) ) ; 

-- CREATE INDEX seqfeaturerel_trm ON seqfeature_relationship ( term_id ); 
-- CREATE INDEX seqfeaturerel_child ON seqfeature_relationship ( subject_seqfeature_id ); 
-- CONFIG: you may want to add this if you can't get the optimizer to
-- use the composite index for the initial keys 
--CREATE INDEX seqfeaturerel_parent ON seqfeature_relationship(object_seqfeature_id); 

-- for deep (depth > 1) seqfeature relationship trees we need a transitive 
-- closure table too 
CREATE TABLE seqfeature_path ( 
	 object_seqfeature_id INTEGER NOT NULL , 
	 subject_seqfeature_id INTEGER NOT NULL , 
	 term_id INTEGER NOT NULL , 
	 distance INTEGER ) ;
-- 	 UNIQUE ( object_seqfeature_id , subject_seqfeature_id , term_id , distance ) ) ; 

-- CREATE INDEX seqfeaturepath_trm ON seqfeature_path ( term_id ); 
-- CREATE INDEX seqfeaturepath_child ON seqfeature_path ( subject_seqfeature_id );
-- CONFIG: you may want to add this if you can't get the optimizer to
-- use the composite index for the initial keys 
--CREATE INDEX seqfeaturerel_parent ON seqfeature_path(object_seqfeature_id); 

-- tag/value associations - or ontology annotations 
CREATE TABLE seqfeature_qualifier_value ( 
	 seqfeature_id INTEGER NOT NULL , 
	 term_id INTEGER NOT NULL , 
	 rank INTEGER NOT NULL DEFAULT 0 , 
	 value TEXT NOT NULL ) ; 
-- 	 PRIMARY KEY ( seqfeature_id , term_id , rank ) ) ; 

-- CREATE INDEX seqfeaturequal_trm ON seqfeature_qualifier_value ( term_id ); 

-- DBXrefs for features. This is necessary for genome oriented viewpoints, 
-- where you have a few have long sequences (contigs, or chromosomes) with many
-- features on them. In that case the features are the semantic scope for 
-- their annotation bundles, not the bioentry they are attached to. 
CREATE TABLE seqfeature_dbxref ( 
	 seqfeature_id INTEGER NOT NULL , 
	 dbxref_id INTEGER NOT NULL , 
	 rank INTEGER ) ;
-- 	 PRIMARY KEY ( seqfeature_id , dbxref_id ) ) ; 

-- CREATE INDEX feadblink_dbx ON seqfeature_dbxref ( dbxref_id ); 

-- basically we model everything as potentially having 
-- any number of locations, ie, a split location. SimpleLocations 
-- just have one location. We need to have a location id for the qualifier 
-- associations of fuzzy locations. 
--
-- please do not try to model complex assemblies with this thing. It wont 
-- work. Check out the ensembl schema for this. 
--
-- we allow nulls for start/end - this is useful for fuzzies as 
-- standard range queries will not be included 
--
-- for remote locations, the join to make is to DBXref 
--
-- the FK to term is a possibility to store the type of the 
-- location for determining in one hit whether it's a fuzzy or not 
CREATE SEQUENCE location_pk_seq;
CREATE TABLE location ( 
	 location_id INTEGER DEFAULT nextval ( 'location_pk_seq' ) NOT NULL , 
	 seqfeature_id INTEGER NOT NULL , 
	 dbxref_id INTEGER , 
	 term_id INTEGER , 
	 start_pos INTEGER , 
	 end_pos INTEGER , 
	 strand INTEGER NOT NULL DEFAULT 0 , 
	 rank INTEGER NOT NULL DEFAULT 0 ) ; 
-- 	 PRIMARY KEY ( location_id ) , 
-- 	 UNIQUE ( seqfeature_id , rank ) ) ; 

-- CREATE INDEX seqfeatureloc_start ON location ( start_pos, end_pos ); 
-- CREATE INDEX seqfeatureloc_dbx ON location ( dbxref_id ); 
-- CREATE INDEX seqfeatureloc_trm ON location ( term_id ); 

-- location qualifiers - mainly intended for fuzzies but anything 
-- can go in here 
-- some controlled vocab terms have slots; 
-- fuzzies could be modeled as min_start(5), max_start(5) 
--  
-- there is no restriction on extending the fuzzy ontology 
-- for your own nefarious aims, although the bio* apis will 
-- most likely ignore these 
CREATE TABLE location_qualifier_value ( 
	 location_id INTEGER NOT NULL , 
	 term_id INTEGER NOT NULL , 
	 value VARCHAR ( 255 ) NOT NULL , 
	 int_value INTEGER ) ;
-- 	 PRIMARY KEY ( location_id , term_id ) ) ; 

-- CREATE INDEX locationqual_trm ON location_qualifier_value ( term_id ); 

--
-- This is to solve a problem arising from how transactions are implemented
-- in Postgres as opposed to, e.g., Oracle and InnoDB (MySQL). In short, the
-- difference is that in the latter RDBMSs' implementation, if a particular
-- statement within a transaction fails, the preceding (and possibly
-- subsequent) statements are still valid. On commit, all succeeded statements
-- are committed. In Postgres, the failure of a statement invalidates all
-- preceding statements within the same transaction as well as all subsequent,
-- if any.
--
-- This leads to a problem if you program SQL insert and update statements
-- such that presence of the record you attempt to insert is indicated by
-- failure of the statement due to a unique key constraint violation. Even
-- if your code is prepared to handle the failure by e.g. looking up the
-- record, in the case of Postgres this approach cannot work unless you
-- commit every single statement.
--
-- The bioperl-db adaptor code uses the aforementioned approach and is
-- currently dependent on the following support code. If you are not going
-- to use bioperl-db to populate the database, you may comment out all
-- rules, as then they might add another look-up to one already done on the
-- code that you use and hence add unnecessary overhead.
--

-- CREATE RULE rule_bioentry_i1
--        AS ON INSERT TO bioentry
--        WHERE (
--              SELECT oid FROM bioentry 
--              WHERE identifier     = new.identifier
--              AND   biodatabase_id = new.biodatabase_id
--              ) 
-- 	     IS NOT NULL
--        DO INSTEAD NOTHING
-- ;
-- CREATE RULE rule_bioentry_i2
--        AS ON INSERT TO bioentry
--        WHERE (
--        	     SELECT oid FROM bioentry 
-- 	     WHERE accession      = new.accession
-- 	     AND   biodatabase_id = new.biodatabase_id
-- 	     AND   version	  = new.version
-- 	     ) 
-- 	     IS NOT NULL
--        DO INSTEAD NOTHING
-- ;

-- CREATE RULE rule_biodatabase_i
--        AS ON INSERT TO biodatabase
--        WHERE (SELECT oid FROM biodatabase WHERE name = new.name)
--        	     IS NOT NULL
--        DO INSTEAD NOTHING
-- ;

-- CREATE RULE rule_bioentry_dbxref_i
--        AS ON INSERT TO bioentry_dbxref
--        WHERE (
--        	     SELECT oid FROM bioentry_dbxref 
-- 	     WHERE bioentry_id = new.bioentry_id
-- 	     AND   dbxref_id   = new.dbxref_id
-- 	     ) 
-- 	     IS NOT NULL
--        DO INSTEAD NOTHING
-- ;

-- CREATE RULE rule_bioentry_path_i
--        AS ON INSERT TO bioentry_path
--        WHERE (
--        	     SELECT oid FROM bioentry_relationship 
-- 	     WHERE object_bioentry_id = new.object_bioentry_id
-- 	     AND   subject_bioentry_id= new.subject_bioentry_id
-- 	     AND   term_id	      = new.term_id
-- 	     )
-- 	     IS NOT NULL
--        DO INSTEAD NOTHING
-- ;

-- CREATE RULE rule_bioentry_qualifier_value_i
--        AS ON INSERT TO bioentry_qualifier_value
--        WHERE (
--        	     SELECT oid FROM bioentry_qualifier_value
-- 	     WHERE bioentry_id = new.bioentry_id
-- 	     AND   term_id     = new.term_id
-- 	     AND   rank	       = new.rank
-- 	     )
-- 	     IS NOT NULL
--        DO INSTEAD NOTHING
-- ;

-- CREATE RULE rule_bioentry_reference_i
--        AS ON INSERT TO bioentry_reference
--        WHERE (
--        	     SELECT oid FROM bioentry_reference 
-- 	     WHERE bioentry_id  = new.bioentry_id
-- 	     AND   reference_id = new.reference_id
-- 	     AND   rank		= new.rank
-- 	     )
-- 	     IS NOT NULL
--        DO INSTEAD NOTHING
-- ;

-- CREATE RULE rule_bioentry_relationship_i
--        AS ON INSERT TO bioentry_relationship
--        WHERE (
--        	     SELECT oid FROM bioentry_relationship 
-- 	     WHERE object_bioentry_id = new.object_bioentry_id
-- 	     AND   subject_bioentry_id= new.subject_bioentry_id
-- 	     AND   term_id	      = new.term_id
-- 	     )
-- 	     IS NOT NULL
--        DO INSTEAD NOTHING
-- ;

-- CREATE RULE rule_biosequence_i
--        AS ON INSERT TO biosequence
--        WHERE (SELECT oid FROM biosequence WHERE bioentry_id = new.bioentry_id)
--        	     IS NOT NULL
--        DO INSTEAD NOTHING
-- ;

-- CREATE RULE rule_comment_i
--        AS ON INSERT TO comment
--        WHERE (
--        	     SELECT oid FROM comment
-- 	     WHERE bioentry_id = new.bioentry_id
-- 	     AND   rank	       = new.rank
-- 	     )
-- 	     IS NOT NULL
--        DO INSTEAD NOTHING
-- ;

-- CREATE RULE rule_dbxref_i
--        AS ON INSERT TO dbxref
--        WHERE (
--        	     SELECT oid FROM dbxref
-- 	     WHERE accession = new.accession
-- 	     AND   dbname    = new.dbname
-- 	     AND   version   = new.version
-- 	     )
-- 	     IS NOT NULL
--        DO INSTEAD NOTHING
-- ;

-- CREATE RULE rule_dbxref_qualifier_value_i
--        AS ON INSERT TO dbxref_qualifier_value
--        WHERE (
--        	     SELECT oid FROM dbxref_qualifier_value
-- 	     WHERE dbxref_id = new.dbxref_id
-- 	     AND   term_id   = new.term_id
-- 	     AND   rank	     = new.rank
-- 	     )
-- 	     IS NOT NULL
--        DO INSTEAD NOTHING
-- ;

-- CREATE RULE rule_location_i
--        AS ON INSERT TO location
--        WHERE (
--        	     SELECT oid FROM location
-- 	     WHERE seqfeature_id = new.seqfeature_id
-- 	     AND   rank		 = new.rank
-- 	     )
-- 	     IS NOT NULL
--        DO INSTEAD NOTHING
-- ;

-- CREATE RULE rule_location_qualifier_value_i
--        AS ON INSERT TO location_qualifier_value
--        WHERE (
--        	     SELECT oid FROM location_qualifier_value
-- 	     WHERE location_id = new.location_id
-- 	     AND   term_id     = new.term_id
-- 	     ) 
--        	     IS NOT NULL
--        DO INSTEAD NOTHING
-- ;

-- CREATE RULE rule_ontology_i
--        AS ON INSERT TO ontology
--        WHERE (SELECT oid FROM ontology WHERE name = new.name) 
--        	     IS NOT NULL
--        DO INSTEAD NOTHING
-- ;

-- CREATE RULE rule_reference_i1
--        AS ON INSERT TO reference
--        WHERE (SELECT oid FROM reference WHERE crc = new.crc) 
--        	     IS NOT NULL
--        DO INSTEAD NOTHING
-- ;
-- CREATE RULE rule_reference_i2
--        AS ON INSERT TO reference
--        WHERE (SELECT oid FROM reference WHERE dbxref_id = new.dbxref_id)
--        	     IS NOT NULL
--        DO INSTEAD NOTHING
-- ;

-- CREATE RULE rule_seqfeature_i
--        AS ON INSERT TO seqfeature
--        WHERE (
--        	     SELECT oid FROM seqfeature 
-- 	     WHERE bioentry_id    = new.bioentry_id
-- 	     AND   type_term_id   = new.type_term_id
-- 	     AND   source_term_id = new.source_term_id
-- 	     AND   rank		  = new.rank
-- 	     )
-- 	     IS NOT NULL
--        DO INSTEAD NOTHING
-- ;

-- CREATE RULE rule_seqfeature_dbxref_i
--        AS ON INSERT TO seqfeature_dbxref
--        WHERE (	    
--        	     SELECT oid FROM seqfeature_dbxref
-- 	     WHERE seqfeature_id = new.seqfeature_id
-- 	     AND   dbxref_id	 = new.dbxref_id
-- 	     )
-- 	     IS NOT NULL
--        DO INSTEAD NOTHING
-- ;

-- CREATE RULE rule_seqfeature_path_i
--        AS ON INSERT TO seqfeature_path
--        WHERE (
--        	     SELECT oid FROM seqfeature_path
-- 	     WHERE object_seqfeature_id = new.object_seqfeature_id
-- 	     AND   subject_seqfeature_id= new.subject_seqfeature_id
-- 	     AND   term_id		= new.term_id
-- 	     )
-- 	     IS NOT NULL
--        DO INSTEAD NOTHING
-- ;

-- CREATE RULE rule_seqfeature_qualifier_value_i
--        AS ON INSERT TO seqfeature_qualifier_value
--        WHERE (
--        	     SELECT oid FROM seqfeature_qualifier_value
-- 	     WHERE seqfeature_id = new.seqfeature_id
-- 	     AND   term_id	 = new.term_id
-- 	     AND   rank		 = new.rank
-- 	     )
-- 	     IS NOT NULL
--        DO INSTEAD NOTHING
-- ;

-- CREATE RULE rule_seqfeature_relationship_i
--        AS ON INSERT TO seqfeature_relationship
--        WHERE (
--        	     SELECT oid FROM seqfeature_relationship
-- 	     WHERE object_seqfeature_id = new.object_seqfeature_id
-- 	     AND   subject_seqfeature_id= new.subject_seqfeature_id
-- 	     AND   term_id		= new.term_id
-- 	     )
-- 	     IS NOT NULL
--        DO INSTEAD NOTHING
-- ;

CREATE RULE rule_taxon_i
       AS ON INSERT TO taxon
       WHERE (SELECT oid FROM taxon WHERE ncbi_taxon_id = new.ncbi_taxon_id)
       	     IS NOT NULL
       DO INSTEAD NOTHING
;

CREATE RULE rule_taxon_name_i
       AS ON INSERT TO taxon_name
       WHERE (
       	     SELECT oid FROM taxon_name
	     WHERE taxon_id   = new.taxon_id
	     AND   name	      = new.name
	     AND   name_class = new.name_class
	     ) 
	     IS NOT NULL
       DO INSTEAD NOTHING
;

-- CREATE RULE rule_term_i1
--        AS ON INSERT TO term
--        WHERE (SELECT oid FROM term WHERE identifier = new.identifier)
--        	     IS NOT NULL
--        DO INSTEAD NOTHING
-- ;
-- CREATE RULE rule_term_i2
--        AS ON INSERT TO term
--        WHERE (
--        	     SELECT oid FROM term 
-- 	     WHERE name        = new.name
-- 	     AND   ontology_id = new.ontology_id
--              AND   is_obsolete = new.is_obsolete
-- 	     )
--        	     IS NOT NULL
--        DO INSTEAD NOTHING
-- ;

-- CREATE RULE rule_term_dbxref_i
--        AS ON INSERT TO term_dbxref
--        WHERE (
--        	     SELECT oid FROM term_dbxref
-- 	     WHERE dbxref_id = new.dbxref_id
-- 	     AND   term_id   = new.term_id
-- 	     )
-- 	     IS NOT NULL
--        DO INSTEAD NOTHING
-- ;

-- CREATE RULE rule_term_path_i
--        AS ON INSERT TO term_path
--        WHERE (
--        	     SELECT oid FROM term_path
-- 	     WHERE subject_term_id   = new.subject_term_id
-- 	     AND   predicate_term_id = new.predicate_term_id
-- 	     AND   object_term_id    = new.object_term_id
-- 	     AND   ontology_id	     = new.ontology_id
-- 	     AND   distance	     = new.distance
-- 	     )
-- 	     IS NOT NULL
--        DO INSTEAD NOTHING
-- ;

-- CREATE RULE rule_term_relationship_i
--        AS ON INSERT TO term_relationship
--        WHERE (
--        	     SELECT oid FROM term_relationship
-- 	     WHERE subject_term_id   = new.subject_term_id
-- 	     AND   predicate_term_id = new.predicate_term_id
-- 	     AND   object_term_id    = new.object_term_id
-- 	     AND   ontology_id	     = new.ontology_id
-- 	     )
-- 	     IS NOT NULL
--        DO INSTEAD NOTHING
-- ;

-- CREATE RULE rule_term_relationship_term_i1
--        AS ON INSERT TO term_relationship_term
--        WHERE (
--        	     SELECT oid FROM term_relationship_term
-- 	     WHERE term_relationship_id   = new.term_relationship_id
-- 	     )
-- 	     IS NOT NULL
--        DO INSTEAD NOTHING
-- ;

-- CREATE RULE rule_term_relationship_term_i2
--        AS ON INSERT TO term_relationship_term
--        WHERE (
--        	     SELECT oid FROM term_relationship_term
-- 	     WHERE term_id   = new.term_id
-- 	     )
-- 	     IS NOT NULL
--        DO INSTEAD NOTHING
-- ;

-- CREATE RULE rule_term_synonym_i
--        AS ON INSERT TO term_synonym
--        WHERE (
--        	     SELECT oid FROM term_synonym
-- 	     WHERE synonym = new.synonym
-- 	     AND   term_id = new.term_id
-- 	     )
-- 	     IS NOT NULL
--        DO INSTEAD NOTHING
-- ;

-- --
-- -- Functions that may be used as an API by applications, e.g. load scripts etc.
-- -- 

-- this is used by load_ncbi_taxonomy.pl to speed up loading into the taxon
-- table by 1 to 2 orders of magnitude
CREATE OR REPLACE FUNCTION unconstrain_taxon ()
RETURNS INTEGER
AS
'
DROP RULE rule_taxon_i ON taxon;
SELECT 1;
'
LANGUAGE SQL
VOLATILE STRICT SECURITY DEFINER
;

-- this function re-establishes what unconstrain_taxon() removed temporarily
CREATE OR REPLACE FUNCTION constrain_taxon ()
RETURNS INTEGER
AS
'
CREATE RULE rule_taxon_i
       AS ON INSERT TO taxon
       WHERE (SELECT oid FROM taxon WHERE ncbi_taxon_id = new.ncbi_taxon_id)
       	     IS NOT NULL
       DO INSTEAD NOTHING
;
SELECT 1;
'
LANGUAGE SQL
VOLATILE STRICT SECURITY DEFINER
;
