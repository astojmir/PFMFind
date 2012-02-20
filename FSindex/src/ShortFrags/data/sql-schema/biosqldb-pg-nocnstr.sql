CREATE SEQUENCE biodatabase_pk_seq;
CREATE TABLE biodatabase (
	 biodatabase_id INTEGER DEFAULT nextval ( 'biodatabase_pk_seq' ) NOT NULL ,
	 name VARCHAR ( 128 ) NOT NULL ,
	 authority VARCHAR ( 128 ) ,
	 description TEXT ) ;


CREATE SEQUENCE taxon_pk_seq;
CREATE TABLE taxon (
	 taxon_id INTEGER DEFAULT nextval ( 'taxon_pk_seq' ) NOT NULL ,
	 ncbi_taxon_id INTEGER ,
	 parent_taxon_id INTEGER ,
	 node_rank VARCHAR ( 32 ) ,
	 genetic_code SMALLINT ,
	 mito_genetic_code SMALLINT ,
	 left_value INTEGER ,
	 right_value INTEGER );


CREATE TABLE taxon_name (
	 taxon_id INTEGER NOT NULL ,
	 name VARCHAR ( 255 ) NOT NULL ,
	 name_class VARCHAR ( 32 ) NOT NULL ) ;


CREATE SEQUENCE ontology_pk_seq;
CREATE TABLE ontology (
	 ontology_id INTEGER DEFAULT nextval ( 'ontology_pk_seq' ) NOT NULL ,
	 name VARCHAR ( 32 ) NOT NULL ,
	 definition TEXT );


CREATE SEQUENCE term_pk_seq;
CREATE TABLE term (
	 term_id INTEGER DEFAULT nextval ( 'term_pk_seq' ) NOT NULL ,
	 name VARCHAR ( 255 ) NOT NULL ,
	 definition TEXT ,
	 identifier VARCHAR ( 40 ) ,
	 is_obsolete CHAR ( 1 ) ,
	 ontology_id INTEGER NOT NULL ) ;


CREATE TABLE term_synonym (
	 synonym VARCHAR(255) NOT NULL,
	 term_id INTEGER NOT NULL );


CREATE TABLE term_dbxref (
	 term_id INTEGER NOT NULL ,
	 dbxref_id INTEGER NOT NULL ,
	 rank INTEGER );


CREATE SEQUENCE term_relationship_pk_seq;
CREATE TABLE term_relationship (
	 term_relationship_id INTEGER DEFAULT nextval ( 'term_relationship_pk_seq' ) NOT NULL ,
	 subject_term_id INTEGER NOT NULL ,
	 predicate_term_id INTEGER NOT NULL ,
	 object_term_id INTEGER NOT NULL ,
	 ontology_id INTEGER NOT NULL );


CREATE TABLE term_relationship_term (
        term_relationship_id INTEGER NOT NULL,
        term_id              INTEGER NOT NULL,
        PRIMARY KEY ( term_relationship_id ),
        UNIQUE ( term_id ) );


CREATE SEQUENCE term_path_pk_seq;
CREATE TABLE term_path (
         term_path_id INTEGER DEFAULT nextval ( 'term_path_pk_seq' ) NOT NULL ,
	 subject_term_id INTEGER NOT NULL ,
	 predicate_term_id INTEGER NOT NULL ,
	 object_term_id INTEGER NOT NULL ,
	 ontology_id INTEGER NOT NULL ,
	 distance INTEGER );


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


CREATE SEQUENCE bioentry_relationship_pk_seq;
CREATE TABLE bioentry_relationship (
	 bioentry_relationship_id INTEGER DEFAULT nextval ( 'bioentry_relationship_pk_seq' ) NOT NULL ,
	 object_bioentry_id INTEGER NOT NULL ,
	 subject_bioentry_id INTEGER NOT NULL ,
	 term_id INTEGER NOT NULL ,
	 rank INTEGER );


CREATE TABLE bioentry_path (
	 object_bioentry_id INTEGER NOT NULL ,
	 subject_bioentry_id INTEGER NOT NULL ,
	 term_id INTEGER NOT NULL ,
	 distance INTEGER );


CREATE TABLE biosequence (
	 bioentry_id INTEGER NOT NULL ,
	 version INTEGER ,
	 length INTEGER ,
	 alphabet VARCHAR ( 10 ) ,
	 seq TEXT ) ;


CREATE SEQUENCE dbxref_pk_seq;
CREATE TABLE dbxref (
	 dbxref_id INTEGER DEFAULT nextval ( 'dbxref_pk_seq' ) NOT NULL ,
	 dbname VARCHAR ( 40 ) NOT NULL ,
	 accession VARCHAR ( 40 ) NOT NULL ,
	 version INTEGER NOT NULL );


CREATE TABLE dbxref_qualifier_value (
	 dbxref_id INTEGER NOT NULL ,
	 term_id INTEGER NOT NULL ,
	 rank INTEGER NOT NULL DEFAULT 0 ,
	 value TEXT );


CREATE TABLE bioentry_dbxref (
	 bioentry_id INTEGER NOT NULL ,
	 dbxref_id INTEGER NOT NULL ,
	 rank INTEGER );


CREATE SEQUENCE reference_pk_seq;
CREATE TABLE reference (
	 reference_id INTEGER DEFAULT nextval ( 'reference_pk_seq' ) NOT NULL ,
	 dbxref_id INTEGER ,
	 location TEXT NOT NULL ,
	 title TEXT ,
	 authors TEXT ,
	 crc VARCHAR ( 32 ));


CREATE TABLE bioentry_reference (
	 bioentry_id INTEGER NOT NULL ,
	 reference_id INTEGER NOT NULL ,
	 start_pos INTEGER ,
	 end_pos INTEGER ,
	 rank INTEGER NOT NULL DEFAULT 0 );


CREATE SEQUENCE comment_pk_seq;
CREATE TABLE comment (
	 comment_id INTEGER DEFAULT nextval ( 'comment_pk_seq' ) NOT NULL ,
	 bioentry_id INTEGER NOT NULL ,
	 comment_text TEXT NOT NULL ,
	 rank INTEGER NOT NULL DEFAULT 0);


CREATE TABLE bioentry_qualifier_value (
	 bioentry_id INTEGER NOT NULL ,
	 term_id INTEGER NOT NULL ,
	 value TEXT ,
	 rank INTEGER NOT NULL DEFAULT 0 );


CREATE SEQUENCE seqfeature_pk_seq;
CREATE TABLE seqfeature (
	 seqfeature_id INTEGER DEFAULT nextval ( 'seqfeature_pk_seq' ) NOT NULL ,
	 bioentry_id INTEGER NOT NULL ,
	 type_term_id INTEGER NOT NULL ,
	 source_term_id INTEGER NOT NULL ,
	 display_name VARCHAR ( 64 ) ,
	 rank INTEGER NOT NULL DEFAULT 0 ) ;


CREATE SEQUENCE seqfeature_relationship_pk_seq;
CREATE TABLE seqfeature_relationship (
	 seqfeature_relationship_id INTEGER DEFAULT nextval ( 'seqfeature_relationship_pk_seq' ) NOT NULL ,
	 object_seqfeature_id INTEGER NOT NULL ,
	 subject_seqfeature_id INTEGER NOT NULL ,
	 term_id INTEGER NOT NULL ,
	 rank INTEGER );


CREATE TABLE seqfeature_path (
	 object_seqfeature_id INTEGER NOT NULL ,
	 subject_seqfeature_id INTEGER NOT NULL ,
	 term_id INTEGER NOT NULL ,
	 distance INTEGER );


CREATE TABLE seqfeature_qualifier_value (
	 seqfeature_id INTEGER NOT NULL ,
	 term_id INTEGER NOT NULL ,
	 rank INTEGER NOT NULL DEFAULT 0 ,
	 value TEXT NOT NULL ) ;


CREATE TABLE seqfeature_dbxref (
	 seqfeature_id INTEGER NOT NULL ,
	 dbxref_id INTEGER NOT NULL ,
	 rank INTEGER );


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


CREATE TABLE location_qualifier_value (
	 location_id INTEGER NOT NULL ,
	 term_id INTEGER NOT NULL ,
	 value VARCHAR ( 255 ) NOT NULL ,
	 int_value INTEGER);

