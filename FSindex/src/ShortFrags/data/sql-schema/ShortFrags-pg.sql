-- Table of hits. We don't keep the sequence: this should be
-- retrieved from the BioSQL schema. The search_id refers to
-- the searches table.

CREATE TABLE hits (
	experiment_id SMALLINT NOT NULL ,
	len_iter_id SMALLINT NOT NULL ,
	frag_id INTEGER NOT NULL ,
	accession VARCHAR ( 40 ) NOT NULL , 
	start INTEGER NOT NULL ,
	distance SMALLINT ,
	similarity SMALLINT ,
	pvalue REAL ,
	Evalue REAL ,
	PRIMARY KEY ( experiment_id, len_iter_id, frag_id, accession, start ) ); 

CREATE TABLE searches (
	experiment_id SMALLINT NOT NULL ,
	len_iter_id SMALLINT NOT NULL ,
	frag_id INTEGER NOT NULL ,
	query_frag VARCHAR ( 32 ) NOT NULL ,
	score_matrix BYTEA ,
	conv_type SMALLINT ,
	sim_range SMALLINT ,
	dist_range SMALLINT ,
	kNN SMALLINT ,
	num_hits INTEGER NOT NULL ,
	PRIMARY KEY ( experiment_id, len_iter_id, frag_id ) ); 

CREATE SEQUENCE experiments_pk_seq;

CREATE TABLE experiments (
	experiment_id SMALLINT DEFAULT nextval ( 'experiments_pk_seq' ) NOT NULL ,
	name VARCHAR ( 40 ) NOT NULL ,
	description TEXT ,
	query_sequence TEXT NOT NULL,
	query_description TEXT ,
	min_len INTEGER NOT NULL ,
	max_len INTEGER NOT NULL ,
	PRIMARY KEY ( experiment_id ) );

CREATE INDEX hits_acc on hits ( accession ); 

-- ALTER TABLE seqfeature_qualifier_value ADD CONSTRAINT FKseqfeature_featqual
--      FOREIGN KEY ( seqfeature_id ) REFERENCES seqfeature ( seqfeature_id )
--      ON DELETE CASCADE ;
