ALTER TABLE biodatabase ADD PRIMARY KEY ( biodatabase_id ) ;
ALTER TABLE biodatabase ADD UNIQUE ( name ) ;
CREATE INDEX db_auth on biodatabase ( authority );


ALTER TABLE taxon ADD PRIMARY KEY ( taxon_id ) ;
ALTER TABLE taxon ADD CONSTRAINT XAKtaxon_ncbi_taxon_id UNIQUE ( ncbi_taxon_id ) ;
ALTER TABLE taxon ADD CONSTRAINT XAKtaxon_left_value UNIQUE ( left_value ) ;
ALTER TABLE taxon ADD CONSTRAINT XAKtaxon_right_value UNIQUE ( right_value ) ;
CREATE INDEX taxparent ON taxon ( parent_taxon_id );


ALTER TABLE taxon_name ADD UNIQUE ( name , name_class, taxon_id ) ;
CREATE INDEX taxnametaxonid ON taxon_name ( taxon_id );
CREATE INDEX taxnamename ON taxon_name ( name );


ALTER TABLE ontology ADD PRIMARY KEY ( ontology_id ) ;
ALTER TABLE ontology ADD UNIQUE ( name ) ;


ALTER TABLE term ADD PRIMARY KEY ( term_id ) ;
ALTER TABLE term ADD UNIQUE (name, ontology_id) ;
-- ALTER TABLE term ADD UNIQUE (name, ontology_id, is_obsolete);
ALTER TABLE term ADD UNIQUE ( identifier ) ;
CREATE INDEX term_ont ON term ( ontology_id );


ALTER TABLE term_synonym ADD PRIMARY KEY ( term_id , dbxref_id ) ;


ALTER TABLE term_dbxref ADD PRIMARY KEY ( term_id , synonym ) ;
CREATE INDEX trmdbxref_dbxrefid ON term_dbxref ( dbxref_id );


ALTER TABLE term_relationship ADD PRIMARY KEY ( term_relationship_id ) ;
ALTER TABLE term_relationship ADD UNIQUE ( subject_term_id , predicate_term_id , object_term_id , ontology_id ) ;
CREATE INDEX trmrel_predicateid ON term_relationship ( predicate_term_id );
CREATE INDEX trmrel_objectid ON term_relationship ( object_term_id );
CREATE INDEX trmrel_ontid ON term_relationship ( ontology_id );
--CREATE INDEX trmrel_subjectid ON term_relationship(subject_term_id);


ALTER TABLE term_relationship_term ADD PRIMARY KEY ( term_relationship_id );
ALTER TABLE term_relationship_term ADD UNIQUE ( term_id );


ALTER TABLE term_path ADD PRIMARY KEY ( term_path_id );
ALTER TABLE term_path ADD UNIQUE ( subject_term_id , predicate_term_id , object_term_id , ontology_id , distance );
CREATE INDEX trmpath_predicateid ON term_path ( predicate_term_id );
CREATE INDEX trmpath_objectid ON term_path ( object_term_id );
CREATE INDEX trmpath_ontid ON term_path ( ontology_id );
--CREATE INDEX trmpath_subjectid ON term_path(subject_term_id);


ALTER TABLE bioentry ADD PRIMARY KEY ( bioentry_id ) ;
ALTER TABLE bioentry ADD UNIQUE ( accession , biodatabase_id , version ) ;
ALTER TABLE bioentry ADD UNIQUE ( identifier ) ;
-- ALTER TABLE bioentry ADD UNIQUE ( identifier , biodatabase_id );
CREATE INDEX bioentry_name ON bioentry ( name );
CREATE INDEX bioentry_db ON bioentry ( biodatabase_id );
CREATE INDEX bioentry_tax ON bioentry ( taxon_id );


ALTER TABLE bioentry_relationship ADD PRIMARY KEY ( bioentry_relationship_id );
ALTER TABLE bioentry_relationship ADD UNIQUE ( object_bioentry_id , subject_bioentry_id , term_id );
CREATE INDEX bioentryrel_trm ON bioentry_relationship ( term_id );
CREATE INDEX bioentryrel_child ON bioentry_relationship (subject_bioentry_id);
-- CREATE INDEX bioentryrel_parent ON bioentry_relationship(object_bioentry_id);


ALTER TABLE bioentry_path ADD UNIQUE ( object_bioentry_id , subject_bioentry_id , term_id , distance );
CREATE INDEX bioentrypath_trm ON bioentry_path ( term_id );
CREATE INDEX bioentrypath_child ON bioentry_path ( subject_bioentry_id );


ALTER TABLE biosequence ADD PRIMARY KEY ( bioentry_id ) ;
-- ALTER TABLE biosequence ADD COLUMN ( isoelec_pt NUMERIC(4,2) );
-- ALTER TABLE biosequence ADD COLUMN (	mol_wgt DOUBLE PRECISION );
-- ALTER TABLE biosequence ADD COLUMN ( perc_gc DOUBLE PRECISION );


ALTER TABLE dbxref ADD PRIMARY KEY ( dbxref_id );
ALTER TABLE dbxref ADD UNIQUE ( accession , dbname , version );
CREATE INDEX dbxref_db ON dbxref ( dbname );


ALTER TABLE dbxref_qualifier_value ADD PRIMARY KEY ( dbxref_id , term_id , rank );
CREATE INDEX dbxrefqual_dbx ON dbxref_qualifier_value ( dbxref_id );
CREATE INDEX dbxrefqual_trm ON dbxref_qualifier_value ( term_id );


ALTER TABLE bioentry_dbxref ADD PRIMARY KEY ( bioentry_id , dbxref_id );
CREATE INDEX dblink_dbx ON bioentry_dbxref ( dbxref_id );


ALTER TABLE reference ADD PRIMARY KEY ( reference_id );
ALTER TABLE reference ADD UNIQUE ( dbxref_id );
ALTER TABLE reference ADD UNIQUE ( crc );


ALTER TABLE bioentry_reference ADD PRIMARY KEY ( bioentry_id , reference_id , rank );
CREATE INDEX bioentryref_ref ON bioentry_reference ( reference_id );


ALTER TABLE comment ADD PRIMARY KEY ( comment_id );
ALTER TABLE comment ADD UNIQUE ( bioentry_id , rank );


ALTER TABLE bioentry_qualifier_value ADD UNIQUE ( bioentry_id , term_id , rank );
CREATE INDEX bioentryqual_trm ON bioentry_qualifier_value ( term_id );


ALTER TABLE seqfeature ADD PRIMARY KEY ( seqfeature_id );
ALTER TABLE seqfeature ADD UNIQUE ( bioentry_id , type_term_id , source_term_id , rank ) ;
CREATE INDEX seqfeature_trm ON seqfeature ( type_term_id );
CREATE INDEX seqfeature_fsrc ON seqfeature ( source_term_id );
-- CREATE INDEX seqfeature_bioentryid ON seqfeature(bioentry_id);


ALTER TABLE seqfeature_relationship ADD PRIMARY KEY ( seqfeature_relationship_id );
ALTER TABLE seqfeature_relationship ADD UNIQUE ( object_seqfeature_id , subject_seqfeature_id , term_id );
CREATE INDEX seqfeaturerel_trm ON seqfeature_relationship ( term_id );
CREATE INDEX seqfeaturerel_child ON seqfeature_relationship ( subject_seqfeature_id );
-- CREATE INDEX seqfeaturerel_parent ON seqfeature_relationship(object_seqfeature_id);

ALTER TABLE seqfeature_path ADD UNIQUE ( object_seqfeature_id , subject_seqfeature_id , term_id , distance );
CREATE INDEX seqfeaturepath_trm ON seqfeature_path ( term_id );
CREATE INDEX seqfeaturepath_child ON seqfeature_path ( subject_seqfeature_id );
-- CREATE INDEX seqfeaturerel_parent ON seqfeature_path(object_seqfeature_id);


ALTER TABLE seqfeature_qualifier_value ADD PRIMARY KEY ( seqfeature_id , term_id , rank ) ;
CREATE INDEX seqfeaturequal_trm ON seqfeature_qualifier_value ( term_id );


ALTER TABLE seqfeature_dbxref ADD PRIMARY KEY ( seqfeature_id , dbxref_id );
CREATE INDEX feadblink_dbx ON seqfeature_dbxref ( dbxref_id );


ALTER TABLE location ADD PRIMARY KEY ( location_id ) ;
ALTER TABLE location ADD UNIQUE ( seqfeature_id , rank ) ;
CREATE INDEX seqfeatureloc_start ON location ( start_pos, end_pos );
CREATE INDEX seqfeatureloc_dbx ON location ( dbxref_id );
CREATE INDEX seqfeatureloc_trm ON location ( term_id );


ALTER TABLE location_qualifier_value ADD PRIMARY KEY ( location_id , term_id );
CREATE INDEX locationqual_trm ON location_qualifier_value ( term_id );


ANALYZE;
