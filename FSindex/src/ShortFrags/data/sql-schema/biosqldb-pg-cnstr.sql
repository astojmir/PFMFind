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
ALTER TABLE term ADD UNIQUE ( identifier ) ; 
CREATE INDEX term_ont ON term ( ontology_id ); 

ALTER TABLE bioentry ADD PRIMARY KEY ( bioentry_id ) ; 
ALTER TABLE bioentry ADD UNIQUE ( accession , biodatabase_id , version ) ;
ALTER TABLE bioentry ADD UNIQUE ( identifier ) ;
CREATE INDEX bioentry_name ON bioentry ( name ); 
CREATE INDEX bioentry_db ON bioentry ( biodatabase_id ); 
CREATE INDEX bioentry_tax ON bioentry ( taxon_id ); 

ALTER TABLE biosequence ADD PRIMARY KEY ( bioentry_id ) ; 

ALTER TABLE seqfeature ADD PRIMARY KEY ( seqfeature_id );
ALTER TABLE seqfeature ADD UNIQUE ( bioentry_id , type_term_id , source_term_id , rank ) ; 
CREATE INDEX seqfeature_trm ON seqfeature ( type_term_id ); 
CREATE INDEX seqfeature_fsrc ON seqfeature ( source_term_id ); 

ALTER TABLE seqfeature_qualifier_value ADD PRIMARY KEY ( seqfeature_id , term_id , rank ) ; 
CREATE INDEX seqfeaturequal_trm ON seqfeature_qualifier_value ( term_id ); 

ALTER TABLE location ADD PRIMARY KEY ( location_id ) ;
ALTER TABLE location ADD UNIQUE ( seqfeature_id , rank ) ; 
CREATE INDEX seqfeatureloc_start ON location ( start_pos, end_pos ); 
CREATE INDEX seqfeatureloc_dbx ON location ( dbxref_id ); 
CREATE INDEX seqfeatureloc_trm ON location ( term_id ); 

ANALYZE;
