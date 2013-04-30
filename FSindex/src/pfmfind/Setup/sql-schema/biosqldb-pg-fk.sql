-- 
-- Create the foreign key constraints 
-- 

-- ontology term
ALTER TABLE term ADD CONSTRAINT FKont_term
      FOREIGN KEY ( ontology_id ) REFERENCES ontology ( ontology_id ) 
      ON DELETE CASCADE ;

-- term synonyms
ALTER TABLE term_synonym ADD CONSTRAINT FKterm_syn
      FOREIGN KEY ( term_id ) REFERENCES term ( term_id )
      ON DELETE CASCADE ;

-- term_dbxref 
ALTER TABLE term_dbxref ADD CONSTRAINT FKdbxref_trmdbxref
      FOREIGN KEY ( dbxref_id ) REFERENCES dbxref ( dbxref_id )
      ON DELETE CASCADE ;
ALTER TABLE term_dbxref ADD CONSTRAINT FKterm_trmdbxref
      FOREIGN KEY ( term_id ) REFERENCES term ( term_id )
      ON DELETE CASCADE ;

-- term_relationship 
ALTER TABLE term_relationship ADD CONSTRAINT FKtrmsubject_trmrel
      FOREIGN KEY ( subject_term_id ) REFERENCES term ( term_id )
      ON DELETE CASCADE ;
ALTER TABLE term_relationship ADD CONSTRAINT FKtrmpredicate_trmrel
      FOREIGN KEY ( predicate_term_id ) REFERENCES term ( term_id )
      ON DELETE CASCADE ;
ALTER TABLE term_relationship ADD CONSTRAINT FKtrmobject_trmrel
      FOREIGN KEY ( object_term_id ) REFERENCES term ( term_id )
      ON DELETE CASCADE ;
ALTER TABLE term_relationship ADD CONSTRAINT FKontology_trmrel
      FOREIGN KEY ( ontology_id ) REFERENCES ontology ( ontology_id )
      ON DELETE CASCADE ;

-- term_relationship_term
ALTER TABLE term_relationship_term ADD CONSTRAINT FKtrmrel_trmreltrm
      FOREIGN KEY (term_relationship_id) REFERENCES term_relationship(term_relationship_id)
      ON DELETE CASCADE ;
ALTER TABLE term_relationship_term ADD CONSTRAINT FKtrm_trmreltrm
      FOREIGN KEY (term_id) REFERENCES term(term_id)
      ON DELETE CASCADE ;

-- term_path 
ALTER TABLE term_path ADD CONSTRAINT FKtrmsubject_trmpath
      FOREIGN KEY ( subject_term_id ) REFERENCES term ( term_id )
      ON DELETE CASCADE ;
ALTER TABLE term_path ADD CONSTRAINT FKtrmpredicate_trmpath
      FOREIGN KEY ( predicate_term_id ) REFERENCES term ( term_id )
      ON DELETE CASCADE ;
ALTER TABLE term_path ADD CONSTRAINT FKtrmobject_trmpath
      FOREIGN KEY ( object_term_id ) REFERENCES term ( term_id )
      ON DELETE CASCADE ;
ALTER TABLE term_path ADD CONSTRAINT FKontology_trmpath
      FOREIGN KEY ( ontology_id ) REFERENCES ontology ( ontology_id )
      ON DELETE CASCADE ;

-- taxon, taxon_name 
-- unfortunately, we can't constrain parent_taxon_id as it is violated
-- occasionally by the downloads available from NCBI
ALTER TABLE taxon ADD CONSTRAINT FKtaxon_taxon
      FOREIGN KEY ( parent_taxon_id ) REFERENCES taxon ( taxon_id )
      DEFERRABLE;
ALTER TABLE taxon_name ADD CONSTRAINT FKtaxon_taxonname
      FOREIGN KEY ( taxon_id ) REFERENCES taxon ( taxon_id )
      ON DELETE CASCADE ;

-- bioentry 
ALTER TABLE bioentry ADD CONSTRAINT FKtaxon_bioentry
      FOREIGN KEY ( taxon_id ) REFERENCES taxon ( taxon_id ) ; 
ALTER TABLE bioentry ADD CONSTRAINT FKbiodatabase_bioentry
      FOREIGN KEY ( biodatabase_id ) REFERENCES biodatabase ( biodatabase_id ) ; 
-- bioentry_relationship 
ALTER TABLE bioentry_relationship ADD CONSTRAINT FKterm_bioentryrel
      FOREIGN KEY ( term_id ) REFERENCES term ( term_id ) ; 
ALTER TABLE bioentry_relationship ADD CONSTRAINT FKparentent_bioentryrel
      FOREIGN KEY ( object_bioentry_id ) REFERENCES bioentry ( bioentry_id )
      ON DELETE CASCADE ;
ALTER TABLE bioentry_relationship ADD CONSTRAINT FKchildent_bioentryrel
      FOREIGN KEY ( subject_bioentry_id ) REFERENCES bioentry ( bioentry_id )
      ON DELETE CASCADE ;

-- bioentry_path 
ALTER TABLE bioentry_path ADD CONSTRAINT FKterm_bioentrypath
      FOREIGN KEY ( term_id ) REFERENCES term ( term_id ) ; 
ALTER TABLE bioentry_path ADD CONSTRAINT FKparentent_bioentrypath
      FOREIGN KEY ( object_bioentry_id ) REFERENCES bioentry ( bioentry_id )
      ON DELETE CASCADE ;
ALTER TABLE bioentry_path ADD CONSTRAINT FKchildent_bioentrypath
      FOREIGN KEY ( subject_bioentry_id ) REFERENCES bioentry ( bioentry_id )
      ON DELETE CASCADE ;

-- biosequence 
ALTER TABLE biosequence ADD CONSTRAINT FKbioentry_bioseq
      FOREIGN KEY ( bioentry_id ) REFERENCES bioentry ( bioentry_id )
      ON DELETE CASCADE ;

-- comment 
ALTER TABLE comment ADD CONSTRAINT FKbioentry_comment
      FOREIGN KEY ( bioentry_id ) REFERENCES bioentry ( bioentry_id )
      ON DELETE CASCADE ;

-- bioentry_dbxref 
ALTER TABLE bioentry_dbxref ADD CONSTRAINT FKbioentry_dblink
      FOREIGN KEY ( bioentry_id ) REFERENCES bioentry ( bioentry_id )
      ON DELETE CASCADE ;
ALTER TABLE bioentry_dbxref ADD CONSTRAINT FKdbxref_dblink
      FOREIGN KEY ( dbxref_id ) REFERENCES dbxref ( dbxref_id )
      ON DELETE CASCADE ;

-- dbxref_qualifier_value 
ALTER TABLE dbxref_qualifier_value ADD CONSTRAINT FKtrm_dbxrefqual
      FOREIGN KEY ( term_id ) REFERENCES term ( term_id ) ;
ALTER TABLE dbxref_qualifier_value ADD CONSTRAINT FKdbxref_dbxrefqual
      FOREIGN KEY ( dbxref_id ) REFERENCES dbxref ( dbxref_id )
      ON DELETE CASCADE ;

-- bioentry_reference 
ALTER TABLE bioentry_reference ADD CONSTRAINT FKbioentry_entryref
      FOREIGN KEY ( bioentry_id ) REFERENCES bioentry ( bioentry_id )
      ON DELETE CASCADE ;
ALTER TABLE bioentry_reference ADD CONSTRAINT FKreference_entryref
      FOREIGN KEY ( reference_id ) REFERENCES reference ( reference_id )
      ON DELETE CASCADE ;

-- reference 
ALTER TABLE reference ADD CONSTRAINT FKdbxref_reference
      FOREIGN KEY ( dbxref_id ) REFERENCES dbxref ( dbxref_id ) ;

-- bioentry_qualifier_value 
ALTER TABLE bioentry_qualifier_value ADD CONSTRAINT FKbioentry_entqual
      FOREIGN KEY ( bioentry_id ) REFERENCES bioentry ( bioentry_id )
      ON DELETE CASCADE ;
ALTER TABLE bioentry_qualifier_value ADD CONSTRAINT FKterm_entqual
      FOREIGN KEY ( term_id ) REFERENCES term ( term_id ) ;

-- seqfeature 
ALTER TABLE seqfeature ADD CONSTRAINT FKterm_seqfeature
      FOREIGN KEY ( type_term_id ) REFERENCES term ( term_id ) ; 
ALTER TABLE seqfeature ADD CONSTRAINT FKsourceterm_seqfeature
      FOREIGN KEY ( source_term_id ) REFERENCES term ( term_id ) ; 
ALTER TABLE seqfeature ADD CONSTRAINT FKbioentry_seqfeature
      FOREIGN KEY ( bioentry_id ) REFERENCES bioentry ( bioentry_id )
      ON DELETE CASCADE ;

-- seqfeature_relationship 
ALTER TABLE seqfeature_relationship ADD CONSTRAINT FKterm_seqfeatrel
      FOREIGN KEY ( term_id ) REFERENCES term ( term_id ) ;
ALTER TABLE seqfeature_relationship ADD CONSTRAINT FKparentfeat_seqfeatrel
      FOREIGN KEY ( object_seqfeature_id ) REFERENCES seqfeature ( seqfeature_id )
      ON DELETE CASCADE ;
ALTER TABLE seqfeature_relationship ADD CONSTRAINT FKchildfeat_seqfeatrel
      FOREIGN KEY ( subject_seqfeature_id ) REFERENCES seqfeature ( seqfeature_id )
      ON DELETE CASCADE ;

-- seqfeature_path 
ALTER TABLE seqfeature_path ADD CONSTRAINT FKterm_seqfeatpath
      FOREIGN KEY ( term_id ) REFERENCES term ( term_id ) ;
ALTER TABLE seqfeature_path ADD CONSTRAINT FKparentfeat_seqfeatpath
      FOREIGN KEY ( object_seqfeature_id ) REFERENCES seqfeature ( seqfeature_id )
      ON DELETE CASCADE ;
ALTER TABLE seqfeature_path ADD CONSTRAINT FKchildfeat_seqfeatpath
      FOREIGN KEY ( subject_seqfeature_id ) REFERENCES seqfeature ( seqfeature_id )
      ON DELETE CASCADE ;

-- seqfeature_qualifier_value 
ALTER TABLE seqfeature_qualifier_value ADD CONSTRAINT FKterm_featqual
      FOREIGN KEY ( term_id ) REFERENCES term ( term_id ) ;
ALTER TABLE seqfeature_qualifier_value ADD CONSTRAINT FKseqfeature_featqual
      FOREIGN KEY ( seqfeature_id ) REFERENCES seqfeature ( seqfeature_id )
      ON DELETE CASCADE ;

-- seqfeature_dbxref 
ALTER TABLE seqfeature_dbxref ADD CONSTRAINT FKseqfeature_feadblink
      FOREIGN KEY ( seqfeature_id ) REFERENCES seqfeature ( seqfeature_id )
      ON DELETE CASCADE ;
ALTER TABLE seqfeature_dbxref ADD CONSTRAINT FKdbxref_feadblink
      FOREIGN KEY ( dbxref_id ) REFERENCES dbxref ( dbxref_id )
      ON DELETE CASCADE ;

-- location 
ALTER TABLE location ADD CONSTRAINT FKseqfeature_location
      FOREIGN KEY ( seqfeature_id ) REFERENCES seqfeature ( seqfeature_id )
      ON DELETE CASCADE ;
ALTER TABLE location ADD CONSTRAINT FKdbxref_location
      FOREIGN KEY ( dbxref_id ) REFERENCES dbxref ( dbxref_id ) ;
ALTER TABLE location ADD CONSTRAINT FKterm_featloc
      FOREIGN KEY ( term_id ) REFERENCES term ( term_id ) ;

-- location_qualifier_value 
ALTER TABLE location_qualifier_value ADD CONSTRAINT FKfeatloc_locqual
      FOREIGN KEY ( location_id ) REFERENCES location ( location_id )
      ON DELETE CASCADE ;
ALTER TABLE location_qualifier_value ADD CONSTRAINT FKterm_locqual
      FOREIGN KEY ( term_id ) REFERENCES term ( term_id ) ;

