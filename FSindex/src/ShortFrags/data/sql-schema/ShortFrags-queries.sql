

-- Getting SwissProt Features: qualifier should be 'By similarity', 'Potential' etc.
-- InterPro will not have seqfeature_qualifier_value_entries

SELECT c.bioentry_id, c.type_term_id, c.source_term_id, c.start_pos, c.end_pos, d.value, d.qualifier FROM (SELECT a.seqfeature_id, a.bioentry_id, a.type_term_id, a.source_term_id, b.start_pos, b.end_pos FROM seqfeature AS a JOIN location AS b ON a.seqfeature_id = b.seqfeature_id) AS c JOIN (SELECT e.seqfeature_id, e.term_id, e.value, f.term_id AS qualifier FROM seqfeature_qualifier_value AS e FULL OUTER JOIN seqfeature_qualifier_value AS f ON e.seqfeature_id = f.seqfeature_id WHERE e.term_id=22 AND f.term_id=21) AS d ON c.seqfeature_id = d.seqfeature_id



---- 1. SwissProt keywords (optionally filtered by category) with qualifiers
---- 2. Uniref clusters filtered by type (e.g. Uniref90, Uniref50)
---- 3. SwissProt features with descriptions, qualifiers, filtered by type
---- 4. Interpro domains (could be later filtered by type - GO)

---- Perhaps best to keep all in seqfeature.


---- Typing query (term_path): give me a list of all terms belonging to a given category,
---- then join with seqfeature on type_term_id

--- Then:

---- 1. get qualifiers (if any)
---- 2. OK
---- 3. get location, value, qualifiers
---- 4. get location

---- Then filter by location

---- Join with bioentries, sort




---- ************* NEW ******************
---- KEYWORDS - get ontology_id
SELECT bioentry.accession, b.name FROM bioentry JOIN (SELECT seqfeature.bioentry_id, a.name FROM seqfeature JOIN (SELECT term_id, name FROM term WHERE ontology_id=2) AS a ON seqfeature.type_term_id=a.term_id) AS b ON bioentry.bioentry_id=b.bioentry_id;

---- UNIPROT FEATURES
SELECT g.*, h.name AS qualifier FROM (SELECT c.*, d.start_pos, d.end_pos FROM (SELECT bioentry.accession, b.name, b.seqfeature_id FROM bioentry JOIN (SELECT s.bioentry_id, s.seqfeature_id, a.name FROM seqfeature AS s JOIN (SELECT term_id, name FROM term WHERE ontology_id=3) AS a ON s.type_term_id=a.term_id) AS b ON bioentry.bioentry_id=b.bioentry_id) AS c JOIN location AS d ON c.seqfeature_id=d.seqfeature_id) AS g LEFT OUTER JOIN (SELECT e.seqfeature_id, f.name FROM seqfeature_qualifier_value AS e JOIN term AS f ON e.term_id=f.term_id) AS h ON g.seqfeature_id=h.seqfeature_id;

---- UNIREF CLUSTER
SELECT bioentry.accession, b.name, b.definition FROM bioentry JOIN (SELECT seqfeature.bioentry_id, a.name, a.definition FROM seqfeature JOIN (SELECT term_id, name, definition FROM term WHERE ontology_id=7) AS a ON seqfeature.type_term_id=a.term_id) AS b ON bioentry.bioentry_id=b.bioentry_id ORDER BY name;

---- INTERPRO DOMAIN
SELECT c.*, d.start_pos, d.end_pos FROM (SELECT bioentry.accession, b.name, b.seqfeature_id FROM bioentry JOIN (SELECT s.bioentry_id, s.seqfeature_id, a.name FROM seqfeature AS s JOIN (SELECT term_id, name FROM term WHERE ontology_id=8) AS a ON s.type_term_id=a.term_id) AS b ON bioentry.bioentry_id=b.bioentry_id) AS c JOIN location AS d ON c.seqfeature_id=d.seqfeature_id;
