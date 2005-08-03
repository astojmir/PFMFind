import xml.parsers.expat
from cStringIO import StringIO
      
def default_handler(name, identifier, accessions):
    pass

class UnirefParser(object):
    def __init__(self, handler=default_handler):
        self.handler = handler
        self.definition = None
        self.name = None
        self.accessions = []
        self.total_clusters = 0
        self.total_accessions = 0
        self._def_flag = False
      
    def _start_element(self, name, attrs):
        if name == 'entry':
            self.name = str(attrs['id'])
        elif name == 'name':
            self._def_flag = True
            self._fs = StringIO()
        elif name == 'property' and \
                 (attrs['type'] == "UniProt accession" or \
                  attrs['type'] == "UniProtKB accession"): 
            self.accessions.append(str(attrs['value']))

    def _end_element(self, name):
        if name == 'entry':
            if len(self.accessions):
                self.total_clusters += 1
                self.total_accessions += len(self.accessions)
                self.handler(self.name, self.definition, self.accessions)
            self.definition = None
            self.name = None
            self.accessions = []
        elif name == 'name':
            self._def_flag = False
            self.definition = self._fs.getvalue()

    def _char_data(self, data):
        if self._def_flag:
            self._fs.write(str(data))
         
    def parse(self, fp):
        p = xml.parsers.expat.ParserCreate()

        p.StartElementHandler = self._start_element
        p.EndElementHandler = self._end_element
        p.CharacterDataHandler = self._char_data
        
        p.ParseFile(fp) 


class UnirefLoader(object):
    def __init__(self, adaptor, dbid, uniref_name):
        self.adaptor = adaptor
        self.dbid = dbid

        self.ontology_id = self._get_ontology_id(uniref_name) 
        self.source_term_id = self._get_source_term_id(uniref_name)
        self.term = self._get_terms()
        self.bioentry = self._get_bioentries()

        self.total_clusters = 0
        self.total_accessions = 0

    def _get_ontology_id(self, name):
        oids = self.adaptor.execute_and_fetch_col0(
            "SELECT ontology_id FROM ontology WHERE name = %s",
            (name,))
        if oids:
            return oids[0]
        self.adaptor.execute(
            "INSERT INTO ontology(name, definition) VALUES (%s, %s)",
            (name, None))
        return self.adaptor.last_id("ontology")

    def _get_source_term_id(self, uniref_name):
        source_cat_id = self._get_ontology_id('Source Tags')

        oids = self.adaptor.execute_and_fetch_col0(
            "SELECT term_id FROM term WHERE name = %s AND ontology_id=%s", 
            (uniref_name, source_cat_id))
        if oids:
            return oids[0]
        
        sql = r"INSERT INTO term (name, ontology_id)" \
              r" VALUES (%s, %s)"
        self.adaptor.execute(sql, (uniref_name, source_cat_id)) 
        return self.adaptor.last_id("term")

    def _get_terms(self):
        terms = self.adaptor.execute_and_fetchall(
            "SELECT name, term_id FROM term WHERE ontology_id=%s",
            (self.ontology_id, ))
        return dict(terms)
        
    def _get_bioentries(self):
        sql = """SELECT a.accession, a.bioentry_id, b.maxrank FROM
        bioentry AS a LEFT OUTER JOIN (SELECT bioentry_id, max(rank)
        AS maxrank FROM seqfeature GROUP by bioentry_id) AS b ON
        (a.bioentry_id=b.bioentry_id) WHERE a.biodatabase_id=%s"""

        cols = self.adaptor.execute_and_fetchall(sql, (self.dbid, ))
        bioentries = {}
        for a, b, r in cols:
            if r == None:
                r = 1
            else:
                r += 1
            bioentries[a] = [b,r]
            
        return bioentries
 
    def _get_cluster_id(self, name, definition=None): 

        if name not in self.term:
            sql = r"INSERT INTO term (name," \
                  r" definition, ontology_id)" \
                  r" VALUES (%s, %s, %s)"
            self.adaptor.execute(sql, (name, definition,
                                       self.ontology_id)) 
            self.term[name] = self.adaptor.last_id("term")
            
        return self.term[name]

    def _load_cluster_member(self, term_id, accession):

        if accession not in self.bioentry: return
        
        bioentry_id, rank = self.bioentry[accession]
        
        sql = r"INSERT INTO seqfeature (bioentry_id, type_term_id, " \
              r"source_term_id, rank) VALUES (%s, %s, %s, %s)"
        self.adaptor.execute(sql, (bioentry_id, term_id,
                                   self.source_term_id, rank))

        self.bioentry[accession][1] += 1
        self.total_accessions += 1
    
    def _check_cluster(self, accessions):
        for a in accessions:
            if a in self.bioentry:
                return True
        return False
    
    def load_cluster(self, name, definition, accessions):

        if not self._check_cluster(accessions): return

        term_id = self._get_cluster_id(name, definition)
        for a in accessions:
            self._load_cluster_member(term_id, a)
        self.total_clusters += 1

if __name__ == "__main__":
    uniref_file = '/home/aleksand/data/bio/databases/UniProt/uniref-Feb2005/uniref50/uniref50.xml'
    commit = True


    from BioSQL import BioSeqDatabase

    server = BioSeqDatabase.open_database(driver = "psycopg", user = "aleksand", db = "test3")
## server.adaptor.execute("SET search_path TO %s;" % schema)
    db = server["m1"]



    UL = UnirefLoader(db.adaptor, db.dbid, 'Uniref50')
    UP = UnirefParser(UL.load_cluster)
    fp = file(uniref_file)
    UP.parse(fp)
    fp.close()
    if commit:
        server.adaptor.commit()
    server.adaptor.close()
    print "Parsed:", UP.total_clusters, UP.total_accessions
    print "Loaded:", UL.total_clusters, UL.total_accessions




