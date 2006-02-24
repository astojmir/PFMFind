#
# Copyright (C) 2005-2006 Victoria University of Wellington
#
# This file is part of the PFMFind module.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2,
# or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
#


from string import split

def default_handler(ac, start, end, name, ipid, desc):
    pass

class Protein2InterProParser(object):
    def __init__(self, handler=default_handler):
        self.handler = handler
        self.total_lines = 0

    def parse(self, fp):
        for line in fp:
            ls = split(line, None, 5)
            if len(ls) != 6: continue
            ac, ipid, start, end, name, desc = ls
            desc = desc.rstrip()
            try:
                start = int(start)
                end = int(end)
            except ValueError:
                continue
            self.handler(ac, start, end, name, ipid, desc)
            self.total_lines += 1

# Apparently a new format was introduced, tab-separated:
# protein_ac      IPR_ac  IPR_name        method_ac       pos_from        pos_to

class Protein2IprParser(object):
    def __init__(self, handler=default_handler):
        self.handler = handler
        self.total_lines = 0

    def parse(self, fp):
        for line in fp:
            ls = line.split('\t', 5)
            if len(ls) != 6: continue
            ac, ipid, desc, name, start, end = ls 
            desc = desc.rstrip()
            try:
                start = int(start)
                end = int(end)
            except ValueError:
                continue
            self.handler(ac, start, end, name, ipid, desc)
            self.total_lines += 1





class InterProLoader(object):
    def __init__(self, adaptor, dbid):
        self.adaptor = adaptor
        self.dbid = dbid

        self.ontology_id = self._get_ontology_id('InterPro domain') 
        self.source_term_id = self._get_source_term_id('InterPro')
        self.term = self._get_terms()
        self.bioentry = self._get_bioentries()

        self.total_domains = 0

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
 
    def _get_domain_id(self, name, definition): 

        if name not in self.term:
            sql = r"INSERT INTO term (name," \
                  r" definition, ontology_id)" \
                  r" VALUES (%s, %s, %s)"
            self.adaptor.execute(sql, (name, definition,
                                       self.ontology_id))    
            self.term[name] = self.adaptor.last_id("term")
            
        return self.term[name]

    def load_domain(self, ac, start, end, name, ipid, desc):

        if ac not in self.bioentry: return
        bioentry_id, rank = self.bioentry[ac]

        term_id = self._get_domain_id(name, desc)
        
        sql = r"INSERT INTO seqfeature (bioentry_id, type_term_id, " \
              r"source_term_id, rank) VALUES (%s, %s, %s, %s)"
        self.adaptor.execute(sql, (bioentry_id, term_id,
                                   self.source_term_id, rank))
        feature_id = self.adaptor.last_id("seqfeature")
        self.bioentry[ac][1] += 1

        sql = r"INSERT INTO location (seqfeature_id, " \
              r"start_pos, end_pos, strand, rank) " \
              r"VALUES (%s, %s, %s, %s, %s)"
        self.adaptor.execute(sql, (feature_id, start, end,0, 1))
        self.total_domains += 1
