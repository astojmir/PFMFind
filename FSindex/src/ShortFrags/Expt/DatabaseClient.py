#
# Copyright (C) 2004-2006 Victoria University of Wellington
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


import cPickle, threading
from itertools import ifilterfalse
from ShortFrags.Expt.hit_list import HitList
from BioSQL import BioSeqDatabase


_hit_view = "(SELECT b.bioentry_id, ht.fragment, ht.start" \
            " FROM bioentry AS b JOIN" \
            " (SELECT h.fragment, h.accession, h.start FROM hits AS h" \
            " WHERE h.experiment_id=%d  AND h.length=%d AND" \
            " h.iteration=%d AND fragment BETWEEN %d AND %d)" \
            " AS ht USING (accession))"

_kw_view = "SELECT bt.bioentry_id,bt.fragment,f.term_id FROM" \
           " %s AS bt JOIN (SELECT s.bioentry_id," \
           " t.term_id FROM seqfeature AS s JOIN" \
           " (SELECT term_id FROM term WHERE" \
           " ontology_id=%d) AS t" \
           " ON s.type_term_id=t.term_id) AS f" \
           " USING (bioentry_id)"

_dom_view = "SELECT f.bioentry_id, f.fragment, f.term_id" \
            " FROM (SELECT bt.bioentry_id, bt.fragment, bt.start," \
            " e.term_id, e.seqfeature_id FROM %s AS bt JOIN" \
            " (SELECT s.bioentry_id, s.seqfeature_id, t.term_id" \
            " FROM seqfeature AS s JOIN (SELECT term_id FROM" \
            " term WHERE ontology_id=%d) AS t ON" \
            " s.type_term_id=t.term_id) AS e USING (bioentry_id))" \
            " AS f JOIN location AS l USING (seqfeature_id)" \
            " WHERE ((f.start+1) BETWEEN l.start_pos AND" \
            " l.end_pos) OR ((f.start+%d) BETWEEN l.start_pos"\
            " AND l.end_pos)" 


# This is becoming most difficult so I'll split it into components

_qualifiers = "SELECT sqv1.q1, sqv2.q2, sqv1.seqfeature_id FROM" \
              " (SELECT sqv.term_id AS q1, sqv.seqfeature_id FROM" \
              " seqfeature_qualifier_value AS sqv JOIN" \
              " (SELECT term_id FROM term WHERE ontology_id=%d)" \
              " AS tt USING (term_id)) AS sqv1 LEFT OUTER JOIN" \
              " (SELECT sqv.term_id AS q2, sqv.seqfeature_id FROM" \
              " seqfeature_qualifier_value AS sqv JOIN" \
              " (SELECT term_id FROM term WHERE ontology_id=%d)" \
              " AS tt USING (term_id)) AS sqv2 USING (seqfeature_id)"

_features = "SELECT s.bioentry_id, s.seqfeature_id, t.term_id" \
            " FROM seqfeature AS s JOIN (SELECT term_id FROM" \
            " term WHERE ontology_id=%d) AS t ON" \
            " s.type_term_id=t.term_id"

_full_qualifiers = "SELECT sft.seqfeature_id, sft.bioentry_id," \
                   " sft.term_id, sqv.q1, sqv.q2 FROM (%s) AS sft" \
                   " LEFT OUTER JOIN (%s) AS sqv USING" \
                   " (seqfeature_id)" % (_features, _qualifiers)


_feature_view = "SELECT f.bioentry_id, f.fragment, f.term_id," \
                " f.q1, f.q2, l.start_pos, l.end_pos" \
                " FROM (SELECT bt.bioentry_id, bt.fragment, bt.start," \
                " e.term_id, e.q1, e.q2, e.seqfeature_id FROM %%s" \
                " AS bt JOIN (%s) AS e USING (bioentry_id))" \
                " AS f JOIN location AS l USING (seqfeature_id)" \
                " WHERE ((f.start+1) BETWEEN l.start_pos AND" \
                " l.end_pos) OR ((f.start+%%d) BETWEEN l.start_pos"\
                " AND l.end_pos)" % (_full_qualifiers)


def _format_with_length(hit_table, oid, length):
    return (hit_table, oid, length)

def _format_without_length(hit_table, oid, length):
    return (hit_table, oid)


class DatabaseClient(object):
    _expt_attrs = ['expt_name',
                   'expt_description',
                   'query_sequence',
                   'query_description', 
                   'min_len', 
                   'max_len',
                   ]

    def __init__(self):

        self.driver = None
        self.dbargs = {}
        self.db_schema = None
        self.PFMF_schema = None

        self.conn = None
        self.module = None
        self.reset_current_experiment()

        self._views = {}
         
       # Only one thread at the time should be able to write to the
        # database. 
        self.db_write_lock = threading.Lock()
       
    def __del__(self):
        self.close()

    def open(self, driver="pgdb", **kwargs):
        """
        Opens a connection to a relational database.
        """

        self.driver = driver
        self.dbargs = kwargs
        self.server = BioSeqDatabase.open_database(driver=driver, **kwargs)
        self.conn = self.server.adaptor.conn
        self.crs = self.server.adaptor.cursor
               
    def close(self):
        """
        Cleaning up.
        """
        
        if self.conn:
            self.conn.commit()
            self.conn.close()
            self.conn = None
            self.module = None
        self.reset_current_experiment()

        self.driver = None
        self.dbargs = {}
        self.db_schema = None
        self.PFMF_schema = None

        self._views = {}

    # ************************************************************
    # ******* Schema *********************************************
    # ************************************************************

    def set_schema(self, db_schema='public', PFMF_schema=None):
        """
        Sets the database search path to consist of PFMF_schema,
        db_schema. The db_schema defaults to 'public' while
        PFMF_schema defaults to the database user (created if non
        existent).
        """

        if not PFMF_schema:
            PFMF_schema = self.dbargs['user']

        # PostgreSQL is case-insensitive so we need to covert schema
        # names to lower case
        db_schema = db_schema.lower()
        PFMF_schema = PFMF_schema.lower()

        sql = "SELECT schema_name FROM information_schema.schemata"
        self.crs.execute(sql)
        schemata = [s[0] for s in self.crs.fetchall()]

        # Validate db schema
        if db_schema not in schemata:
            self.close()
            raise RuntimeError("Missing db_schema '%s'." % db_schema)

        sql = "SELECT table_name FROM information_schema.tables"\
              " WHERE table_schema='%s'" % db_schema
        self.crs.execute(sql)
        BioSQL_tables = [s[0] for s in self.crs.fetchall()]

        if 'bioentry' not in BioSQL_tables:
            self.close()
            raise RuntimeError("Invalid db_schema '%s'." % db_schema)

        # Create PFMF_schema if necessary
        if PFMF_schema not in schemata:
            self.crs.execute("CREATE SCHEMA %s" % PFMF_schema)
            self.conn.commit()

        search_path = "%s, %s" % (PFMF_schema, db_schema)
        self.crs.execute("SET search_path TO %s" % search_path)

        if PFMF_schema not in schemata:
            self.load_PFMF_schema()

        self.db_schema = db_schema
        self.PFMF_schema = PFMF_schema


        # Get available ontologies
        self.crs.execute('SELECT name, ontology_id FROM ontology')
        self._ontologies = dict(self.crs.fetchall())

        # Get those ontologies that point to features
        def filter_tags(tag):
            return tag in ['Source Tags',
                           'Uniprot Feature Qualifiers',
                           'Annotation Tags',
                           'Non Experimental Qualifiers',
                           ]
        self.feature_types = [tag for tag in \
            ifilterfalse(filter_tags, self._ontologies.keys())]   

        # Assign functions to retrieve features

        for k in self.feature_types:
            if k == 'Uniprot Keywords':
                self._views[k] = (_kw_view, "name",
                                  _format_without_length)
            elif k == 'Uniprot Feature Keys':
                self._views[k] = (_dom_view, "name",
                                  _format_with_length )
            elif k == 'InterPro domain':
                self._views[k] = (_dom_view,
                                  "definition || ' (' || name || ')'",
                                  _format_with_length)
            else: #uniref cluster
                self._views[k] = (_kw_view,
                                  "definition || ' (' || name || ')'",
                                  _format_without_length)  


    def load_PFMF_schema(self):
        """
        Creates the tables to store PFMF hits, searches and experiments.
        """

        self.db_write_lock.acquire()
        sql = """
        CREATE TABLE hits (
            experiment_id SMALLINT NOT NULL ,
            length SMALLINT NOT NULL ,
            iteration SMALLINT NOT NULL ,
            fragment INTEGER NOT NULL ,
            accession VARCHAR ( 40 ) NOT NULL , 
            start INTEGER NOT NULL ,
            distance SMALLINT ,
            similarity SMALLINT ,
            pvalue REAL ,
            Evalue REAL ,
            PRIMARY KEY ( experiment_id, length, iteration, fragment, accession, start ) ); 

        CREATE TABLE searches (
            experiment_id SMALLINT NOT NULL ,
            length SMALLINT NOT NULL ,
            iteration SMALLINT NOT NULL ,
            fragment INTEGER NOT NULL ,
            query_frag VARCHAR ( 32 ) NOT NULL ,
            matrix_name VARCHAR ( 40 ) ,
            score_matrix TEXT ,
            conv_type SMALLINT ,
            sim_range SMALLINT ,
            dist_range SMALLINT ,
            kNN SMALLINT ,
            num_hits INTEGER NOT NULL ,
            PRIMARY KEY ( experiment_id, length, iteration, fragment ) ); 

        CREATE SEQUENCE experiments_pk_seq;

        CREATE TABLE experiments (
            experiment_id SMALLINT DEFAULT nextval ( 'experiments_pk_seq' ) NOT NULL ,
            name VARCHAR ( 40 ) UNIQUE NOT NULL ,
            description TEXT ,
            query_sequence TEXT NOT NULL,
            query_description TEXT ,
            min_len INTEGER NOT NULL ,
            max_len INTEGER NOT NULL ,
            PRIMARY KEY ( experiment_id ) );
        """
        self.crs.execute(sql)
        self.conn.commit()
        self.db_write_lock.release()

    # ************************************************************
    # ******* Experiments ****************************************
    # ************************************************************

    def get_experiment_names(self):
        """
        Retrieves a list of pairs (experiment_id, name)
        denoting the experiments stored in the database.
        """

        sql = "SELECT experiment_id, name FROM experiments ORDER BY" \
              " experiment_id DESC" 
        cursor = self.conn.cursor()
        cursor.execute(sql)
        results = cursor.fetchall()
        cursor.close()
        return results

    def get_experiment_data(self, experiment_id):
        """
        Retrieves full experiment data from the database.
        """
        
        sql = """SELECT name, description, query_sequence,
                   query_description, min_len, max_len
                   FROM experiments WHERE experiment_id=%s"""        
        self.crs.execute(sql, (experiment_id,))
        return self.crs.fetchone()

    def get_experiment_id(self, name):
        """
        Retrieves experiment_id from a name.
        """
        
        sql = "SELECT experiment_id FROM experiments WHERE name=%s"
        cursor = self.conn.cursor()
        cursor.execute(sql, (name,))
        res = cursor.fetchall()
        cursor.close()
        if len(res) == 0:
            return None
        else:
            return res[0][0]

    def create_experiment(self, name, description, query_sequence, 
                          query_description, min_len, max_len):
        """
        Creates a new experiment in the database.
        """

        self.db_write_lock.acquire()
        sql = """
        INSERT INTO experiments (name, description, query_sequence,
        query_description, min_len, max_len)
        VALUES (%s, %s, %s, %s, %s, %s)"""
        self.crs = self.conn.cursor()
        self.crs.execute(sql, (name, description, query_sequence,
                               query_description, min_len, max_len))
        self.conn.commit()
        self.db_write_lock.release()
        return self.get_experiment_id(name)

    def update_experiment(self, experiment_id, description=None):
        """
        Updates the experiment description. Other columns may not be
        updated.   
        """
        
        self.db_write_lock.acquire()
        if description:
            sql = "UPDATE experiments SET description=%s WHERE" \
                  " experiment_id=%s"
            self.crs.execute(sql, (description, experiment_id))
            self.conn.commit()
        self.db_write_lock.release()
        
    def set_current_experiment(self, experiment_id):
        """
        Sets the experiment all searches refer to.
        """
        self.cur_expt = experiment_id
        vals = self.get_experiment_data(experiment_id)

        for k,v in zip(self._expt_attrs, vals):
            setattr(self, k, v)
        self.get_max_iters()

    def reset_current_experiment(self):
        """
        Sets all current experiment data to None.
        """
        self.cur_expt = None
        for k in self._expt_attrs:
            setattr(self, k, None)
        self.max_iters = {}

    def delete_experiment(self, experiment_id):
        """
        Deletes all searches and hits associated with
        a given experiment.
        
        Deleting the current experiment is not allowed. 
        """

        if experiment_id == self.cur_expt:
            raise ValueError("Deleting the current experiment is not allowed.")

        self.db_write_lock.acquire()
            
        sql = """ DELETE FROM hits WHERE experiment_id = %s """
        self.crs.execute(sql, (experiment_id,))
        sql = """ DELETE FROM searches WHERE experiment_id = %s """
        self.crs.execute(sql, (experiment_id,))
        sql = """ DELETE FROM experiments WHERE experiment_id = %s """
        self.crs.execute(sql, (experiment_id,))

        self.conn.commit()
        self.db_write_lock.release()

    # ************************************************************
    # ******* Searches *******************************************
    # ************************************************************

    def get_max_iters(self):
        sql = "SELECT fragment, length, max(iteration) FROM searches" \
              " WHERE experiment_id=%s GROUP BY fragment, length" \
              " ORDER BY fragment, length"
        cursor = self.conn.cursor()
        cursor.execute(sql, (self.cur_expt,))
        vals = cursor.fetchall()
        cursor.close()
        self.max_iters = {}
        for f, l, i in vals:
            if f in self.max_iters:
                self.max_iters[f][l] = i
            else:
                self.max_iters[f] = {l : i}
        
    def get_num_hits(self, length, iteration, fragment):
        """
        Returns the number of hits for given parameters or
        None if no search was performed. 
        """
        
        sql = """
        SELECT num_hits FROM searches WHERE
          experiment_id=%s AND length=%s AND iteration=%s AND fragment=%s"""
        cursor = self.conn.cursor()
        cursor.execute(sql, (self.cur_expt, length,
                                  iteration, fragment))
        res = cursor.fetchall()
        cursor.close()
        if len(res) == 0:
            return 0
        else:
            return res[0][0]
        
    def insert_search(self, length, iteration, fragment, HL, lock=True):
        """
        Inserts a hit list into relational database. Scoring matrix can
        be optionally inserted as well.
        """

        if lock:
            self.db_write_lock.acquire()

        matrix_name = getattr(HL, 'matrix_name')
        matrix = getattr(HL, 'matrix')

        # Pack matrix
        if matrix:
            matrix = cPickle.dumps(matrix, 0)

        sql = """ INSERT INTO searches (
          experiment_id, length, iteration, fragment,
          query_frag, matrix_name, score_matrix, conv_type,
          sim_range, dist_range, kNN, num_hits)
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)"""
        data = (self.cur_expt, length, iteration, fragment,
                HL.query_seq, matrix_name, matrix, HL.conv_type,
                HL.sim_range, HL.dist_range, HL.kNN, len(HL))
        self.crs.execute(sql, data)

        sql = """ INSERT INTO hits (
          experiment_id, length, iteration, fragment,
          accession, start, distance, similarity, pvalue, Evalue)
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)"""
        for ht in HL:
            data = (self.cur_expt, length, iteration, fragment,
                    ht.accession, ht.seq_from, ht.dist,
                    ht.sim, ht.pvalue, ht.Evalue)
            self.crs.execute(sql, data)
        self.conn.commit()
        if lock:
            self.db_write_lock.release()

    def delete_search(self, length, iteration, fragment):
        """
        Deletes the search associated with (length, iteration, fragment)
        """

        self.db_write_lock.acquire()

        sql = """ DELETE FROM hits WHERE experiment_id = %s AND
                    length=%s AND iteration = %s AND fragment = %s """
        self.crs.execute(sql, (self.cur_expt, length,
                                  iteration, fragment))

        sql = """ DELETE FROM searches WHERE experiment_id = %s AND
                    length=%s AND iteration = %s AND fragment = %s """
        self.crs.execute(sql, (self.cur_expt, length,
                                  iteration, fragment))

        self.conn.commit()
        self.db_write_lock.release()

    def replace_search(self, length, iteration, fragment, HL):
        """
        First deletes all hits having the same length and fragment and
        iteration greater than or equal to given iteration, then inserts
        a new hit list and matrix (optional)
        """

        self.db_write_lock.acquire()

        sql = """ DELETE FROM hits WHERE experiment_id = %s AND
                    length=%s AND iteration >= %s AND fragment = %s """
        self.crs.execute(sql, (self.cur_expt, length,
                                  iteration, fragment))

        sql = """ DELETE FROM searches WHERE experiment_id = %s AND
                    length=%s AND iteration >= %s AND fragment = %s """
        self.crs.execute(sql, (self.cur_expt, length,
                                  iteration, fragment))

        self.insert_search(length, iteration, fragment, HL, lock=False)
        self.db_write_lock.release()
    
    def select_lif_search(self, length, iteration, fragment,
                          get_matrix=False, features=[]):
        """
        Retrieves the list of hits for a given search from the database
        """

        cursor = self.conn.cursor()
        # Get the search parameters first
        HL_dict = {'query_def': self.query_description}
        attribs1 = ['query_seq', 'conv_type', 'sim_range',
                    'dist_range', 'kNN', 'matrix_name']
        attribs2 = ['_bioentry_id', '_taxon_id', 'accession',
                    'seq_from', 'dist', 'sim', 'pvalue', 'Evalue',
                    'defline', 'sequence'] 


        if get_matrix:
            sql = """
            SELECT query_frag, conv_type, sim_range, dist_range, kNN,
              matrix_name, score_matrix FROM searches WHERE experiment_id=%s
              AND length=%s AND iteration=%s AND fragment=%s"""
        else:
            sql = """
            SELECT query_frag, conv_type, sim_range, dist_range, kNN,
              matrix_name FROM searches WHERE experiment_id=%s AND
              length=%s AND iteration=%s AND fragment=%s"""

        cursor.execute(sql, (self.cur_expt, length,
                                  iteration, fragment))
        res = cursor.fetchone()
        if res == None: return None

        # Unpack the matrix
        if get_matrix:
            if res[-1] != None:
                HL_dict['matrix'] = cPickle.loads(res[-1])
            else:
                HL_dict['matrix'] = None
            res = res[:-1]
            
        for key, val in zip(attribs1,res):
            HL_dict[key] = val

            
        sql = """SELECT a.*, substring(s.seq from a.start+1 for %s ),
        %s AS fragment FROM
        (SELECT b.bioentry_id, b.taxon_id, h1.*,
          b.accession || '|' || b.description FROM
          (SELECT h.accession, h.start, h.distance, h.similarity,
          h.pvalue, h.Evalue FROM hits AS h WHERE
          h.experiment_id=%s AND h.length=%s 
          AND h.iteration=%s AND h.fragment=%s)
          AS h1 JOIN bioentry AS b ON h1.accession=b.accession)
        AS a JOIN biosequence AS s ON a.bioentry_id=s.bioentry_id
        ORDER BY a.distance, a.accession, a.start;"""
        cursor.execute(sql, (length, fragment, self.cur_expt,
                             length, iteration, fragment))

        HL_dict['hits'] = []
        while 1: 
            res =cursor.fetchone()
            if res == None: break
            hd = {}

            for key, val in zip(attribs2, res[:-1]):
                hd[key] = val
            hd['seq_to'] = hd['seq_from'] + length
            if hd['defline'] == None:
                hd['defline'] = hd['accession']

            HL_dict['hits'].append(hd)

        # Now get the species
        taxons = """SELECT b.taxon_id FROM
         (SELECT h.accession FROM hits AS h WHERE
         h.experiment_id=%s AND h.length=%s 
         AND h.iteration=%s AND h.fragment=%s) AS h1
         JOIN bioentry AS b ON h1.accession=b.accession
         GROUP BY taxon_id """
        

        sql = "SELECT i.taxon_id, n.name FROM"\
              " (%s) AS i JOIN"\
              " (SELECT taxon_id, name FROM taxon_name"\
              " WHERE name_class='scientific name') AS n"\
              " USING (taxon_id)" % taxons

        cursor.execute(sql, (self.cur_expt, length, iteration, fragment))
        res = cursor.fetchall()
        taxon_dict = dict(res)
        for hd in HL_dict['hits']:
            if hd['_taxon_id'] and hd['_taxon_id'] in taxon_dict:
                hd['taxon'] = taxon_dict[hd['_taxon_id']]
            hd['features'] = {}


        # It seems the full SQL query is faster than temporary table
        hit_table = _hit_view % (self.cur_expt, length, iteration,
                                 fragment, fragment)
        for ft in features:
            bt_dict = self.get_bioentry_features(hit_table,
                                                 ft,
                                                 length)
            for hd in HL_dict['hits']:
                if hd['_bioentry_id'] in bt_dict:
                    hd['features'][ft] = bt_dict[hd['_bioentry_id']]

        cursor.close()
        return HitList(HL_dict)

    # ************************************************************
    # ******* Features/Views *************************************
    # ************************************************************

    def select_li_features(self, length, iteration, feature_type,
                           frag_range): 
        """
        Returns the keywords sorted for KeywordView GUI.
        """

        if feature_type not in self.feature_types: return []

        cursor = self.conn.cursor()

        oid = self._ontologies[feature_type]
        hit_table = _hit_view % (self.cur_expt, length, iteration,
                                 frag_range[0], frag_range[1]-1)
        feature_table = self._views[feature_type][0] % \
                        self._views[feature_type][2](hit_table,
                                                     oid, length)

        sql = "SELECT term_id, fragment, bioentry_id FROM (%s) AS fttbl" \
              " ORDER BY term_id, fragment" % feature_table
        cursor.execute(sql)

        view_list = []
        old_tid = None
        while 1:
            res = cursor.fetchone()
            if res == None: break
            tid = int(res[0])
            frag = int(res[1])
            acc = res[2]
            
            if tid == old_tid:
                if frag == old_frag:
                    view_list[-1][1][-1][1].append(acc)
                else:
                    view_list[-1][1].append((frag, [acc]))
            else:
                view_list.append([tid, [(frag, [acc])]])

            old_tid = tid
            old_frag = frag
            
        for i in xrange(len(view_list)):
            tid = view_list[i][0]
            cursor.execute("SELECT %s FROM term WHERE term_id=%d" \
                           % (self._views[feature_type][1], tid))
            view_list[i][0] = cursor.fetchone()[0]

        cursor.close()
        return view_list

    def get_full_uniprot_features(self, bioentry_table, length):
        """
        Retrieve full description (not only the keys) of Uniprot
        features of bioentries in bioentry_table from the database.
        """
        
        cursor = self.conn.cursor()

        oid_ft_key = self._ontologies['Uniprot Feature Keys']
        oid_qf = self._ontologies['Uniprot Feature Qualifiers']
        oid_nexp = self._ontologies['Non Experimental Qualifiers']

        feature_table = _feature_view % \
            (bioentry_table,
             self._ontologies['Uniprot Feature Keys'],
             self._ontologies['Uniprot Feature Qualifiers'],
             self._ontologies['Non Experimental Qualifiers'],
             length)

        sql = "SELECT bioentry_id,term_id,q1,q2,start_pos,end_pos" \
              " FROM (%s) AS fttbl ORDER BY bioentry_id" % \
              feature_table 
        cursor.execute(sql)

        bioentry_dict ={}
        term_dict = {}
        qf_dict = {}
        nexp_dict = {}
        while 1:
            res = cursor.fetchone()
            if res == None: break
            bid = res[0]
            ft_data = res[1:]
            if bid in bioentry_dict:
                bioentry_dict[bid].append(ft_data)
            else:
                bioentry_dict[bid] = [ft_data]
            term_dict[ft_data[0]] = None
            if ft_data[1]:
                qf_dict[ft_data[1]] = None
            if ft_data[2]:
                nexp_dict[ft_data[2]] = None
            
        for tid in term_dict.iterkeys():
            cursor.execute("SELECT name FROM term WHERE term_id=%d" \
                           % tid)
            term_dict[tid] = cursor.fetchone()[0]
        for tid in qf_dict.iterkeys():
            cursor.execute("SELECT name FROM term WHERE term_id=%d" \
                           % tid)
            qf_dict[tid] = cursor.fetchone()[0]
        for tid in nexp_dict.iterkeys():
            cursor.execute("SELECT name FROM term WHERE term_id=%d" \
                           % tid)
            nexp_dict[tid] = cursor.fetchone()[0]


        for tids in bioentry_dict.itervalues():
            for i in xrange(len(tids)):
                desc = "%s (%d, %d)" % (term_dict[tids[i][0]],
                                        tids[i][3]-1,
                                        tids[i][4])
                
                if tids[i][1]:
                    desc += ": %s" % qf_dict[tids[i][1]]
                if tids[i][2]:
                    if tids[i][1]:
                        desc += " (%s)" % nexp_dict[tids[i][2]]
                    else:
                        desc += ": (%s)" % nexp_dict[tids[i][2]]
                tids[i] = desc
            tids.sort()

        cursor.close()
        return bioentry_dict

    def get_bioentry_features(self, bioentry_table,
                              feature_type, length): 
        """
        Retrieve descriptions of features of type feature_type of
        bioentries given in bioentry_table from the database.

        If feature type is 'Uniprot Feature Keys', the description is
        given by get_full_uniprot_features.
        """

        if feature_type not in self.feature_types:
            return []
        elif feature_type == 'Uniprot Feature Keys':
            return self.get_full_uniprot_features(bioentry_table,
                                                  length) 
            
        cursor = self.conn.cursor()

        oid = self._ontologies[feature_type]
        feature_table = self._views[feature_type][0] % \
                        self._views[feature_type][2](bioentry_table,
                                                     oid, length)

        sql = "SELECT bioentry_id, term_id FROM (%s) AS fttbl" \
              " ORDER BY bioentry_id" % feature_table
        cursor.execute(sql)

        bioentry_dict ={}
        term_dict = {}
        while 1:
            res = cursor.fetchone()
            if res == None: break
            bid, tid = res
            if bid in bioentry_dict:
                bioentry_dict[bid].append(tid)
            else:
                bioentry_dict[bid] = [tid]
            term_dict[tid] = None
            
        for tid in term_dict.iterkeys():
            cursor.execute("SELECT %s FROM term WHERE term_id=%d" \
                           % (self._views[feature_type][1], tid))
            term_dict[tid] = cursor.fetchone()[0]
        
        for tids in bioentry_dict.itervalues():
            for i in xrange(len(tids)):
                tids[i] = term_dict[tids[i]]
            tids.sort()

        cursor.close()
        return bioentry_dict

