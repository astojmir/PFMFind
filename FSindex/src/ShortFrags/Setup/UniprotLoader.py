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

import re
from Bio import SeqIO
import Bio.SeqIO.SwissIO
from BioSQL.Loader import DatabaseLoader
from BioSQL import BioSeqDatabase


# We must resort to monkey-patching Bio.SeqIO.SwissIO in order to
# avoid code duplication. It is not possible to subclass because
# Bio.SeqIO.SwissIO._make_seqfeature is a module function rather
# than method
_old_make_seqfeature = Bio.SeqIO.SwissIO._make_seqfeature
_nonexp = re.compile(r"(.*?) ?(?:\(?(Potential|Probable|By similarity)\)?)?\.?$")

def _new_make_seqfeature(name, from_res, to_res, description, ft_id):
    feature = _old_make_seqfeature(name, from_res, to_res, description, ft_id)
    # Must do this because the qualifiers for these
    # feature types are too long to fit as term names
    if type == 'CONFLICT' or type == 'VARIANT' or \
       type == 'VARSPLIC' or type == 'MUTAGEN':
        qualifiers = {'value': [description.rstrip('. ')]}
    else:
        qualifiers = {}
        m = _nonexp.match(description)
        if m is not None:
            g = m.groups()
            if len(g[0]):
                qualifiers['term'] = [g[0]]
            if g[1] is not None:
                qualifiers['non_exp'] = [g[1]]

    feature.qualifiers.update(qualifiers)
    return feature

Bio.SeqIO.SwissIO._make_seqfeature = _new_make_seqfeature



class UniprotLoader(DatabaseLoader):
    def __init__(self, *args, **kwargs):
        DatabaseLoader.__init__(self, *args, **kwargs)

        self.src_id = self._get_source_term_id()
        self.kw_tag_id = self._get_ontology_id('Uniprot Keywords')
        self.ft_tag_id = self._get_ontology_id('Uniprot Feature Keys')
        onts = 'Uniprot Feature Qualifiers'
        self.qf_tag_id = self._get_ontology_id(onts)
        onts = 'Non Experimental Qualifiers'
        self.nx_tag_id = self._get_ontology_id(onts)
        self.term = self._get_terms()

        ann_tag_id = self._get_ontology_id('Annotation Tags')
        self.val_tag = self._get_term_id('Change indicator value',
                                         ann_tag_id)

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

    def _get_source_term_id(self):
        source_cat_id = self._get_ontology_id('Source Tags')

        oids = self.adaptor.execute_and_fetch_col0(
            "SELECT term_id FROM term WHERE name = %s AND ontology_id=%s",
            ("Uniprot", source_cat_id))
        if oids:
            return oids[0]

        sql = r"INSERT INTO term (name, ontology_id)" \
              r" VALUES (%s, %s)"
        self.adaptor.execute(sql, ("Uniprot", source_cat_id))
        return self.adaptor.last_id("term")

    def _get_terms(self):
        terms = self.adaptor.execute_and_fetchall(
            "SELECT name, ontology_id, term_id FROM term")

        return dict([[(t[0],t[1]),t[2]] for t in terms])

    def _get_term_id(self, name, ontology_id,
                     definition=None, identifier=None):

        # Sometimes name is too long
        name = name[:255]
        if (name, ontology_id) not in self.term:
            sql = r"INSERT INTO term (name, definition," \
                  r" identifier, ontology_id)" \
                  r" VALUES (%s, %s, %s, %s)"
            self.adaptor.execute(sql, (name, definition,
                                       identifier, ontology_id))
            self.term[(name, ontology_id)] = self.adaptor.last_id("term")

        return self.term[(name, ontology_id)]

    def load_seqrecord(self, record):
        """Load a Biopython SeqRecord into the database.
        """
        self.current_rank = 1
        bioentry_id = self._load_bioentry_table(record)
        self._load_biosequence(record, bioentry_id)
        self._load_keywords(record, bioentry_id)
        for seq_feature_num in range(len(record.features)):
            seq_feature = record.features[seq_feature_num]
            if seq_feature.location.start.position is None or \
               seq_feature.location.end.position is None:
                continue
            self._load_seqfeature(seq_feature,  self.current_rank, bioentry_id)
            self.current_rank += 1

    def _get_taxon_id(self, record):
        """Get the id corresponding to a taxon.

        We assume the NCBI taxon data is preloaded and that
        ncbi_taxon_id = taxon_id.
        """

        ncbi_taxon_id = record.annotations.get("ncbi_taxid")
        if not ncbi_taxon_id:
            return None
        else:
            return int(ncbi_taxon_id[0])

    def _load_bioentry_table(self, record):
        """Fill the bioentry table with sequence information.
        """
        # get the pertinent info and insert it

        if record.id.find('.') >= 0: # try to get a version from the id
            accession, version = record.id.split('.')
            version = int(version)
        else: # otherwise just use a version of 0
            accession = record.id
            version = 0

        taxon_id = self._get_taxon_id(record)
        identifier = record.annotations.get('gi')
        description = getattr(record, 'description', None)
        division = record.annotations.get("data_file_division", "UNK")

        sql = """
        INSERT INTO bioentry (
         biodatabase_id, taxon_id, name, accession, identifier, division,
         description, version) VALUES (%s, %s, %s, %s, %s, %s, %s, %s)"""
        self.adaptor.execute(sql, (self.dbid,
                                   taxon_id,
                                   record.name,
                                   accession,
                                   identifier,
                                   division,
                                   description,
                                   version))
        # now retrieve the id for the bioentry
        bioentry_id = self.adaptor.last_id('bioentry')

        return bioentry_id

    def _load_keywords(self, record, bioentry_id):
        """Adds keywords as bioentry qualifiers.
        """
        keywords = record.annotations.get("keywords", [])

        sql = r"INSERT INTO seqfeature" \
              r" (bioentry_id, type_term_id, source_term_id, rank)" \
              r" VALUES (%s, %s, %s, %s)"

        for keyword in keywords:
            kw_id = self._get_term_id(keyword, self.kw_tag_id)
            self.adaptor.execute(sql, (bioentry_id, kw_id,
                                       self.src_id, self.current_rank))
            self.current_rank += 1

    def _load_seqfeature_basic(self, feature_type, feature_rank,
                               bioentry_id):
        """Load the first tables of a seqfeature and returns the id.

        This loads the "key" of the seqfeature (ie. CDS, gene) and
        the basic seqfeature table itself.
        """

        seqfeature_key_id = self._get_term_id(feature_type,
                                              self.ft_tag_id)

        sql = r"INSERT INTO seqfeature (bioentry_id, type_term_id, " \
              r"source_term_id, rank) VALUES (%s, %s, %s, %s)"
        self.adaptor.execute(sql, (bioentry_id, seqfeature_key_id,
                                   self.src_id, feature_rank + 1))
        seqfeature_id = self.adaptor.last_id('seqfeature')

        return seqfeature_id


    def _load_seqfeature_qualifiers(self, qualifiers, seqfeature_id):

        qval = ''
        for qk in qualifiers.keys():
            if qk == 'value':
                qk_id = self.val_tag
                qval = qualifiers[qk][0]
            elif qk == 'term':
                qk_id = self._get_term_id(qualifiers[qk][0],
                                          self.qf_tag_id)
            elif qk == 'non_exp':
                qk_id = self._get_term_id(qualifiers[qk][0],
                                          self.nx_tag_id)
            else:
                continue

            # Assume only one qualifier value per term_id
            sql = r"INSERT INTO seqfeature_qualifier_value VALUES" \
                  r" (%s, %s, %s, %s)"
            self.adaptor.execute(sql, (seqfeature_id, qk_id, 1, qval))


def load_Uniprot(adaptor, dbid, filename):
    """Load Uniprot records from filename into the BioSQL database.

    Returns the number of records loaded.
    """

    db_loader = UniprotLoader(adaptor, dbid)
    fp = file(filename)
    record_iterator = SeqIO.parse(fp, "swiss")

    num_records = 0
    while 1:
        try:
            cur_record = record_iterator.next()
        except StopIteration:
            break
        if cur_record is None:
            break
        num_records += 1
        db_loader.load_seqrecord(cur_record)
        if (num_records+1) % 10000 == 0 :
            print '      %d records loaded.' % (num_records+1)

    fp.close()
    return num_records
