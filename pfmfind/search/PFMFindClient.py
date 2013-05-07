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


"""
PFMFind client code.
"""

import os
import os.path
import imp
import xml.parsers.expat
from cStringIO import StringIO

import pfmfind
from pfmfind.search.matrix import SCORE
from pfmfind.search.matrix import POSITIONAL
from pfmfind.search.matrix import ScoreMatrix
from pfmfind.search.matrix import ProfileMatrix
from pfmfind.search.SearchClient import SearchClient
from pfmfind.search.DatabaseClient import DatabaseClient
from pfmfind.search.SearchServer import REL_SRCH

_MAX_SUBMITTED_SEARCHES = 20

def _write_xml_tag(fp, tag, indent=0):
    s = '  ' * indent
    fp.write("%s%s" % (s, tag))


class PFMFindClient(SearchClient, DatabaseClient):
    _default_plugin_dir = \
        os.path.join(pfmfind.__path__[0], 'plugins')

    def __init__(self, plugin_dir=None):
        SearchClient.__init__(self)
        DatabaseClient.__init__(self)
        self.search_lock = False
        self.init_plugins(plugin_dir)

    def init_plugins(self, plugin_dir=None):
        """
        Goes through plugin directories and gets all plugin modules.
        """

        self.plugin_dir = plugin_dir
        self.start_plugins = {}
        self.iteration_plugins = {}

        if plugin_dir:
            pdirs = [plugin_dir]
        else:
            pdirs = []
        pdirs.append(self._default_plugin_dir)

        suffixes = [s[0] for s in imp.get_suffixes()]

        for d in pdirs:
            dir_list = os.listdir(d)
            for filename in dir_list:
                for s in suffixes:
                    if filename[-len(s):] != s:
                        continue
                    name = os.path.split(filename[:-len(s)])[1]
                    if name in self.iteration_plugins or \
                       name in self.start_plugins or \
                       name == '__init__':
                        break

                    fp, pathname, description = \
                        imp.find_module(name, [d])
                    try:
                        M = imp.load_module(name, fp,
                                            pathname,
                                            description)
                    finally:
                        if fp: fp.close()

                    if hasattr(M, 'iteration') and \
                       hasattr(M, 'arg_list') and \
                       hasattr(M, 'get_matrix'):
                        print_info = getattr(M, 'print_info', None)
                        if M.iteration:
                            self.iteration_plugins[M.__name__] = \
                            (M.get_matrix, print_info, M.arg_list)
                        else:
                            self.start_plugins[M.__name__] = \
                            (M.get_matrix, print_info, M.arg_list)
                    break

    def plugin_matrix_str(self, PM, matrix_type, ctype):

        if PM == None:
            return "INVALID MATRIX"

        if matrix_type == SCORE:
            M = ScoreMatrix(PM)
        elif  matrix_type == POSITIONAL:
            M = ProfileMatrix(PM.pssm)
        else:
            return ""

        if ctype:
            M.conv_type = ctype
            M = M.matrix_conv()

        return M.__str__()

    def is_valid_fragment(self, qseq):
        """
        Checks that a given sequence is non-empty and consists only of
        letters from the standard amino acid alphabet.
        """

        valid_alphabet = 'ACDEFGHIKLMNPQRSTVWY'
        if not qseq or qseq.upper().translate('#'*256, valid_alphabet):
            return False
        return True

    def _set_matrices_to_hit_lists(self, srch_args, results):
        """
        Assigns a query matrix to each hit list.
        """
        # NOTE: Only the name is assigned at this stage because
        # storage of PSSMs is broken.

        for i, HL in enumerate(results):
            if HL == None: continue
            PM = srch_args[i][2]
            results[i].matrix_name = getattr(PM, 'name', None)
            results[i].matrix = None

    def _insert_results_single(self, jobs, srch_args, results,
                               job_index):
        """
        Inserts results into the database one by one.
        """

        self._set_matrices_to_hit_lists(srch_args, results)

        for i, HL in enumerate(results):
            key = job_index[i]
            length, fragment = key
            iteration, search_type, cutoff, plugin, plugin_args \
                       = jobs[key]

            if HL == None: continue

            if fragment in self.max_iters and\
                   iteration <= self.max_iters[fragment].get(length, -1):
                self.replace_search(length, iteration, fragment, HL)
            else:
                self.insert_search(length, iteration, fragment, HL)


    def _insert_results_batch(self, jobs, srch_args, results,
                              job_index):
        """
        Inserts results into database in batches, using psycopg
        routine copy_to.
        """

        results_batch = []
        self._set_matrices_to_hit_lists(srch_args, results)

        for i, HL in enumerate(results):
            key = job_index[i]
            length, fragment = key
            iteration, search_type, cutoff, plugin, plugin_args \
                       = jobs[key]


            if HL == None: continue

            if fragment in self.max_iters and\
                   iteration <= self.max_iters[fragment].get(length, -1):
                self.delete_search2(length, iteration, fragment)

            results_batch.append((length, iteration, fragment, HL))

        self.copy_searches(results_batch)

    def _insert_results(self, jobs, srch_args, results, job_index):
        """
        Inserts results into the database.
        """

        if self.has_copy:
            self._insert_results_batch(jobs, srch_args, results,
                                       job_index)
        else:
            self._insert_results_single(jobs, srch_args, results,
                                        job_index)


    def run_search_jobs(self, jobs):
        """
        Runs a collection of search jobs and stores the results into a
        database.
        """

        if not self.attached:
            raise RuntimeError, "Not attached to a search client."

        if self.search_lock:
            jobs = {}
            raise RuntimeError, "Search in progress"
        self.search_lock = True

        srch_args = []
        job_index = []
        j = 0
        k = 0
        for key, val in jobs.iteritems():

            length, fragment = key
            iteration, search_type, cutoff, plugin, plugin_args \
                       = val
            qseq = self.query_sequence[fragment: fragment + length]

            if not self.is_valid_fragment(qseq): # skip if invalid
                continue

            if iteration == 0:
                plugin_func = self.start_plugins[plugin][0]
            else:
                plugin_func = self.iteration_plugins[plugin][0]

            HL = self.select_lif_search(length, iteration-1, fragment)
            if HL == None and iteration: # skip if no search was done
                continue

            PM, matrix_type, conv_type = plugin_func(HL, *plugin_args)
            if PM == None: # skip if plugin rejects the arguments
                continue

            srch_args.append([search_type, qseq, PM, matrix_type,
                              cutoff, conv_type])
            job_index.append(key)
            j += 1
            k += 1
            # Submit searches in small batches
            if j == _MAX_SUBMITTED_SEARCHES:
                results = self.search(srch_args)
                self._insert_results(jobs, srch_args,
                                    results, job_index)
                srch_args = []
                job_index = []
                j = 0


        # Submit the remaining searches
        results = self.search(srch_args)
        self._insert_results(jobs, srch_args,
                             results, job_index)

        # Delete all jobs, even those that were rejected
        for key in jobs.keys():
            del(jobs[key])

        self.search_lock = False
        self.get_max_iters()
        return k

    def write_config(self, fp):
        """
        Writes the current configuration to a file in XML format.
        """

        fp.write("""<?xml version="1.0" encoding="UTF-8"?>\n""")
        fp.write("""<!DOCTYPE PFMF_config SYSTEM "%s">\n""" %\
             os.path.join(pfmfind.CONFIG_DATA_DIR, 'PFMFcf.dtd' ))
        i = 0
        _write_xml_tag(fp, "<PFMF_config>\n", i)
        i += 1
        _write_xml_tag(fp, "<Database_config>\n", i)
        i += 1

        # Database options
        fs = StringIO()
        fs.write('driver="%s"' % self.driver)
        for k in ['db', 'user', 'password', 'host', 'port']:
            if k in self.dbargs:
                fs.write(' %s="%s"' % (k, self.dbargs[k]))
        _write_xml_tag(fp, '<Database %s/>\n' % fs.getvalue(), i)

        # Schemas
        if self.db_schema:
            _write_xml_tag(fp, '<DbSchema name="%s"/>\n' %\
                self.db_schema, i)
        if self.PFMF_schema:
            _write_xml_tag(fp, '<PFMFSchema name="%s"/>\n' %\
                self.PFMF_schema, i)
        i -= 1
        _write_xml_tag(fp, "</Database_config>\n", i)

        # Index options
        if self.host and self.port:
            _write_xml_tag(fp, "<Index_config>\n", i)
            i += 1
            _write_xml_tag(fp, '<Connection host="%s" port="%s"/>\n'\
                % (self.host, self.port), i)
            i -= 1
            _write_xml_tag(fp, "</Index_config>\n", i)

        # Plugin options
        if self.plugin_dir:
            _write_xml_tag(fp, "<Plugin_config>\n", i)
            i += 1
            _write_xml_tag(fp, '<Plugin_dir>%s</Plugin_dir>\n'\
                % (self.plugin_dir), i)
            i -= 1
            _write_xml_tag(fp, "<<Plugin_config>\n", i)

        i -= 1
        _write_xml_tag(fp, "</PFMF_config>\n", i)


    def read_config(self, fp):
        """
        Reads the configuration from an XML file and attempts to bring
        the client to the same state.
        """

        if self.conn or self.attached:
            raise RuntimeError("Client must be in the initial state"\
                " (no open connections) before reading configuration.")

        p = xml.parsers.expat.ParserCreate()
        p.returns_unicode = 0
        p.StartElementHandler = self._start_element
        p.EndElementHandler = self._end_element
        p.CharacterDataHandler = self._char_data

        for var in ['_sflag', '_fs', '_tmp_dbargs', '_tmp_db_schema',
                    '_tmp_PFMF_schema', '_tmp_ixargs',
                    '_tmp_plugin_dir']:
            setattr(self, var, None)

        p.ParseFile(fp)

        if 'port' in self._tmp_dbargs:
            self._tmp_dbargs['port'] = int(self._tmp_dbargs['port'])

        # Try to open connections - quietly leave if failing
        try:
            try:
                # Try to connect to database - will throw exception if
                # something goes wrong
                if self._tmp_dbargs:
                    self.open(**self._tmp_dbargs)
                else:
                    raise RuntimeError("Missing arguments.")

                # Set database schema search path
                schema_args = dict()
                if self._tmp_db_schema:
                    schema_args['db_schema'] = self._tmp_db_schema
                if hasattr(self, '_tmp_PFMF_schema'):
                    schema_args['PFMF_schema'] = self._tmp_PFMF_schema
                self.set_schema(**schema_args)

                # Attach to index
                if self._tmp_ixargs:
                    self._tmp_ixargs['port'] =\
                        int(self._tmp_ixargs['port'])
                    self.attach(**self._tmp_ixargs)

                # Initialise plugins
                if self._tmp_plugin_dir:
                    self.init_plugins(self._tmp_plugin_dir)

            finally:
                # Clean up
                for var in ['_sflag', '_fs', '_tmp_dbargs',
                            '_tmp_db_schema', '_tmp_PFMF_schema',
                            '_tmp_ixargs', '_tmp_plugin_dir']:
                    delattr(self, var)
        except:
            raise
            pass

    def _start_element(self, name, attrs):
        if name == 'Database':
            self._tmp_dbargs = attrs
        elif name == 'DbSchema':
            self._tmp_db_schema = attrs['name']
        elif name == 'PFMFSchema':
            self._tmp_PFMF_schema = attrs['name']
        elif name == 'Connection':
            self._tmp_ixargs = attrs
        elif name == 'Plugin_dir':
            self._sflag = True
            self._fs = StringIO()

    def _end_element(self, name):
        if name == 'Sql_dir':
            self._sflag = False
            self._tmp_plugin_dir = self._fs.getvalue()

    def _char_data(self, data):
        if self._sflag:
            self._fs.write(data)



    def iterated_rel_search(self, experiment_name,
                            experiment_description,
                            sequence, description,
                            frag_len,
                            start_plugin, start_params,
                            iter_plugin, iter_params,
                            cutoff_Evalues):
        """
        Runs automated iterated search for all fragments of a
        certain length of a given sequence.

        """

        if not self.attached:
            raise RuntimeError, "Not attached to a search client."

        tot_frags = len(sequence) + 1 - frag_len

        # Create experiment if it doesn't already exist
        eid = self.get_experiment_id(experiment_name)

        if eid == None:
            eid = self.create_experiment(experiment_name,
                                         experiment_description,
                                         sequence,
                                         description,
                                         frag_len,
                                         frag_len+1)

        # Set experiment
        self.set_current_experiment(eid)


        # ** Prepare the first iteration **
        # Get the matrix first (it is the same for all searches)
        plugin_func = self.start_plugins[start_plugin][0]
        PM, matrix_type, conv_type = plugin_func(None, *start_params)


        fragments = []
        srch_args = []
        Eval = cutoff_Evalues[0]
        for f in xrange(tot_frags):
            qseq = self.query_sequence[f: f + frag_len]
            if not self.is_valid_fragment(qseq): continue
            fragments.append(f)
            srch_args.append([REL_SRCH, qseq, PM, matrix_type,
                              Eval, conv_type])

        plugin_func = self.iteration_plugins[iter_plugin][0]

        # *** Run all iterations ***
        final_results = []
        step = _MAX_SUBMITTED_SEARCHES
        for i in xrange(len(cutoff_Evalues)):
            profiles = []
            new_fragments = []
            for j in xrange(0, len(srch_args), step):
                # Run batch search
                results = self.search(srch_args[j:j+step])

                # Immediately process the results to obtain profiles
                if i < len(cutoff_Evalues) - 1:
                    for k,HL in enumerate(results):
                        if HL == None: continue
                        PM_data = plugin_func(HL, *iter_params)
                        if PM_data[0] != None:
                            profiles.append(PM_data)
                            new_fragments.append(fragments[j+k])

                final_results.extend(results)
                del(results)

            # Put results into database - we need to generate
            # standard 'jobs' etc.
            job_index = [(frag_len, f) for f in fragments]
            jobs = {}.fromkeys(job_index,
                               (i, REL_SRCH, cutoff_Evalues[i],
                                iter_plugin, iter_params))
            self._insert_results(jobs, srch_args,
                                 final_results, job_index)

            final_len = len(final_results)
            final_results = []

            if i < len(cutoff_Evalues) - 1:
                # Create arguments for next iteration
                fragments = new_fragments
                srch_args = [[REL_SRCH,
                              self.query_sequence[f: f + frag_len],
                              profiles[k][0],
                              profiles[k][1],
                              cutoff_Evalues[i+1],
                              profiles[k][2]] \
                             for k,f in enumerate(fragments) ]

        self.reset_current_experiment()
        return eid, final_len
