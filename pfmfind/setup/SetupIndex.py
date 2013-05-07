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


import sys
import os
import os.path
import time
import xml.parsers.expat
from cStringIO import StringIO
from BioSQL import BioSeqDatabase
from pfmfind.search.index import FSIndex


_Create_Index_Script = \
"""
import os
from pfmfind.search.index import FSIndex

index_dir = '%s'
fasta_name = '%s'
index_name = '%s'
pttn = %s

os.chdir(index_dir)
I = FSIndex(fasta_name, pttn, 0, 0)
I.save(index_name)
"""

_COL_WIDTH = 60

class PFMF_IndexCreator(object):
    def __init__(self):
        self._setup_flag = False
        self._sflag = False

        self.dbargs = None
        self.index_dir = './'
        self.dataset = []

    def _start_element(self, name, attrs):

        if name == 'PFMF_index_setup':
            self._setup_flag = True

        if not self._setup_flag: return

        if name == 'Database':
            self.dbargs = attrs
        elif name == 'Index_dir' or name == 'Partition':
            self._sflag = True
            self._fs = StringIO()
        elif name == 'Dataset':
            self.dataset.append(([],attrs['name'],
                                 attrs.get('schema', None),
                                 attrs.get('namespace', None),
                                 int(attrs.get('max_residues', '2147483648'))))
        elif name == 'Index':
            self.dataset[-1][0].append(([], int(attrs['length'])))

    def _end_element(self, name):

        if not self._setup_flag: return

        if name == 'PFMF_db_setup':
            self._setup_flag = False
        elif name == 'Index_dir':
            self._sflag = False
            self.index_dir = self._fs.getvalue()
        elif name == 'Schema':
            self.schema_name = attrs.get('name', None)
        elif name == 'Partition':
            self._sflag = False
            pttn = self._fs.getvalue()
            if len(self.dataset[-1][0][-1][0]) < self.dataset[-1][0][-1][1]:
                self.dataset[-1][0][-1][0].append(pttn)
        elif name == 'Index':
            pttn = self.dataset[-1][0][-1][0][-1]
            n = self.dataset[-1][0][-1][1] - len(self.dataset[-1][0][-1][0])
            for i in xrange(n):
                self.dataset[-1][0][-1][0].append(pttn)

    def _char_data(self, data):

        if self._setup_flag and self._sflag:
            self._fs.write(data)


    def parse_config(self, fp):
        """
        Parses the XML configuration file.
        """

        p = xml.parsers.expat.ParserCreate()
        p.returns_unicode = 0
        p.StartElementHandler = self._start_element
        p.EndElementHandler = self._end_element
        p.CharacterDataHandler = self._char_data

        p.ParseFile(fp)


    def create_fasta_files(self):
        """
        Creates FASTA files from database sequences.
        """

        print "Creating FASTA datasets."
        print "  Opening database."
        server = BioSeqDatabase.open_database(**self.dbargs)
        cur = server.adaptor.cursor

        for indexes, name, schema, namespace, max_res in self.dataset:

            print "  Creating dataset %s." % name

            if schema:
                print "    Setting schema %s." % schema
                cur.execute("SET search_path TO %s" % schema)
                server.adaptor.commit()

            # Skip the whole FASTA file creation if the first file
            # already exists
            fasta_name = name + '%3.3d.fas' % 0
            if fasta_name in os.listdir(self.index_dir):
                print "    File %s present - skipping all." % fasta_name
                continue

            # Execute the SQL query (specific namespace or all sequences)
            if namespace:
                dbid = server[namespace].dbid
                sql = """SELECT e.name || ' (' || e.accession || ') '
                         || e.description AS header, s.seq AS residues
                         FROM bioentry e, biosequence s
                         WHERE e.bioentry_id = s.bioentry_id AND
                         e.biodatabase_id = %s ORDER BY e.name;"""
                cur.execute(sql, (dbid,))
            else:
                sql = """SELECT e.name || ' (' || e.accession || ') '
                         || e.description AS header, s.seq AS residues
                         FROM bioentry e, biosequence s
                         WHERE e.bioentry_id = s.bioentry_id
                         ORDER BY e.name;"""
                cur.execute(sql)


            # Number of residues in the current file set so that new file is opened
            num_residues = max_res + 1
            # Counter and file are dummy
            file_counter = -1
            fp = StringIO()

            while 1:
                res = cur.fetchone()
                if not res: break
                title = res[0][:_COL_WIDTH-1]
                sequence = res[1]

                if num_residues + len(sequence) > max_res:
                    num_residues = 0
                    file_counter += 1
                    fp.close()
                    fasta_name = name + '%3.3d.fas' % file_counter
                    print "    Creating file %s." % fasta_name
                    fasta_path = os.path.join(self.index_dir,
                                              fasta_name)
                    # Open file in binary mode:
                    # We write in UNIX format with line separator '\n'
                    fp = file(fasta_path, 'wb')

                # Now write the sequence
                fp.write('>%s\n' % title)
                i = 0
                while i < len(sequence):
                    fp.write('%s\n' % sequence[i:i+_COL_WIDTH])
                    i += _COL_WIDTH

                num_residues += len(sequence)
            fp.close()

        server.adaptor.close()

    def create_indexes(self):

        def now():
            return time.ctime(time.time())

        print "Creating indexes."
        old_path = os.getcwd()
        os.chdir(self.index_dir)

        # Change resource limits (memory) - only if POSIX platform
	if os.name == 'posix':
            print "  Changing resource limits."
            import resource
            if 'RLIMIT_DATA' in dir(resource):
                lim = resource.getrlimit(resource.RLIMIT_DATA)
                resource.setrlimit(resource.RLIMIT_DATA, (lim[1],lim[1]))
            if 'RLIMIT_RSS' in dir(resource):
                lim = resource.getrlimit(resource.RLIMIT_RSS)
                resource.setrlimit(resource.RLIMIT_RSS, (lim[1],lim[1]))
            if 'RLIMIT_MEMLOCK' in dir(resource):
                lim = resource.getrlimit(resource.RLIMIT_MEMLOCK)
                resource.setrlimit(resource.RLIMIT_MEMLOCK, (lim[1],lim[1]))

        for indexes, name, schema, namespace, max_res in self.dataset:
            print "  Dataset %s." % name

            for pttn, m in indexes:
                j = 0;  # FASTA file counter
                while 1:
                    fasta_name = name + '%3.3d.fas' % j
                    index_name = name + '%2.2d%3.3d.ix' % (m,j)
                    if fasta_name not in os.listdir(os.getcwd()):
                        break
                    if index_name in os.listdir(os.getcwd()):
                        j += 1
                        continue
                    print "    Creating index %s (%s)." % (index_name, now())

                    # It appears there is a memory leak (Python related)
                    # if more than one index is created in the same
                    # process. Therefore, we create each index in a
                    # separate process.

                    # For UNIX, use the standard fork() way:
                    if os.name == 'posix':
                        pid = os.fork()
                        if pid == 0:
                            I = FSIndex(fasta_name, pttn, 0, 0)
                            I.save(index_name)
                            os._exit(0)
                        os.wait()

                    # For other platforms (Windows), open a pipe to a
                    # new python interpreter and send the code and
                    # variables as stdin.
                    else:
                        args = (self.index_dir, fasta_name, index_name, pttn.__repr__())
                        pipe_in, pipe_out = os.popen4(sys.executable)
                        pipe_in.write(_Create_Index_Script % args)
                        pipe_in.close()

                        # Print the output in case there are errors.
                        for line in pipe_out:
                            print line
                        pipe_out.close()

                    j += 1

        os.chdir(old_path)
        print "Finished creating indexes at %s." % now()

