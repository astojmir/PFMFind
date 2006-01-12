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


from ftplib import FTP
from bisect import bisect_left, bisect_right
import os, os.path, sys, tarfile 


# Stage constants
DOWNLOAD = 0
CREATE_TABLES = 1
COPY_TABLES = 2

class TaxonLoader(object):
    """
    Replaces the load_ncbi_taxonomy.pl script from BioSQL
    distribution. PostgreSQL specific. Deletes all entries from the taxon
    and taxon_name tables and re-creates them - this assumes no foreign
    key constraints. 
    """

    def __init__(self, adaptor, taxon_dir, stage=DOWNLOAD):

        self.adaptor = adaptor
        self.taxon_dir = taxon_dir
        self.stage = stage

    def load_ncbi_taxonomy(self):

        # Change to taxon_dir
        print 'Entering %s directory.' % self.taxon_dir
        try:
            os.mkdir(self.taxon_dir)
        except OSError:
            pass
        self.old_dir = os.getcwd()
        os.chdir(self.taxon_dir)

        # Download data if needed
        if self.stage <= DOWNLOAD:
            print 'Downloading taxdump.tar.gz from NCBI ftp site.'
            self._download_taxondb()

        if self.stage <= CREATE_TABLES:
            # Create taxon table
            print 'Creating taxon.tbl file.'
            self._get_nodes()

            # Create taxon_name table
            print 'Creating taxon_name.tbl file.'
            self._get_names()

        if self.stage <= COPY_TABLES:
            # Load  taxon_name table
            print 'Loading taxon table.'
            self._load_nodes()

            # Load  taxon_name table
            print 'Loading taxon_name table.'
            self._load_names()

        os.chdir(self.old_dir)
       

    def _download_taxondb(self):
        # download
        ftp = FTP('ftp.ncbi.nlm.nih.gov')
        ftp.login()
        ftp.cwd('/pub/taxonomy')
        ftp.retrbinary('RETR taxdump.tar.gz', open('taxdump.tar.gz', 'wb').write)
        ftp.quit()

        # extract taxdump.tar.gz
        tar = tarfile.open("taxdump.tar.gz", "r:gz")
        for tarinfo in tar:
            tar.extract(tarinfo)
        tar.close()
        os.remove("taxdump.tar.gz")

    def _get_names(self):
        fpin = file('names.dmp')
        fpout = file('taxon_name.tbl', 'w')

        for line in fpin:
            cols = line.split('\t|')
            new_line = ''.join((cols[0], cols[1], cols[3]))
            fpout.write(new_line)
            fpout.write('\n')

        fpin.close()
        fpout.close()

    def _get_nodes(self):
        self._nodes = list()
        # First import all nodes into memory
        print '  Importing nodes into memory.'
        fpin = file('nodes.dmp')
        for line in fpin:
            cols = line.split('\t|\t')
            # Creating node row - note that first col is parent_id 
            node = [cols[1], cols[0], cols[2],
                    cols[6], cols[8], None, None]
            self._nodes.append(node)
        fpin.close()
        print '  Sorting nodes.'
        self._nodes.sort()
        self._parents = [node[0] for node in self._nodes]

        # Get nested values.
        print '  Calculating nested values.'
        self._ctr = 0
        # First node is a root - everything is under it
        self._traverse_subtree(self._nodes[0])

        # Now write to file
        print '  Writing to file.'
        fpout = file('taxon.tbl', 'w')
        for node in self._nodes:
            new_line = '\t'.join((node[1], node[1], node[0], node[2],
                                  node[3], node[4], node[5], node[6])) 
            fpout.write(new_line)
            fpout.write('\n')
        fpout.close()

    def _traverse_subtree(self, node):
        if node[5] != None:
            return
        
        self._ctr += 1
        node[5] = str(self._ctr)  # left value
        C = self._get_children(node)
        for child in C:
            if child[1] != node[1]:
               self._traverse_subtree(child)
        self._ctr += 1
        node[6] = str(self._ctr)  # right value

    def _get_children(self, node):
        p = node[1]
        a = bisect_left(self._parents, p)
        b = bisect_right(self._parents, p)
        return self._nodes[a:b]
    
    def _load_names(self):
        cur = self.adaptor.cursor

        # First delete all old entries
        cur.execute("DELETE FROM taxon_name")
        self.adaptor.commit()
        self.adaptor.conn.autocommit(1)
        cur.execute("VACUUM FULL taxon_name")
        self.adaptor.conn.autocommit(0)

        # Now copy our new table
        path = os.path.join(os.getcwd(), 'taxon_name.tbl')
        cur.execute("COPY taxon_name FROM %s", (path,))
        cur.execute("ANALYZE taxon_name")
        self.adaptor.commit()

    def _load_nodes(self):
        cur = self.adaptor.cursor

        # First delete all old entries
        cur.execute("DELETE FROM taxon")
        self.adaptor.commit()
        self.adaptor.conn.autocommit(1)
        cur.execute("VACUUM FULL taxon")
        self.adaptor.conn.autocommit(0)

        # Now copy our new table
        path = os.path.join(os.getcwd(), 'taxon.tbl')
        cur.execute("COPY taxon FROM %s", (path,))
        cur.execute("ANALYZE taxon")
        self.adaptor.commit()

