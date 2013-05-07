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


import FS

class db(object):
    def __init__(self, db_name, new=True, own=True):
        self.thisown = 0
        if new == True:
            self.this = FS.new_db(db_name)
            self.thisown = 1
        elif own == True:
            # SWIG string - thisown = 1
            self.this = db_name
            self.thisown = 1
        else:
            # SWIG string - thisown = 0
            self.this = db_name
            self.thisown = 0
            
    def __del__(self):
        if self.thisown:
            FS.delete_db(self) 

    def __str__(self):
      return 'FASTA SEQUENCE DATABASE\nFile: %s\nLength: %d\nSequences: %d\n' % (self.db_name, self.length , self.no_seq)    

    def get_seq(self, i):
        return FS.db_get_seq(self, i)

    def get_def(self, i):
        return FS.db_get_def(self, i)

    def get_frag(self, i, a, b):
        return FS.db_get_frag(self, i, a, b)

    # Sequence iterator
    def sequences(self):
        for i in xrange(self.no_seq):
            yield self.get_seq(i)
            
    # Fragment iterator
    def fragments(self, length):
        for i in xrange(self.no_seq):
            seq = self.get_seq(i)
            slen = len(seq)
            for j in xrange(slen-length+1):
                yield seq[j:j+length]
        
    db_name = property(FS.db_db_name_get)
    length = property(FS.db_length_get)
    no_seq = property(FS.db_no_seq_get)
