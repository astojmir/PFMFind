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
Wrapper around FSindex library.
"""

import os, os.path
from pfmfind.search import FS
from pfmfind.search.db import db
from pfmfind.search.hit_list import HitList

FS_BINS = FS.FS_BINS
SUFFIX_ARRAY = FS.SUFFIX_ARRAY
SEQ_SCAN = FS.SEQ_SCAN

SARRAY = FS.SARRAY
DUPS_ONLY = FS.DUPS_ONLY
FULL_SCAN = FS.FULL_SCAN

def _get_db(I):
    return db(FS.Index_s_db_get(I), new=False, own=False)

def _get_ix_data(I):
    return FS.Index_get_data(I)

class FSIndex(object):
    def __init__(self, filename, sepn=None, use_sa=1, print_flag=0):
        self.thisown = 0
        if sepn == None:
            sepn = []

        # Need to change directory because the C code works only on UNIX
        old_workdir = os.getcwd()
        new_workdir, new_filename = os.path.split(filename)
        if len(new_workdir):
            os.chdir(new_workdir)

        self.this = FS.new_Index(new_filename, sepn, use_sa, print_flag)

        # Return to old working directory
        if len(new_workdir):
            os.chdir(old_workdir)

        self.thisown = 1
        self.__dict__.update(FS.Index_get_data(self))

    def __del__(self):
        if self.thisown:
            FS.delete_Index(self)

    def save(self, filename):
        return FS.Index_save(self, filename)

    def __str__(self):
        return FS.Index___str__(self)

    def seq2bin(self, seq):
        return FS.Index_seq2bin(self, seq)

    def print_bin(self, bin, options=1):
        return FS.Index_print_bin(self, bin, options)

    def print_stats(self, options=3):
        return FS.Index_print_stats(self, options)

    def get_bin_size(self, bin):
        return FS.Index_get_bin_size(self, bin)

    def get_unique_bin_size(self, bin):
        return FS.Index_get_unique_bin_size(self, bin)

    def rng_srch(self, qseq, M, rng, stype=FS_BINS,
                 ptype=SARRAY, qdef=""):


        hits_dict = FS.Index_rng_srch(self, qseq, M, rng, M.conv_type,
                                      stype, ptype, qdef)
        return HitList(hits_dict)

    def kNN_srch(self, qseq, M, kNN, stype=FS_BINS,
                 ptype=SARRAY, qdef=""):
        hits_dict = FS.Index_kNN_srch(self, qseq, M, kNN,
                                      stype, ptype, qdef)
        return HitList(hits_dict)

    def threaded_search(self, srch_args):
        results = FS.Index_threaded_search(self, srch_args)
        return [HitList(r) for r in results]

    db = property(_get_db)
    ix_data = property(_get_ix_data)
