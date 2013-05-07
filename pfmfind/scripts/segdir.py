#!/usr/bin/env python

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
Runs SEG program to mask low complexity regions of all FASTA files
(identified by extension .fas) in the current directory. All FASTA
files are replaced by masked files.

"""

__version__ = "$Revision: 1.0 $"
__date__ = "$Date: 2006-03-08 16:50:09 -0500 (Wed, 08 Mar 2006) $"
__author__ = "Aleksandar Stojmirovic"


import os, os.path

fasta_list = [f for f in os.listdir(os.getcwd()) if os.path.splitext(f)[1] == '.fas']
command = 'seg %s -x > temp.fas; mv temp.fas %s'

for f in fasta_list:
    print "Running SEG on %s" % f
    os.system(command % (f,f))

