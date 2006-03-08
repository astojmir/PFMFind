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


# Module ShortFrags
__all__ = ['GUI', 'Expt', 'Setup', 'plugins']


import sys
from os.path import join

#CONFIG_DATA_DIR = join(sys.exec_prefix, 'share/ShortFrags/setup_config')
#SQL_DATA_DIR =  join(sys.exec_prefix, 'share/ShortFrags/sql-schema')
CONFIG_DATA_DIR = join('share', 'ShortFrags', 'setup_config')
SQL_DATA_DIR = join('share', 'ShortFrags', 'sql-schema')
