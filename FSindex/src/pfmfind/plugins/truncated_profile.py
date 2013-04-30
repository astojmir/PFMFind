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

from pfmfind.search.DirichletInfo import NAMES
from default_profile import get_matrix as default_get_matrix
from default_profile import print_info as default_print_info

iteration = True
arg_list = [('Scale', 'real', 2.0),
            ('Weighting', ['None', 'Henikoff'], 'None'),
            ('Regulariser', NAMES, 'recode3.20comp'),
            ('Minimum Hits', 'numeric', 30),
            ]

def get_matrix(HL, scale, weight_type, dirichlet_type, min_hits):
    """
    Returns default profile if the number of hits is not less than
    min_hits. 
    """

    if len(HL) < min_hits:
        return None, 0, 0

    return default_get_matrix(HL, scale, weight_type, dirichlet_type) 

def print_info(HL, scale, weight_type, dirichlet_type, min_hits):
    """
    Print the same information as default_profile.
    """

    if len(HL) < min_hits:
        return ""
    
    return default_print_info(HL, scale, weight_type, dirichlet_type) 


