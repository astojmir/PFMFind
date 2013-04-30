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



from Bio.SubsMat import MatrixInfo
from pfmfind.Expt.matrix import QUASI, MAX, AVG, SCORE
from pfmfind.Expt.matrix import SubstitutionMatrix

_MATRIX_CTYPE = {'None': 0, 'Quasi': QUASI, 'Avg': AVG, 'Max': MAX} 


iteration = False
arg_list = [('Matrix Name', MatrixInfo.available_matrices, 'blosum62'),
            ('Conversion', _MATRIX_CTYPE.keys(), 'None'),
            ]

_std_alphabet_map = {}.fromkeys(list("ACDEFGHIKLMNPQRSTVWY"))


def _filter_non_standard_letters(S):
    for a, b in S.keys():
        if a not in _std_alphabet_map or \
           b not in _std_alphabet_map:
            del(S[(a,b)])

def get_matrix(HL, matrix_name, conv_type):

    S = SubstitutionMatrix()
    S.update(MatrixInfo.__dict__[matrix_name])
    S.name = matrix_name
    
##     S = MatrixInfo.__dict__[matrix_name]
    _filter_non_standard_letters(S)
    matrix_type = SCORE
    ctype = _MATRIX_CTYPE[conv_type]
    return S, matrix_type, ctype

