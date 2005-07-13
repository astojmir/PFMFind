
from Bio.SubsMat import MatrixInfo
from ShortFrags.Expt.matrix import QUASI, MAX, AVG, SCORE 

_MATRIX_CTYPE = {'None': 0, 'Quasi': QUASI, 'Avg': AVG, 'Max': MAX} 


iteration = False
arg_list = [('Matrix Name', MatrixInfo.available_matrices, 'blosum62'),
            ('Conversion', _MATRIX_CTYPE.keys(), 'None'),
            ]


def _filter_non_standard_letters(S):
    for a, b in S.keys():
        if a not in "ACDEFGHIKLMNPQRSTVWY" or \
           b not in "ACDEFGHIKLMNPQRSTVWY":
            del(S[(a,b)])

def get_matrix(HL, matrix_name, conv_type):

    S = MatrixInfo.__dict__[matrix_name]
    _filter_non_standard_letters(S)
    matrix_type = SCORE
    ctype = _MATRIX_CTYPE[conv_type]
    return S, matrix_type, ctype

