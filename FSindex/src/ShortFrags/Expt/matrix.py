import FS
import sys
import types

SCORE = FS.SCORE
POSITIONAL = FS.POSITIONAL

SIMILARITY = FS.SIMILARITY
DISTANCE = FS.DISTANCE

POS = FS.POS
QUASI = FS.QUASI
MAX = FS.MAX
AVG = FS.AVG


class Smatrix(object):
    def __init__(self, inarg, Stype=FS.SIMILARITY, new=True):
        if type(inarg) is types.StringType:
            if new == True:
                # Load from file 
                self.this = FS.Smatrix_load(inarg)
                self.thisown = 1
            else:
                # SWIG string
                self.this = inarg
                self.thisown = 1
                
            # Dictionary - Biopython matrix 
        elif isinstance(inarg, dict):
            self.this = FS.new_Smatrix(inarg, Stype)
            self.thisown = 1
            # List - Biopython PSSM
        elif type(inarg) is types.ListType:
            self.this = FS.Smatrix_pssm(inarg)
            self.thisown = 1
            # Instance of Smatrix
        elif isinstance(inarg, Smatrix):
            self.this = FS.Smatrix_copy(inarg)
            self.thisown = 1
        else:
            raise TypeError

    def __del__(self):
        FS.delete_Smatrix(self)

    def __str__(self):
        return FS.Smatrix___str__(self)

    def fprint(fp=sys.stdout):
        return FS.Smatrix_fprint(self, fp)

    def eval_score(self, seq1, seq2):
        return FS.Smatrix_eval_score(self, seq1, seq2)

    def range_conv(self, seq, r):
        return FS.Smatrix_range_conv(self, seq, r)

    def item_conv(self, c, i, j):
        return FS.Smatrix_item_conv(self, c, i, j)

    def matrix_conv(self, seq=None):
        return Smatrix(FS.Smatrix_matrix_conv(self, seq),
                       new=False)

    Mtype = property(FS.Smatrix_Mtype_get)
    
    Stype = property(FS.Smatrix_Stype_get)

    len = property(FS.Smatrix_len_get)

    qseq = property(FS.Smatrix_qseq_get)

    alphabet = property(FS.Smatrix_alphabet_get)

    conv_type = property(FS.Smatrix_conv_type_get, FS.Smatrix_conv_type_set)
