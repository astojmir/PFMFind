import FS
import sys
import types
from DirichletMix import BKGRND_PROBS as bg_dict

# Mtype values
SCORE = FS.SCORE
POSITIONAL = FS.POSITIONAL

# Stype values
SIMILARITY = FS.SIMILARITY
DISTANCE = FS.DISTANCE

# Conversion values
POS = FS.POS
QUASI = FS.QUASI
MAX = FS.MAX
AVG = FS.AVG


class MatrixColumn(object):
    def __init__(self, col, this, alphabet):
        self._col = col
        self.this = this
        self.thisown = 0
        self._abet = list(alphabet)
    
    def __len__(self):
        return len(self._abet)

    def __getitem__(self, key):
        if type(key) is not types.StringType and len(key) != 1:
            raise TypeError, 'Invalid key type.'
        elif key not in self._abet:
            raise IndexError, 'Invalid letter.'
        return self.get_func(self, self._col, key)
        
    def __setitem__(self, key, value):
        if type(key) is not types.StringType and len(key) != 1:
            raise TypeError, 'Invalid key type.'
        elif key not in self._abet:
            raise IndexError, 'Invalid letter.'
        self.set_func(self, self._col, key, value)

    def __iter__(self):
        return self.iteritems()

    def has_key(self, key):
        return key in self._abet

    def items(self):
        return zip(self.keys(), self.values())

    def keys(self):
        return self._abet

    def values(self):
        return [self.__getitem__(a) for a in self._abet]
    
    def get(self, key, x=None):
        if self.has_key(key):
            return self.__getitem__(key)
        else:
            return x
    
    def iteritems(self):
        return iter(self.items)

    def iterkeys(self):
        return iter(self.keys)

    def itervalues(self):
        return iter(self.values)
    
class ScoreColumn(MatrixColumn):
    def __init__(self, col, this, alphabet):
        MatrixColumn.__init__(self, col, this, alphabet)
        self.get_func = FS.Smatrix_score_get_item
        self.set_func = FS.Smatrix_score_set_item

class ProfileColumn(MatrixColumn):
    def __init__(self, col, this, alphabet):
        MatrixColumn.__init__(self, col, this, alphabet)
        self.get_func = FS.Smatrix_profile_get_item
        self.set_func = FS.Smatrix_profile_set_item


        

class Smatrix(object): # Base class for matrices
    def __del__(self):
        FS.delete_Smatrix(self)

    def __str__(self):
        return FS.Smatrix___str__(self)

    def __iter__(self):
        return self.iteritems()

    def has_key(self, key):
        return key in self.keys()

    def items(self):
        return zip(self.keys(), self.values())

    def get(self, key, x=None):
        if self.has_key(key):
            return self.__getitem__(key)
        else:
            return x
    
    def iteritems(self):
        return iter(self.items)

    def iterkeys(self):
        return iter(self.keys)

    def itervalues(self):
        return iter(self.values)

    def fprint(fp=sys.stdout):
        return FS.Smatrix_fprint(self, fp)

    def score2rel(self, score, length=1):
        return score / (length * self.mean)

    def rel2score(self, rel, length=1):
        return rel * length * self.mean


    Mtype = property(FS.Smatrix_Mtype_get)
    
    Stype = property(FS.Smatrix_Stype_get)

    conv_type = property(FS.Smatrix_conv_type_get, FS.Smatrix_conv_type_set)


class ScoreMatrix(Smatrix):
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
            # Instance of ScoreMatrix
        elif isinstance(inarg, ScoreMatrix):
            self.this = FS.Smatrix_copy(inarg)
            self.thisown = 1
        else:
            raise TypeError, "Constructor only accepts"
        " a filename, a dictionary"
        " and an instance of ScoreMatrix." 

        self.mean = self._calc_mean()
        self.alphabet = FS.Smatrix_alphabet_get(self)

    def __len__(self):
        return len(self.alphabet)

    def __getitem__(self, key):
        if type(key) is not types.StringType and len(key) != 1:
            raise TypeError, 'Invalid key type.'
        elif key not in self.alphabet:
            raise IndexError, 'Invalid letter.'
        return ScoreColumn(key, self.this, self.alphabet)

    def keys(self):
        return list(self.alphabet)

    def values(self):
        return [self.__getitem__(a) for a in self.alphabet]

    def eval_score(self, seq1, seq2):
        return FS.Smatrix_eval_score(self, seq1, seq2)

    def range_conv(self, seq, r):
        return FS.Smatrix_range_conv(self, seq, r)

    def matrix_conv(self, seq=""):
        if (self.conv_type & 1) == POS:
            return ProfileMatrix(FS.Smatrix_matrix_conv(self, seq),
                                 new=False)
        else:
            return ScoreMatrix(FS.Smatrix_matrix_conv(self, seq),
                               new=False)
    def _calc_mean(self):
        sum = 0
        for i in self.alphabet:
            for j in self.alphabet;
            sum += bg_dict[i] * bg_dict[j] * self[i][j] 
        return sum 



class ProfileMatrix(Smatrix):
    def __init__(self, inarg, Stype=FS.SIMILARITY, new=True):
        if type(inarg) is types.StringType:
            # SWIG string
            self.this = inarg
            self.thisown = 1
        elif type(inarg) is types.ListType:
            # List - Biopython PSSM
            self.this = FS.Smatrix_pssm(inarg)
            self.thisown = 1
        elif isinstance(inarg, ProfileMatrix):
            # Instance of ProfileMatrix
            self.this = FS.Smatrix_copy(inarg)
            self.thisown = 1
        else:
            raise TypeError, "Constructor only accepts"
        " a Biopython PSSM"
        " and an instance of ProfileMatrix." 

        self.mean = self._calc_mean()
        self.alphabet = FS.Smatrix_alphabet_get(self)
        self._len = FS.Smatrix_len_get(self)
        self.qseq = FS.Smatrix_qseq_get(self)
    def __len__(self):
        return self._len

    def __getitem__(self, key):
        if type(key) is not types.IntType:
            raise TypeError, 'Invalid key type.'
        elif key not in self.keys():
            raise IndexError, 'Invalid column.'
        return ProfileColumn(key, self.this, self.alphabet)

    def keys(self):
        return range(self._len)

    def values(self):
        return [self.__getitem__(i) for i in range(self._len)]

    def eval_score(self, seq1):
        seq2 = 'X'*len(seq1)
        return FS.Smatrix_eval_score(self, seq1, seq2)

    def range_conv(self, r):
        seq = 'X'*self._len
        return FS.Smatrix_range_conv(self, seq, r)

    def matrix_conv(self):
        return ProfileMatrix(FS.Smatrix_matrix_conv(self, ''),
                             new=False)





