from string import split, join, strip
import sys, os, copy
from transcendental import gamma, log, array, product, beta
from Bio.Align.AlignInfo import PSSM
from cStringIO import StringIO

# Background amino acid frequencies - as dictionary
## BKGRND_PROBS = {'A': 7.5610, 'C': 1.6504, 'D': 5.0605,
##                 'E': 6.1824, 'F': 4.0928, 'G': 6.9696,
##                 'H': 2.2818, 'I': 5.7749, 'K': 5.6397,
##                 'L': 9.5874, 'M': 2.3668, 'N': 4.3354,
##                 'P': 5.0832, 'Q': 3.9400, 'R': 5.3245,
##                 'S': 7.3014, 'T': 5.7451, 'V': 6.5107,
##                 'W': 1.3716, 'Y': 3.1464}

BKGRND_PROBS = {'A': 0.075666295724018673, 'C': 0.016516288118360062,
                'E': 0.061870031303289649, 'D': 0.050642678152545494,
                'G': 0.069747892431969383, 'F': 0.040958473103989368,
                'I': 0.05779199724595098,  'H': 0.022834989231988598,
                'K': 0.056438990609013105, 'M': 0.023685622102844514,
                'L': 0.095945383365223735, 'N': 0.043386279391867545,
                'Q': 0.039429335425556614, 'P': 0.050869847166291719,
                'S': 0.073068362861969299, 'R': 0.053284643774968568,
                'T': 0.057493775368874442, 'W': 0.013726212301952652,
                'V': 0.065155475673901384, 'Y': 0.031487426645424192}

def _print_vector(v, indent=0, floats_per_row=10, fp=sys.stdout,
                  indent_first_row=True):
    if indent_first_row: indent1 = indent
    else: indent1 = 0
    fp.write("%*s[" % (indent1, ""))
    k = range(len(v))
    for i in range(0, len(v), floats_per_row):
        for j in k[i: i+floats_per_row]:
            fp.write("%.7f, " % v[j])
        fp.write("\n%*s" % (indent+1, ""))
    fp.write("]")

def freq_counts(seqs, alphabet, weights=None):
    if weights == None: weights = [1.0]*len(seqs)
    counts = [[0.0 for a in alphabet] for s in seqs[0]]
    order = {}
    for i, a in enumerate(alphabet):
        order[a] = i
    for i in range(len(seqs[0])):
        for j in range(len(seqs)):
            a = order[seqs[j][i]]
            counts[i][a] += weights[j]
    return counts

def henikoff_weights(seqs, alphabet, counts):
    order = {}
    for i, a in enumerate(alphabet):
        order[a] = i
    weights = [0.0]*len(seqs)
    w_sum = 0.0
    for i in range(len(counts)):
        r = 0
        for j in range(len(counts[i])):
            if counts[i][j] > 0.0: r += 1
        for k in range(len(seqs)):
            w = 1 / (r * counts[i][order[seqs[k][i]]])
            weights[k] += w
            w_sum += w
    c = len(seqs) / w_sum
    return [w*c for w in weights]


class DirichletMix:
    def __init__(self, dic=None):
        if dic == None:
            self.mixture = []
            self.alpha = []
            self.order = {}
        else:
            self.__dict__.update(dic)
            self.order = {}
            for i, a in enumerate(self.alphabet):
                self.order[a] = i
            self.num_distr = len(self.mixture)
           
    def read(self, filename):
        fp = file(filename, 'r')
        actions = {'Name': "self.name = atoms[2]",
                   'Name=': "self.name = atoms[1]",
                   'Order': "self.alphabet = join(atoms[2:],'')",
                   'Order=': "self.alphabet = join(atoms[1:],'')",
                   'Mixture=': "self.mixture.append(float(atoms[1]))",
                   'Alpha=': "self.alpha.append(map(float, atoms[1:]))",
                   }

        for line in fp:
            atoms = split(strip(line))
            if not len(atoms): continue
            if atoms[0] == 'EndClassName': break
            if atoms[0] not in actions: continue
            else:
                exec(actions[atoms[0]])
        fp.close()
        for i, a in enumerate(self.alphabet):
            self.order[a] = i
        self.num_distr = len(self.mixture)

    def print_as_PyDict(self, indent=0, floats_per_row=10, fp=sys.stdout,
                        indent_first_row=True):
        if indent_first_row: indent1 = indent
        else: indent1 = 0
        fp.write("%*s{'name': '%s',\n" % (indent1, "", self.name))
        fp.write("%*s 'alphabet': '%s',\n" % (indent, "", self.alphabet))
        fp.write("%*s 'mixture': " % (indent, ""))
        _print_vector(self.mixture, indent+12, floats_per_row, fp, False)
        fp.write(",\n")
        fp.write("%*s 'alpha': [" % (indent, ""))
        _print_vector(self.alpha[0], indent+11, floats_per_row, fp, False)
        fp.write(",\n")
        for i in range(1,len(self.alpha)):
            _print_vector(self.alpha[i], indent+11, floats_per_row, fp, True)
            fp.write(",\n")
        fp.write("%*s           ]}" % (indent, ""))

    def _coeffs(self, k, n, sum_n):
        alpha = array(self.alpha[k][1:])
        sum_alpha = self.alpha[k][0]
        q = self.mixture[k]
        P = 1.0
        for i in range(len(alpha)):
            if n[i]: P *= n[i] * beta(n[i], alpha[i])

        
        ## P = product(gamma(n+alpha) / (gamma(alpha) * gamma(n+1)))
        B = sum_n * beta(sum_n, sum_alpha)
##         print "sum_n=%d" % int(sum_n)
##         print "sum_alpha=%.4e" % sum_alpha
##         print "beta(sum_n, sum_alpha)=%.4e" % beta(sum_n, sum_alpha)
##         print "sum_n * beta(sum_n, sum_alpha)=%.4e" % (sum_n * beta(sum_n, sum_alpha))
##         print "q=%.4e" % q
##         print "P=%.4e" % P
##         print "gamma(sum_n+1)=%.4e" % gamma(sum_n+1)
##         print "gamma(sum_alpha)=%.4e" % gamma(sum_alpha)
##         print "gamma(sum_n + sum_alpha)=%.4e" % gamma(sum_n + sum_alpha)
##        return q * (gamma(sum_n+1) * gamma(sum_alpha) / gamma(sum_n + sum_alpha)) * P
        return q * B / P

    def _probs(self, i, coeff, n):
        a = array([self.alpha[k][i+1] for k in range(self.num_distr)])
        return sum(a * coeff) + n[i] * sum(coeff)

    def aa_vector(self, aa_dict):
        return [aa_dict[a] for a in self.alphabet]

    def block_counts(self, seqs, weights=None):
        return freq_counts(seqs, self.alphabet, weights)

        
    def pos_probs(self, counts):
        n = array(counts)
        sum_n = sum(n)

        # Calculate coefficients - reused for all amino acids,
        # depend only on k
##         for k in range(self.num_distr):
##             print "k=%d, %.4e" % (k, self._coeffs(k, n, sum_n)) 
        
        coeff = array([self._coeffs(k, n, sum_n) for k in range(self.num_distr)])
##         print "coeff="
##         print coeff
##         print "sum_coeff=", sum(coeff)
##         print "self.alpha[k][0] + sum_n=", (self.alpha[k][0] + sum_n)
        coeff = coeff / (sum(coeff) * (self.alpha[k][0] + sum_n))
        #print "*********"
        #print "n=", n
        #print "coeff=", coeff
        return [self._probs(i, coeff, n) for i in range(len(counts))]

    def pos_log_odds(self, pos_probs, bkgrnd, scale=1.0):
        n = len(pos_probs)
        return [scale * log(pos_probs[i]/bkgrnd[i]) / log(2.0) for i in range(n)]

    def block_probs(self, block_counts):
        return [self.pos_probs(counts) for counts in block_counts]

    def block_log_odds(self, block_probs, bkgrnd, scale=1.0):
        return [self.pos_log_odds(probs, bkgrnd, scale) for probs in block_probs]
    
    def block2pssm(self, block_data, seq):
        pssm_info  = []
        for i in range(len(block_data)):
            score_dict = {}
            for a in self.alphabet:
                score_dict[a] = block_data[i][self.order[a]]
            pssm_info.append((seq[i], score_dict))
        return PSSM(pssm_info)

    def print_block_data(self, block_data, field_width=4, precision=4, dtype='int'):
        file_str = StringIO()
        if dtype == 'int':
            pdata = lambda x: "%*d " % (field_width, x)
        else:
            pdata = lambda x: "%*.*f " % (field_width, precision, x)
            
        file_str.write("   ")
        for i in range(len(block_data)):
            file_str.write("%*d " % (field_width, i))
        file_str.write("\n")
        for j, a in enumerate(self.alphabet):
            file_str.write(" %c " % a)
            for i in range(len(block_data)):
                file_str.write(pdata(block_data[i][j]))
            file_str.write("\n")
        return file_str.getvalue()
                      
def Comp2PyScript(path, filename=None):
    """
    Converts all Dirichlet mixtures files (ending in comp)
    in the path to Python dictionaries and stores them
    in filename if provided (otherwise prints to stdout).
    """

    if filename == None:
        fp = sys.stdout
    else:
        fp = file(filename, 'w')

    names = os.listdir(path)
    ffunc = lambda s: s[-4:] == 'comp'
    names = filter(ffunc ,names)
    names.sort()

    fp.write("import DirichletMix\n\n")
    fp.write("get_mix = lambda name: DirichletMix.DirichletMix(DIR_MIX[name])\n") 
    fp.write("get_names = lambda : NAMES\n\n")

    k = range(len(names))
    fp.write("NAMES = [")
    for i in range(0, len(k), 3):
        for j in k[i: i+3]:
            fp.write("'%s', " % names[j])
        fp.write("\n%*s" % (9+1, ""))
    fp.write("]\n\n")
    
    fp.write("DIR_MIX = {\n")
    fp.write("    '%s': " % names[0])
    DM = DirichletMix()
    DM.read(names[0])
    DM.print_as_PyDict(4+len(names[0])+4, 4, fp, False)
    fp.write(",\n")
    
    for s in names[1:]:
        fp.write("%*s'%s': " % (4, "", s))
        DM = DirichletMix()
        DM.read(s)
        DM.print_as_PyDict(4+len(s)+4, 4, fp, False)                
        fp.write(",\n")
    fp.write("           }\n")
    
    if filename != None:
        fp.close()

