from Numeric import zeros, ones, Float, Complex
from ShortFrags.Expt.DirichletMix import BKGRND_PROBS as BG
import math
from FFT import real_fft, inverse_real_fft
from cStringIO import StringIO
from ShortFrags.Expt.matrix import POSITIONAL, SCORE, SIMILARITY, DISTANCE 
from ShortFrags.Expt.matrix import ScoreMatrix, ProfileMatrix
import cPickle

# Use Biopython matrices

class ScoreDistr(object):
    def __init__(self, PM, matrix_type, conv_type, seq=None):

        if matrix_type == SCORE:
            M = ScoreMatrix(PM)
        else: # matrix_type == POSITIONAL:
            M = ProfileMatrix(PM.pssm)
            seq = M.qseq
        
        if conv_type:
            M.conv_type = conv_type
            M = M.matrix_conv()

        # Noise level
        self.noise = min(BG.values())**len(seq)

        a = 100000000.0
        b = -100000000.0

        if M.Mtype == SCORE:
            I = seq
        else:
            I = xrange(len(seq))

        for c0 in I:
            a = min(a ,min(M[c0].values()))
            b = max(b ,max(M[c0].values()))

        self.a = a
        self.b = b
        self.start = len(seq)*self.a



        N = 2*(len(seq)+1)*(b-a+1)+1
        N1 = int(math.ldexp(1,int(math.ceil(math.log(N,2)))))
        
        ft = ones(N1//2+1, Complex)
        if M.Mtype == SCORE:
            I = seq
        else:
            I = xrange(len(seq))

            
        self.min_score = 0
        self.max_score = 0

        for c0 in I:
            dtmp = zeros(N1, Float)
            self.min_score += min(M[c0].values())
            self.max_score += max(M[c0].values())
            
            for c1 in BG.keys():
                dtmp[M[c0][c1]-self.a] += BG[c1]
            ft *= real_fft(dtmp)
        self.density = inverse_real_fft(ft, N1)
        self.F = zeros(N, Float)

        for i in range(len(self.density)):
            if self.density[i] < self.noise:
                self.density[i] = 0.0 

        if M.Stype == SIMILARITY:
            self.decr = 1
            I = xrange(N-1,-1,-1)
        else:
            self.decr = 0
            I = xrange(N)

        FF = 0.0
        for i in I:
            FF += self.density[i]
            self.F[i] = FF


    def __str__(self):
        file_str = StringIO()
        file_str.write('data range [%d, %d]\n' % (self.min_score, self.max_score))
        file_str.write('density sum: %8.6f\n' % sum(self.density))
        for x in xrange(self.min_score-self.start, self.max_score-self.start+1):
            file_str.write('%4d %4d %8.6f %8.2e ' % (x,
                                                     x+self.start,
                                                     self.density[x],
                                                     self.F[x]))
            k = int(self.density[x]*100)
            for i in range(k):
                file_str.write('*')
            file_str.write('\n')
        return file_str.getvalue()

    def pvalue(self, score):
        if score < self.start:
            score = self.start
        elif score > len(self.F) + self.start - 1:
            score = len(self.F) + self.start - 1

        return self.F[score-self.start]

    def cutoff(self, pvalue):
        a = self.min_score-self.start
        b = self.max_score-self.start+1
        if self.decr:
            I = xrange(b-1,a-1,-1)
            j = b
        else:
            I = xrange(a,b)
            j = a

        for i in I:
            if self.F[i] > pvalue:
                break
            else:
                j = i

        return j+self.start


## from Bio.SubsMat import MatrixInfo

## M = MatrixInfo.__dict__['blosum62']
## Sdistr = ScoreDistr(M, SCORE, 0, 'YPFPGP')
## print Sdistr

## pickled = cPickle.dumps(Sdistr, 2)
## print "Pickled size", len(pickled)
## print "double array size", len(cPickle.dumps(zeros(2048, Float)))
    


