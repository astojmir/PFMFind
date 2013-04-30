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


import math
import cPickle
from cStringIO import StringIO
import numpy as np
import numpy.fft as fft
from pfmfind.Expt.matrix import POSITIONAL
from pfmfind.Expt.matrix import SCORE
from pfmfind.Expt.matrix import SIMILARITY
from pfmfind.Expt.matrix import DISTANCE
from pfmfind.Expt.matrix import ScoreMatrix
from pfmfind.Expt.matrix import ProfileMatrix
from pfmfind.Expt.DirichletMix import BKGRND_PROBS as BG

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

        ft = np.ones(N1//2+1, np.complex128)
        if M.Mtype == SCORE:
            I = seq
        else:
            I = xrange(len(seq))


        self.min_score = 0
        self.max_score = 0

        for c0 in I:
            dtmp = np.zeros(N1, np.float64)
            self.min_score += min(M[c0].values())
            self.max_score += max(M[c0].values())

            for c1 in BG.keys():
                dtmp[M[c0][c1]-self.a] += BG[c1]
            ft *= fft.rfft(dtmp)
        self.density = fft.irfft(ft, N1)
        self.F = np.zeros(N, np.float64)

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
