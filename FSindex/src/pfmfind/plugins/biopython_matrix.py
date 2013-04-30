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



from pfmfind.search.DirichletMix import freq_counts, \
     henikoff_weights, BKGRND_PROBS as bg_dict
from Bio import SubsMat
from Bio.SubsMat.FreqTable import FreqTable, FREQ
from Bio.Alphabet import IUPAC
from Bio.Align import Generic, AlignInfo
from pfmfind.search.matrix import SCORE 

iteration = True
arg_list = [('Scale', 'real', 2.0),
            ('Weighting', ['None', 'Henikoff'], 'Henikoff'),
            ]

def get_matrix(HL, scale, weight_type):

    seqs = HL.get_seqs()

    bcounts = freq_counts(seqs, "ACDEFGHIKLMNPQRSTVWY")
    if weight_type == 'None':
        weights = [1.0]*len(seqs)
    elif weight_type == 'Henikoff':
        weights = henikoff_weights(seqs, "ACDEFGHIKLMNPQRSTVWY",
                                   bcounts) 

    align = Generic.Alignment(IUPAC.protein)
    for i in range(len(seqs)):
        align.add_sequence("Seq #%d" % i, seqs[i], weight=weights[i])
    summary_align = AlignInfo.SummaryInfo(align)

    # Must get expected frequencies from our own (i.e. whole
    # database) background frequencies. Otherwise they would be
    # derived from the alignment which wouldn't be good if we
    # have a small sample
    ftab = FreqTable(bg_dict, FREQ)
    arm = SubsMat.SeqMat(summary_align.replacement_dictionary())
    lom = SubsMat.make_log_odds_matrix(arm, ftab, factor=scale,
                                       round_digit=0, keep_nd=0)

    PM = lom
    matrix_type = SCORE
    return PM, matrix_type, 0
