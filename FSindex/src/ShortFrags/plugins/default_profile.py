

from cStringIO import StringIO

from ShortFrags.Expt.DirichletMix import DirichletMix, freq_counts, \
     henikoff_weights, BKGRND_PROBS as bg_dict
from ShortFrags.Expt.DirichletInfo import get_mix, NAMES
from ShortFrags.Expt.matrix import POSITIONAL 




iteration = True
arg_list = [('Scale', 'real', 2.0),
            ('Weighting', ['None', 'Henikoff'], 'Henikoff'),
            ('Regulariser', NAMES, 'recode3.20comp'),
            ]
    
def _get_matrix_counts(HL, scale, weight_type, dirichlet_type): 

    seqs = HL.get_seqs()
        
    # Calculate sequence weights
    DM = get_mix(dirichlet_type)
    bcounts = DM.block_counts(seqs)

    if weight_type == 'None':
        weights = [1.0]*len(seqs)
        wcounts = bcounts
    elif weight_type == 'Henikoff':
        weights = henikoff_weights(seqs, DM.alphabet, bcounts)  
        wcounts = DM.block_counts(seqs, weights)
        
    wprobs = DM.block_probs(wcounts)
    bkgrnd = DM.aa_vector(bg_dict)

    PM = DM.block2pssm(DM.block_log_odds(wprobs, bkgrnd, scale),
                       HL.query_seq) 
    PM.name = 'PSSM'
    PM.module = __name__
    matrix_type = POSITIONAL
    ctype = 0

    return PM, matrix_type, ctype, bcounts, weights, wcounts, wprobs 


def get_matrix(HL, scale, weight_type, dirichlet_type):

    return _get_matrix_counts(HL, scale, weight_type,
                              dirichlet_type)[0:3]


def print_info(HL, scale, weight_type, dirichlet_type):

    if weight_type == None: return ""
    seqs = HL.get_seqs()
    deflines = HL.get_deflines()

    PM, matrix_type, ctype, bcounts, weights, wcounts, wprobs = \
        _get_matrix_counts(HL, scale, weight_type, dirichlet_type)

    DM = get_mix(dirichlet_type)
    file_str = StringIO()
    file_str.write('***** ALIGNMENT *****\n')
    for i in range(len(seqs)):
        file_str.write('%8.4f %s %s\n' % (weights[i], seqs[i],
                                          deflines[i])) 
    file_str.write('\n***** COUNTS *****\n')
    file_str.write(DM.print_block_data(bcounts))
    file_str.write('\n***** WEIGHTED COUNTS *****\n')
    file_str.write(DM.print_block_data(wcounts, 5, 1, 'float'))
    file_str.write('\n***** DIRICHLET MIXTURE PROBABILITIES *****\n')
    bprobs = DM.block_probs(wcounts)
    file_str.write(DM.print_block_data(bprobs, 6, 4, 'float'))
    # file_str.write('\n***** INFORMATION CONTENT *****\n')
    # file_str.write(DM.block2pssm(wcounts, seqs[0]).__str__())
    file_str.write("\n"+ PM.__str__())
    return file_str.getvalue()
