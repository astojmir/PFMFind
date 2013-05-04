Beyond GUI
==========

While the graphical interface is very useful for visualising search results and performing indivudual searches, it also has some limitations. For example, the selection of search matrices is limited and large iterated batch searches are not easy to perform. To perform these advanced tasks, it may be necessary to write your own Python code that interfaces with PFMFind as a library. PFMFind allows custom user plugins that provide additional search matrices, and it is relatively easy to write short scripts to automate some tasks.

Search matrix plugins
---------------------

.. py:module:: search_plugin

Writing plugins for PFMFind is quite easy as far as the interface is concerned. A PFMFind search matrix plugin is a Python module that defines two global variables, :py:data:`iteration` and :py:data:`arg_list`, as well as two functions: :py:func:`get_matrix` (mandatory) and :py:func:`print_info` (optional).

.. py:data:: iteration

A boolean constant (either ``True`` or ``False``). If it is set to ``True``, the plugin can be used in the second or subsequent iterations, otherwise it is only used for the first iteration.

.. py:data:: arg_list

A list specifying the arguments of the functions :py:func:`get_matrix` and :py:func:`print_info`. Its elements are triplets of the form ``(name, type, default_value)``, where ``name`` is a string identifying the variable (it is displayed on the GUI), the ``type`` is either a string or a list of strings. In the former case the GUI sets up a ``Pmw.EntryField`` widget whose value type is given by the string (Please refer to Pmw documentation). In the latter case, a ``Pmw.OptionMenu`` is setup with options being the members of lists. In both cases the given ``default_value`` is preselected.

.. py:function:: get_matrix(HL, *args)

   Construct the scoring matrix for similarity search.

   :param HL: a hit list from the previous iteration
   :type HL: :py:class:`pfmfind.search.hit_list.HitList` instance
   :param args: additional positional arguments in the order they appear in :py:data:`arg_list`
   :return:  A tuple of the form ``(M, matrix_type, ctype)``, where ``M`` is a Biopython-style score matrix or PSSM, ``matrix_type`` is 0 if the matrix is a score matrix and 1 if it is a PSSM, while ``ctype`` should be set to 0 if the matrix contains similarity scores (the other values are for distance based matrices used by FSIndex).

.. py:function:: print_info(HL, *args)

   Construct a printable representation of the matrix and the method used to obtain it. Can be omitted, in which case the default printout is produced. It takes the same arguments as :py:func:`get_matrix`.

   :param HL: a hit list from the previous iteration
   :type HL: :py:class:`pfmfind.search.hit_list.HitList` instance
   :param args: additional positional arguments in the order they appear in :py:data:`arg_list`
   :return: A string showing (in a human-readable way) the matrix produced by :py:func:`get_matrix`.


Examples
^^^^^^^^

The listing below shows the code of the default first iteration search plugin as an example. It extracts available amino acid scoring matrices from Biopython and removes non-standard letters from them before returning them::

    from Bio.SubsMat import MatrixInfo
    from pfmfind.search.matrix import QUASI, MAX, AVG, SCORE
    from pfmfind.search.matrix import SubstitutionMatrix

    _MATRIX_CTYPE = {'None': 0, 'Quasi': QUASI, 'Avg': AVG, 'Max': MAX}

    iteration = False
    arg_list = [('Matrix Name', MatrixInfo.available_matrices, 'blosum62'),
                ('Conversion', _MATRIX_CTYPE.keys(), 'None'),
                ]

    _std_alphabet_map = {}.fromkeys(list("ACDEFGHIKLMNPQRSTVWY"))


    def _filter_non_standard_letters(S):
        for a, b in S.keys():
            if a not in _std_alphabet_map or b not in _std_alphabet_map:
                del(S[(a,b)])


    def get_matrix(HL, matrix_name, conv_type):

        S = SubstitutionMatrix()
        S.update(getattr(MatrixInfo, matrix_name))
        S.name = matrix_name

        _filter_non_standard_letters(S)
        matrix_type = SCORE
        ctype = _MATRIX_CTYPE[conv_type]
        return S, matrix_type, ctype

The default profile plugin is more complicated::

    from cStringIO import StringIO

    from pfmfind.search.DirichletMix import DirichletMix
    from pfmfind.search.DirichletMix import freq_counts
    from pfmfind.search.DirichletMix import henikoff_weights
    from pfmfind.search.DirichletMix import BKGRND_PROBS as bg_dict
    from pfmfind.search.DirichletInfo import get_mix
    from pfmfind.search.DirichletInfo import NAMES
    from pfmfind.search.matrix import POSITIONAL


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

        if not len(HL):
            return None, 0, 0
        return _get_matrix_counts(HL, scale, weight_type,
                                  dirichlet_type)[0:3]


    def print_info(HL, scale, weight_type, dirichlet_type):

        if not len(HL):
            return "Too few hits to construct PSSM"

        if weight_type is None:
            return ""

        seqs = HL.get_seqs()
        deflines = HL.get_deflines()

        PM, matrix_type, ctype, bcounts, weights, wcounts, wprobs = \
            _get_matrix_counts(HL, scale, weight_type, dirichlet_type)

        DM = get_mix(dirichlet_type)
        file_str = StringIO()
        file_str.write('***** ALIGNMENT *****\n')
        for i in range(len(seqs)):
            file_str.write('%8.4f %s %s\n' % (weights[i], seqs[i], deflines[i]))
        file_str.write('\n***** COUNTS *****\n')
        file_str.write(DM.print_block_data(bcounts))
        file_str.write('\n***** WEIGHTED COUNTS *****\n')
        file_str.write(DM.print_block_data(wcounts, 5, 1, 'float'))
        file_str.write('\n***** DIRICHLET MIXTURE PROBABILITIES *****\n')
        bprobs = DM.block_probs(wcounts)
        file_str.write(DM.print_block_data(bprobs, 6, 4, 'float'))
        file_str.write("\n"+ str(PM))
        return file_str.getvalue()
