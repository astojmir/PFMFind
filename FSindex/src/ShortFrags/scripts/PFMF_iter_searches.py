#! /usr/bin/env python

"""
Runs iterated Evalue-based searches for all fragments of given length
in a given set of sequences.

SYNOPSIS:

    PFMF_iter_searches.py host port tag fasta_file frag_len first_seq last_seq

"""
## ARGUMENTS:

## client_config_file: XML file with database and index configuration
##                     information (same format as used by PFMFind.py
##                     GUI). 
## fasta_file:         File (FASTA format) with query sequences. 
## frag_len:           Fragment length.
## skip:               Number of sequences to skip.
## """

__version__ = "$Rev: 1.0 $"
__date__ = "$Date: 2006-03-02 15:48:06 -0500 (Thu, 02 Mar 2006) $"
__author__ = "Aleksandar Stojmirovic"


from ShortFrags.Expt.SearchClient import SearchClient
from ShortFrags.Expt.fragexpt import PFMFindClient
from ShortFrags.Expt.db import db 
from ShortFrags.Expt.SearchServer import REL_SRCH
import sys, time, os

_RUNS_PER_FORK = 500
_MAX_SUBMITTED_SEARCHES = 20

def now():
    return time.ctime(time.time())


def write_experiment(efp, experiment_id,
                     name,
                     description,
                     query_sequence,
                     query_description,
                     min_len,
                     max_len):
    efp.write('%d\t' % experiment_id)
    efp.write('%s\t' % name)
    efp.write('%s\t' % description)
    efp.write('%s\t' % query_sequence)
    efp.write('%s\t' % query_description)
    efp.write('%d\t' % min_len)
    efp.write('%d\n' % max_len)


def write_search(sfp, hfp, experiment_id,
                 length, iteration, fragment,
                 HL):

    sfp.write('%d\t' % experiment_id)
    sfp.write('%d\t' % length)
    sfp.write('%d\t' % iteration)
    sfp.write('%d\t' % fragment)
    
    matrix_name = getattr(HL, 'matrix_name')
    matrix = '\N' # We don't write matrix at this stage

    sfp.write('%s\t' % HL.query_seq)
    sfp.write('%s\t' % HL.matrix_name)
    sfp.write(matrix)
    sfp.write('\t')
    sfp.write('%d\t' % HL.conv_type)
    sfp.write('%d\t' % HL.sim_range)
    sfp.write('%d\t' % HL.dist_range)
    sfp.write('%d\t' % HL.kNN)
    sfp.write('%d\n' % len(HL))

    # Hits next              

    for ht in HL:
        hfp.write('%d\t' % experiment_id)
        hfp.write('%d\t' % length)
        hfp.write('%d\t' % iteration)
        hfp.write('%d\t' % fragment)
        hfp.write('%s\t' % ht.accession)
        hfp.write('%d\t' % ht.seq_from)
        hfp.write('%d\t' % ht.dist)
        hfp.write('%s\t' % ht.sim)
        hfp.write('%e\t' % ht.pvalue)
        hfp.write('%e\n' % ht.Evalue)

def is_valid_fragment(qseq):
    """
    Checks that a given sequence is non-empty and consists only of
    letters from the standard amino acid alphabet.
    """
        
    valid_alphabet = 'ACDEFGHIKLMNPQRSTVWY'
    if not qseq or qseq.upper().translate('#'*256, valid_alphabet):
        return False
    return True


def iterated_rel_search(search_client,
                        efp, sfp, hfp,
                        experiment_id,
                        experiment_name,
                        experiment_description,
                        sequence, description,
                        frag_len,
                        start_func, start_params,
                        iter_func, iter_params,
                        cutoff_Evalues):
    """
    Runs automated iterated search for all fragments of a
    certain length of a given sequence.

    """

    tot_frags = len(sequence) + 1 - frag_len

    write_experiment(efp, experiment_id,
                     experiment_name,
                     experiment_description,
                     sequence,
                     description,
                     6, 18)

    # ** Prepare the first iteration **
    # Get the matrix first (it is the same for all searches)
    plugin_func = start_func
    PM, matrix_type, conv_type = plugin_func(None, *start_params)

    fragments = []
    srch_args = []
    Eval = cutoff_Evalues[0]
    for f in xrange(tot_frags):
        qseq = sequence[f: f + frag_len]
        if not is_valid_fragment(qseq): continue
        fragments.append(f)
        srch_args.append([REL_SRCH, qseq, PM, matrix_type,
                          Eval, conv_type])

    plugin_func = iter_func
        
    # *** Run all iterations ***
    final_results = []
    step = _MAX_SUBMITTED_SEARCHES
    for i in xrange(len(cutoff_Evalues)):
        profiles = []
        new_fragments = []
        for j in xrange(0, len(srch_args), step): 
            # Run batch search
            results = search_client.search(srch_args[j:j+step])

            # Immediately process the results to obtain profiles,
            for k,HL in enumerate(results):
                HL.matrix_name = getattr(srch_args[j+k][2], 'name', None)
                write_search(sfp, hfp, experiment_id, frag_len, i,
                             fragments[j+k], HL)  

            if i < (len(cutoff_Evalues) - 1):
                for k,HL in enumerate(results):
                    if HL == None: continue
                    PM_data = plugin_func(HL, *iter_params)
                    if PM_data[0] != None:
                        profiles.append(PM_data)
                        new_fragments.append(fragments[j+k])

                        
            # I don't know whether this would help !?
            for HL in results:
                del(HL)
            del(results)

        # Create arguments for next iteration
        if i < (len(cutoff_Evalues) - 1):
            fragments = new_fragments
            srch_args = [[REL_SRCH,
                          sequence[f: f + frag_len],
                          profiles[k][0],
                          profiles[k][1],
                          cutoff_Evalues[i+1],
                          profiles[k][2]] \
                         for k,f in enumerate(fragments) ] 

    return len(fragments)




def one_search_fork(server_host, server_port, tag,
                    fasta_file, frag_len, first_seq,
                    last_seq):

    finished_all = False

    # This is used only in order to get the plugin data
    PFMF_client = PFMFindClient()
    start_func = PFMF_client.start_plugins[start_plugin][0]
    iter_func = PFMF_client.iteration_plugins[iter_plugin][0]

    # Load fasta file
    fasta_db = db(fasta_file)
    n = fasta_db.no_seq
    if first_seq >= n: return True
    if last_seq >= n:
        last_seq = n-1
        finished_all = True
        
    # Connect to search server
    search_client = SearchClient()
    search_client.attach(server_host, server_port)

    # Open output files
    efp = file('%s_experiments.tbl' % tag,'a+')
    sfp = file('%s_searches.tbl' % tag, 'a+')
    hfp = file('%s_hits.tbl' % tag, 'a+')

    # Iterate over all sequences

    m = last_seq - first_seq + 1
    print "*** PFMF_iter_searches.py: starting %d runs at %s. ***" % (m, now())
    sys.stdout.flush()

    for i in xrange(first_seq, last_seq+1):
        sequence = fasta_db.get_seq(i)
        description = fasta_db.get_def(i)
        
        exp_name = 'Expt%4.4d %s' % (i, description.split(None,1)[0]) 


        final_motifs = iterated_rel_search(search_client,
                                           efp, sfp, hfp,
                                           i+1,
                                           exp_name,
                                           experiment_description,
                                           sequence, description,
                                           frag_len,
                                           start_func, start_params,
                                           iter_func, iter_params,
                                           cutoff_Evalues)


        print "Experiment %s finished at %s (%d profiles)" %\
              (exp_name, now(), final_motifs)
        sys.stdout.flush()

    efp.close()
    sfp.close()
    hfp.close()

    return finished_all


def copy_results(tag):

    import tarfile
    tar_filename = "%s_tables.tar.bz2" % tag
    tar = tarfile.open(tar_filename, "w:bz2")
    for name in ['%s_experiments.tbl' % tag,
                 '%s_searches.tbl' % tag,
                 '%s_hits.tbl' % tag]:
        tar.add(name)
    tar.close()

    os.system('scp %s 137.122.49.191:data/' % tar_filename)
    os.system('rm %s_*.tbl' % tag)

if __name__=='__main__':

    experiment_description = 'Iterated E-value searches'

    start_plugin = 'default_matrix'
    start_params = ('blosum62', 'None')

    iter_plugin = 'truncated_profile'
    iter_params = (2.0, 'Henikoff', 'recode3.20comp', 30)

    cutoff_Evalues = [1.0, 0.1, 0.01]

    skip = 0
    if len(sys.argv) < 8:
        print __doc__
        sys.exit(1)

    server_host = sys.argv[1]
    server_port = int(sys.argv[2])
    tag = sys.argv[3]
    fasta_file = sys.argv[4]
    frag_len = int(sys.argv[5])
    first_seq = int(sys.argv[6])
    last_seq = int(sys.argv[7])

    j = 0
    for first_seq1 in xrange(first_seq, last_seq+1, _RUNS_PER_FORK):
        last_seq1 = first_seq1 + _RUNS_PER_FORK - 1
        if last_seq1 > last_seq:
            last_seq1 = last_seq
        tag1 = '%s%4.4d' % (tag, j)

        pid = os.fork()
        if pid == 0:
            one_search_fork(server_host, server_port, tag1,
                            fasta_file, frag_len, first_seq1,
                            last_seq1)
            print "** Copying results to king **"
            sys.stdout.flush()
            copy_results(tag1)            
            os._exit(0)
        os.wait()
        j += 1
