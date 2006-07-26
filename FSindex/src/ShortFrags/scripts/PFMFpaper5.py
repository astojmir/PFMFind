

import sys, psycopg2
from sets import Set

def true_results(filename):
    """
    Iterates over all the domains associated with the yeast sequence from a
    tab-separated file.

    Format: sequence name, domain id, start, end
    """

    fp = open(filename, 'r')
    coverage = Set([])

    # Add first line
    line = fp.next()
    ls = line.split('\t')
    yeast_seq = ls[0]
    
    coverage.union_update(Set(range(int(ls[2])-1,int(ls[3]))))

    # All other lines
    for line in fp:
        ls = line.split('\t')
        if ls[0] != yeast_seq:
            yield (yeast_seq, coverage)
            coverage = Set([])
            yeast_seq = ls[0]

        coverage.union_update(Set(range(int(ls[2])-1,int(ls[3]))))
            
    # Yield last yeast sequence
    fp.close()
    yield (yeast_seq, coverage)


def PFMF_results(filename, length, cutoff):
    """
    Iterates over all the PFMF hits from a tab-separated file.
    """

    fp = open(filename, 'r')


    # Add first line
    line = fp.next()
    ls = line.split('\t')

    true_hits = Set([])
    false_hits = Set([])
    yeast_seq = ls[0]

    score = float(ls[6])
    frag = int(ls[2])
    
    if score >= cutoff:
        true_hits.union_update(Set(range(frag,frag+length)))
    else:
        false_hits.union_update(Set(range(frag,frag+length)))

    # All other lines
    for line in fp:
        ls = line.split('\t')
        if ls[0] != yeast_seq:
            true_hits.difference_update(false_hits)
            yield (yeast_seq, true_hits, false_hits)
            true_hits = Set([])
            false_hits = Set([])
            yeast_seq = ls[0]

        score = float(ls[6])
        frag = int(ls[2])
        if score >= cutoff:
            true_hits.union_update(Set(range(frag,frag+length)))
        else:
            false_hits.union_update(Set(range(frag,frag+length)))
            
    # Yield last yeast sequence
    fp.close()
    true_hits.difference_update(false_hits)
    yield (yeast_seq, true_hits, false_hits)



def get_experiments(dbcur):
    sql_base = "SELECT name, experiment_id, length(query_sequence) FROM experiments WHERE %s"

    sql_verified = "query_description ~ 'Verified'"
    sql_dubious = "query_description ~ 'Dubious'"
    sql_unknown = "query_description ~ 'Uncharacterized'"
    sql_others = "query_description !~ 'Verified' AND query_description !~ 'Dubious' AND query_description !~ 'Uncharacterized'"

    sqls = [sql_verified, sql_dubious, sql_unknown, sql_others]

    exp_dict = {}
    for i,s in enumerate(sqls):
        dbcur.execute(sql_base % s)
        exp_list = dbcur.fetchall()
        tmp = {}
        for name, eid, length in exp_list:
            tmp[name.split(' ')[1]] = (eid, length, i)
        exp_dict.update(tmp)
        
    return exp_dict

DBHOST = 'localhost'
DBDATABASE = 'PFMFind'
DBUSER = 'aleksand'

if __name__=='__main__':
    interpro_filename = sys.argv[1]
    PFMF_filename = sys.argv[2]
    length = int(sys.argv[3])
    cutoff = float(sys.argv[4])
#    output_filename = sys.argv[3]
#    ofp = open(output_filename, 'w')

    # Connect to the database and get the mapping of yeast sequence
    # names to experiment ids
    dbcon = psycopg2.connect('host=%s dbname=%s user=%s' %(DBHOST, DBDATABASE, DBUSER))
    dbcur = dbcon.cursor()
    exp_dict = get_experiments(dbcur)


    IP_iter = true_results(interpro_filename)
    c = 0
    c1 = 0
    T = [0]*4
    F = [0]*4
    O = [0]*4
    H = [0]*4
    CVR = [0]*4
    
    for yeast_seq, true_hits, false_hits in PFMF_results(PFMF_filename, length, cutoff):
        c1 += 1
        while 1:
            yeast_seq1, coverage = IP_iter.next()
            eid, l, cat = exp_dict[yeast_seq1]
            CVR[cat] += len(coverage)
#            print yeast_seq1, len(coverage)
            if yeast_seq1 == yeast_seq: break
        c += 1

        t = len(coverage & true_hits)
        f = len(coverage & false_hits)
        o = t + f
        h = len((true_hits | false_hits) - coverage)

        T[cat] += t
        F[cat] += f
        O[cat] += o
        H[cat] += h
        
    print CVR
    print "TOTAL COUNT: %d. TOTAL YEAST_SEQUENCES: %d\n\n" % (c, c1)
    

    print "%10s %10s %10s %10s %10s" % ('', 'VER', 'DUB', 'UNK', 'OTH')
    print "%10s %10d %10d %10d %10d" % tuple(['OVERLAP'] + O)
    print "%10s %10d %10d %10d %10d" % tuple(['TRUE'] + T)
    print "%10s %10d %10d %10d %10d" % tuple(['FALSE'] + F)
    print "%10s %10d %10d %10d %10d" % tuple(['HANGING'] + H)
    
