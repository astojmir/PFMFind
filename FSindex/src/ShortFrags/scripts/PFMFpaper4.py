

import sys, psycopg2


def coverage(filename, exp_dict):
    fp = open(filename, 'r')
    
    # Add first line
    line = fp.next()
    ls = line.split('\t')
    yeast_seq, a, b = ls[0], int(ls[2]), int(ls[3])
    eid, l, cat = exp_dict[yeast_seq]
    cvr = [0]*l
    cvr[a-1:b] = [1]*(b-a+1)

    # All other lines
    for line in fp:
        ls = line.split('\t')
        if ls[0] != yeast_seq:
            yield (yeast_seq, cat, l, sum(cvr), 100.0*sum(cvr)/l)
            yeast_seq = ls[0]
            eid, l, cat = exp_dict[yeast_seq]
            cvr = [0]*l
        a, b = int(ls[2]), int(ls[3])
        cvr[a-1:b] = [1]*(b-a+1)

    # Yield last yeast sequence
    fp.close()
    yield (yeast_seq, cat, l, sum(cvr), 100.0*sum(cvr)/l)

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
    input_filename = sys.argv[1]

    # Connect to the database and get the mapping of yeast sequence
    # names to experiment ids
    dbcon = psycopg2.connect('host=%s dbname=%s user=%s' %(DBHOST, DBDATABASE, DBUSER))
    dbcur = dbcon.cursor()
    exp_dict = get_experiments(dbcur)

    tot_cvr = [0]*4
    for res in coverage(input_filename, exp_dict):
        tot_cvr[res[1]] += res[3]
        #print "%s\t%d\t%d\t%d\t%5.1f" % res
        
    print "\n\n", tot_cvr, sum(tot_cvr)
