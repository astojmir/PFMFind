"""
Gets the coverage of each yeast sequence by our hits.
"""

import psycopg2 
DBHOST = 'localhost'
DBDATABASE = 'PFMFind'
DBUSER = 'aleksand'


def get_experiments(dbcur):
    sql_base = "SELECT name, experiment_id, length(query_sequence) FROM experiments WHERE %s"

    sql_verified = "query_description ~ 'Verified'"
    sql_dubious = "query_description ~ 'Dubious'"
    sql_unknown = "query_description ~ 'Uncharacterized'"
    sql_others = "query_description !~ 'Verified' AND query_description !~ 'Dubious' AND query_description !~ 'Uncharacterized'"

    sqls = [sql_verified, sql_dubious, sql_unknown, sql_others]

    exp_dicts = []
    for s in sqls:
        dbcur.execute(sql_base % s)
        exp_list = dbcur.fetchall()
        exp_dict = {}
        for name, eid, length in exp_list:
            exp_dict[name.split(' ')[1]] = (eid, length)
        exp_dicts.append(exp_dict)
        
    return exp_dicts


if __name__=='__main__':
    dbcon = psycopg2.connect('host=%s dbname=%s user=%s' %(DBHOST, DBDATABASE, DBUSER))
    dbcur = dbcon.cursor()

    exp_dicts = get_experiments(dbcur)

    cat = 0;
    for seq_category in exp_dicts:
        for name,v in seq_category.iteritems():
            eid = v[0]
            l = v[1]
            
            # Length 9 
            sql_command = "SELECT fragment, num_hits FROM searches WHERE"\
                          " experiment_id=%d AND length=9 AND"\
                          " iteration=2" % eid
            dbcur.execute(sql_command)
            search_list = dbcur.fetchall()
            hits9 = len(search_list)

            cvr = [0]*l
            for f in search_list:
                cvr[f[0]:f[0]+9] = [1]*9
            cvr9 = sum(cvr)
            pcvr9 = 100.00 * cvr9/l
            
            # Length 12
            sql_command = "SELECT fragment, num_hits FROM searches WHERE"\
                          " experiment_id=%d AND length=12 AND"\
                          " iteration=2" % eid
            dbcur.execute(sql_command)
            search_list = dbcur.fetchall()
            hits12 = len(search_list)

            cvr = [0]*l
            for f in search_list:
                cvr[f[0]:f[0]+12] = [1]*12
            cvr12 = sum(cvr)
            pcvr12 = 100.00 * cvr12/l


            print '\t'.join((name, str(cat), str(l), str(hits9), str(cvr9), str(int(pcvr9)), str(hits12), str(cvr12), str(int(pcvr12))))
        cat += 1
