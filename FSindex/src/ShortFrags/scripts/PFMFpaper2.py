

import sys, psycopg2


def true_results(filename):
    """
    Iterates over all the domains associated with the yeast sequence from a
    tab-separated file.

    Format: sequence name, domain id, start, end
    """

    fp = open(filename, 'r')
    domains = []

    # Add first line
    line = fp.next()
    ls = line.split('\t')
    yeast_seq = ls[0]
    domains.append((int(ls[2]), int(ls[3]), ls[1]))

    # All other lines
    for line in fp:
        ls = line.split('\t')
        if ls[0] != yeast_seq:
            yield (yeast_seq, domains)
            domains = []
            yeast_seq = ls[0]

        domains.append((int(ls[2]), int(ls[3]), ls[1]))
            
    # Yield last yeast sequence
    fp.close()
    yield (yeast_seq, domains)


def get_matching_domains(domains, length, frag):
    matching_domains = {}

    fa = frag+1
    fb = frag+length
    for a, b, dom in domains:
        if (fa >= a and fa <= b) or\
           (fb >= a and fb <= b):
            matching_domains[dom] = None
    return matching_domains.keys()


def get_domain_histogram(eid, length, frag):

    # THIS IS NOT GENERAL ENOUGH (ontology_id = 9)
    sql_command = "SELECT f.name, count(f.name) FROM (SELECT bt.bioentry_id, bt.fragment, bt.start, e.name, e.seqfeature_id FROM (SELECT b.bioentry_id, ht.fragment, ht.start FROM bioentry AS b JOIN  (SELECT h.fragment, h.accession, h.start FROM hits AS h WHERE h.experiment_id=%d  AND h.length=%d AND h.iteration=2 AND fragment=%d) AS ht USING (accession)) AS bt JOIN (SELECT s.bioentry_id, s.seqfeature_id, t.name FROM seqfeature AS s JOIN (SELECT term_id, name FROM term WHERE ontology_id=9) AS t ON s.type_term_id=t.term_id) AS e USING (bioentry_id)) AS f JOIN location AS l USING (seqfeature_id) WHERE ((f.start+1) BETWEEN l.start_pos AND l.end_pos) OR ((f.start+%d) BETWEEN l.start_pos AND l.end_pos) GROUP by f.name" % (eid, length, frag, length)
    dbcur.execute(sql_command)
    hist_list = dbcur.fetchall()
    total_keywords = 0
    hist_dict = {}
    for name, count in hist_list:
        total_keywords += count
        hist_dict[name] = count

    return hist_dict, total_keywords

def get_experiments(dbcur):
    dbcur.execute("SELECT name, experiment_id FROM experiments")
    exp_list = dbcur.fetchall()
    exp_dict = {}
    for name, eid in exp_list:
        exp_dict[name.split(' ')[1]] = eid
    return exp_dict


DBHOST = 'localhost'
DBDATABASE = 'PFMFind'
DBUSER = 'aleksand'

if __name__=='__main__':
    input_filename = sys.argv[1]
    length = int(sys.argv[2])

    # Connect to the database and get the mapping of yeast sequence
    # names to experiment ids
    dbcon = psycopg2.connect('host=%s dbname=%s user=%s' %(DBHOST, DBDATABASE, DBUSER))
    dbcur = dbcon.cursor()
    exp_dict = get_experiments(dbcur)


    # Iterate through known sequences
    c = 0
    for yeast_seq, domains in true_results(input_filename):
        eid = exp_dict[yeast_seq]

        # Get all successfull searches
        sql_command = "SELECT fragment, num_hits FROM searches WHERE"\
                      " experiment_id=%d AND length=%d AND"\
                      " iteration=2" % (eid, length)
        dbcur.execute(sql_command)
        search_list = dbcur.fetchall()
        c += len(search_list)

    print "TOTAL TO BE SCANNED: %d" % c

        


