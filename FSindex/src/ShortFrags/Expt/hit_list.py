from cStringIO import StringIO

class Hit:
    Idata = None
    attr = {'defline': 'self.Idata.get_def(self.seq_id)',
            'seq': 'self.Idata.get_frag(self.seq_id, self.seq_from, self.seq_to)',
            'cluster': 'self.Idata.seq2cluster(self.seq_id)',
            'accession': 'self.Idata.get_accession(self.seq_id)',
            }
    def __init__(self, dict):
        self.__dict__ = dict

    def __getattr__(self, name):
        try:
            val = eval(self.attr[name])
        except KeyError:
            raise AttributeError
        return val


# Functions for printing hits (as dictionaries) and
# their lists

DEFLINE_MAX = 57

def summary(hits,
	    show_rank=0,
	    defline_width=DEFLINE_MAX,
            rank_offset=0):
    """
    Prints the summary of a list of hits
    """
	
    file_str = StringIO()

    if show_rank:
        line_func = lambda i,  w, defline, dist, sim:\
                        '%4d. %-*.*s %4d %4d' % (i, w, w, defline, dist, sim)   
        if len(hits) > 1:
            file_str.write('%4.4s  %-*.*s %4.4s %4.4s\n' %\
                           ('Rank', defline_width, defline_width,
                            'Description', 'Dist', 'Sim'))
    else:
        line_func = lambda i,  w, defline, dist, sim:\
                    '%-*.*s %4d %4d' % (w, w, defline, dist, sim)   
        if len(hits) > 1:
            file_str.write('%-*.*s %4.4s %4.4s\n' %\
                           (defline_width, defline_width,
                            'Description', 'Dist', 'Sim'))

    for i, ht in enumerate(hits):
        defline = ht.defline
        if len(defline) > defline_width:
            defline = defline[:defline_width-3] + '...'
            file_str.write(line_func(i+rank_offset, defline_width,
                                     defline, ht.dist, ht.sim))
            if len(hits) > 1: file_str.write('\n')
    return file_str.getvalue()

def description(hits, query_seq, qs=0):
    """
    Prints full details for a list of hits.
    """
    file_str = StringIO()

    for ht in hits:
        defline = ht.defline
        a = qs
        b = qs + len(query_seq)
        file_str.write('%s\n' % defline)
        file_str.write('dist = %d  sim = %d\n\n' % (ht.dist,
                                                    ht.sim))
        file_str.write('Query: %5d %s %5d\n' % (a, query_seq, b))
        file_str.write('Sbjct: %5d %s %5d\n' % (ht.seq_from,
                                                ht.seq,
                                                ht.seq_to,
                                                ))
        file_str.write('\n')
        
    return file_str.getvalue()

class HitList:
    def __init__(self, dict):
        self.__dict__ = dict
        tmphits = [Hit(ht) for ht in self.hits]
        self.hits = tmphits

    def header_str(self):
        """
        Prints the query details.
        """

        file_str = StringIO()
        file_str.write("***** Query Parameters *****\n")
        file_str.write("Query fragment: %s\n" % self.query_seq)
        file_str.write("Query decription: %s\n" % self.query_def)
        file_str.write("Score matrix: %s\n" % self.matrix_name)
        file_str.write("Matrix conversion type: %s\n" % self.conv_type)
        file_str.write("Distance range: %d\n" % self.dist_range)
        file_str.write("Similarity score range: %d\n" % self.sim_range)
        file_str.write("Nearest neighbours cutoff: %d\n" % self.kNN)
        file_str.write("Actual neighbours: %d\n\n" % len(self.hits))
        return file_str.getvalue()

    def summary_str(self):
        """
        Prints the full summary.
        """

        file_str = StringIO()
        file_str.write("***** Summary *****\n")
        file_str.write(summary(self.hits, 1))
        file_str.write('\n')
        return file_str.getvalue()

    def full_str(self, qs=0):
        """
        Prints full details for all hits.
        """

        file_str = StringIO()
        file_str.write("***** Full Details *****\n")
        file_str.write(description(self.hits, self.query_seq, qs))
        file_str.write('\n')
        return file_str.getvalue()

    def perf_str(self, Idata=None):
        """
        Index performance statistics as string which can be printed.
        """

        file_str = StringIO()
        if Idata == None or Idata.I == None:
            perc_func = lambda val, item: "\n"
        else:
            index_data = Idata.I.get_data()
            file_str.write(Idata.print_str())
            perc_func = lambda val, item: " (%.2f %%)\n" \
                        % (100.0 * val / index_data[item])

        file_str.write("***** Index Performance *****\n")
        file_str.write("Number of checked bins : %d" % self.bins_visited)
        file_str.write(perc_func(self.bins_visited, 'bins'))
        file_str.write("Number of accepted bins : %d" % self.bins_hit)
        file_str.write(perc_func(self.bins_hit, 'bins'))
        file_str.write("Accepted out of generated bins:  %.2f %%\n\n" \
                       % (100.0 * self.bins_hit / self.bins_visited))
        
        file_str.write("Number of checked fragments : %d" % self.frags_visited)
        file_str.write(perc_func(self.frags_visited, 'fragments'))
        file_str.write("Number of accepted fragments : %d" % self.frags_hit)
        file_str.write(perc_func(self.frags_hit, 'fragments'))
        file_str.write("Accepted out of generated fragments:  %.2f %%\n\n" \
                       % (100.0 * self.frags_hit / self.frags_visited))

        if self.unique_frags_visited:
             file_str.write("Number of checked distinct fragments : %d" %\
                            self.unique_frags_visited)
             file_str.write(perc_func(self.unique_frags_visited,
                                      'unique_fragments'))
             file_str.write("Number of accepted distinct fragments : %d" %\
                            self.unique_frags_hit)
             file_str.write(perc_func(self.unique_frags_hit, 'unique_fragments'))
             file_str.write("Accepted out of generated distinct fragments:  %.2f %%\n\n" \
                            % (100.0 * self.unique_frags_hit / self.unique_frags_visited))

        file_str.write("Search time: %.2f sec.\n" % self.search_time)
        return file_str.getvalue()

    
    def print_str(self, Idata=None, qs=0):
        file_str = StringIO()
        file_str.write(self.header_str())
        file_str.write(self.summary_str(defline_func))
        file_str.write(self.full_str(defline_func, qs))
        file_str.write(self.perf_str(Idata))
        return file_str.getvalue()
        
    def __str__(self):
        return self.print_str()
    
    def _sort_hits(self, incr, attribs):
        """
        General hits sort - pass attributes in order
        of priority.
        """

        tmp = [[eval('ht.%s' % k) for k in attribs] + [i] for i, ht in enumerate(self.hits)]
        tmp.sort()
        if not incr:
            tmp.reverse()
        n = len(attribs)
        htmp = [self.hits[a[n]] for a in tmp]
        self.hits = htmp
        
    def sort_by_similarity(self, incr=True):
        """
        Sort hits by similarity score.
        """

        attribs = ['sim',
                   'seq_id',
                   'seq_from',
                   'seq_to',
                   ]
        self._sort_hits(incr, attribs)
        
    def sort_by_distance(self, incr=True):
        """
        Sort hits by distance.
        """
        
        attribs = ['dist',
                   'seq_id',
                   'seq_from',
                   'seq_to',
                   ]
        self._sort_hits(incr, attribs)

    def sort_by_seq(self, incr=True):
        """
        Sort hits by sequence (alphabetical).
        """

        attribs = ['seq',
                   'seq_id',
                   'seq_from',
                   'seq_to',
                   ]
        self._sort_hits(incr, attribs)

    def sort_by_seqid(self, incr=True):
        """
        Sort hits by sequence id.
        """
        
        attribs = ['seq_id',
                   'seq_from',
                   'seq_to',
                   ]
        self._sort_hits(incr, attribs)

    def sort_by_orthologs(self, incr=True):
        """
        Sort hits by orthologous cluster.
        """

        attribs = ['cluster',
                   'seq_id',
                   'seq_from',
                   'seq_to',
                   ]
        self._sort_hits(incr, attribs)

    def get_seqs(self):
        seqs = [ht.seq for ht in self.hits]
        return seqs

    def get_deflines(self):
        deflines = [ht.defline for ht in self.hits]
        return deflines
