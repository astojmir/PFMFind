"""
Provides classes to define a hit and construct a 
string to print its information.

Hits are fragments from the searched database.

Idata is assigned dynamically in B{index.py}.

B{Exceptions}:
    - KeyError

B{Classes}:
    - Hit  
        - Defines a dictionary for each hit [a fragment of the sequence].
    - HitList  
        - Provides the list of hits and the functions to manipulate them.

B{Functions}:
    - description(hits, query_seq, qs=0):
        - Return a printable string providing full details for a list 
	of hits.
    - summary(hits, show_rank=0, defline_width=57, rank_offset=0)
        - A list of hit class instances to be displayed.
    
"""


from cStringIO import StringIO

class Hit(object):

    """
    Defines a dictionary for each hit [a fragment of the sequence].
    
    For each hit a dictionary is defined.  Each class 
    contains a list of dictionaries for individual hits. 
    The dictionary for each hit provides the following 
    attributes:
    
        - B{seq_id} 
	     - Sequence description.
        - B{seq_to} 
	     - End of sequence. 
        - B{seq_from}
	     - Starting point for sequence.
        - B{dist}
	     - Distance.
        - B{sim} 
             - Similarity score.   

    Therefore, an instance of Hit will contain all of the above 
    attributes plus the following ones which are assigned 
    dynamically and for which an instance of IndexedDb must 
    exist and have a sequence database loaded.
    
        - B{defline}
            - Description of the hit.
        - B{accession}
	    - Key to uniquely identify protein.
        - B{seq}  
	    - The sequence of the hit fragment.
 
    B{Exceptions:}
        - KeyError
    
    """     


    
    Idata = None
    attr = {'defline': 'self.Idata.get_def(self.seq_id)',
            'seq': 'self.Idata.get_frag(self.seq_id, self.seq_from, self.seq_to)',
            'accession': 'self.Idata.get_accession(self.seq_id)',
            }
	    
    def __init__(self, dict):
        """
	Constructor, assigns all non dynamic attributes.
	See class description for details on which attributes
	are dynamic.
	"""
        self.__dict__ = dict

    def __getattr__(self, name):
  
        try:
            val = eval(self.attr[name])
        except KeyError:
            raise AttributeError
        return val


# Functions for printing hits (as dictionaries) and
# their lists

def summary(hits,
	    show_rank=0,
	    defline_width=57,
            rank_offset=0):
    """
    Construct a printable string from a list of 
    hit class instances.
    
        @param hits: The hit to be displayed.
        @param show_rank: Display the rank of the hits if true.
	Default value 0.
	@param defline_width: Width of the description.
	@param rank_offset: Control how much space is between
	the rank and the hit.  Default value 0.
	
        @return: String providing summary.
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

    for j, ht in enumerate(hits):
        defline = ht.defline
        if len(defline) > defline_width:
            defline = defline[:defline_width-3] + '...'
        file_str.write(line_func(j+rank_offset, defline_width,
                                 defline, ht.dist, ht.sim))
        if len(hits) > 1: file_str.write('\n')
    return file_str.getvalue()

def description(hits, query_seq, qs=0):
    """
    Return a printable string providing full details for a 
    list of hits.
    
    @param hits: A list of hits to be displayed.
    @param query_seq: A query sequence.
    @param qs:  Starting offset of query fragment in B{query_seq}. The 
    length of the query fragment is the same as the length of the hit.
    
    @return: String providing details of hits.
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

class HitList(list):

    """
    Provides the list of hits.
    
    The class may be used to print query details, print a full summary, 
    print full details for all hits, and to print performance statistics.  
    The class also provides the functions necessary to sort and print 
    the results by priority, similarity score, distance, 
    sequence(alphabetical and sequence id.
    
    This class uses Idata.  B{Idata} is assigned dynamically in B{index.py}
    and provides the index information.  This includes the following 
    attributes:
    
        - B{I}
            - The index.
	    - Includes the index name, alphabet, and partitions.
        - B{get_frag}
            - The sequence fragment.
        - B{deflines}
            - Description of the hit.
        - B{accessions}
            - Accession numbers.
        - B{sdb}
            - Pointer to FASTA database.
    """
    
    def __init__(self, dict):
        """
	Constructor, assigns the dictionary providing the search 
	related information.  For each hit in the hits dictionary 
	B{Hit} is called to create a list of B{Hit} instances. The 
	information provided in the hits dictionary includes the 
	following attributes:
	    - B{query_seq}
	        - The query fragment.
	    - B{query_def}
	        - The query description.
	    - B{conv_type}
	        - Matrix conversion type.
	    - B{sim_range}
	        - Similarity score range.
	    - B{dist_range}
	        - Distance range.
	    - B{kNN}
	        - Nearest neighbor cutoff.
	    - B{bins_visited}
	        - Bins checked.
	    - B{bins_hit}
	        - Accepted bins.
	    - B{frags_visited}
	        - Number of checked fragments.
	    - B{frags_hit} 
	        - Number of accepted fragments.
	    - B{search_time}
	        - Time to complete the search.
	    - B{hits}
	        - A list holding individual dictionaries for each hit.  
		For the attributes provided in the dictionary of 
		each hit please refer to the help for class B{Hit}.
	"""
        HL = [Hit(ht) for ht in dict['hits']] 
        del dict['hits']
        list.__init__(self, HL)
        self.__dict__ = dict
	
    def __str__(self):
        """
        Call B{print_str} and return the string generated.
        """
        return self.print_str()

    def header_str(self):
        """
        Return a printable string with header information.
	
	@return: String with header information.
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
        file_str.write("Actual neighbours: %d\n\n" % len(self))
        return file_str.getvalue()

    def summary_str(self):
        """
        Return a printable string with full summary.
	
	@return: String with summary information.
        """

        file_str = StringIO()
        file_str.write("***** Summary *****\n")
        file_str.write(summary(self, 1))
        file_str.write('\n')
        return file_str.getvalue()

    def full_str(self, qs=0):
        """
        Return a printable string with full details for all hits.
	
	@param qs:  Starting offset of query fragment in B{query_seq}. The 
        length of the query fragment is the same as the length of the hit.
	@return: String with full details for all hits.
        """

        file_str = StringIO()
        file_str.write("***** Full Details *****\n")
        file_str.write(description(self, self.query_seq, qs))
        file_str.write('\n')
        return file_str.getvalue()

    def perf_str(self, Idata=None):
        """
        Return a printable string with performance statistics.
	
	@param Idata: Index data.  Please see the HitList 
	class help for more information on what attributes 
	are provided by Idata.
	@return: String providing performance statistics.
        """

        file_str = StringIO()
        if Idata == None or Idata.I == None:
            perc_func = lambda val, item: "\n"
        else:
            file_str.write(Idata.print_str())
            perc_func = lambda val, item: " (%.2f %%)\n" \
                        % (100.0 * val / item)

        file_str.write("***** Index Performance *****\n")
        file_str.write("Number of checked bins : %d" % self.bins_visited)
        file_str.write(perc_func(self.bins_visited, Idata.I.bins))
        file_str.write("Number of accepted bins : %d" % self.bins_hit)
        file_str.write(perc_func(self.bins_hit, Idata.I.bins))
        file_str.write("Accepted out of generated bins:  %.2f %%\n\n" \
                       % (100.0 * self.bins_hit / self.bins_visited))
        
        file_str.write("Number of checked fragments : %d" % self.frags_visited)
        file_str.write(perc_func(self.frags_visited, Idata.I.fragments))
        file_str.write("Number of accepted fragments : %d" % self.frags_hit)
        file_str.write(perc_func(self.frags_hit, Idata.I.fragments))
        file_str.write("Accepted out of generated fragments:  %.2f %%\n\n" \
                       % (100.0 * self.frags_hit / self.frags_visited))

#         if self.unique_frags_visited:
#              file_str.write("Number of checked distinct fragments : %d" %\
#                             self.unique_frags_visited)
#              file_str.write(perc_func(self.unique_frags_visited,
#                                       Idata.I.unique_fragments))
#              file_str.write("Number of accepted distinct fragments : %d" %\
#                             self.unique_frags_hit)
#              file_str.write(perc_func(self.unique_frags_hit, Idata.I.unique_fragments))
#              file_str.write("Accepted out of generated distinct fragments:  %.2f %%\n\n" \
#                             % (100.0 * self.unique_frags_hit / self.unique_frags_visited))

        file_str.write("Search time: %.2f sec.\n" % self.search_time)
        return file_str.getvalue()
    
    def print_str(self, Idata=None, qs=0):
        """
	Return a printable string with the header, summary, full details, 
	and performance statistics.
	
	@param Idata: Index data.  Please see the HitList 
	class help for more information on what attributes 
	are provided by Idata.
	@param qs: Starting offset of query fragment in B{query_seq}. The 
        length of the query fragment is the same as the length of the hit.
	
	@return:  String providing header, summary, full details and 
	performance statistics.
	"""
        file_str = StringIO()
        file_str.write(self.header_str())
        file_str.write(self.summary_str(defline_func))
        file_str.write(self.full_str(defline_func, qs))
        file_str.write(self.perf_str(Idata))
        return file_str.getvalue()
        
    
    def _sort_hits(self, incr, attribs):
        """
        General hits sort - pass attributes in order
        of priority.
        """

        tmp = [[eval('ht.%s' % k) for k in attribs] + [i] for i, ht in enumerate(self)]
        tmp.sort()
        if not incr:
            tmp.reverse()
        n = len(attribs)
        htmp = [self[a[n]] for a in tmp]
        self[:] = htmp
        
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

    def get_seqs(self):
        """
	Extract sequences of all hits from the dictionary.
	
	@return: Sequences from hits.
	"""
        seqs = [ht.seq for ht in self]
        return seqs

    def get_deflines(self):
        """
	Extract descriptions of all hits from the dictionary.
	
	@return: Descriptions from hits.
	"""
        deflines = [ht.defline for ht in self]
        return deflines
