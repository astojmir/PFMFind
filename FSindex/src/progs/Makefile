SHELL=/bin/sh

FSLIB=../lib/libFSindex.a

all: $(FSLIB) \
     FSindex_search \
     FSindex_cluster \
     FSindex_bin_size \
     FSindex_self_sim \
     FSindex_create \
     FSindex_stats \
     SS_search \
     count_letters \
     Sscan_eps_net \
     print_matrices \
     QD_ddist \
     split_db \
     check_db \
     print_doublets \
     FSindex_delta \
     FSindex_rand_perf \
     FSindex_eps_net1 \
     fragdist
#     cluster_test \
#     FSindex_rand_perf \
#     FSindex_benchmark \
#     FSindex_eps_net1 \

check_db:

split_db:

cluster_test: 

FSindex_rand_perf:

FSindex_search:

FSindex_self_sim:  

QD_ddist:   

print_matrices:   

FSindex_create:  

FSindex_eps_net1:  

FSindex_delta:  

FSindex_stats:  

SS_search:  

Sscan_eps_net:  

count_letters:  

fragdist:  

clean:
	for x in *; do [ -x $$x -a -f $$x.c ] && echo $$x; done\
	 | xargs rm -f; \
	rm -f *.o; \
	rm -f core;

srcclean: clean
	rm -f *~;


