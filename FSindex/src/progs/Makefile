SHELL=/bin/sh

FSLIB=../lib/libFSindex.a

all: $(FSLIB) \
     FSindex_create \
     FSindex_stats \
     FSindex_search \
     FSindex_bin_size \
     SS_search \
     count_letters \
     Sscan_eps_net \
     print_matrices \
     QD_ddist \
     Zdist \
     split_db \
     check_db \
     FSindex_show_bin \
     fragdist \
     KWindex_create
#     print_doublets \
#     FSindex_self_sim \
#     FSindex_cluster \
#     cluster_test \
#     FSindex_rand_perf \
#     FSindex_benchmark \
#     FSindex_eps_net1 \
#     distances \
#    FSindex_delta \
#    FSindex_delta1 \
#    FSindex_rand_perf \
#    FSindex_rand_perf1 \
#    FSindex_rand_perf2 \
#    FSindex_eps_net1 \

distances:

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

FSindex_show_bin:  

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


