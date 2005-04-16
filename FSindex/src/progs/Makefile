SHELL=/bin/sh

FSLIB=libFSindex.a libsarray.a
VPATH=../lib:../sarray
all: $(FSLIB) \
     test_smatrix \
     count_frags \
     FSindex_create \
     FSindex_rng_perf \
     FSindex_rng_perf_thr \
     FSindex_kNN_perf \
     FSindex_ball_size \
     FSindex_delta2 \
     fastadb_test
#     FSindex_stats 
#     partition_test \




distances:

check_db:

split_db:

cluster_test: 

FSindex_rand_perf:

FSindex_rand_queries: 

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

fastadb_sample: 

clean:
	for x in *; do [ -x $$x -a -f $$x.c ] && echo $$x; done\
	 | xargs rm -f; \
	rm -f *.o; \
	rm -f core;

srcclean: clean
	rm -f *~;


