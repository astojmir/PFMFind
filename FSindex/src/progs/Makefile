SHELL=/bin/sh

FSLIB=../lib/libFSindex.a

all: FSindex_search \
     FSindex_self_sim \
     FSindex_create \
     FSindex_stats \
     SS_search \
     count_letters \
     Sscan_eps_net \
     print_matrices \
     QD_ddist \
#     FSindex_delta \
#     cluster_test \
#     FSindex_rand_perf \
#     FSindex_benchmark \
#     FSindex_eps_net1 \

cluster_test: $(FSLIB)

FSindex_rand_perf: $(FSLIB)

FSindex_search: $(FSLIB)

FSindex_self_sim: $(FSLIB)

QD_ddist: $(FSLIB) 

print_matrices: $(FSLIB) 

FSindex_create: $(FSLIB)

FSindex_eps_net1: $(FSLIB)

FSindex_delta: $(FSLIB)

FSindex_stats: $(FSLIB)

SS_search: $(FSLIB)

Sscan_eps_net: $(FSLIB)

count_letters: $(FSLIB)

clean:
	for x in *; do [ -x $$x -a -f $$x.c ] && echo $$x; done\
	 | xargs rm -f; \
	rm -f *.o; \
	rm -f core; \

srcclean: clean
	rm -f *~; \


