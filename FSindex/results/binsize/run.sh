#! /bin/sh

bin/FSindex_stats db/sprotfull_8_6.ind > sprot_8_6.txt
bin/FSindex_stats db/sprotfull_10_4.ind > sprot_10_4.txt
bin/FSindex_stats db/sprotfull_10_5.ind > sprot_10_5.txt
bin/FSindex_stats db/sprot_F_10_4.ind > sprot_F_10_4.txt
bin/FSindex_stats db/sprot_F_10_5.ind > sprot_F_10_5.txt
bin/FSindex_stats db/nr_10_5.ind > nr_10_5.txt
bin/FSindex_stats db/nr_10_4.ind > nr_10_4.txt
bin/FSindex_stats db/nr_8_6.ind > nr_8_6.txt


bin/FSindex_bin_size -p STAN#ILVM#KR#DEQ#WFYH#GPC -f data/sprot.count \
-n 35289647 -m 8 -M 100 -o sprot_F_8_6_T.txt -b;
bin/FSindex_bin_size -p STAN#ILVM#KR#DEQ#WFYH#GPC -f data/sprot.count \
-n 231749566 -m 8 -M 100 -o nr_8_6_T.txt -b;
bin/FSindex_bin_size -p STAN#ILVM#KR#DEQ#WFYH#GPC -f data/sprot.count \
-n 38970685 -m 8 -M 100 -o sprot_8_6_T.txt -b;
bin/FSindex_bin_size -p STAN#ILVM#KRDEQ#WFYH#GPC -f data/sprot.count \
-n 34833919 -m 10 -M 100 -o sprot_F_10_5_T.txt -b;
bin/FSindex_bin_size -p STAN#ILVM#KRDEQ#WFYH#GPC -f data/sprot.count \
-n 38751330 -m 10 -M 100 -o sprot_10_5_T.txt -b;
bin/FSindex_bin_size -p STAN#ILVM#KRDEQ#WFYH#GPC -f data/sprot.count \
-n 219486720 -m 10 -M 100 -o nr_10_5_T.txt -b;
bin/FSindex_bin_size -p STAN#ILVM#KRDEQ#WFYHGPC -f data/sprot.count \
-n 34833919 -m 10 -M 100 -o sprot_F_10_4_T.txt -b;
bin/FSindex_bin_size -p STAN#ILVM#KRDEQ#WFYHGPC -f data/sprot.count \
-n 38751330 -m 10 -M 100 -o sprot_10_4_T.txt -b;
bin/FSindex_bin_size -p STAN#ILVM#KRDEQ#WFYHGPC -f data/sprot.count \
-n 219486720 -m 10 -M 100 -o nr_10_4_T.txt -b;
