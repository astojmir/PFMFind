#!/bin/sh
#
#
# Perhaps we should check if we really have two params
#
printf "file <- \"%s\"\npsfilename <- \"%s\"\n" $1 $2 \
| cat - QD_ddist.plot.R | R --no-save --slave
