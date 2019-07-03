#!/usr/bin/env bash
#
# This wrapper script is now obsolete, it now simply invokes
# medal_ceremony.pl with the specified -genome parameter.
#
# Params:
# $1 = Genome (Ex: GRCh37-lite)
# $2,... = Any additional params that need to be passed to medal_ceremony.pl

# Show usage if insufficient number of params provided
if [ "$#" -eq 0 ]; then about.sh $0; exit 1; fi

medal_ceremony.pl -genome $@
