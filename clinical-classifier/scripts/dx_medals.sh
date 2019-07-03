#!/bin/bash
# run medal ceremony on dnanexus
# MNE 3/2017
#
# parameters:
# 1. genome: must be GRCh37-lite or GRCh38 (placeholder, GRCh38 not impl yet)
# 2. data file to process
# 
# TO DO: cloud output file cleanup!!

[ -z "$1" ] && (echo "specify genome (GRCh37-lite or GRCh38)"; exit 1)
[ -z "$2" ] && (echo "specify file of variants (VCF or tab-delimited)"; exit 1)

GENOME=$1
VARIANTS_FILE=$2

OUTFILE=`basename $VARIANTS_FILE`.medals.tab.medals.tab

variants_file_id=$(dx upload "$VARIANTS_FILE" --brief)
mc_job_id=$(dx run app-stjude_medal_ceremony -iinfile=$variants_file_id -y --brief)
dx wait $mc_job_id

dx download $mc_job_id:output_file -o $OUTFILE -f

dx rm $variants_file_id

# TO DO: how to delete the cloud versions of the output files??
