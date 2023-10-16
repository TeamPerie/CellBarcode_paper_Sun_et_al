#!/usr/bin/env bash


# sample 1
f="../../data/A1007/A1007R099.fastq.gz"
f_out="../../data/A1007/sample1_filtered.fastq"
zgrep -B 1 -A 2 -E 'ACGGAATG' $f | grep -v "^--$" > $f_out
PATH=$PATH:../../lib/bartender-1.1
bartender_extractor_com -f $f_out -o out -u0,15 -m 0 -p ATCAG[20]ATGTG

bartender_single_com -f out_barcode.txt -o 1M_barcode -d 1

## sample 2
f="../../data/A1007/A1007R100.fastq.gz"
grep -A 3 -E 'PATTERN' $f | grep -v "^--$" > output.fastq
PATH=$PATH:../../lib/bartender-1.1
bartender_extractor_com -f $f -o out_non_umi -m 0 -p TCAAG[0-60]CCGTA

bartender_single_com -f out_non_umi_barcode.txt -o 1M_barcode_non_umi -d 1
