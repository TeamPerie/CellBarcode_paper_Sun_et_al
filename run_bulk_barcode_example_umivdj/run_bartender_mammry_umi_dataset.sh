#!/usr/bin/env bash

f="../../data/VDJ_mammry_gland_meghan_UMI/UMI_m9250_P0_Lum_Rep2/L345R46.R1.fastq"
PATH=$PATH:../../lib/bartender-1.1
bartender_extractor_com -f $f -o out -u0,15 -m 0 -p TCAAG[0-60]CCGTA

bartender_single_com -f out_barcode.txt -o 1M_barcode -d 1

f="../../data/VDJ_mammry_gland_meghan_UMI/noUMI_m9250_P0_Lum_Rep2/L412R11.R1.fastq"
PATH=$PATH:../../lib/bartender-1.1
bartender_extractor_com -f $f -o out_non_umi -m 0 -p TCAAG[0-60]CCGTA

bartender_single_com -f out_non_umi_barcode.txt -o 1M_barcode_non_umi -d 1
