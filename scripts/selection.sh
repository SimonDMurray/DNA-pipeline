#!/bin/bash
set -euo pipefail

grep -w '^#\|^#CHROM\|chr[1-9]\|chr[1-2][0-9]\|chr[X,Y,M]' initial-selected-gatk.vcf > selected-gatk.vcf
