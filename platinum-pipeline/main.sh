#!/bin/bash
set -euo pipefail

##original path = /opt/conda/bin:/opt/conda/bin:/opt/conda/condabin:/usr/local/go/bin:/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/lib/rstudio-server/bin
scripts/get-fastq.sh
PATH="/usr/lib/:/opt/conda/condabin:/usr/local/go/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/lib/rstudio-server/bin"
scripts/qc.sh
PATH="/opt/conda/bin:/opt/conda/condabin:/usr/local/go/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/lib/rstudio-server/bin"
scripts/alignment.sh
scripts/refinement.sh
scripts/calling-variants.sh
scripts/annotating-variants.sh
