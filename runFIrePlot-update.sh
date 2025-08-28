#!/bin/bash -l
set -euo pipefail

REMOTE=derecho.hpc.ucar.edu
WORKDIR=/glade/u/home/zyang/FIREX-AQ-2

date
echo "Host is $(hostname)"

ssh "$REMOTE" 'bash -lc "
  set -euo pipefail
  module purge
  module load ncarenv
  module load conda
  source \$(conda info --base)/etc/profile.d/conda.sh
  conda activate npl

  cd '"$WORKDIR"'
  python3 firex-plot.py year=today species=isopr,O3_SFC,PM2_5_DRY_SFC,co_fire,nh3,co_anth,c2h6,PM2_5_DRY source=aq-watch domain=d01 height=-1.0
"'

date
