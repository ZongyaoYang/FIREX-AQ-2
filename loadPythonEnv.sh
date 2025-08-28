#!/bin/bash -l
set -euo pipefail

module purge
module load ncarenv
module load conda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate npl
