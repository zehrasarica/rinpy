#!/usr/bin/env bash

set -e
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
if conda env list | grep -q "rinpy_env"; then
  conda activate "rinpy_env"
else
  exit 1
fi

python rinpy_kras_test.py
