#!/usr/bin/bash

PY_ENV="/home/workspace/environment/pbmc_flow_py_env"

# Create conda environment with Python and R packages
conda create -y -p $PY_ENV -c conda-forge \
    python=3.9 \
    ipykernel jupyterlab \
    numpy==1.24.4 pandas \
    scipy scikit-learn \
    matplotlib seaborn plotly \
    tqdm requests openpyxl h5py \
    dill

# Activate the environment
conda activate $PY_ENV

# Install the Jupyter kernel
python -m ipykernel install --user --name=pbmc_flow_py_env --display-name="Python (PBMC Flow)"