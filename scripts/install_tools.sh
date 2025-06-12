#!/bin/bash
set -e

# Check if Conda is installed
if ! command -v conda &> /dev/null; then
    echo "Conda not found. Installing Miniconda..."
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh
    else
        echo "Unsupported OS. Please install Miniconda manually."
        exit 1
    fi
    bash miniconda.sh -b -p $HOME/miniconda
    export PATH="$HOME/miniconda/bin:$PATH"
    source $HOME/miniconda/etc/profile.d/conda.sh
    conda init
else
    echo "Conda is already installed."
fi

# Create and activate Conda environment
echo "Creating Conda environment 'genome-annotation'..."
conda env create -f environment.yml
conda activate genome-annotation

# Install Python package
echo "Installing genome-annotation-pipeline..."
pip install .

# Verify installation
echo "Verifying installation..."
if command -v genome-annotation &> /dev/null; then
    echo "Installation successful! Run 'genome-annotation --help' to see usage."
else
    echo "Installation failed. Check the logs above."
    exit 1
fi
