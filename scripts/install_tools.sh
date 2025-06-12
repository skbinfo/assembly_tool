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
    eval "$($HOME/miniconda/bin/conda shell.bash hook)"
    conda init bash
else
    echo "Conda is already installed."
    eval "$(conda shell.bash hook)"
fi

# Ensure Miniconda is prioritized over Anaconda
export PATH="$HOME/miniconda/bin:$PATH"

# Create and activate Conda environment
echo "Creating Conda environment 'genome-annotation'..."
conda env create -f environment.yml || {
    echo "Failed to create Conda environment. Check environment.yml and Conda logs."
    exit 1
}
conda activate genome-annotation || {
    echo "Failed to activate genome-annotation environment."
    exit 1
}

# Install Python package
echo "Installing genome-annotation-pipeline..."
pip install . || {
    echo "Failed to install Python package."
    exit 1
}

# Verify installation
echo "Verifying installation..."
if command -v genome-annotation &> /dev/null; then
    echo "Installation successful! Run 'genome-annotation --help' to see usage."
else
    echo "Installation failed. Check the logs above."
    exit 1
fi
