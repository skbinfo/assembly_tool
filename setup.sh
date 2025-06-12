#!/bin/bash

# setup.sh
# Script to automate the setup of the Assembly Tool environment
# Creates Conda environment, activates it, and installs required tools

set -e  # Exit on any error

echo "Starting setup for Assembly Tool..."

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Check if Conda is installed
if ! command_exists conda; then
    echo "Error: Conda is not installed. Please install Miniconda or Anaconda first."
    echo "Download from: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

# Initialize Conda for the shell
echo "Initializing Conda..."
# Use eval to handle different Conda installations
CONDA_BASE=$(conda info --base)
if [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
    source "$CONDA_BASE/etc/profile.d/conda.sh"
else
    echo "Error: Cannot find Conda initialization script. Please ensure Conda is properly installed."
    exit 1
fi

# Check if environment.yml exists
if [ ! -f "environment.yml" ]; then
    echo "Error: environment.yml not found in the current directory."
    echo "Please ensure you are in the assembly_tool repository directory."
    exit 1
fi

# Create Conda environment
echo "Creating Conda environment 'assembly_tool_env'..."
conda env create -f environment.yml || {
    echo "Error: Failed to create Conda environment. Check environment.yml for errors."
    exit 1
}

# Activate the environment
echo "Activating environment 'assembly_tool_env'..."
conda activate assembly_tool_env || {
    echo "Error: Failed to activate Conda environment."
    exit 1
}

# Add Bioconda channel if not already added
echo "Configuring Bioconda channels..."
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# Install additional Bioconda tools
echo "Installing ncbi-genome-download, quast, and prokka..."
conda install -y -c bioconda ncbi-genome-download quast prokka || {
    echo "Error: Failed to install Bioconda packages."
    exit 1
}

# Verify installations
echo "Verifying installed tools..."
for tool in ncbi-genome-download quast.py prokka; do
    if command_exists $tool; then
        echo "$tool is installed successfully."
    else
        echo "Error: $tool is not installed correctly."
        exit 1
    fi
done

echo "Setup completed successfully!"
echo "To use the tool, activate the environment with:"
echo "  conda activate assembly_tool_env"
echo "Then run:"
echo "  python assembly_tool.py -c config.yaml -o output"

# Keep the environment activated in the current session
echo "Environment 'assembly_tool_env' is now active in this session."
echo "You can immediately run: python assembly_tool.py -c config.yaml -o output"
