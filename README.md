# Genome Annotation Pipeline

A Python-based pipeline for downloading, processing, and annotating genomic assemblies using QUAST and Prokka.

## Features
- Downloads assemblies from NCBI based on organism and location.
- Performs quality assessment with QUAST.
- Annotates genomes with Prokka.
- Supports custom parameters via a YAML configuration file.
- Logs detailed tool outputs to a file, with high-level steps shown interactively.

## Prerequisites
- Linux or macOS (Windows not fully supported due to Conda dependencies).
- Internet connection for downloading dependencies and genomes.
- At least 8GB RAM and 20GB disk space for processing.

## Installation

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/yourusername/genome-annotation-pipeline.git
   cd genome-annotation-pipeline
