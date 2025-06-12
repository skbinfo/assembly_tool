#!/bin/bash
echo "🔧 Setting up Genome Tool..."

# Step 1: Create a virtual environment (optional)
python3 -m venv venv
source venv/bin/activate

# Step 2: Install requirements
pip install --upgrade pip
pip install -r requirements.txt

# Step 3: Make output directories (optional)
mkdir -p output/genomes output/assemblies output/quast_output output/prokka_output

# Step 4: Print usage
echo "✅ Setup complete."
echo "👉 To run: python test.py -c config.yaml -o output/"

