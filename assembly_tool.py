import subprocess
import pandas as pd
import re
import os
from Bio import Entrez
from termcolor import colored
import yaml
import argparse
from tqdm import tqdm
import logging
import sys
import time

# Custom filter for console to show only step messages
class StepFilter(logging.Filter):
    def filter(self, record):
        return record.msg.startswith("[Step ")

# Set up logging
def setup_logging(output_dir):
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    
    # File handler for all logs
    log_file = os.path.join(output_dir, "pipeline.log")
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)
    file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(file_formatter)
    
    # Console handler for step messages only
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.addFilter(StepFilter())
    
    # Clear existing handlers to avoid duplicates
    logger.handlers.clear()
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    
    return logger

def load_config(config_file):
    """Load and validate configuration from YAML file."""
    try:
        with open(config_file, 'r') as file:
            if config_file.endswith('.yaml') or config_file.endswith('.yml'):
                config = yaml.safe_load(file)
            else:
                raise ValueError("Config file must be .yaml or .yml")
        
        # Required fields
        required_fields = [
            'organism_query', 'email', 'organism', 'location', 'reference_genome_url',
            'gff_url', 'group', 'organism_type', 'quast_threads', 'prokka_kingdom', 'prokka_env'
        ]
        missing_fields = [field for field in required_fields if field not in config]
        if missing_fields:
            raise ValueError(f"Missing required fields in config: {', '.join(missing_fields)}")
        
        # Validate prokka_kingdom
        valid_kingdoms = ["Bacteria", "Archaea", "Viruses", "Mitochondria"]
        if config['prokka_kingdom'] not in valid_kingdoms:
            raise ValueError(f"Invalid prokka_kingdom. Must be one of: {', '.join(valid_kingdoms)}")
        
        # Validate organism_type
        if config['organism_type'] not in ['PK', 'EK']:
            raise ValueError("organism_type must be 'PK' or 'EK'")
        
        # Initialize optional params dictionaries if not present
        config.setdefault('quast_params', {})
        config.setdefault('prokka_params', {})
        
        return config
    except Exception as e:
        logger.error(f"Error loading configuration file: {e}")
        raise

def format_params(params_dict):
    """Convert a dictionary of parameters to command-line arguments."""
    params = []
    for key, value in params_dict.items():
        param = f"--{key.replace('_', '-')}"
        if isinstance(value, bool) and value:
            params.append(param)
        elif not isinstance(value, bool):
            params.append(f"{param} {value}")
    return " ".join(params)

def download_reference_genome(url, output_dir1, logger):
    """Download reference genome or GFF file with a spinner."""
    filename = os.path.join(output_dir1, os.path.basename(url))
    if not os.path.exists(filename):
        logger.info(f"Downloading {url} to {output_dir1}")
        with tqdm(total=100, desc="Downloading", bar_format="{desc}: {percentage:3.0f}%|{bar}|") as pbar:
            for _ in range(10):
                time.sleep(0.1)
                pbar.update(10)
            result = subprocess.run(
                f"wget -P {output_dir1} {url}",
                shell=True,
                check=True,
                capture_output=True,
                text=True
            )
            logger.info(f"wget output: {result.stdout}")
            if result.stderr:
                logger.error(f"wget error: {result.stderr}")
    return filename[:-3]

def fetch_assembly_data_from_config(config, logger):
    """Fetch assembly data from NCBI using Entrez."""
    try:
        organism = config.get('organism')
        if not organism:
            raise ValueError("The 'organism' is missing in the config file.")
        
        command = f"""
        esearch -db assembly -query "{organism}[Organism]" |
        esummary |
        xtract -pattern DocumentSummary -element AssemblyAccession,LastMajorReleaseAccession,AssemblyName,AssemblyStatus,Organism,SpeciesName,LastUpdateDate,SubmitterOrganization,FtpPath_Assembly_rpt,FtpPath_Stats_rpt
        """
        
        logger.info(f"Fetching assembly data for {organism}")
        result = subprocess.run(
            command,
            shell=True,
            text=True,
            capture_output=True,
            check=True
        )
        logger.info(f"esearch output: {result.stdout}")
        if result.stderr:
            logger.error(f"esearch error: {result.stderr}")
        
        if not result.stdout.strip():
            raise ValueError("No data returned from the query.")
        
        output = result.stdout.strip()
        data = [line.split("\t") for line in output.split("\n")]
        columns = [
            "AssemblyAccession", "LastMajorReleaseAccession", "AssemblyName", "AssemblyStatus",
            "Organism", "SpeciesName", "LastUpdateDate", "SubmitterOrganization",
            "FtpPath_Assembly_rpt", "FtpPath_Stats_rpt"
        ]
        df = pd.DataFrame(data, columns=columns)
        return df
    except Exception as e:
        logger.error(f"Error fetching assembly data: {e}")
        raise

def filter_rows_by_cities(df, location, num_assemblies=None, logger=None):
    """Filter rows by location and limit to num_assemblies."""
    indian_cities = [
        "IND", "Indian", "india", "India", "Agartala", "Ahmedabad", "Aizawl", "Ajmer", "Allahabad", "Amritsar",
        "Anand", "Avikanagar", "Aurangabad", "Amravati", "Bangalore", "Bareilly", "Batinda", "Belgavi", "Banglore", "Bengaluru",
        "Bhopal", "Berhampur", "Bhagalpur", "Bhilai", "Bhubaneswar", "Bhubaneshwar", "Bilaspur", "Bibinagar", "Calicut",
        "Chandigarh", "Chennai", "Cochin", "Dehradun", "Deoghar", "Delhi", "Dindigul", "Dhanbad",
        "Durgapur", "Faridabad", "Gangtok", "Gandhinagar", "Goa", "Greater Noida", "Gorakhpur", "Gurugram", "Guwahati",
        "Hyderabad", "Telangana", "Imphal", "Indore", "Itanagar", "Jaipur", "Jabalpur", "Jammu", "Jamshedpur",
        "Jalandhar", "Jodhpur", "Jorhat", "Kalaburagi", "Kalpakkam", "Kerala", "Kalyani", "Kanpur", "Karnal", "Kangra", "Kannur",
        "Kashipur", "Kochi", "Kolkata", "Kolkata, Hindupur", "Leh", "Lucknow", "Malda", "Manesar",
        "Manipur", "Mangalagiri", "Mandi", "Mumbai", "Mysuru", "Nagpur", "Nadiad", "Nellore", "Navi Mumbai", "New Delhi",
        "Noida", "Patiala", "Pilani", "Puducherry", "Pune", "Prayagraj", "Bareli", "Raichur", "Raipur", "Rajendranagar",
        "Rajkot", "Rishikesh", "Rohtak", "Roorkee", "Ropar", "Rourkela", "Sambalpur", "Secunderabad", "Sikkim", "Shibpur",
        "Shillong", "Shimla", "Silchar", "Sonepat", "Sri City", "Srinagar", "Surat", "Tadepalligudem", "Tezpur",
        "Thiruvananthapuram", "Thiruvananthapura", "Tiruchirappalli", "Tirupati", "Trichy", "Tuljapur", "Udaipur", "Una",
        "Vadodara", "Varanasi", "Vijayawada", "Visakhapatnam", "Warangal", "Yupia",
        "RCB", "ICMR", "IMTECH", "CSIR", "SPPU", "IIT", "NIT", "IISC", "IIIT", "Central of University", "IGIB", "ICGEB", "CDRI", "SRM", "AIIMS"
    ]
    
    if location.lower() == "india":
        city_pattern = re.compile(r"\b(?:" + "|".join([re.escape(city) for city in indian_cities]) + r")\b", re.IGNORECASE)
        filtered_df = df[df["SubmitterOrganization"].str.contains(city_pattern, na=False)]
    else:
        filtered_df = df[df["SubmitterOrganization"].str.lower().str.contains(location.lower(), na=False)]
    
    if num_assemblies is not None and len(filtered_df) > num_assemblies:
        filtered_df = filtered_df.head(num_assemblies)
        logger.info(f"Limited to {num_assemblies} assemblies for location: {location}")
    else:
        logger.info(f"Filtered {len(filtered_df)} assemblies for location: {location}")
        if num_assemblies is not None and len(filtered_df) < num_assemblies:
            logger.warning(f"Requested {num_assemblies} assemblies, but only {len(filtered_df)} are available")
    
    return filtered_df

def save_filtered_data(df, output_file, logger):
    """Save filtered DataFrame to a TSV file."""
    try:
        df.to_csv(output_file, sep="\t", index=False)
        logger.info(f"Filtered data saved to {output_file}")
    except Exception as e:
        logger.error(f"Error saving file {output_file}: {e}")
        raise

def rename_files(directory, successful_accessions, logger):
    """Rename files in the directory to <accession>.fna format."""
    accession_pattern = re.compile(r"(GC[AF]_\d{9}\.\d)")
    for filename in os.listdir(directory):
        if filename.endswith('.fna'):
            match = accession_pattern.search(filename)
            if match:
                accession = match.group(1)
                if accession in successful_accessions:
                    old_path = os.path.join(directory, filename)
                    new_path = os.path.join(directory, f"{accession}.fna")
                    if old_path != new_path:
                        os.rename(old_path, new_path)
                        logger.info(f"Renamed {old_path} to {new_path}")
                else:
                    logger.warning(f"File {filename} contains accession {accession} not in successful_accessions")
            else:
                logger.warning(f"No valid accession found in filename {filename}")

def process_accessions(assembly_ids_list, output_dir2, group, logger):
    """Download assemblies for GCF and GCA accessions."""
    gcf_accessions = [id for id in assembly_ids_list if id.startswith("GCF")]
    gca_accessions = [id for id in assembly_ids_list if id.startswith("GCA")]
    successful_accessions = []
    
    logger.info(f"Processing {len(gcf_accessions)} GCF and {len(gca_accessions)} GCA accessions")
    
    for accession in tqdm(gcf_accessions + gca_accessions, desc="Downloading assemblies"):
        source = "refseq" if accession.startswith("GCF") else "genbank"
        command = f"ncbi-genome-download -s {source} -l all -F fasta -o {output_dir2} {group} --assembly-accessions {accession}"
        try:
            result = subprocess.run(
                command,
                shell=True,
                check=True,
                capture_output=True,
                text=True
            )
            logger.info(f"ncbi-genome-download output for {accession}: {result.stdout}")
            if result.stderr:
                logger.error(f"ncbi-genome-download error for {accession}: {result.stderr}")
            successful_accessions.append(accession)
        except subprocess.CalledProcessError as e:
            logger.warning(f"Failed to download {accession}: {e.stderr}")
    
    return successful_accessions

def decompress_and_rename(output_dir1, output_dir2, successful_accessions, logger):
    """Decompress .gz files in place, move and rename .fna files."""
    # Decompress .fna.gz files in their original locations
    gunzip_command = f"find {output_dir2} -type f -name '*.fna.gz' -exec gunzip -f {{}} \\;"
    try:
        with tqdm(total=100, desc="Decompressing assembly files", bar_format="{desc}: {percentage:3.0f}%|{bar}|") as pbar:
            for _ in range(10):
                time.sleep(0.1)
                pbar.update(10)
            result = subprocess.run(
                gunzip_command,
                shell=True,
                check=True,
                capture_output=True,
                text=True
            )
            logger.info(f"gunzip output for assemblies: {result.stdout}")
            if result.stderr:
                logger.error(f"gunzip error for assemblies: {result.stderr}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to decompress assembly files: {e.stderr}")
        raise

    # Decompress reference genome and GFF in output_dir1
    if os.path.exists(output_dir1):
        ref_gunzip_command = f"gunzip -f {output_dir1}/*.gz"
        try:
            result = subprocess.run(
                ref_gunzip_command,
                shell=True,
                check=True,
                capture_output=True,
                text=True
            )
            logger.info(f"Reference gunzip output: {result.stdout}")
            if result.stderr:
                logger.error(f"Reference gunzip error: {result.stderr}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to decompress reference files: {e.stderr}")
            raise

    # Move .fna files to output_dir2 and rename them
    for accession in successful_accessions:
        find_command = f"find {output_dir2} -type f -name '*.fna'"
        try:
            result = subprocess.run(
                find_command,
                shell=True,
                check=True,
                capture_output=True,
                text=True
            )
            fna_files = result.stdout.strip().split('\n')
            for fna_file in fna_files:
                if fna_file:
                    new_path = os.path.join(output_dir2, f"{accession}.fna")
                    os.rename(fna_file, new_path)
                    logger.info(f"Moved and renamed {fna_file} to {new_path}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to find or move .fna files for {accession}: {e.stderr}")
            raise

    # Clean up empty directories
    cleanup_command = f"find {output_dir2} -type d -empty -delete"
    try:
        result = subprocess.run(
            cleanup_command,
            shell=True,
            check=True,
            capture_output=True,
            text=True
        )
        logger.info(f"Cleanup output: {result.stdout}")
        if result.stderr:
            logger.error(f"Cleanup error: {result.stderr}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to clean up empty directories: {e.stderr}")

def run_quast(successful_accessions, output_dir2, output_dir3, reference_genome, gff_genome, threads, quast_params, logger):
    """Run QUAST on downloaded genomes with specified threads and additional parameters."""
    quast_params_str = format_params(quast_params)
    for accession in tqdm(successful_accessions, desc="Running QUAST"):
        assembly_path = os.path.join(output_dir2, f"{accession}.fna")
        if not os.path.exists(assembly_path):
            logger.error(f"Assembly file {assembly_path} not found for QUAST")
            continue
        os.makedirs(f"{output_dir3}/{accession}", exist_ok=True)
        quast_command = (
            f"quast.py {assembly_path} -r {reference_genome} -g {gff_genome} "
            f"-o {output_dir3}/{accession} --threads {threads} {quast_params_str}"
        )
        try:
            result = subprocess.run(
                quast_command,
                shell=True,
                check=True,
                capture_output=True,
                text=True
            )
            logger.info(f"QUAST output for {accession}: {result.stdout}")
            if result.stderr:
                logger.error(f"QUAST error for {accession}: {result.stderr}")
            logger.info(f"QUAST analysis completed for {accession}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to run QUAST for {accession}: {e.stderr}")

def run_prokka_with_env(accession, input_file, output_dir4, env_name, kingdom, prokka_params, logger):
    """Run Prokka annotation within a specified conda environment with additional parameters."""
    if not os.path.exists(input_file):
        logger.error(f"Input file {input_file} not found for Prokka")
        return
    prokka_params_str = format_params(prokka_params)
    
    # Get Conda base path
    try:
        conda_base = subprocess.run(
            "conda info --base",
            shell=True,
            check=True,
            capture_output=True,
            text=True
        ).stdout.strip()
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to get Conda base path: {e.stderr}")
        raise

    # Construct Prokka command with Conda initialization
    prokka_command = (
        f"source {conda_base}/etc/profile.d/conda.sh && "
        f"conda activate {env_name} && "
        f"prokka --setupdb && "
        f"prokka --outdir {output_dir4}/{accession} --prefix {accession} "
        f"{input_file} --force --kingdom {kingdom} {prokka_params_str}"
    )
    try:
        with tqdm(total=100, desc=f"Running Prokka for {accession}", bar_format="{desc}: {percentage:3.0f}%|{bar}|") as pbar:
            for _ in range(10):
                time.sleep(0.1)
                pbar.update(10)
            result = subprocess.run(
                prokka_command,
                shell=True,
                check=True,
                capture_output=True,
                text=True,
                executable="/bin/bash"
            )
            logger.info(f"Prokka output for {accession}: {result.stdout}")
            if result.stderr:
                logger.error(f"Prokka error for {accession}: {result.stderr}")
            logger.info(f"Prokka ran successfully for {accession}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running Prokka for {accession}: {e.stderr}")
        raise

def main():
    parser = argparse.ArgumentParser(description="Genomic Assembly and Annotation Pipeline")
    parser.add_argument(
        "-c", "--config",
        required=True,
        help="Path to the configuration file (YAML)"
    )
    parser.add_argument(
        "-o", "--output-dir",
        required=True,
        help="Output directory for all files"
    )
    parser.add_argument(
        "--quast-threads",
        type=int,
        help="Number of threads for QUAST (overrides config file if specified)"
    )
    parser.add_argument(
        "--num-assemblies",
        type=int,
        help="Maximum number of assemblies to fetch and process (default: all)"
    )
    
    args = parser.parse_args()
    
    # Validate num-assemblies
    if args.num_assemblies is not None and args.num_assemblies <= 0:
        logger.error("Error: --num-assemblies must be a positive integer")
        sys.exit(1)
    
    # Set up logging
    logger = setup_logging(args.output_dir)
    
    logger.info(f"Using output directory: {args.output_dir}")
    
    # Load configuration
    config = load_config(args.config)
    
    # Set Entrez email
    Entrez.email = config.get('email')
    
    # Create output directories
    parent_dir = args.output_dir
    os.makedirs(parent_dir, exist_ok=True)
    output_dir1 = os.path.join(parent_dir, "genomes")
    output_dir2 = os.path.join(parent_dir, "assemblies")
    output_dir3 = os.path.join(parent_dir, "quast_output")
    output_dir4 = os.path.join(parent_dir, "prokka_output")
    for dir_path in [output_dir1, output_dir2, output_dir3, output_dir4]:
        os.makedirs(dir_path, exist_ok=True)
    
    # Get QUAST threads
    quast_threads = args.quast_threads if args.quast_threads is not None else config.get('quast_threads')
    
    # Step 1: Download reference genome and GFF
    logger.info(colored("[Step 1/6] Downloading reference genome and annotation file...", "green"))
    reference_genome = download_reference_genome(config.get('reference_genome_url'), output_dir1, logger)
    gff_genome = download_reference_genome(config.get('gff_url'), output_dir1, logger)
    
    # Step 2: Fetch and filter assemblies
    logger.info(colored("[Step 2/6] Fetching and filtering assemblies...", "green"))
    df = fetch_assembly_data_from_config(config, logger)
    filtered_df = filter_rows_by_cities(df, config.get('location'), args.num_assemblies, logger)
    if not filtered_df.empty and not df.empty:
        save_filtered_data(df, os.path.join(parent_dir, "all_assemblies.tsv"), logger)
        save_filtered_data(filtered_df, os.path.join(parent_dir, "filtered_assemblies.tsv"), logger)
    else:
        logger.warning("No data matched the filter.")
        return
    
    # Step 3: Download filtered assemblies
    logger.info(colored("[Step 3/6] Downloading filtered assemblies...", "green"))
    assembly_ids_list = filtered_df["AssemblyAccession"]
    successful_accessions = process_accessions(assembly_ids_list, output_dir2, config.get('group'), logger)
    
    # Step 4: Decompress and rename files
    logger.info(colored("[Step 4/6] Decompressing and renaming files...", "green"))
    decompress_and_rename(output_dir1, output_dir2, successful_accessions, logger)
    
    # Step 5: Run QUAST
    logger.info(colored("[Step 5/6] Running QUAST for assemblies...", "green"))
    run_quast(successful_accessions, output_dir2, output_dir3, reference_genome, gff_genome, quast_threads, config.get('quast_params'), logger)
    
    # Step 6: Run Prokka
    logger.info(colored("[Step 6/6] Running Prokka for assemblies...", "green"))
    if config.get('organism_type') == 'PK':
        for accession in successful_accessions:
            input_file = os.path.join(output_dir2, f"{accession}.fna")
            run_prokka_with_env(
                accession,
                input_file,
                output_dir4,
                config.get('prokka_env'),
                config.get('prokka_kingdom'),
                config.get('prokka_params'),
                logger
            )
    elif config.get('organism_type') == 'EK':
        logger.error(colored("Eukaryotic annotation not implemented yet.", "red"))
    else:
        logger.error("Invalid organism type! Use 'PK' for prokaryotes or 'EK' for eukaryotes")

if __name__ == "__main__":
    main()
