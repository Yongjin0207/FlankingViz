# FlankingViz
An Automated Tool for T-DNA Insertion Site Analysis and Visualization in Rice, Maize, and Soybean

## Overview

FlankingViz is a command-line tool for mapping flanking sequences to plant reference genomes and identifying nearby annotated genes based on genomic coordinates. It supports species-specific configuration through YAML profiles and generates standalone interactive HTML output for insertion site inspection.

The current release includes default profiles for rice, soybean, and maize.

Before running FlankingViz for the first time, users should initialize the local reference resources with `bash setup.sh`.

## Key Features

- BLAST-based mapping of flanking sequences to reference genomes
- Identification of nearby annotated genes from genome annotation files
- Species-specific settings through simple YAML profiles
- Standalone interactive HTML output
- Useful for genome walking, flanking PCR, and insertion site inspection

## Supported Species

Default profiles are provided for:

- Rice (*Oryza sativa*, IRGSP 1.0)
- Soybean (*Glycine max*, Williams 82 v4.0)
- Maize (*Zea mays*, B104)

Additional species can be added by creating new profile files.

## Project Structure

```text
FlankingViz/
├─ flankingviz.py              # main pipeline script
├─ flankingviz.html            # HTML visualization template
├─ setup.sh                    # reference data setup script
├─ profiles/                   # species configuration files
│  ├─ rice.yaml
│  ├─ soybean.yaml
│  └─ maize.yaml
├─ Rice/                       # genome and BLAST database files for rice
├─ Soybean/                    # genome and BLAST database files for soybean
├─ Maize/                      # genome and BLAST database files for maize
├─ seqs/                       # input flanking sequence files
└─ out/                        # output HTML files
```

## Requirements

- Python 3.8 or higher
- BLAST+ (`makeblastdb`, `blastn`)
- Python package: `pyyaml`
- Standard command-line tools: `wget`, `curl`, `gunzip`, `tar`, `awk`

Install the Python dependency with:

```bash
pip install pyyaml
```

On Ubuntu or WSL, BLAST+ can be installed with:

```bash
sudo apt update
sudo apt install -y ncbi-blast+ wget curl gawk tar gzip
```

## Reference Data Setup

FlankingViz uses a setup script to download reference genome and annotation files from the original public sources and to build local BLAST databases.

Run:

```bash
bash setup.sh
```

This script will:

1. Download reference genome and annotation files for rice, soybean, and maize
2. Extract and organize the files into species-specific directories
3. Apply species-specific preprocessing where needed
4. Build local BLAST databases

After setup is complete, FlankingViz can be run directly with the provided profiles.

## Input

Place flanking sequence files in the `seqs/` directory.

Accepted extensions:

- `.seq`
- `.fa`
- `.fasta`
- `.txt`

Each file should contain a single flanking sequence. Multi-record FASTA files are not supported and should be split into separate files before use.

Example:

```fasta
>line1_LB
ATGCTAGCTAGCTAGCTAGC...
```

## Usage

Run the pipeline with a species profile:

```bash
python3 flankingviz.py --profile rice
```

To switch species:

```bash
python3 flankingviz.py --profile soybean
python3 flankingviz.py --profile maize
```

To use a custom profile file:

```bash
python3 flankingviz.py --profile profiles/custom_profile.yaml
```

To specify an output file:

```bash
python3 flankingviz.py --profile rice --out out/rice_flankingviz.html
```

## Profile System

Species-specific settings are defined in YAML profile files. This allows users to switch species or add new genome resources without modifying the main script.

Example (`profiles/rice.yaml`):

```yaml
name: rice
species: rice
reference: irgsp_rap_msu

genome_fasta: Rice/IRGSP1.0.fa
blast_db: Rice/IRGSP1.0

gff_files:
  - Rice/rapdb.gff3
  - Rice/rapdb_predicted.gff3
  - Rice/msu_chr_genes.gff3

gene_feature: gene
attr_priority: [Name, ID]
seqid_normalization: rice
```

## Output

The output is a standalone HTML file that can be opened in any modern web browser.

The HTML report is intended for:

- insertion site inspection
- annotation review
- figure preparation
- record keeping and sharing

## Typical Use Cases

- Validation of transgene insertion sites
- Genome walking and flanking PCR analysis
- Identification of genic and intergenic insertions
- Comparative inspection across multiple plant species

## Adding a New Species

To add a new species:

1. Prepare the reference genome FASTA file
2. Prepare one or more GFF annotation files
3. Build a BLAST database
4. Create a new YAML profile in the `profiles/` directory
5. Run FlankingViz with the new profile

Example:

```bash
python3 flankingviz.py --profile new_species
```

## Notes on Third-Party Data

FlankingViz distributes code and configuration only.

Reference genome, annotation, and archive files are downloaded by each user from the original providers during setup. FlankingViz does not redistribute third-party source files or prebuilt BLAST databases.

Users are responsible for complying with the terms of use, licensing conditions, citation requirements, and access policies of each source database.

Current default profiles use resources from:

- RAP-DB
- Rice Genome Annotation Project / UGA
- NCBI RefSeq
- SoyBase
- MaizeGDB

Some reference resources require minor species-specific preprocessing during setup so that BLAST output and annotation coordinates remain internally consistent.

## Citation of Source Databases

If you use FlankingViz in academic work, please also cite the original genome and annotation resources used in your analysis, according to the policies of the corresponding source databases.

## Citation

If you use FlankingViz in academic work, please cite the corresponding publication when available. Until then, please cite or link to this repository.

## License

MIT License
