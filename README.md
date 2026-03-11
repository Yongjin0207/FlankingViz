# FlankingViz
An Automated Tool for T-DNA Insertion Site Analysis and Visualization in Rice, Maize, and Soybean

Key Features
------------

- BLAST-based mapping of flanking sequences to reference genomes
- Identification of nearby annotated genes from genome annotation files
- Species-specific settings through simple YAML profiles
- Standalone interactive HTML output
- Useful for genome walking, flanking PCR, and insertion site inspection

Supported Species
-----------------

Default profiles are provided for:

- Rice (*Oryza sativa*, IRGSP 1.0)
- Soybean (*Glycine max*, Williams 82)
- Maize (*Zea mays*, B104)

Additional species can be added by creating new profile files.

Project Structure
-----------------

```text
FlankingViz/
├─ flankingviz.py              # main pipeline script
├─ flankingviz.html            # HTML visualization template
├─ profiles/                   # species configuration files
│  ├─ rice.yaml
│  ├─ soybean.yaml
│  └─ maize.yaml
├─ Rice/                       # genome and BLAST database files for rice
├─ Soybean/                    # genome and BLAST database files for soybean
├─ Maize/                      # genome and BLAST database files for maize
├─ seqs/                       # input flanking FASTA files
└─ out/                        # output HTML files
````

## Installation

Requirements:

* Python 3.8 or higher
* BLAST+ (`makeblastdb`, `blastn`)
* `pyyaml`

Install the Python dependency with:

```bash
pip install pyyaml
```

## Input

Place flanking sequence FASTA files in the `seqs/` directory.

Each FASTA file should contain a single flanking sequence. Multi-record FASTA files are not supported and should be split into separate files before use.

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

This will:

1. Map flanking sequences to the selected reference genome using BLAST
2. Identify the nearest annotated gene
3. Generate an interactive HTML visualization

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

Species-specific settings are defined in YAML profile files. This makes it easy to switch species or add new genome resources without changing the main script.

Example (`profiles/rice.yaml`):

```yaml
name: rice
blast_db: Rice/IRGSP1.0
gff_files:
  - Rice/rapdb.gff3
  - Rice/rapdb_predicted.gff3
  - Rice/msu_chr_genes.gff3
```

## Output

The output is a standalone HTML file that can be opened in any modern web browser. It is intended for inspection, figure preparation, and sharing.

## Typical Use Cases

* Validation of transgene insertion sites
* Genome walking and flanking PCR analysis
* Identification of genic and intergenic insertions
* Comparative inspection across multiple plant species

## Adding a New Species

To add a new species:

1. Prepare the reference genome FASTA file and GFF annotation file
2. Build a BLAST database
3. Create a new YAML profile in the `profiles/` directory
4. Run FlankingViz with the new profile

Example:

```bash
python3 flankingviz.py --profile new_species
```

## Citation

If you use FlankingViz in academic work, please cite the corresponding publication when available. Until then, please cite or link to this repository.

## License

MIT License

```
```
