[![DOI](https://zenodo.org/badge/1178350495.svg)](https://doi.org/10.5281/zenodo.19062013)

# FlankingViz

FlankingViz for interpreting and visualizing flanking sequences in rice maize and soybean

## Version Update

**Version: v1.1.2**

This version adds optional vector-assisted junction interpretation, BLAST-based boundary estimation without a vector, alternative mapping candidates with alignment scores, and explicit uncertainty reporting. The complete analysis runs through a single Python entry point, `flankingviz.py`.

The viewer provides candidate selection, editable figures, CSV and SVG exports, and ZIP export through an externally loaded JSZip library. No `vendor/` directory is required. **ZIP export requires an internet connection when opening the HTML report.** Analysis, viewing, CSV export, and individual SVG export do not require JSZip.

The DOI badge links to the existing archive. Check the version recorded in that archive before citing it as this version; a repository update does not automatically update an archived release.

## Overview

FlankingViz maps flanking sequencing reads to plant reference genomes using BLAST, identifies overlapping and neighboring annotated genes, and generates interactive HTML and structured result files. It supports genome walking, flanking PCR, and insertion-site inspection through species-specific YAML profiles.

Vector sequence is optional. Reads may retain vector-derived sequence and may be supplied in either orientation. A single flank can be analyzed independently; paired-flank information is not required.

Before the first analysis, install the dependencies and initialize the reference resources with `bash setup.sh`.

## Key Features

- BLAST-based mapping with automatic handling of both query orientations
- Optional full-vector FASTA or GenBank input to support host–vector boundary interpretation
- Boundary estimation without a vector when a qualifying unaligned read end is present
- Gene-level intragenic/intergenic classification and neighboring-gene distances
- Alternative mapping candidates with bit score, E-value, identity, and alignment length
- Separate reporting of mapping ambiguity and junction uncertainty
- Interactive candidate selection, figure editing, and CSV/SVG/ZIP export
- Species profiles and recorded analysis parameters for reproducibility

## Supported Species

| Profile | Species | Reference resources |
|---|---|---|
| `rice` | Rice (*Oryza sativa*) | IRGSP 1.0; RAP annotated, RAP predicted, and MSU/UGA gene models |
| `soybean` | Soybean (*Glycine max*) | Williams 82, `Wm82.gnm6.S97D`; annotation `gnm6.ann1.PKSW` |
| `maize_b104` | Maize (*Zea mays*) | `Zm-B104-REFERENCE-CORTEVA-1.0` |
| `maize_b73` | Maize (*Zea mays*) | `Zm-B73-REFERENCE-NAM-5.0` |

Exact resource paths are defined in `profiles/` and `setup.sh`. Additional references can be configured through custom profiles.

## Project Structure

```text
FlankingViz/
├─ flankingviz.py              # complete analysis and command-line entry point
├─ flankingviz.html            # interactive viewer template
├─ setup.sh                    # reference download and BLAST database setup
├─ README.md
├─ LICENSE
├─ .gitignore
├─ profiles/
│  ├─ rice.yaml
│  ├─ soybean.yaml
│  ├─ maize_b104.yaml
│  └─ maize_b73.yaml
├─ benchmark/                  # rice example and current expected outputs
│  ├─ benchmark_100_seqs/
│  ├─ benchmark_100.csv
│  ├─ benchmark_100.candidates.csv
│  ├─ benchmark_100.html
│  └─ benchmark_100_truth.csv
├─ validation/                 # synthetic inputs, truth coordinates and checked results
└─ seqs/
   └─ example.seq
```

`Rice/`, `Soybean/`, and `Maize/` are created locally by setup. Analysis output is written to the selected output directory. Downloaded references, private inputs, local validation records, and development scripts are not part of this minimal distribution.

## Requirements

- Python 3 with PyYAML
- NCBI BLAST+ (`blastn` and `makeblastdb`)
- Bash, `wget`, `curl`, `awk`, `tar`, and `gzip`
- A modern web browser for the HTML viewer
- Internet access for initial reference downloads and the viewer's ZIP export library

The analysis was tested on Ubuntu 24.04 under WSL2 with Python 3.12.3, BLAST+ 2.12.0+, and PyYAML. These are the verified environment versions, not minimum-version compatibility claims.

On Ubuntu or WSL:

```bash
sudo apt update
sudo apt install -y python3 python3-yaml ncbi-blast+ wget curl gawk tar gzip
```

Alternatively, install PyYAML into your active Python environment:

```bash
python3 -m pip install pyyaml
```

On shared servers or HPC systems, dependencies can be supplied through system modules, Conda environments, or user-level installations.

## Reference Data Setup

Run from the project directory:

```bash
bash setup.sh
```

The script downloads genome and annotation resources, extracts and organizes them, applies the required species-specific preprocessing, and builds local BLAST databases. It does not install system packages; the required programs must already be available.

Setup prepares all four supplied profiles. Existing resources are reused where the script's checks permit. A completed setup does not need to be repeated before each analysis. Downloaded `.gff3.gz` files are decompressed to the `.gff3` paths used by the profiles; compression itself does not change annotation content.

## Input

Place sequence files in `seqs/` or another input directory. Accepted extensions are `.seq`, `.fa`, `.fasta`, and `.txt`.

Each file must contain **one raw sequence or one FASTA record**. Split multi-record FASTA files before analysis. File stems identify the input lines and must be unique within the input directory. Standard IUPAC ambiguity bases are retained. The default minimum input length is 80 bp; shorter inputs and reads without a qualifying hit are reported with status information.

Reads can include vector sequence; trimming it away is not required. Input sequence text is expected, not chromatogram files. FlankingViz does not perform chromatogram base calling or quality trimming.

## Usage

### Without a vector

```bash
python3 flankingviz.py \
  --profile rice \
  --seqs_dir seqs \
  --out out/rice_flankingviz.html
```

Intragenic/intergenic classification and neighboring-gene distances remain available without vector input. Their coordinate basis depends on the boundary evidence described below.

### With an optional vector

```bash
python3 flankingviz.py \
  --profile rice \
  --seqs_dir seqs \
  --vector vectors/my_vector.fa \
  --out out/rice_vector.html
```

Supply the actual vector sequence, preferably including the T-DNA and backbone. Only a matching fragment needs to be present in a read. Vector FASTA may contain multiple records with unique identifiers; GenBank sequence in the ORIGIN section is also accepted. Export Clone Manager `.cm5` files to FASTA or GenBank first.

An optional `vector_fasta` entry can also be added to a profile. The command-line `--vector` takes precedence.

### Other profiles and options

```bash
python3 flankingviz.py --profile soybean --seqs_dir seqs --out out/soybean.html
python3 flankingviz.py --profile maize_b104 --seqs_dir seqs --out out/maize_b104.html
python3 flankingviz.py --profile maize_b73 --seqs_dir seqs --out out/maize_b73.html
python3 flankingviz.py --profile profiles/custom_profile.yaml --out out/custom.html
```

Use `--threads 4` to request four BLAST threads. The default is one. Choose a new output prefix for each run, or add `--overwrite` to explicitly replace existing outputs. Run `python3 flankingviz.py --help` for all options. Existing legacy rice options (`--db`, `--rap_annot`, `--rap_pred`, and `--msu`) remain available.

Optional paired-flank screening accepts `--pairs pairs.csv` with columns `pair_id,line_a,line_b`. It checks coordinate relationships and flags inconsistencies; it does not reconstruct insertion structure. Ordinary single-flank analysis does not require this option.

## Interpreting Boundaries and Gene Distances

FlankingViz distinguishes three main situations:

| Evidence | Interpretation and annotation basis |
|---|---|
| Adjacent qualifying vector and host alignments | A vector-supported boundary at the adjacent aligned host base |
| A qualifying single unaligned read end without a usable vector boundary | A labelled BLAST-based boundary estimate |
| Boundary unresolved, including ambiguous alignment context | Gene overlap and distances based on the complete aligned host interval |

An unaligned read end does not by itself prove the presence of vector DNA. The automatic no-vector estimate requires a qualifying tail of at least 20 bases, at least 20 unambiguous A/C/G/T bases, no more than 10% ambiguous bases, an opposite-end overhang of at most 10 bases, and no competing host alignment covering at least 20 bases of that tail. Overlapping or separated host/vector alignments are reported with their uncertainty rather than forced into an exact junction.

Coordinates are 1-based. When a point is available, intragenic status means that point falls within an annotated gene interval; otherwise it is intergenic relative to the loaded gene models. Distances use coordinate differences to neighboring gene boundaries. For an unresolved host interval, any gene overlap yields interval-based intragenic context, and distances are measured from the corresponding interval edge. These interval-based distances should not be presented as confirmed junction distances.

Left and right refer to increasing reference coordinates, not transcription direction. `NA` in a neighboring-gene field means there is no gene on that side in the loaded annotation, including reference-end contexts. Gene identifiers and classifications depend on the selected reference and annotation set.

The marker label **Insertion** is a display label. Check the analysis details and fields such as `junction_status`, `coordinate_basis`, and `annotation_basis` before interpreting its location.

## Mapping Candidates

The highest-scoring reported locus is selected initially. Alternative loci, including those on different chromosomes, remain available with their own annotations and alignment scores.

- The viewer displays **Top 10** candidates by default. Choose **All reported** to inspect the full retained set.
- This display setting does not filter the underlying results. A selected candidate outside the top 10 remains accessible when returning to the compact list.
- The default 0.95 bit-score ratio is used to flag similar top scores, **not to discard all candidates below 95% of the best score**.
- The viewer adds no separate bit-score or E-value cutoff. BLAST search settings and pipeline eligibility criteria still apply.
- Default host anchors require at least 40 aligned bases and 80% identity; vector anchors require at least 20 bases and 90% identity. The BLAST E-value threshold is 0.1.
- Default search/reporting limits are 100 target sequences, 100 HSPs per target, and 100 candidate loci per input. Reached limits are flagged; reported candidates should not be assumed exhaustive.

A highest-scoring candidate is a provisional interpretation, not proof of a unique insertion site.

## Output and Interactive Viewer

For `--out out/rice_flankingviz.html`, the pipeline generates:

| File | Contents |
|---|---|
| `rice_flankingviz.html` | Interactive report that opens directly in a browser |
| `rice_flankingviz.csv` | Initial selected candidate or status for each input |
| `rice_flankingviz.candidates.csv` | All reported mapping candidates |
| `rice_flankingviz.json` | Structured analysis results |
| `rice_flankingviz.run.json` | Parameters, versions, hashes, and run provenance |
| `rice_flankingviz.pairs.csv` | Pair-screening results, only when requested |

Open the HTML file, choose an input in **Line selection**, and inspect its **Mapping candidate** entries. Scores and analysis details help distinguish alternative interpretations. Figure labels, arrow direction, and layout can be edited for presentation; these edits do not change the analytical coordinates or inferred evidence. Figure notes can be shown or hidden.

CSV/TSV tables and tables copied from Excel can be loaded into the viewer. Import a `.candidates.csv` file to retain alternatives grouped by species and line. For manually prepared candidate tables without candidate IDs, select the candidate-grouping import mode. Direct `.xlsx` import is not supported. A primary CSV cannot restore alternative candidates that it does not contain.

**Export chosen results (CSV)** exports the current candidate choice for every input line, not just the highlighted line. Individual SVG export saves the current figure. **Export ALL (ZIP)** packages the figures together. Exported SVGs retain analysis metadata.

JSZip 3.10.1 is loaded from jsDelivr when the HTML opens. If that service cannot be reached, ZIP export is unavailable and the viewer displays an explanation when requested. Reload the HTML after restoring connectivity. Core viewing, candidate selection, CSV export, and individual SVG export remain available without that library. No local web server is required.

## Benchmark Example

The original 100 rice reference-derived input sequences are retained in `benchmark/`. Expected CSV and HTML outputs have been regenerated with the current implementation:

```bash
python3 flankingviz.py \
  --profile rice \
  --seqs_dir benchmark/benchmark_100_seqs \
  --out out/benchmark_100.html
```

This processes 100 pseudo-flanking sequences without an optional vector reference. The repository provides these expected outputs:

- `benchmark_100.csv`: one initial selected result or status per input (100 rows).
- `benchmark_100.candidates.csv`: all reported mapping candidates, including alternatives.
- `benchmark_100.html`: the current interactive viewer with the same results and candidate data.

Compare newly generated results with these files using the same input sequences, rice reference resources, and default analysis settings. The expected files were regenerated on 2026-11-14 using BLAST+ 2.12.0+ and four BLAST threads. Different reference resources or software settings may change the results. ZIP export from the HTML requires internet access.

The original `benchmark_100_truth.csv` is preserved as the reference-derived source table; it has not been rewritten from the pipeline output. Its original coordinate and annotation conventions should be distinguished from inferred junction coordinates. These examples support reproducibility and regression checks, not independent experimental junction validation.

## Synthetic Validation

The [`validation/`](validation/) directory contains synthetic junction reads made from the actual **pCAMBIA1300 sequence ([GenBank AF234296.1](https://www.ncbi.nlm.nih.gov/nuccore/AF234296.1))** and host reference sequence. These replace the previous examples containing a randomly generated vector segment; they are not experimental reads from transformed plants.

Rice, soybean, maize B104 and maize B73 each have **40 reads**, representing 20 reference positions in two read orientations. Each read contains 80 bp from an LB- or RB-adjacent T-DNA region and 200 bp of host sequence. The complete 8,958-bp vector reference is supplied in `vector.fa`. Reverse-complement reads test read orientation; they are not paired flanks of one insertion.

Each profile provides:

- `seqs/`: 40 input sequences.
- `vector.fa`: the matching vector reference.
- `truth.json`: coordinates recorded when constructing the reads.
- `result.html`, `result.csv` and `result.candidates.csv`: the report, initial selected results and all reported candidates.
- `assessment.json`: comparison with the construction coordinates, including unresolved and ambiguous outcomes.
- `expected.csv`: a copy of the selected-result CSV for reproducibility, not independent ground truth.

The vector accession and fragment coordinates are documented in `validation/source/`.

After installing dependencies and running `bash setup.sh`, run from the repository root:

```bash
mkdir -p out/validation
for profile in rice soybean maize_b104 maize_b73; do
  python3 flankingviz.py \
    --profile "$profile" \
    --seqs_dir "validation/$profile/seqs" \
    --vector "validation/$profile/vector.fa" \
    --threads 4 \
    --out "out/validation/$profile.html" || break
done
```

To analyze the same inputs without a vector reference, omit `--vector` and choose a different output name. Add `--overwrite` only when intentionally replacing earlier outputs. Download `result.html` and open it in a browser to inspect the stored results. ZIP export requires internet access.

Exact points, boundary intervals and ambiguous or search-limited mappings are reported separately. These examples evaluate computational behavior at known construction coordinates; they do not establish independent biological junction accuracy.

## Profile System and Additional Species

Profiles allow reference resources to be changed without editing the analysis code. For example, `profiles/rice.yaml` contains:

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

To add another reference, prepare its genome FASTA and compatible GFF gene annotations, build a BLAST database, and create a corresponding YAML profile. Ensure that chromosome identifiers and coordinate systems agree between the genome and annotations. Place the profile in `profiles/` and use its name or supply its path with `--profile`.

## Interpretation Limits

- The retained `insertion_strand` field and the alignment arrow describe query alignment orientation relative to the reference. They do **not** establish T-DNA orientation. `alignment_strand` records this explicitly.
- Annotation is at the gene-interval level. Exon, intron, UTR, and promoter classifications are not assigned, even when the source GFF contains these features.
- One flank supports interpretation of the observed side. Complete insertion architecture, copy number, rearrangements, absence of backbone DNA, and effects on gene function require additional evidence.
- Similar loci, repetitive sequences, low-quality bases, and truncated borders can prevent a unique mapping or exact boundary assignment. Review uncertainty fields and alternative candidates.
- Large datasets with many candidate loci can generate large HTML reports and require substantial browser memory. The Top 10 display does not reduce the stored candidate data.

## Notes on Third-Party Data

Reference genomes, annotation archives, and prebuilt BLAST databases are not bundled with this repository. Setup downloads resources from the original providers, including RAP-DB, the Rice Genome Annotation Project/UGA, NCBI, SoyBase, and MaizeGDB as specified in `setup.sh`. The repository includes the original example sequences and reference-derived source table, together with regenerated expected outputs.

Observe the providers' terms of use and cite the genome and annotation resources used in your analysis. Species-specific preprocessing in setup aligns resource identifiers with the supplied profiles.

## Citation

If you use FlankingViz in academic work, cite the corresponding publication when available and identify the software version and reference resources used. The unchanged analysis code is archived at [Zenodo](https://doi.org/10.5281/zenodo.24214700). The pCAMBIA1300 validation update postdates that archive; archive the updated repository release to preserve this dataset snapshot.

## License

MIT License
