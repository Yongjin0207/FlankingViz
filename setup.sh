#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

check_program() {
    local prog="$1"
    command -v "$prog" >/dev/null 2>&1 || {
        echo "ERROR: '$prog' not found"
        exit 1
    }
}

download_if_missing() {
    local url="$1"
    local out="$2"

    if [[ -f "$out" && -s "$out" ]]; then
        echo "  - Exists: $out"
        return 0
    fi

    echo "  - Downloading: $url"

    if [[ "$url" == *"data.soybase.org"* ]]; then
        curl --http1.1 -L --fail \
             --retry 50 \
             --retry-delay 5 \
             --retry-all-errors \
             --connect-timeout 30 \
             --speed-time 30 \
             --speed-limit 1024 \
             -C - \
             -o "$out" \
             "$url"
    else
        wget -c --tries=20 --waitretry=5 -O "$out" "$url"
    fi
}

gunzip_to_if_missing() {
    local gz="$1"
    local out="$2"

    if [[ -f "$out" && -s "$out" ]]; then
        echo "  - Exists: $out"
        return 0
    fi

    echo "  - Uncompressing: $gz -> $out"
    gunzip -c "$gz" > "$out"
}

test_gzip_or_fail() {
    local gz="$1"
    echo "  - Testing gzip integrity: $gz"
    gunzip -t "$gz"
}

build_blastdb_if_missing() {
    local fasta="$1"
    local prefix="$2"

    if [[ -f "${prefix}.nin" || -f "${prefix}.nhr" || -f "${prefix}.nsq" ]]; then
        echo "  - BLAST DB already exists: $prefix"
        return 0
    fi

    echo "  - Building BLAST DB: $prefix"
    makeblastdb -in "$fasta" -dbtype nucl -parse_seqids -out "$prefix"
}

build_blastdb_no_parse_if_missing() {
    local fasta="$1"
    local prefix="$2"

    if [[ -f "${prefix}.nin" || -f "${prefix}.nhr" || -f "${prefix}.nsq" ]]; then
        echo "  - BLAST DB already exists: $prefix"
        return 0
    fi

    echo "  - Building BLAST DB without parse_seqids: $prefix"
    makeblastdb -in "$fasta" -dbtype nucl -out "$prefix"
}

write_note() {
    local out="$1"
    shift
    {
        echo "setup_date=$(date '+%Y-%m-%d %H:%M:%S')"
        for line in "$@"; do
            echo "$line"
        done
    } > "$out"
}

convert_msu_gff() {
    local in_file="$1"
    local out_file="$2"

    if [[ -f "$out_file" && -s "$out_file" ]]; then
        echo "  - Exists: $out_file"
        return 0
    fi

    echo "  - Converting MSU GFF: $in_file -> $out_file"
    awk -F'\t' '
    BEGIN{OFS="\t"}
    /^#/ {print; next}
    $3!="gene" {next}
    {
      if ($1 ~ /^Chr[0-9]+$/) {
        n = substr($1,4) + 0
        $1 = sprintf("chr%02d", n)
      }
      print
    }' "$in_file" > "$out_file"
}

normalize_rice_fasta_headers() {
    local in_fa="$1"
    local out_fa="$2"

    if [[ -f "$out_fa" && -s "$out_fa" ]]; then
        echo "  - Exists: $out_fa"
        return 0
    fi

    echo "  - Normalizing rice FASTA headers: $in_fa -> $out_fa"
    awk '
    /^>/ {
        hdr = substr($0, 2)
        if (hdr ~ /^Chr[0-9]+$/) {
            n = substr(hdr, 4) + 0
            printf(">chr%02d\n", n)
        } else if (hdr ~ /^chr[0-9]+$/) {
            n = substr(hdr, 4) + 0
            printf(">chr%02d\n", n)
        } else if (hdr ~ /^[0-9]+$/) {
            n = hdr + 0
            printf(">chr%02d\n", n)
        } else {
            print $0
        }
        next
    }
    { print }
    ' "$in_fa" > "$out_fa"
}

extract_archive_auto() {
    local archive="$1"
    local dest_parent="$2"

    echo "  - Extracting: $archive"
    case "$archive" in
        *.tar.gz|*.tgz)
            tar -xzf "$archive" -C "$dest_parent"
            ;;
        *.tar)
            tar -xf "$archive" -C "$dest_parent"
            ;;
        *)
            echo "ERROR: Unsupported archive format: $archive"
            exit 1
            ;;
    esac
}

setup_rice() {
    local DIR="${ROOT_DIR}/Rice"
    mkdir -p "$DIR"

    echo "========================================"
    echo "Setting up Rice"
    echo "========================================"

    local RICE_GENOME_URL="https://rice.uga.edu/osa1r7_download/osa1_r7.asm.fa.gz"
    local RICE_MSU_RAW_URL="https://rice.uga.edu/osa1r7_download/osa1_r7.all_models.gff3.gz"

    local RAP_REP_URL="https://rapdb.dna.naro.go.jp/download/archive/irgsp1/IRGSP-1.0_representative_2026-02-05.tar.gz"
    local RAP_PRED_URL="https://rapdb.dna.naro.go.jp/download/archive/irgsp1/IRGSP-1.0_predicted_2026-02-05.tar.gz"

    download_if_missing "$RICE_GENOME_URL" "${DIR}/IRGSP1.0.raw.fa.gz"
    test_gzip_or_fail "${DIR}/IRGSP1.0.raw.fa.gz"
    gunzip_to_if_missing "${DIR}/IRGSP1.0.raw.fa.gz" "${DIR}/IRGSP1.0.raw.fa"
    normalize_rice_fasta_headers "${DIR}/IRGSP1.0.raw.fa" "${DIR}/IRGSP1.0.fa"

    download_if_missing "$RICE_MSU_RAW_URL" "${DIR}/msu_raw.gff3.gz"
    test_gzip_or_fail "${DIR}/msu_raw.gff3.gz"
    gunzip_to_if_missing "${DIR}/msu_raw.gff3.gz" "${DIR}/msu_raw.gff3"
    convert_msu_gff "${DIR}/msu_raw.gff3" "${DIR}/msu_chr_genes.gff3"

    download_if_missing "$RAP_REP_URL" "${DIR}/IRGSP-1.0_representative_2026-02-05.tar.gz"
    test_gzip_or_fail "${DIR}/IRGSP-1.0_representative_2026-02-05.tar.gz"

    download_if_missing "$RAP_PRED_URL" "${DIR}/IRGSP-1.0_predicted_2026-02-05.tar.gz"
    test_gzip_or_fail "${DIR}/IRGSP-1.0_predicted_2026-02-05.tar.gz"

    if [[ ! -d "${DIR}/IRGSP-1.0_representative" ]]; then
        extract_archive_auto "${DIR}/IRGSP-1.0_representative_2026-02-05.tar.gz" "$DIR"
    else
        echo "  - Already extracted: ${DIR}/IRGSP-1.0_representative"
    fi

    if [[ ! -d "${DIR}/IRGSP-1.0_predicted" ]]; then
        extract_archive_auto "${DIR}/IRGSP-1.0_predicted_2026-02-05.tar.gz" "$DIR"
    else
        echo "  - Already extracted: ${DIR}/IRGSP-1.0_predicted"
    fi

    if [[ ! -f "${DIR}/rapdb.gff3" || ! -s "${DIR}/rapdb.gff3" ]]; then
        cp "${DIR}/IRGSP-1.0_representative/locus.gff" "${DIR}/rapdb.gff3"
    else
        echo "  - Exists: ${DIR}/rapdb.gff3"
    fi

    if [[ ! -f "${DIR}/rapdb_predicted.gff3" || ! -s "${DIR}/rapdb_predicted.gff3" ]]; then
        cp "${DIR}/IRGSP-1.0_predicted/locus.gff" "${DIR}/rapdb_predicted.gff3"
    else
        echo "  - Exists: ${DIR}/rapdb_predicted.gff3"
    fi

    build_blastdb_if_missing "${DIR}/IRGSP1.0.fa" "${DIR}/IRGSP1.0"

    write_note "${DIR}/download_sources.txt" \
        "genome_url=${RICE_GENOME_URL}" \
        "msu_raw_url=${RICE_MSU_RAW_URL}" \
        "rap_representative_archive_url=${RAP_REP_URL}" \
        "rap_predicted_archive_url=${RAP_PRED_URL}" \
        "rapdb_used_file=IRGSP-1.0_representative/locus.gff -> rapdb.gff3" \
        "rapdb_predicted_used_file=IRGSP-1.0_predicted/locus.gff -> rapdb_predicted.gff3" \
        "msu_conversion=gene_only + ChrN_to_chrNN" \
        "rice_fasta_header_normalization=ChrN_to_chrNN"
}

setup_soybean() {
    local DIR="${ROOT_DIR}/Soybean"
    mkdir -p "$DIR"

    echo "========================================"
    echo "Setting up Soybean"
    echo "========================================"

    local SOYBASE_GENOME_URL="https://data.soybase.org/Glycine/max/genomes/Wm82.gnm6.S97D/glyma.Wm82.gnm6.S97D.genome_main.fna.gz"
    local SOYBASE_ANN_URL="https://data.soybase.org/Glycine/max/annotations/Wm82.gnm6.ann1.PKSW/glyma.Wm82.gnm6.ann1.PKSW.gene_models_main.gff3.gz"

    download_if_missing \
        "${SOYBASE_GENOME_URL}" \
        "${DIR}/glyma.Wm82.gnm6.S97D.genome_main.fna.gz"
    test_gzip_or_fail "${DIR}/glyma.Wm82.gnm6.S97D.genome_main.fna.gz"
    gunzip_to_if_missing \
        "${DIR}/glyma.Wm82.gnm6.S97D.genome_main.fna.gz" \
        "${DIR}/glyma.Wm82.gnm6.S97D.genome_main.fna"

    download_if_missing \
        "${SOYBASE_ANN_URL}" \
        "${DIR}/glyma.Wm82.gnm6.ann1.PKSW.gene_models_main.gff3.gz"
    test_gzip_or_fail "${DIR}/glyma.Wm82.gnm6.ann1.PKSW.gene_models_main.gff3.gz"
    gunzip_to_if_missing \
        "${DIR}/glyma.Wm82.gnm6.ann1.PKSW.gene_models_main.gff3.gz" \
        "${DIR}/glyma.Wm82.gnm6.ann1.PKSW.gene_models_main.gff3"

    build_blastdb_no_parse_if_missing \
        "${DIR}/glyma.Wm82.gnm6.S97D.genome_main.fna" \
        "${DIR}/soybean_soybase"

    write_note "${DIR}/download_sources.txt" \
        "genome_url=${SOYBASE_GENOME_URL}" \
        "soybase_gff_url=${SOYBASE_ANN_URL}" \
        "blastdb_mode=no_parse_seqids"
}
setup_maize() {
    local DIR="${ROOT_DIR}/Maize"
    mkdir -p "$DIR"

    echo "========================================"
    echo "Setting up Maize"
    echo "========================================"

    local MAIZE_BASE="https://download.maizegdb.org/Zm-B104-REFERENCE-CORTEVA-1.0"

    download_if_missing \
        "${MAIZE_BASE}/Zm-B104-REFERENCE-CORTEVA-1.0.fa.gz" \
        "${DIR}/Zm-B104-REFERENCE-CORTEVA-1.0.fa.gz"
    test_gzip_or_fail "${DIR}/Zm-B104-REFERENCE-CORTEVA-1.0.fa.gz"
    gunzip_to_if_missing \
        "${DIR}/Zm-B104-REFERENCE-CORTEVA-1.0.fa.gz" \
        "${DIR}/Zm-B104-REFERENCE-CORTEVA-1.0.fa"

    download_if_missing \
        "${MAIZE_BASE}/Zm-B104-REFERENCE-CORTEVA-1.0_Zm00057aa.1.gff3.gz" \
        "${DIR}/Zm-B104-REFERENCE-CORTEVA-1.0_Zm00057aa.1.gff3.gz"
    test_gzip_or_fail "${DIR}/Zm-B104-REFERENCE-CORTEVA-1.0_Zm00057aa.1.gff3.gz"
    gunzip_to_if_missing \
        "${DIR}/Zm-B104-REFERENCE-CORTEVA-1.0_Zm00057aa.1.gff3.gz" \
        "${DIR}/Zm-B104-REFERENCE-CORTEVA-1.0_Zm00057aa.1.gff3"

    build_blastdb_no_parse_if_missing \
        "${DIR}/Zm-B104-REFERENCE-CORTEVA-1.0.fa" \
        "${DIR}/Zm-B104"

    write_note "${DIR}/download_sources.txt" \
        "genome_url=${MAIZE_BASE}/Zm-B104-REFERENCE-CORTEVA-1.0.fa.gz" \
        "gff_url=${MAIZE_BASE}/Zm-B104-REFERENCE-CORTEVA-1.0_Zm00057aa.1.gff3.gz" \
        "blastdb_mode=no_parse_seqids"
}

main() {
    check_program wget
    check_program curl
    check_program gunzip
    check_program awk
    check_program tar
    check_program makeblastdb

    echo "This script downloads third-party reference files from official sources."
    echo "FlankingViz does not redistribute these files."
    echo

    setup_rice
    setup_soybean
    setup_maize

    echo "========================================"
    echo "All setup completed successfully."
    echo "========================================"
}

main "$@"