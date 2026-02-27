#!/bin/bash
# ==================================================================================================
# Max Planck Institute for Evolutionary Anthropology 
# Lesage Louison | louison_lesage[at]eva.mpg.de
# Last update: 27.02.2026
# ==================================================================================================
# This script is designed to generate a cross-mappability filter for X given genomes.
# Each fasta are splitted into k-mers and cross-mappability is computed against the target.
# This allows the identification of area of the target where other k-mers align.
# ==================================================================================================
# This script is the BWA-version of the pipeline. 
# It allows indels (-o) and consumes less memory than the GenMap one.
# It however consumes a lot ot disk storage and is slower than the GenMap one.
# It requires seqkit, parallel, bwa, samtools and bedtools.
# ==================================================================================================

usage() {
    echo "Usage: $0 -i=<input_fasta_directory> -t=<input_target> -k=<kmer_length> [options]"
    echo
    echo "Required Arguments:"
    echo "  -i, --input_fasta_directory    Path to directory containing FASTA files to split into k-mers."
    echo "  -t, --input_target             Path to the target reference FASTA (must be located in the input_fasta_directory)."
    echo "  -k, --kmer_length              Length of k-mers to produce."
    echo
    echo "Optional Arguments:"
    echo "  -bn, --bwa_missing_prob_err_rate    Missing probability error rate, bwa aln -n (default: 0.01)."
    echo "  -bo, --bwa_max_gap_opens            Maximum number or fraction of gap opens, bwa aln -o (default: 2)."
    echo "  -bl, --bwa_seed_length              Seed length, bwa aln -l (default: 16500)."
    echo "  -s,  --offset_step                  Offset step for k-mers sliding (default: 1)."
    echo "  -rc, --cross_stringency             Minimum ratio (0.0 to 1.0) of overlapping shared k-mers required to mask a base during cross-mappability. A value of 0.0 masks the target as soon as a single k-mer from another genome aligns to it (highly conservative). (default: 0.0)"
    echo "  -u,  --unique_only                  If set, only k-mers from 'other' mapping exactly ONCE to the target are considered for cross-mappability masking."
    echo "  -o,  --output_prefix                Output directory/prefix (default: ./output/)."
    echo "  -j,  --n_threads                    Number of threads for parallelisation (default: 1)."
    echo 
    echo "Description:"
    echo "  This script generates a mappability filter excluding the target and the region where k-mers from other FASTA align."
    echo "  It also produces also a bedgraph file per FASTA to allows futher analyses with a genome browser." 
    echo "  It requires samtools, bedtools and bwa."
    exit 0
}

set -euo pipefail
output_prefix="./output/"
n_threads=1
offset_step=1
cross_stringency=0.0
unique_only=0

bwa_missing_prob_err_rate=0.01
bwa_max_gap_opens=2
bwa_seed_length=16500

if [[ "$*" == *"-h"* || "$*" == *"--help"* || $# -eq 0 ]]; then
    usage
fi

for arg in "$@"; do
    case $arg in
        -i=*|--input_fasta_directory=*)
            input_fasta_directory="${arg#*=}"
            ;;
        -t=*|--input_target=*)
            input_target="${arg#*=}"
            ;;
        -k=*|--kmer_length=*)
            kmer_length="${arg#*=}"
            ;;
        -s=*|--offset_step=*)
            offset_step="${arg#*=}"
            ;;
        -bn=*|--bwa_missing_prob_err_rate=*)
            bwa_missing_prob_err_rate="${arg#*=}"
            ;;
        -bo=*|--bwa_max_gap_opens=*)
            bwa_max_gap_opens="${arg#*=}"
            ;;
        -bl=*|--bwa_seed_length=*)
            bwa_seed_length="${arg#*=}"
            ;;
        -rc=*|--cross_stringency=*)
            cross_stringency="${arg#*=}"
            ;;
        -u|--unique_only)
            unique_only=1
            ;;
        -o=*|--output_prefix=*)
            output_prefix="${arg#*=}"
            ;;
        -j=*|--n_threads=*)
            n_threads="${arg#*=}"
            ;;
        *)
            echo "Unknown option $arg"
            exit 1
            ;;
    esac
done

if [[ -z "$input_fasta_directory" || -z "$input_target" || -z "$kmer_length" ]]; then
    echo "Error: Missing required arguments."
    usage
fi

for tool in bwa samtools bedtools seqkit parallel; do
    if ! command -v "$tool" &>/dev/null; then
        echo "Error: '$tool' not found in PATH."
        exit 1
    fi
done

# Validate kmer_length
if ! [[ "$kmer_length" =~ ^[0-9]+$ ]] || [[ "$kmer_length" -lt 1 ]]; then
    echo "Error: kmer_length must be a positive integer."
    exit 1
fi

# Validate cross_stringency
if ! awk -v v="$cross_stringency" 'BEGIN { exit !(v >= 0.0 && v <= 1.0) }'; then
    echo "Error: cross_stringency must be between 0.0 and 1.0."
    exit 1
fi

# Validate n_threads
if ! [[ "$n_threads" =~ ^[0-9]+$ ]] || [[ "$n_threads" -lt 1 ]]; then
    echo "Error: n_threads must be a positive integer."
    exit 1
fi

LOCAL_SCRATCH="${TMPDIR:-/tmp}/${USER}_crossmap_$$"
trap 'kill $(jobs -p) 2>/dev/null || true; rm -rf "${LOCAL_SCRATCH:-/tmp/placeholder}"' EXIT TERM INT
echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Pipeline temporary location: $LOCAL_SCRATCH/"
mkdir -p "$LOCAL_SCRATCH"
mkdir -p "$output_prefix"

target_basename=$(basename "$input_target")
target_filename="${target_basename%.*}"
overlap_dir="${output_prefix}tracks"
mkdir -p "$overlap_dir"

# Get each mapping on the target instead of only a random one
# Generate awk script based on the chosen stringency mode
if [[ "$unique_only" -eq 1 ]]; then
    echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Mode: STRICT ORTHOLOGY (-u). Keeping only strictly unique mappings."
    cat << 'EOF' > "$LOCAL_SCRATCH/parse_xa.awk"
        BEGIN { OFS="\t" }
        function get_ref_len(cigar, def_k,   len, tmp_num, char, i) {
            len = 0; tmp_num = "";
            for (i = 1; i <= length(cigar); i++) {
                char = substr(cigar, i, 1);
                if (char ~ /[0-9]/) { tmp_num = tmp_num char; } 
                else {
                    if (char ~ /[MDN=X]/) len += (tmp_num + 0);
                    tmp_num = "";
                }
            }
            return len > 0 ? len : def_k;
        }
        /^@/ { next }
        # $4 > 0 : mapped | $5 > 0 : MAPQ > 0 (prevents BWA 100k hit limit trap)
        $4 > 0 && $5 > 0 {
            is_unique = 1;
            for (i=12; i<=NF; i++) {
                if ($i ~ /^XA:Z:/ && length($i) > 5) { is_unique = 0; break; }
            }
            if (is_unique == 1) {
                c_len = get_ref_len($6, k);
                print $3, $4 - 1, $4 - 1 + c_len;
            }
        }
EOF
else
    echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Mode: ALL HOMOLOGIES. Keeping primary and XA alternative hits."
    cat << 'EOF' > "$LOCAL_SCRATCH/parse_xa.awk"
        BEGIN { OFS="\t" }
        function get_ref_len(cigar, def_k,   len, tmp_num, char, i) {
            len = 0; tmp_num = "";
            for (i = 1; i <= length(cigar); i++) {
                char = substr(cigar, i, 1);
                if (char ~ /[0-9]/) { tmp_num = tmp_num char; } 
                else {
                    if (char ~ /[MDN=X]/) len += (tmp_num + 0);
                    tmp_num = "";
                }
            }
            return len > 0 ? len : def_k;
        }
        /^@/ { next }
        $4 > 0 {
            c_len = get_ref_len($6, k);
            print $3, $4 - 1, $4 - 1 + c_len;
            
            for (i=12; i<=NF; i++) {
                if ($i ~ /^XA:Z:/) {
                    split(substr($i, 6), hits, ";");
                    for (j=1; j<=length(hits); j++) {
                        if (hits[j] != "") {
                            split(hits[j], fields, ",");
                            chrom = fields[1];
                            pos = int(fields[2]);
                            if (pos < 0) pos = -pos;
                            xa_cigar = fields[3];
                            xa_len = get_ref_len(xa_cigar, k);
                            print chrom, pos - 1, pos - 1 + xa_len;
                        }
                    }
                }
            }
        }
EOF
fi

thresh_cross=$(awk -v k="$kmer_length" -v s="$offset_step" -v r="$cross_stringency" '
function ceil(x) { return (x == int(x)) ? x : int(x) + 1 }
BEGIN {
    max_cov = k / s;
    t = ceil(max_cov) * r;
    t_int = ceil(t);
    print (t_int <= 1 ? 1 : t_int);
}')

if [[ ! -f "${input_target}.bwt" ]]; then
    echo "[$(date +"%Y.%m.%d-%H:%M:%S")] BWA target index not found. Indexing..."
    bwa index "$input_target"
else
    echo "[$(date +"%Y.%m.%d-%H:%M:%S")] BWA target index found. Skipping this step."
fi 

if [[ ! -f "${output_prefix}target.genome.sizes" ]]; then
    echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Genome sizes not found. Computing..."
    samtools faidx "$input_target"
    cut -f1,2 "${input_target}.fai" > "${output_prefix}target.genome.sizes"
fi

# ==================================================================================================
# Compute cross-mappability with the target (to identify shared k-mers)
# ==================================================================================================

echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Processing input one by one for cross-mappability..."

shopt -s nullglob 
for other_path in "$input_fasta_directory"/*.{fa,fasta}; do
    other_base=$(basename "$other_path")
    other_file="${other_base%.*}"

    if [[ "$other_file" == "$target_filename" ]]; then continue; fi

    if [[ ! -f "$overlap_dir/overlap_${other_file}.bed" ]]; then

        echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Comparing $target_filename vs $other_file"
        
        tmp_local="$LOCAL_SCRATCH/splits_${other_file}"
        mkdir -p "$tmp_local"
        
        echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Generating $kmer_length-mers from $other_file... "
        seqkit sliding -W "$kmer_length" -s "$offset_step" "$other_path" | \
        seqkit rmdup -s | \
        seqkit split2 -s 1000000 -O "$tmp_local" --extension .fasta
        
        echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Aligning k-mers from $other_file on $target_filename and extracting all hits..."

        find "$tmp_local/" -name "*.fasta" | parallel -j "$n_threads" \
            "bwa aln -t 1 -n $bwa_missing_prob_err_rate -o $bwa_max_gap_opens -l $bwa_seed_length '$input_target' {} > {.}.sai && \
            bwa samse -n 100000 '$input_target' {.}.sai {} | awk -v k=$kmer_length -f '$LOCAL_SCRATCH/parse_xa.awk' > {.}.bed && \
            rm -f {} {.}.sai"

        echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Sorting coordinates and merging overlaps..."
        
        # Use -u to remove repetitions on a same position and to make -rc reliable
        cat /dev/null "$tmp_local/"*.bed | LC_ALL=C sort -T "$tmp_local" -u -k1,1 -k2,2n > "$tmp_local/overlap_${other_file}_sorted.bed"
        
        bedtools genomecov -i "$tmp_local/overlap_${other_file}_sorted.bed" -g "${output_prefix}target.genome.sizes" -bg | \
        awk -v t="$thresh_cross" '$4 >= t' | \
        bedtools merge -i stdin > "$overlap_dir/overlap_${other_file}.bed"

        rm -rf "$tmp_local"
    fi
done
shopt -u nullglob 

echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Merging all overlaps into final mask..."

shopt -s nullglob
overlap_files=("$overlap_dir"/overlap_*.bed)
shopt -u nullglob

if [ ${#overlap_files[@]} -gt 0 ]; then
    cat "${overlap_files[@]}" | sort -k1,1 -k2,2n | bedtools merge > "${output_prefix}all_overlaps_mask.bed"
else
    echo "[$(date +"%Y.%m.%d-%H:%M:%S")] No overlap found."
fi

rm -rf "$LOCAL_SCRATCH"

echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Cross-mappability done!"