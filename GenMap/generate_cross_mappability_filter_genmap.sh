#!/bin/bash
# ==================================================================================================
# Max Planck Institute for Evolutionary Anthropology 
# Lesage Louison | louison_lesage[at]eva.mpg.de
# Last update: 27.02.2026
# ==================================================================================================
# This script is designed to generate a cross-mappability filter for X given genomes.
# The mappability of the target genome is first computed in order to identify inner repetitions.
# Then, each other fasta are splitted into k-mers and cross-mappability is computed against the target.
# This allows the identification of area of the target where other k-mers align.
# Finally, a mask is produced to allows the identification of the target region which are both unique 
# into the target and not found into the other fasta.
# ==================================================================================================
# This script is the GenMap-version of the pipeline. 
# It does not allows indels but mimic the -n BWA-behavior. It consumes higher RAM than the BWA-version.
# It however consumes less disk storage and is quite quicker.
# It requires Python, GenMap and seqkit.
# ==================================================================================================
# Note: it's also possible to save all positions using --csv or --exclude-pseudo 
# It however consumes A LOT of RAM, >>512 Gb for ~ 2+0.5 Gb genomes
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
    echo "  -n, --missing_prob_err_rate    Missing probability error rate, same behavior as bwa aln -n (default: 0.04)."
    echo "  -o, --output_prefix            Output directory/prefix (default: ./output/)."
    echo "  -j, --n_threads                Number of threads for GenMap (default: 1)."
    echo "  -s, --sampling_rate            Sampling rate for GenMap indexing, between 1 and 64. The lower the value, the faster the mapping, but the higher RAM consumption. (default: 10)."
    echo "  -p, --input_split_chunks       If >1, split the FASTA between sequences, reducing RAM consumption but slowing the cross-mappability process as each chunk requires a new indexion. (default: 1)."
    echo "  -rs, --self_stringency         Minimum ratio (0.0 to 1.0) of overlapping unique k-mers required to validate a base in the target genome. For example, 0.5 requires 50% of the covering k-mers to be strictly unique. (default: 0.5)"
    echo "  -rc, --cross_stringency        Minimum ratio (0.0 to 1.0) of overlapping shared k-mers required to mask a base during cross-mappability. A value of 0.0 masks the target as soon as a single k-mer from another genome aligns to it (highly conservative). (default: 0.0)"
    echo "  -mf, --min_frequency           Minimum frequency of the k-mers in other to be masked (default: 1)"
    echo
    echo "Description:"
    echo "  This script generates a highly specific mappability mask for the target genome. It retains only the strictly unique regions of the target by excluding both internal repetitions (self-mappability) and regions where k-mers from other provided genomes can align (cross-mappability)."
    echo "  It also produces also a bedgraph file per FASTA to allows futher analyses with a genome browser." 
    echo "  It requires Python, bedtools and GenMap."
    exit 0
}

set -euo pipefail
LOCAL_SCRATCH="${TMPDIR:-/tmp}/${USER}_crossmap_$$"
trap 'kill $(jobs -p) 2>/dev/null || true; rm -rf "${LOCAL_SCRATCH:-/tmp/placeholder}"' EXIT TERM INT
echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Pipeline temporary location: $LOCAL_SCRATCH/"
mkdir -p "$LOCAL_SCRATCH"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
missing_prob_err_rate=0.04
input_split_chunks=1
sampling_rate=10
output_prefix="./output/"
n_threads=1

self_stringency=0.5
cross_stringency=0.0
min_frequency=1

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
        -n=*|--missing_prob_err_rate=*)
            missing_prob_err_rate="${arg#*=}"
            ;;
        -o=*|--output_prefix=*)
            output_prefix="${arg#*=}"
            ;;
        -s=*|--sampling_rate=*)
            sampling_rate="${arg#*=}"
            ;;
        -p=*|--input_split_chunks=*)
            input_split_chunks="${arg#*=}"
            ;;
        -rs=*|--self_stringency=*)
            self_stringency="${arg#*=}"
            ;;
        -rc=*|--cross_stringency=*)
            cross_stringency="${arg#*=}"
            ;;
        -mf=*|--min_frequency=*)
            min_frequency="${arg#*=}"
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

# Validate kmer_length
if ! [[ "$kmer_length" =~ ^[0-9]+$ ]] || [[ "$kmer_length" -lt 1 ]]; then
    echo "Error: kmer_length must be a positive integer."
    exit 1
fi

# Validate sampling_rate
if ! [[ "$sampling_rate" =~ ^[0-9]+$ ]] || [[ "$sampling_rate" -lt 1 ]] || [[ "$sampling_rate" -gt 64 ]]; then
    echo "Error: sampling_rate must be an integer between 1 and 64."
    exit 1
fi

# Validate self_stringency and cross_stringency
for val_name in self_stringency cross_stringency; do
    val="${!val_name}"
    if ! awk -v v="$val" 'BEGIN { exit !(v >= 0.0 && v <= 1.0) }'; then
        echo "Error: $val_name must be between 0.0 and 1.0."
        exit 1
    fi
done

# Validate n_threads
if ! [[ "$n_threads" =~ ^[0-9]+$ ]] || [[ "$n_threads" -lt 1 ]]; then
    echo "Error: n_threads must be a positive integer."
    exit 1
fi

for tool in genmap seqkit bedtools python3; do
    if ! command -v "$tool" &>/dev/null; then
        echo "Error: '$tool' not found in PATH."
        exit 1
    fi
done

if [[ ! -f "$SCRIPT_DIR/max_diff.py" ]]; then
    echo "Error: max_diff.py not found in $SCRIPT_DIR."
    exit 1
fi

# Compute max errors as BWA do in bwa_cal_maxdiff (bwtaln.c)
max_e=$(python3 "$SCRIPT_DIR/max_diff.py" "$kmer_length" 0.02 "$missing_prob_err_rate")
echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Computing mappability filter for target (allowing up to $max_e errors)..."

# Sweep line filter to compute stringency quickly
# As GenMap unfortunately do not save each k-mers position and depth, it's mandatory to compute it again after mapping
# This code allows such computation 
# Instead of calculating and saving each position, it calculates the depth variation using the second derivative
apply_sweep_line_filter() {
    local in_bg="$1"; local out_bed="$2"; local k="$3"; local T="$4"; 
    local mode="$5"; local chrom_sizes="$6"; local mf="$7"
    local tmp_events="${LOCAL_SCRATCH}/events_${mode}_$$.tmp"

    awk -v k="$k" -v mode="$mode" -v mf="$mf" 'BEGIN{OFS="\t"} 
    {
        bad = 0;
        if (mode == "self" && $4 > 1) bad = 1; # Non-unique = risque
        else if (mode == "cross" && $4 >= mf) bad = 1; # Partagé = risque
        
        if (bad) {
            print $1, $2,      1;  # Start segment
            print $1, $3,     -1;  # End source
            print $1, $2 + k, -1;  # Transition
            print $1, $3 + k,  1;  # End segment
        }
    }' "$in_bg" | LC_ALL=C sort -k1,1 -k2,2n > "$tmp_events"

    awk -v T="$T" -v chrom_sizes_file="$chrom_sizes" '
        BEGIN { OFS="\t"; while ((getline < chrom_sizes_file) > 0) chrom_len[$1] = $2; }
        {
            chrom=$1; pos=$2; delta=$3;
            if (chrom != cur_c) {
                if (in_r) print cur_c, start_r, (cur_c in chrom_len ? chrom_len[cur_c] : cur_p);
                cur_c=chrom; cur_p=pos; slope=0; depth=0; in_r=0;
            }
            if (pos > cur_p) {
                steps = pos - cur_p;
                if (!in_r) {
                    if (slope > 0 && (T - depth)/slope < steps) { start_r = cur_p + int((T - depth)/slope + 0.999); in_r = 1; }
                    else if (depth > T) { start_r = cur_p; in_r = 1; }
                }
                if (in_r && slope < 0 && (depth - T)/(-slope) < steps) { print cur_c, start_r, cur_p + int((depth - T)/(-slope) + 0.001); in_r = 0; }
                depth += slope * steps; cur_p = pos;
            }
            slope += delta;
        }
        END { if (in_r) print cur_c, start_r, (cur_c in chrom_len ? chrom_len[cur_c] : cur_p); }
    ' "$tmp_events" | bedtools merge > "$out_bed"
}

target_basename=$(basename "$input_target")
target_filename="${target_basename%.*}"
target_idx="${output_prefix}target_index"
overlap_dir="${output_prefix}tracks"
mkdir -p "$overlap_dir"

echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Generating chromosome sizes for bedtools..."
seqkit fx2tab -n -l "$input_target" > "${output_prefix}target.chrom.sizes"

# Get absolute threshold (depth = k * stringence)
thresh_self=$(awk -v k="$kmer_length" -v r="$self_stringency" 'BEGIN { print k * (1 - r) }')
# If stringency is 0, threshold=1 (1 shared k-mer is enough)
thresh_cross=$(awk -v k="$kmer_length" -v r="$cross_stringency" 'BEGIN { print k * r }')

echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Self-mappability depth threshold: $thresh_self"
echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Cross-mappability depth threshold: $thresh_cross"

# ==================================================================================================
# Compute self-mappability for the target (to identify unique k-mers)
# ==================================================================================================

if [[ ! -d $target_idx ]]; then
    echo "[$(date +"%Y.%m.%d-%H:%M:%S")] GenMap target index not found. Indexing..."
    genmap index -F "$input_target" -I "$target_idx" -S "$sampling_rate"
else
    echo "[$(date +"%Y.%m.%d-%H:%M:%S")] GenMap target index found. Skipping this step."
fi 

if [[ ! -f "${output_prefix}target.bedgraph" ]]; then
    echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Computing self-mappability for target..."
    genmap map -I "$target_idx" -K "$kmer_length" -E "$max_e" -bg -fl -O "${output_prefix}target" -T "$n_threads"
    sort -k1,1 -k2,2n "${output_prefix}target.bedgraph" > "${output_prefix}target_sorted.bedgraph"
    mv "${output_prefix}target_sorted.bedgraph" "${output_prefix}target.bedgraph"
else
    echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Self-mappability bedgraph found. Skipping this step."
fi 

if [[ ! -f "${output_prefix}target_unique.bed" ]]; then
    echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Applying stringency filter to target..."
    apply_sweep_line_filter "${output_prefix}target.bedgraph" "${output_prefix}target_non_unique.bed" "$kmer_length" "$thresh_self" "self" "${output_prefix}target.chrom.sizes" 0
    awk '{print $1"\t0\t"$2}' "${output_prefix}target.chrom.sizes" > "${output_prefix}full_genome.bed"
    bedtools subtract -a "${output_prefix}full_genome.bed" -b "${output_prefix}target_non_unique.bed" > "${output_prefix}target_unique.bed"
else
    echo "[$(date +"%Y.%m.%d-%H:%M:%S")] target_unique.bed found. Skipping stringency calculation."
fi

echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Self-mappability done!"

# ==================================================================================================
# Compute cross-mappability with the target (to identify shared k-mers)
# ==================================================================================================

shopt -s nullglob
echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Processing input one by one for cross-mappability..."
for other_path in "$input_fasta_directory"/*.{fa,fasta}; do
    other_basename=$(basename "$other_path")
    other_file="${other_basename%.*}"
    
    if [[ "$other_file" == "$target_filename" ]]; then continue; fi

    if [[ ! -f "$overlap_dir/overlap_${other_file}.bed" ]]; then
        echo "[$(date +"%Y.%m.%d-%H:%M:%S")] > Comparing $target_filename vs $other_file"
    
        tmp_duo="${output_prefix}tmp_$other_file"
        mkdir -p "$tmp_duo"

        # To reduce the risk of OOM, divide the task if required (only useful with --exclude-pseudo or --csv or on tiny machine)
        if [[ "$input_split_chunks" -gt 1 ]]; then
            echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Dividing $other_file in $input_split_chunks chunks..."
            mkdir -p "$tmp_duo/chunks" "$tmp_duo/bedgraphs"

            seqkit split2 "$other_path" -p "$input_split_chunks" -O "$tmp_duo/chunks" --extension .fa

            shopt -s nullglob 
            for chunk_fasta in "$tmp_duo/chunks"/*.part_*.fa; do
                chunk_base=$(basename "$chunk_fasta" .fa)
                echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Processing chunk $chunk_base"
                
                chunk_dir="$tmp_duo/run_$chunk_base"
                mkdir -p "$chunk_dir/fasta" "$chunk_dir/map"
                
                ln -sf "$(realpath "$input_target")" "$chunk_dir/fasta/"
                
                # To avoid a name conflict, add a prefix to the other sequences
                seqkit replace -p "^" -r "${other_file}_" "$chunk_fasta" > "$chunk_dir/fasta/$(basename "$chunk_fasta")"
                
                echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Indexing..."
                genmap index -FD "$chunk_dir/fasta" -I "$chunk_dir/idx" -S "$sampling_rate"

                echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Computing cross-mappability..."
                genmap map -I "$chunk_dir/idx" -K "$kmer_length" -E "$max_e" -bg -fl -O "$chunk_dir/map" -T "$n_threads" 

                sort -k1,1 -k2,2n "$chunk_dir/map/${target_filename}.genmap.bedgraph" > "$tmp_duo/bedgraphs/${chunk_base}_on_target.bedgraph"

                rm -rf "$chunk_dir"
            done
            shopt -u nullglob 

            echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Merging chunk frequencies and computing hits..."
            
            bedtools unionbedg -i "${output_prefix}target.bedgraph" "$tmp_duo/bedgraphs"/*_on_target.bedgraph > "$tmp_duo/all_chunks_comparison.bg"

            # Compute other k-mers count as (Target+Other)-Target
            awk 'BEGIN{OFS="\t"} {
                self_hits = $4;
                total_other_hits = 0;
                for (i=5; i<=NF; i++) {
                    chunk_hits = $i - self_hits;
                    if (chunk_hits > 0) {         
                        total_other_hits += chunk_hits;
                    }
                }
                print $1, $2, $3, total_other_hits
            }' "$tmp_duo/all_chunks_comparison.bg" > "$overlap_dir/${other_file}_hits_on_target.bedgraph"

        else
            mkdir -p "$tmp_duo/fasta" "$tmp_duo/map"

            ln -sf "$(realpath "$input_target")" "$tmp_duo/fasta/"

            # To avoid a name conflict, add a prefix to the other sequences
            seqkit replace -p "^" -r "${other_file}_" "$other_path" > "$tmp_duo/fasta/$(basename "$other_path")"

            if [[ ! -d "$tmp_duo/idx" ]]; then
                echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Indexing..."
                genmap index -FD "$tmp_duo/fasta" -I "$tmp_duo/idx" -S "$sampling_rate"
            else
                echo "[$(date +"%Y.%m.%d-%H:%M:%S")] GenMap cross-index found. Skipping."
            fi 

            echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Computing cross-mappability..."
            genmap map -I "$tmp_duo/idx" -K "$kmer_length" -E "$max_e" -bg -fl -O "$tmp_duo/map" -T "$n_threads" 

            # It produces "$tmp_duo/map/<input>.genmap.bedgraph   

            sort -k1,1 -k2,2n "$tmp_duo/map/${target_filename}.genmap.bedgraph" > "$tmp_duo/${other_file}_on_target.bedgraph"
            #sort -k1,1 -k2,2n "$tmp_duo/map/${other_file}.genmap.bedgraph" > "$overlap_dir/target_on_${other_file}.bedgraph"

            # Make the two results comparable
            bedtools unionbedg -i "${output_prefix}target.bedgraph" "$tmp_duo/${other_file}_on_target.bedgraph" > "$tmp_duo/freq_comparison.bg"
            
            # Compute other k-mers count as (Target+Other)-Target
            awk 'BEGIN{OFS="\t"} { 
                input_hits = $5 - $4;
                if (input_hits < 0) input_hits = 0; 
                print $1, $2, $3, input_hits 
            }' "$tmp_duo/freq_comparison.bg" > "$overlap_dir/${other_file}_hits_on_target.bedgraph"
            
        fi

        # Get the homology between target and other (but remove the repetitive target area)
        echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Generating conserved unique tracks for $other_file..."
        bedtools intersect -a "$overlap_dir/${other_file}_hits_on_target.bedgraph" \
                            -b "${output_prefix}target_unique.bed" \
                            > "$overlap_dir/${other_file}_conserved_unique.bg"

        echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Applying cross-stringency filter..."
        apply_sweep_line_filter "$overlap_dir/${other_file}_hits_on_target.bedgraph" "$overlap_dir/overlap_${other_file}.bed" "$kmer_length" "$thresh_cross" "cross" "${output_prefix}target.chrom.sizes" "$min_frequency"

        rm -rf "$tmp_duo"

    else
        echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Cross-mappability overlap found for $other_file. Skipping this step."
    fi 
done
shopt -u nullglob

echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Merging all overlaps into final mask..."
shopt -s nullglob
overlap_files=("$overlap_dir"/overlap_*.bed)
shopt -u nullglob

if [ ${#overlap_files[@]} -gt 0 ]; then
    cat "${overlap_files[@]}" | sort -k1,1 -k2,2n | bedtools merge > "${output_prefix}all_overlaps.bed"
    
    # Mask containing all overlaps on the target but removing the repetitions within the target
    bedtools subtract -a "${output_prefix}target_unique.bed" -b "${output_prefix}all_overlaps.bed" > "${output_prefix}final_mask.bed"
    
    # Same but keeping the repetitions - is it really useful ?
    # awk '{print $1"\t0\t"$2}' "${output_prefix}target.chrom.sizes" > "${output_prefix}full_target_genome.bed"
    # bedtools subtract -a "${output_prefix}full_target_genome.bed" -b "${output_prefix}all_overlaps.bed" > "${output_prefix}final_mask_with_repetitions.bed"
    # rm -f "${output_prefix}full_target_genome.bed"

else
    echo "[$(date +"%Y.%m.%d-%H:%M:%S")] No overlap found: final_mask is the self-mappability mask."
    cp "${output_prefix}target_unique.bed" "${output_prefix}final_mask.bed"
fi

echo "[$(date +"%Y.%m.%d-%H:%M:%S")] Cross-mappability done!"