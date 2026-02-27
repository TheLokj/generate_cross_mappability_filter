# Generate Cross-Mappability Filter (BWA)

This directory contains a script to generate a cross-mappability filter excluding regions where different fasta files (*e.g.* faunal genomes) map on a specific target (*e.g.* human genome).

Each *other* input fasta from `--input_fasta_directory` is divided into `-k`-mers with a sliding offset of `-s`. K-mers are then mapped on the `--input_target`, allowing up to `e` mismatches where `e` is computed thanks to the `-n` argument, which is the one used by `bwa aln -n`. 

## Dependencies

This script requires BWA, samtools, bedtools, seqkit and GNU parallel.

## Output

This script produces different files:

- `overlap_${other_file}.bed`, which contains the ambiguous region of the target where the `-k`-mers from `${other_file}` can map, after filtration based on the input parameters.

## How it works

### Script steps

- For each *other* input, generate a cross-mappability filter:
    - Split the *other* input into `-k`-mers (and remove duplicated k-mers)
    - Map the k-mers on the target using `bwa aln` 
    - Apply the cross-mappability stringency filter described below

### Cross-mappability stringency filter

At any single base $P$ in the target, there are exactly `-k` different k-mers that overlap that position. Even though duplicates are removed from the *other* genome, these `-k` overlapping k-mers have different (shifted) sequences. The depth at position $P$ represents how many of these *other* k-mers are overlapping. The script applies the cross-mappability stringency (`-rc`) filter using the following threshold: `thresh_cross`= (`-k`/`-s`) * `-rc`. 

Then, if `-rc=0.0`, a base is highlighted as problematic if at least 1 overlapping *other* k-mer matches the position (highly conservative). If `-rc=0.99`, a base is highlighted only if nearly all overlapping k-mers match the target position (highly permissive/specific). 

If the `-u` parameter is specified, only unique overlapping k-mers are considered in order to align the behavior with the traditional self-mappability definition. 

## Performances

To make this script quicker, you can only increment the number of threads used for mapping `-j`.
You can also increment sliding offset `-s`, but this will significantly reduce the relevance of the results.

*Last update: 27.02.2026*