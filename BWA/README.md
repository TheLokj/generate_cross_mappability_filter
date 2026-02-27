# Generate Cross-Mappability Filter (BWA)

This directory contains a script to generate a cross-mappability filter excluding regions where different fasta files (*e.g.* faunal genomes) map on a specific target (*e.g.* human genome).

Each *other* input fasta from `--input_fasta_directory` is divided into `-k`-mers with a sliding offset of `-s`. K-mers are then mapped on the `--input_target`, allowing up to `e` mismatches where `e` is computed thanks to the `-n` argument, which is the one used by `bwa aln -n`. 

## Dependencies

This script requires BWA, samtools, bedtools, seqkit and GNU parallel.

## How it works

### Script steps

For each *other* input, generates a cross-mappability filter:
    - Split the *other* input into `-k`-mers (and remove duplicated k-mers)
    - Map the k-mers on the target using `bwa aln` 
    - Apply the cross-mappability stringency filter described below

### Output

This script produces different files:

- `overlap_${other_file}.bed`, which contains the ambiguous region of the target where the `-k`-mers from `${other_file}` can map, after filtration based on the input parameters.

### Cross-mappability stringency filter

At any single base $P$ in the target, there are exactly `-k` different k-mers that overlap that position. Even though duplicates are removed from the *other* genome, these `-k` overlapping k-mers have different (shifted) sequences. The depth at position $P$ represents how many of these *other* k-mers are overlapping.

The script applies the cross-stringency (`-rc`) filter using the following threshold:
`thresh_cross`= (`-k`/`-s`) * `-rc`. Therefore, if `-rc=0.0`, a base is masked if at least 1 overlapping *other* k-mer matches the target genome (highly conservative). If `-rc=0.99`, a base is masked only if nearly all overlapping k-mers match the target genome (highly permissive/specific).

**Note that this current definition will probably evolve.**

*Last update: 27.02.2026*