# Generate Cross-Mappability Filter (BWA)

This directory contains a script to generate a cross-mappability filter excluding regions where different fasta files (*e.g.* faunal genomes) map on a specific target (*e.g.* human genome).

Each *other* input fasta from `--input_fasta_directory` is divided into `-k`-mers with a sliding offset of `-s`. K-mers are then mapped on the `--input_target`, allowing up to `e` mismatches where `e` is computed thanks to the `-n` argument, which is the one used by `bwa aln -n`. 

## How it works

### Script steps

For each *other* input, generates a cross-mappability filter:
    - Split the *other* input into `-k`-mers (and duplicated k-mers)
    - Map the k-mers on the target using `bwa aln`

*Last update: 26.02.2026*