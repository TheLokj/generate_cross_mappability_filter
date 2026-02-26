# Generate Cross-Mappability Filter (BWA)

This directory contains a script to generates a cross-mappability filter excluding regions where different fasta files (*e.g.* faunal genomes) map on a specific target (*e.g.* human genome). 

Each *other* input fasta from `--input_fasta_directory` is divided into `-k`-mers with a sliding offset of `-s`. K-mers are then mapped on the `--input_target`, allowing up to `e` mismatches where `e` is computed thanks to the `-n` argument, which is the one used by `bwa aln -n`. A mask is generated for each *other* fasta, containing regions from the target where no *other* k-mers align.  

## How it works

### Script steps

For each *other* input, it generates a cross-mappability filter:
    - Split the *other* input into `-k`-mers (and duplicated k-mers)
    - Map the k-mers on the target suing `bwa aln`
    - 
    