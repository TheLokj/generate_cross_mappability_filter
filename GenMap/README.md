# Generate Cross-Mappability Filter (GenMap)

This directory contains a script to generate a cross-mappability filter excluding regions where different fasta files (*e.g.* faunal genomes) map on a specific target (*e.g.* human genome).

Each *other* input fasta from `--input_fasta_directory` is divided into `-k`-mers with a sliding offset of 1. K-mers are then mapped on the `--input_target`, allowing up to `e` mismatches where `e` is computed thanks to the `-n` argument, using the same function as `bwa aln -n`. A mask is generated for each *other* fasta, containing regions from the target where no *other* k-mers align.

The script also generates a final mask combining each mask and a self-mappability mask to allow safe filtration. The `-r` parameters (`-rs` for self-mappability, `-rc` for cross-mappability) allow you to specify a stringency, *i.e.*, for a base, the number of overlapping k-mers required to be classified as identified. For example, with `-k`=35 and `-rs`=0.5, a base is only kept if at least 18 overlapping `-k`-mers are unique in the target. For the moment, cross-mappability stringence filter is not available and probably broke. Then, **each base mapped by a k-mer from *other* is labeled**.  

Note that this script **does not** allow indels as it is based on GenMap. If you want to consider indels, please use the BWA alternative to this script. The other version is, however, slower and requires a lot of disk storage.

## Dependencies

This script requires the associated `max_diff.py` Python script, GenMap (v1.3.0), bedtools, and seqkit.

## How it works

To be quicker, this version of the script uses GenMap. It creates a unique index both for the target fasta and the *other* fasta. It then maps the k-mers from the index on each fasta (target and the *other*). By removing the target repeated area, it allows to identify precisely the ambiguous region of the target when compared to the *other*.

As GenMap saves only the frequency of the k-mers, self-mappability stringency is computed by considering each unique k-mer (*i.e.* with a mappability of 1) and by computing, for each base, how many unique overlapping k-mers are associated with this base. 

### Script steps

1. Generates the self-mappability filter of the target:
    - Index the target k-mers
    - Map the k-mers on the target
    - Apply the stringency filter according to the number of overlapping k-mers on a position

2. For each *other* input, generates a cross-mappability filter:
    - Index the target and the *other* k-mers
    - Map the k-mers on the target
    - Subtract the self-mappability mapped k-mers to identify the *other* k-mers mapped on target

3. Creates a final mask excluding both repeated areas of the target and ambiguous regions that can align with the different *other* sequences.

## Performances

To make this pipeline quicker, you can increment the number of threads used for mapping `-j`.
You can also reduce the GenMap sampling rate parameter `-s` (range from 1 to 64). This will make the mapping quicker but will consume more disk storage and RAM.

The `-p` parameter is only useful if you plan to save each position into a csv using the `--csv` parameter of GenMap (or if you activate the `--exclude-pseudo` option), as it can consume more than 512 Gb of RAM.

Using 64 threads, GenMap requires ~24 hours per input genome of 2 Gb.

*Last update: 26.02.2026*