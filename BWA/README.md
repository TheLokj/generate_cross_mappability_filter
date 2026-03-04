# Generate Cross-Mappability Filter (BWA)

This directory contains a script to generate a cross-mappability filter excluding regions where different fasta files (*e.g.* faunal genomes) map on a specific target (*e.g.* human genome).

Each *other* input fasta from `--input_fasta_directory` is divided into `-k`-mers with a sliding offset of `-s`. K-mers are then mapped on the `--input_target`, allowing up to `e` mismatches where `e` is computed thanks to the `-bn` argument, which is the one used by `bwa aln -n`, and to `bo` opening gaps. 

## Dependencies

This script requires BWA, samtools, bedtools, bedGraphToBigWig, seqkit and GNU parallel.

## Output

This script produces different files:

- `cov_map_${other_file}.{bg,bw}`, which contains the depth of `-k`-mers overlapping each position of the genome region of the target.
- `cov_uniq_${other_file}.{bg,bw}`, which contains the depth of unique `-k`-mers overlapping each position of the genome region of the target.
- `overlap_${other_file}_${cross_stringency}.bed`, which contains the ambiguous region of the target where the `-k`-mers from `${other_file}` can map, after filtration based on `-rc`.
- `all_overlaps_mask.bed` a mask combining each `overlap_${other_file}_${cross_stringency}.bed` generated.

## How it works

### Script steps

- For each *other* input, generate a cross-mappability filter:
    - Split the *other* input into `-k`-mers (and remove duplicated k-mers).
    - Map the k-mers on the target using `bwa aln`.
    - For each position, compute the `u`/`n` ratio and apply the stringency filter as described below.
    - Produce a mask referencing the region of the target where the *other* k-mers are ambiguously able to map.

### Stringency filter

For a position $P$ in the target, there are exactly `-k` different possible k-mers that overlap perfectly that position. After mapping, on `n` k-mers overlapping the position (maximum `-k` if there is no mismatches allowed, overwise more), the number `u` is computed as the number of these k-mers which are unique in the target. Then, the position is highlighted as ambiguous if `u`/`n` >= `-rc`, *i.e. if the proportion of unique k-mer overlappping the position comparing to the total number of k-mers overlappping the position is higher than the stringency*. 

For example, considering a situation with one mismatch allowed and with `-rc=0.5`, if a position is overlapped by 10 exact k-mers and by 4 k-mers with 1 mismatch, the position is only retained if 7 of these k-mers are uniques, regardless their exactitude. In a similar case with `-rc=0.99`, all the 14 k-mers must be only associated with this position. 

Unique k-mers are defined as the k-mers without BWA `XA` tag and with `MAPQ>0`. 

This behavior aims to reproduce the [Heng Li's seqbility](https://github.com/lh3/misc/tree/cc0f36a9a19f35765efb9387389d9f3a6756f08f/seq/seqbility) behavior, which is not directly usable in a cross-mappability context. 

## Performances

To make this script quicker, you can only increment the number of threads used for mapping `-j`. You can also increment sliding offset `-s`, but this will significantly reduce the relevance of the results.

*Last update: 04.03.2026*