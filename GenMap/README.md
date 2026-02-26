# Generate Cross-Mappability Filter (GenMap)

This directory contains a script to generates a cross-mappability filter excluding regions where different fasta files (*e.g.* faunal genomes) map on a specific target (*e.g.* human genome). 

Each *other* input fasta from `--input_fasta_directory` is divided into `-k`-mers with a sliding offset of 1. K-mers are then mapped on the `--input_target`, allowing up to `e` mismatches where `e` is computed thanks to the `-n` argument, using the same function than `bwa aln -n`. A mask is generated for each *other* fasta, containing regions from the target where no *other* k-mers align.  

The script also generates a final mask combining each mask and a self-mappability mask in order to allow a safe filtration. The `-r` parameters (`-rs` for self-mappability, `-rc` for cross-mappability) allow to specify a stringency, *i.e*, for a basis, the number of overlapping k-mers required to be classified as identified. For example, with `-k`=35 and `-r`=0.5, a base is only kept if at least 18 overlapping `-k`-mers are unique. This parameter is not available for cross-mappability yet.

Note that this script **DO NOT** allo indels as it is based on GenMap. If you want to consider indels, please use the BWA alternative to this script. The other version is however slower and requires A LOT of disk storage.

## Dependencies

This script is requiring the associated `max_diff.py` Python script, GenMap (v1.3.0), bedtools and seqkit.

## How it works

To be quicker, this version of the script use GenMap. It creates an unique index both for the target fasta and the *other* fasta. It then map the k-mers from the index on each fasta (target and the *other*). By removing the target repeted area, it allows to identify precisely the ambiguous region of the target when compared to the *other*. 

As GenMap save only the frequency of the k-mers, strigency is computed by considering each unique k-mer (*i.e* with a mappability of 1) and by computing, for each basis, how many unique overlapping k-mers are associated to this basis. 

### Script steps

1. Generates the self-mappability filter of the target:
    - Index the target k-mers
    - Map the k-mers on the target
    - Apply the stringence filter according to the number of overlapping k-mers on a position

2. For each *other* input, generate a cross-mappability filter:
    - Index the target and the *other* k-mers
    - Map the k-mers on the target
    - Substract the self-mappability mapped k-mers in order to identify the *other* k-mers mapped on target
    
3. Create a final mask excluding both repeted area of the target and ambiguous regions that can align with the different *other*s sequences.

## Performances

To make this pipeline quicker, you can increment the number of threads used for mapping `-j`. 
You can also reduce the GenMap sampling rate parameter `-s` (range from 1 to 64). This will make the mapping quicker but will coonsumes more disk storage and RAM.

The `-p` parameter is only useful if you plan to save each position into a csv using the `--csv` parameter of GenMap (or if you activate the `--exclude-pseudo` option), as it can consumes more than 512 Gb of RAM.

Using 64 threads, GenMap requires ~24 hours per input genome of 2 Gb.

*Last update: 26.02.2026*