# Generate Cross-Mappability Filter

This repository contains scripts designed to compute cross-mappability filter allowing to identify potential k-mers overlapping regions of a target fasta. 

Currently, two scripts are available:

- `generate_cross_mappability_filter_bwa.sh`, which uses `bwa aln` and allows indels, but which is pretty slow and which costs a lot of computing time and storage. [More information here](./BWA/README.md).

- `generate_cross_mappability_filter_genmap.sh`, which uses GenMap, is quicker, integrates a self-mappability filter but which does not allow the consideration of indels and cross-stringency. To be conceptually compatible with the previous one, this version mimic the behavior of `bwa aln -n`. [More information here](./GenMap/README.md).

Please note that these scripts are still under development, results are not guaranteed.

## What for? 

- `generate_cross_mappability_filter_bwa.sh`:
    - you must considerate indels.
    - you have a lot of time and a lot of storage.
    - you want to study risky k-mers but accept to don't know every target regions they mapped. 

- `generate_cross_mappability_filter_genmap.sh`:
    - you accept to do not considerate indels.
    - you want to be quick and storage-respectful.
    - you want to refer each position of the *other* k-mers on the target, no matter its uniqueness (leading to an ultra-conservative filter, `r=0.0`).
