# Generate Cross-Mappability Filter

This repository contains scripts designed to compute cross-mappability filters that identify potential k-mers overlapping regions of a target fasta.

Currently, two scripts are available:

- `generate_cross_mappability_filter_bwa.sh`, which uses `bwa aln` and allows indels, but is quite slow and requires significant computing time and storage. [More information here](./BWA/README.md).

- `generate_cross_mappability_filter_genmap.sh`, which uses GenMap, is faster, exhaustive, integrates a self-mappability filter but does not allow the consideration of indels and cross-stringency. To be conceptually compatible with the previous one, this version mimics the behavior of `bwa aln -n`. [More information here](./GenMap/README.md).

Please note that these scripts are still under development and results are not guaranteed.

## What for?

- `generate_cross_mappability_filter_bwa.sh`:
    - you need to consider indels.
    - you have plenty of time and storage resources.
    - you want to study risky k-mers but accept not knowing every target region they map to.

- `generate_cross_mappability_filter_genmap.sh`:
    - you do not consider indels.
    - you want to be quick and mindful of storage.
    - you want to reference every position of the *other* k-mers on the target, regardless of their uniqueness (leading to an ultra-conservative filter, `r=0.0` - this will probably evolve in a future version).