# Generate Cross-Mappability Filter

This repository contains scripts designed to compute cross-mappability filters that identify potential k-mers from input that may overlap regions of a target fasta.

Currently, two scripts are available:

- `generate_cross_mappability_filter_bwa.sh`, which uses [BWA](https://www.github.com/lh3/bwa), more precisely `bwa aln`, allows indels, allows a uniqueness-filtering, but is quite slow and requires significant computing time and storage. [More information here](./BWA/README.md).

- `generate_cross_mappability_filter_genmap.sh`, which uses [GenMap](https://www.github.com/cpockrandt/genmap), is faster, exhaustive, integrates a self-mappability filter but does not allow the consideration of indels and cross-stringency. To be conceptually compatible with the previous one, this version mimics the behavior of `bwa aln -n`. [More information here](./GenMap/README.md).

Please note that these scripts are still under development and results are not guaranteed.

## What for?

- `generate_cross_mappability_filter_bwa.sh`:
    - you need to consider indels.
    - you have plenty of time and storage resources, or accept to increment the sliding offset.
    - you want to filter you results using the uniqueness of the overlapping k-mers.
    - you want to study the depth in your favorite genome browser.

- `generate_cross_mappability_filter_genmap.sh`:
    - you do not consider indels.
    - you want to be quick and mindful of storage.
    - you want to reference every position of the *other* k-mers on the target, regardless of their uniqueness in the target genome.
    - you want to study the depth in your favorite genome browser.

Note that the definition of stringency is different between these scripts. The BWA script is computing the stringency according to the uniqueness of the overlapping k-mers, where the GenMap script doesn't. Please refer to each README for complete explaination.  

*Last update: 04.03.2026*