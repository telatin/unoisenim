# unoisenim

Nim library implementing USEARCH algorithms for 16S/ITS amplicon sequencing analysis.

## Algorithms

- **UNOISE3** — error correction and ZOTU clustering (`unoisenim/unoise_algo`)
- **SINTAX** — k-mer bootstrap taxonomic classifier (`unoisenim/sintax_algo`)
- **UCHIME2** — chimera detection (`unoisenim/uchime2_algo`)
- **NBC** — Naive Bayesian Classifier, RDP-style (`unoisenim/nbc_algo`)

## Citations

> Edgar RC. UNOISE2: improved error-correction for Illumina 16S and ITS amplicon sequencing.
> *bioRxiv* 081257; doi: [10.1101/081257](https://doi.org/10.1101/081257)

> Edgar RC. SINTAX: a simple non-Bayesian taxonomy classifier for 16S and ITS sequences.
> *bioRxiv* 074161; doi: [10.1101/074161](https://www.biorxiv.org/content/10.1101/074161v1.abstract)

> Wang Q, Garrity GM, Tiedje JM, Cole JR. Naïve Bayesian Classifier for Rapid Assignment of
> rRNA Sequences into the New Bacterial Taxonomy. *Appl Environ Microbiol.* 2007;
> doi: [10.1128/AEM.00062-07](https://doi.org/10.1128/AEM.00062-07)
