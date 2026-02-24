## unoisenim — Nim library for amplicon sequencing analysis.
##
## Provides efficient implementations of USEARCH algorithms for processing
## 16S/ITS amplicon sequencing data:
##
## * `unoisenim/unoise_algo <unoisenim/unoise_algo.html>`_ —
##   UNOISE3 error correction and ZOTU clustering
## * `unoisenim/uchime2_algo <unoisenim/uchime2_algo.html>`_ —
##   UCHIME2 chimera detection
## * `unoisenim/sintax_algo <unoisenim/sintax_algo.html>`_ —
##   SINTAX taxonomic classification
## * `unoisenim/nbc_algo <unoisenim/nbc_algo.html>`_ —
##   Naive Bayesian Classifier (RDP-style) taxonomy
## * `unoisenim/utils <unoisenim/utils.html>`_ —
##   Utility helpers for FASTA label parsing
##
## **References**
##
## * Edgar RC (2016). UNOISE2: improved error-correction for Illumina 16S
##   and ITS amplicon sequencing. *bioRxiv* doi:10.1101/081257
## * Edgar RC (2016). SINTAX: a simple non-Bayesian taxonomy classifier
##   for 16S and ITS sequences. *bioRxiv* doi:10.1101/074161

import unoisenim/utils
import unoisenim/unoise_algo
import unoisenim/uchime2_algo
import unoisenim/sintax_algo
import unoisenim/nbc_algo

export utils, unoise_algo, uchime2_algo, sintax_algo, nbc_algo
