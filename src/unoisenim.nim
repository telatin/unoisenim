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
## * Edgar RC. UNOISE2: improved error-correction for Illumina 16S and ITS amplicon
##   sequencing. *bioRxiv* 081257; doi: `10.1101/081257 <https://doi.org/10.1101/081257>`_
## * Edgar RC. SINTAX: a simple non-Bayesian taxonomy classifier for 16S and ITS
##   sequences. *bioRxiv* 074161; doi: `10.1101/074161 <https://www.biorxiv.org/content/10.1101/074161v1.abstract>`_
## * Wang Q, Garrity GM, Tiedje JM, Cole JR. Naive Bayesian Classifier for Rapid Assignment
##   of rRNA Sequences into the New Bacterial Taxonomy. *Appl Environ Microbiol.* 2007;
##   doi: `10.1128/AEM.00062-07 <https://doi.org/10.1128/AEM.00062-07>`_

import unoisenim/utils
import unoisenim/unoise_algo
import unoisenim/uchime2_algo
import unoisenim/sintax_algo
import unoisenim/nbc_algo
import unoisenim/remove_phix_utils

export utils, unoise_algo, uchime2_algo, sintax_algo, nbc_algo, remove_phix_utils
