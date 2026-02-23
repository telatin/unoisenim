# unoisenim Library Documentation

`unoisenim` is a high-performance Nim library designed for bioinformatics tasks, specifically focusing on amplicon sequencing data processing.
It provides efficient implementations of key algorithms such as UNOISE for sequence clustering, UCHIME2 for chimera detection, and SINTAX for taxonomic classification.

Citation:

>  UNOISE2: improved error-correction for Illumina 16S and ITS amplicon sequencing
> Robert C. Edgar
> bioRxiv 081257; doi: [10.1101/081257](https://doi.org/10.1101/081257) 

>  SINTAX: a simple non-Bayesian taxonomy classifier for 16S and ITS sequences
> Robert C. Edgar
> bioRxiv 074161; doi: [10.1101/074161](https://www.biorxiv.org/content/10.1101/074161v1.abstract)

## Library Structure

The `unoisenim` library is organized into several modules, primarily located under the `src/unoisenim/` directory:

-   **`unoisenim` (src/unoisenim.nim)**: The main library entry point, intended to re-export core functionalities from its submodules.
    -   `proc add*(x, y: int): int`
        -   Adds two numbers together. (Note: This is currently an example placeholder in the main library file).

-   **`unoisenim/unoise_algo` (src/unoisenim/unoise_algo.nim)**:
    Implements the UNOISE algorithm for clustering unique sequences into Zero-radius OTUs (ZOTUs) based on abundance and edit distance.
    -   `type UnoiseSeq*`: Represents a biological sequence with associated metadata (ID, sequence, size).
    -   `type Centroid*`: Represents a centroid sequence after UNOISE clustering, including the sequence object and total size.
    -   `proc editDistanceLimit*(s1, s2: string, limit: int): int`: Calculates the Levenshtein (edit) distance between two strings, stopping if the distance exceeds a specified `limit`.
    -   `proc extractKmers*(s: string, kmers: var seq[uint16])`: Extracts 8-mer (octamer) words from a DNA/RNA sequence.
    -   `proc unoise*(sequences: seq[UnoiseSeq], alpha: float = 2.0, minsize: int = 8): seq[Centroid]`: Performs the UNOISE algorithm for sequence clustering.

-   **`unoisenim/uchime2_algo` (src/unoisenim/uchime2_algo.nim)**:
    Provides the UCHIME2 algorithm for identifying and filtering out chimeric sequences from a set of candidate ZOTUs.
    -   `type AlignPath*`: Represents an alignment path as a sequence of characters ('M', 'D', 'I').
    -   `proc globalAlign*(q, t: string, path: var AlignPath): int`: Performs a global sequence alignment between two sequences with a banded optimization.
    -   `proc getLeftRight*(q, t: string, path: AlignPath): tuple[...]`: Helps determine positions related to the first and last differences in an alignment path.
    -   `proc uchime*(centroids: seq[Centroid], minAbSkew: float = 16.0): seq[bool]`: Executes the UCHIME2 algorithm to detect chimeric sequences among centroids.

-   **`unoisenim/sintax_algo` (src/unoisenim/sintax_algo.nim)**:
    Contains the SINTAX algorithm for rapid taxonomic classification of DNA/RNA sequences against a reference database.
    -   `type SintaxIndex*`: Represents a k-mer index used for SINTAX classification (k-mer postings, taxonomy strings).
    -   `type SintaxHit*`: Represents the result of a SINTAX classification (rank names, probabilities, strand).
    -   `proc extractTaxRanks*(tax: string): seq[string]`: Extracts taxonomic ranks from a SINTAX-formatted taxonomy string.
    -   `proc buildIndex*(seqs: seq[string], taxStrings: seq[string]): SintaxIndex`: Builds a k-mer index from sequences and their associated taxonomy strings.
    -   `proc rc*(s: string): string`: Computes the reverse complement of a given DNA/RNA sequence.
    -   `proc sintax*(query: string, idx: SintaxIndex, bootSubset: int = 32, bootIters: int = 100): SintaxHit`: Performs SINTAX taxonomic classification on a query sequence.

-   **`unoisenim/utils` (src/unoisenim/utils.nim)**:
    A collection of utility functions.
    -   `proc parseSize*(label: string): int`: Parses the abundance (size) from a FASTA label string (e.g., `size=X;`).

-   **`unoisenim/submodule` (src/unoisenim/submodule.nim)**:
    An example submodule within the library.
    -   `type Submodule*`: Represents an example submodule object.
    -   `proc initSubmodule*(): Submodule`: Initializes a new `Submodule` object.

## Command-Line Tools

The `unoisenim` project also includes two main command-line tools built on top of the library, located in `src/`:

### `sintax` (src/sintax.nim)

A command-line tool for performing SINTAX taxonomic classification.

#### Usage:

```
sintax [options]

USEARCH sintax algorithm implementation in Nim

Options:
  -i, --input        Input query FASTA file
  -d, --database     Reference database FASTA file with taxonomy annotations
  -t, --tabbedout    Output tabbed file
  -c, --cutoff       Confidence cutoff (default 0.8)
  -h, --help         Show this help message
```

#### Example:

```bash
sintax -i query.fasta.gz -d ref_database.fasta -t classification_results.txt -c 0.9
```

### `unoise` (src/unoise.nim)

A command-line tool for running the UNOISE algorithm for sequence clustering followed by UCHIME2 for chimera filtering.

#### Usage:

```
unoise [options] <input>

USEARCH unoise algorithm implementation in Nim

Arguments:
  input              Input FASTA file with size=X; annotations

Options:
  -z, --zotus        Output ZOTUs FASTA file
  -a, --alpha        Alpha parameter for skew calculation (default 2.0)
  -m, --minsize      Minimum abundance for a sequence to be processed (default 8)
  -h, --help         Show this help message
```

#### Example:

```bash
unoise input_sequences.fasta.gz -z final_zotus.fasta -a 3.0 -m 10
```

## Developer Notes

When contributing or extending the `unoisenim` library:

-   Adhere to existing coding style and conventions.
-   Ensure new features or bug fixes are covered by appropriate unit tests in the `tests/` directory.
-   Utilize Nim's `##` syntax for documenting exported types and procedures to facilitate automatic documentation generation.
