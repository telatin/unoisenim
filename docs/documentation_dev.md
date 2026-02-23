# unoisenim Library Documentation for Nim Developers

`unoisenim` is a high-performance Nim library designed for bioinformatics tasks, specifically focusing on amplicon sequencing data processing. It provides efficient native Nim implementations of key algorithms for sequence clustering (UNOISE), chimera detection (UCHIME2), and taxonomic classification (SINTAX).

## Project Structure

The core library modules are located under the `src/unoisenim/` directory. The `src/unoisenim.nim` file serves as the main entry point for the library, ideally re-exporting key functionalities.

```
/Users/telatina/git/unoisenim/
├───src/
│   ├───sintax.nim                  # CLI tool for SINTAX
│   ├───unoise.nim                  # CLI tool for UNOISE/UCHIME2
│   ├───unoisenim.nim               # Main library entry point
│   └───unoisenim/
│       ├───sintax_algo.nim         # SINTAX algorithm implementation
│       ├───submodule.nim           # Example submodule
│       ├───uchime2_algo.nim        # UCHIME2 algorithm implementation
│       ├───unoise_algo.nim         # UNOISE algorithm implementation
│       └───utils.nim               # Utility functions
```

## Library Modules

Below is a detailed breakdown of the `unoisenim` library's modules, including their exported types and procedures.

### `unoisenim` (src/unoisenim.nim)

This is the main library module. It is intended to re-export the primary components from submodules, making them accessible via a single `import unoisenim` statement.

Currently contains:
-   `proc add*(x, y: int): int`
    -   **Description**: Adds two numbers together. *Note: This is an example placeholder and not part of the core bioinformatics functionality.*

### `unoisenim/unoise_algo` (src/unoisenim/unoise_algo.nim)

Implements the UNOISE algorithm for clustering unique sequences into Zero-radius OTUs (ZOTUs) based on abundance and edit distance.

-   `type UnoiseSeq* = object`
    -   **Description**: Represents a biological sequence with associated metadata.
    -   **Fields**:
        -   `id*`: `string` - The identifier (e.g., FASTA header) of the sequence.
        -   `seq*`: `string` - The nucleotide or amino acid sequence string.
        -   `size*`: `int` - The abundance or size of the sequence, typically parsed from the FASTA header.

-   `type Centroid* = object`
    -   **Description**: Represents a centroid sequence after UNOISE clustering.
    -   **Fields**:
        -   `seqObj*`: `UnoiseSeq` - The `UnoiseSeq` object corresponding to this centroid.
        -   `totalSize*`: `int` - The sum of abundances of all sequences clustered into this centroid.

-   `proc editDistanceLimit*(s1, s2: string, limit: int): int`
    -   **Description**: Calculates the Levenshtein (edit) distance between two strings, but stops if the distance exceeds a specified `limit`. Useful for performance when only a maximum allowed distance is of interest.
    -   **Parameters**:
        -   `s1`: `string` - The first string.
        -   `s2`: `string` - The second string.
        -   `limit`: `int` - The maximum edit distance to consider.
    -   **Returns**: `int` - The edit distance between `s1` and `s2`, or `limit + 1` if the distance is greater than `limit`.

-   `proc extractKmers*(s: string, kmers: var seq[uint16])`
    -   **Description**: Extracts 8-mer (octamer) words from a DNA/RNA sequence, encoding each k-mer as a 16-bit unsigned integer. The `kmers` sequence is cleared and populated.
    -   **Parameters**:
        -   `s`: `string` - The input DNA/RNA sequence.
        -   `kmers`: `var seq[uint16]` - The sequence to which extracted k-mers will be appended.

-   `proc unoise*(sequences: seq[UnoiseSeq], alpha: float = 2.0, minsize: int = 8): seq[Centroid]`
    -   **Description**: Implements the UNOISE algorithm. It sorts sequences by abundance and iteratively processes them, either assigning a sequence to an existing centroid or creating a new centroid based on abundance skew and edit distance.
    -   **Parameters**:
        -   `sequences`: `seq[UnoiseSeq]` - Unique sequences with their abundances.
        -   `alpha`: `float` (default: `2.0`) - Alpha parameter for skew calculation.
        -   `minsize`: `int` (default: `8`) - Minimum abundance for a sequence to be processed.
    -   **Returns**: `seq[Centroid]` - UNOISE-derived ZOTUs, sorted by total abundance descending.

### `unoisenim/uchime2_algo` (src/unoisenim/uchime2_algo.nim)

Provides the UCHIME2 algorithm for identifying and filtering out chimeric sequences.

-   `type AlignPath* = seq[char]`
    -   **Description**: Represents an alignment path as a sequence of characters. 'M' for match/mismatch, 'D' for deletion (in query), 'I' for insertion (in query).

-   `proc globalAlign*(q, t: string, path: var AlignPath): int`
    -   **Description**: Performs a global sequence alignment (Needleman-Wunsch-like with banding) between two sequences.
    -   **Parameters**:
        -   `q`: `string` - The query sequence.
        -   `t`: `string` - The target sequence.
        -   `path`: `var AlignPath` - Will store the alignment path.
    -   **Returns**: `int` - The alignment score (number of differences), or `AlignInf` if sequences are too different.

-   `proc getLeftRight*(q, t: string, path: AlignPath): tuple[diffs, posL0, posL1, posR0, posR1: int]`
    -   **Description**: Determines positions related to the first and last differences in an alignment path, mimicking USEARCH's `GetLeftRight`.
    -   **Parameters**:
        -   `q`: `string` - The query sequence.
        -   `t`: `string` - The target sequence.
        -   `path`: `AlignPath` - The alignment path.
    -   **Returns**: `tuple` - Contains total differences (`diffs`), query positions of first/second differences (`posL0`, `posL1`), and last/second-to-last differences (`posR0`, `posR1`).

-   `proc uchime*(centroids: seq[Centroid], minAbSkew: float = 16.0): seq[bool]`
    -   **Description**: Implements the UCHIME2 algorithm to detect chimeric sequences among a set of centroids. Compares each centroid against potential parents based on alignment and abundance skew.
    -   **Parameters**:
        -   `centroids`: `seq[Centroid]` - Centroid objects, usually sorted by abundance.
        -   `minAbSkew`: `float` (default: `16.0`) - Minimum abundance skew ratio for potential parent consideration.
    -   **Returns**: `seq[bool]` - A sequence indicating `true` for chimeric centroids, `false` otherwise.

### `unoisenim/sintax_algo` (src/unoisenim/sintax_algo.nim)

Contains the SINTAX algorithm for rapid taxonomic classification.

-   `type SintaxIndex* = object`
    -   **Description**: Represents a k-mer index used for SINTAX classification.
    -   **Fields**:
        -   `postings*`: `array[65536, seq[int32]]` - K-mer to sequence index mapping.
        -   `taxStrings*`: `seq[seq[string]]` - Parsed taxonomic ranks for database sequences.

-   `type SintaxHit* = object`
    -   **Description**: Represents the result of a SINTAX classification for a single query sequence.
    -   **Fields**:
        -   `rankNames*`: `seq[string]` - Identified taxonomic ranks.
        -   `rankProbs*`: `seq[float]` - Confidence probabilities for ranks.
        -   `strand*`: `char` - Strand ('+' or '-') yielding best classification.

-   `proc extractTaxRanks*(tax: string): seq[string]`
    -   **Description**: Extracts taxonomic ranks from a SINTAX-formatted taxonomy string (e.g., "d:Bacteria,p:Firmicutes").
    -   **Parameters**:
        -   `tax`: `string` - The taxonomy string.
    -   **Returns**: `seq[string]` - A sequence of rank strings.

-   `proc buildIndex*(seqs: seq[string], taxStrings: seq[string]): SintaxIndex`
    -   **Description**: Builds a k-mer index for SINTAX classification from a collection of sequences and their taxonomies.
    -   **Parameters**:
        -   `seqs`: `seq[string]` - DNA/RNA sequences for the database.
        -   `taxStrings`: `seq[string]` - SINTAX-formatted taxonomy strings corresponding to `seqs`.
    -   **Returns**: `SintaxIndex` - The built SINTAX index.

-   `proc rc*(s: string): string`
    -   **Description**: Computes the reverse complement of a given DNA/RNA sequence.
    -   **Parameters**:
        -   `s`: `string` - The input sequence.
    -   **Returns**: `string` - The reverse complement.

-   `proc sintax*(query: string, idx: SintaxIndex, bootSubset: int = 32, bootIters: int = 100): SintaxHit`
    -   **Description**: Performs SINTAX taxonomic classification on a query sequence against a pre-built index using a bootstrap approach. Considers both forward and reverse complement strands.
    -   **Parameters**:
        -   `query`: `string` - The query DNA/RNA sequence.
        -   `idx`: `SintaxIndex` - The pre-built index.
        -   `bootSubset`: `int` (default: `32`) - Number of k-mers sampled per iteration.
        -   `bootIters`: `int` (default: `100`) - Number of bootstrap iterations.
    -   **Returns**: `SintaxHit` - The classification results.

### `unoisenim/utils` (src/unoisenim/utils.nim)

A module for general utility functions used across the library.

-   `proc parseSize*(label: string): int`
    -   **Description**: Parses the abundance (size) from a FASTA label string. Expects a `size=X;` annotation.
    -   **Parameters**:
        -   `label`: `string` - The FASTA sequence label.
    -   **Returns**: `int` - The parsed size, or `0` if not found or unparseable.

### `unoisenim/submodule` (src/unoisenim/submodule.nim)

An example submodule demonstrating how to structure additional modules within the library.

-   `type Submodule* = object`
    -   **Description**: Represents an example submodule object.
    -   **Fields**:
        -   `name*`: `string` - The name associated with this submodule instance.

-   `proc initSubmodule*(): Submodule`
    -   **Description**: Initialises a new `Submodule` object.
    -   **Returns**: `Submodule` - An instance with its name initialized to "Anonymous".

## Getting Started: Basic Usage Examples

### 1. Using UNOISE to Cluster Sequences

```nim
import unoisenim/unoise_algo
import unoisenim/utils # For parseSize if needed

# Example: Assuming you have a list of UnoiseSeq objects
var mySequences: seq[UnoiseSeq] = @[
  UnoiseSeq(id: "seq1;size=100;", seq: "ATGCATGCATGC", size: 100),
  UnoiseSeq(id: "seq2;size=90;", seq: "ATGCATGCATGT", size: 90),
  UnoiseSeq(id: "seq3;size=50;", seq: "GGCCTTAA", size: 50)
]

let centroids = unoise(mySequences, alpha = 2.5, minsize = 5)

echo "Generated Centroids:"
for centroid in centroids:
  echo "  ID: ", centroid.seqObj.id, ", Total Size: ", centroid.totalSize, ", Sequence: ", centroid.seqObj.seq
```

### 2. Performing SINTAX Classification

```nim
import unoisenim/sintax_algo

# Assuming you have database sequences and their taxonomy strings
var dbSeqs = @["ATGCATGC", "GGCCTTAA"]
var dbTaxStrings = @["d:Bacteria,p:Firmicutes", "d:Archaea,p:Euryarchaeota"]

let sintaxIndex = buildIndex(dbSeqs, dbTaxStrings)

let querySeq = "ATGCATGT"
let hit = sintax(querySeq, sintaxIndex)

echo "SINTAX Classification for '", querySeq, "':"
echo "  Strand: ", hit.strand
echo "  Rank Names: ", hit.rankNames
echo "  Rank Probs: ", hit.rankProbs

# Example of how to interpret a hit (simplified)
if hit.rankNames.len > 0:
  echo "  Top classification: ", hit.rankNames[0], " with confidence ", hit.rankProbs[0]
else:
  echo "  No classification found."
```

### 3. Detecting Chimeras with UCHIME2

```nim
import unoisenim/unoise_algo
import unoisenim/uchime2_algo

# Create some dummy centroids (e.g., from a previous UNOISE run)
var candidateCentroids: seq[Centroid] = @[
  Centroid(seqObj: UnoiseSeq(id: "zotu1", seq: "ATGCATGCATGC", size: 200), totalSize: 200),
  Centroid(seqObj: UnoiseSeq(id: "zotu2", seq: "ATGCATGCGGGG", size: 150), totalSize: 150),
  Centroid(seqObj: UnoiseSeq(id: "zotu3", seq: "ATGCATGCATGCGGGG", size: 50), totalSize: 50) # Potential chimera
]

let isChimeraFlags = uchime(candidateCentroids, minAbSkew = 8.0)

echo "UCHIME2 Chimera Detection Results:"
for i, flag in isChimeraFlags:
  echo "  Centroid '", candidateCentroids[i].seqObj.id, "': is chimera = ", flag
```

## Developer Notes

-   **Documentation**: Ensure all exported types (`*`) and procedures (`proc*`) are documented using Nim's `##` syntax. This allows for automatic documentation generation (e.g., with `nim doc`).
-   **Testing**: New features and bug fixes should be accompanied by appropriate unit tests within the `tests/` directory to maintain code quality and correctness.
-   **Memory Management**: Be mindful of Nim's memory management features, especially when working with large biological sequences and data structures (e.g., `seq` and `string` types). Use `newSeqOfCap` where appropriate to pre-allocate memory and avoid reallocations.
-   **Performance**: Nim's performance is a key advantage. Profile critical sections of code and optimize where necessary, leveraging Nim's low-level control and efficient data structures.
-   **External Libraries**: Verify the established usage of any external libraries within the project before introducing new ones.
