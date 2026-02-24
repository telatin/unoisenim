## Utility helpers for parsing USEARCH-format FASTA/FASTQ labels.

import strutils

proc parseSize*(label: string): int =
  ## Parses the abundance value from a USEARCH-format FASTA sequence label.
  ##
  ## Abundance is encoded as ``;size=N`` in the sequence identifier, for example
  ## ``"seq1;size=42"``.  The function splits on ``;``, finds the ``size=``
  ## field, and returns its integer value.
  ##
  ## Returns ``0`` if no valid ``size=`` field is present or if parsing fails.
  ##
  ## Example:
  ##   ```nim
  ##   assert parseSize("seq1;size=100") == 100
  ##   assert parseSize("seq1;size=foo") == 0
  ##   assert parseSize("seq1") == 0
  ##   ```
  let parts = label.split(";")
  for p in parts:
    if p.startsWith("size="):
      try:
        return parseInt(p[5..^1])
      except ValueError:
        return 0
  return 0
