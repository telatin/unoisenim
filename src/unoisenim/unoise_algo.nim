## UNOISE3 error correction and ZOTU clustering algorithm.
##
## Implements the UNOISE3 algorithm for denoising amplicon sequences into
## Zero-radius OTUs (ZOTUs).  Sequences are sorted by abundance and clustered
## using the UNOISE skew formula and banded Levenshtein distance.
##
## **Reference:** Edgar RC (2016). UNOISE2: improved error-correction for
## Illumina 16S and ITS amplicon sequencing. *bioRxiv* doi:10.1101/081257

import math, algorithm

type
  UnoiseSeq* = object
    ## An amplicon sequence with its identifier and abundance count.
    id*: string   ## Sequence identifier (from FASTA header, without ``>``)
    seq*: string  ## DNA sequence string
    size*: int    ## Abundance (copy count) parsed from ``size=`` annotation

  Centroid* = object
    ## A cluster centroid (ZOTU candidate) produced by `unoise`.
    seqObj*: UnoiseSeq  ## The representative sequence for this centroid
    totalSize*: int     ## Accumulated abundance across all merged sequences

proc editDistanceLimitImpl(s1, s2: string, limit: int, v0: var seq[int],
    v1: var seq[int]): int =
  let len1 = s1.len
  let len2 = s2.len
  let limitPlusOne = limit + 1
  if abs(len1 - len2) > limit: return limitPlusOne
  if len1 == 0: return if len2 <= limit: len2 else: limitPlusOne
  if len2 == 0: return if len1 <= limit: len1 else: limitPlusOne
  if limit == 0: return if s1 == s2: 0 else: 1

  if v0.len < len2 + 1:
    v0.setLen(len2 + 1)
    v1.setLen(len2 + 1)
  for j in 0..len2:
    if j <= limit:
      v0[j] = j
    else:
      v0[j] = limitPlusOne

  for i in 1..len1:
    let minJ = max(1, i - limit)
    let maxJ = min(len2, i + limit)
    if minJ > maxJ:
      return limitPlusOne

    v1[0] = if i <= limit: i else: limitPlusOne
    if minJ > 1:
      v1[minJ - 1] = limitPlusOne

    var minVal = limitPlusOne
    for j in minJ..maxJ:
      let cost = if s1[i - 1] == s2[j - 1]: 0 else: 1
      let ins = v1[j - 1] + 1
      let del = v0[j] + 1
      let sub = v0[j - 1] + cost
      var best = if ins < del: ins else: del
      if sub < best:
        best = sub
      v1[j] = best
      if best < minVal:
        minVal = best

    if maxJ < len2:
      v1[maxJ + 1] = limitPlusOne

    if minVal > limit: return limitPlusOne
    swap(v0, v1)

  if v0[len2] <= limit: return v0[len2]
  return limitPlusOne

proc editDistanceLimit*(s1, s2: string, limit: int): int =
  ## Computes the Levenshtein edit distance between ``s1`` and ``s2`` with
  ## early termination when the distance would exceed ``limit``.
  ##
  ## Uses a banded dynamic-programming approach with O(n) space.
  ## Returns ``limit + 1`` as soon as the true distance is known to exceed
  ## ``limit``, making it efficient for threshold-based comparisons.
  ##
  ## **Parameters**
  ## * ``s1``, ``s2`` — DNA sequences to compare
  ## * ``limit`` — maximum edit distance; returns ``limit+1`` if exceeded
  var v0 = newSeq[int](s2.len + 1)
  var v1 = newSeq[int](s2.len + 1)
  return editDistanceLimitImpl(s1, s2, limit, v0, v1)

proc extractKmers*(s: string, kmers: var seq[uint16]) =
  ## Extracts all overlapping 8-mer (octamer) words from sequence ``s``.
  ##
  ## Each 8-mer is encoded as a ``uint16`` using 2 bits per base
  ## (A=0, C=1, G=2, T/U=3).  Ambiguous bases (N, etc.) reset the current
  ## k-mer accumulator.  The ``kmers`` sequence is cleared before filling.
  ## Does nothing if ``s`` is shorter than 8 bases.
  kmers.setLen(0)
  if s.len >= 8:
    var k: uint16 = 0
    for j in 0 ..< 7:
      k = (k shl 2)
      case s[j]:
        of 'C', 'c': k = k or 1
        of 'G', 'g': k = k or 2
        of 'T', 't', 'U', 'u': k = k or 3
        else: discard
    for j in 7 ..< s.len:
      k = (k shl 2)
      case s[j]:
        of 'C', 'c': k = k or 1
        of 'G', 'g': k = k or 2
        of 'T', 't', 'U', 'u': k = k or 3
        else: discard
      kmers.add(k)

proc unoise*(sequences: seq[UnoiseSeq], alpha: float = 2.0,
    minsize: int = 8): seq[Centroid] =
  ## Runs the UNOISE3 clustering algorithm on a set of dereplicated sequences.
  ##
  ## Sequences are processed in descending abundance order.  Each sequence is
  ## compared against existing centroids using edit distance and the UNOISE
  ## skew formula ``skew = targetSize / querySize``.  A sequence is merged into
  ## the best matching centroid when the number of differences satisfies
  ## ``diffs ≤ (log2(skew) - 1) / alpha``; otherwise it becomes a new centroid.
  ##
  ## **Parameters**
  ## * ``sequences`` — input sequences with ``size`` annotations
  ## * ``alpha`` — skew exponent controlling stringency (default: ``2.0``)
  ## * ``minsize`` — sequences with fewer copies are skipped (default: ``8``)
  ##
  ## Returns centroids sorted by total abundance (descending).
  ## These are ZOTU candidates before chimera filtering.
  ##
  ## See also: `uchime <unoisenim/uchime2_algo.html#uchime>`_ for chimera removal.
  # 1. Sort sequences by size, descending
  var sortedSeqs = sequences
  sortedSeqs.sort(proc(a, b: UnoiseSeq): int = cmp(b.size, a.size))

  var centroids = newSeqOfCap[Centroid](1024)
  var row0 = newSeq[int](0)
  var row1 = newSeq[int](0)

  for query in sortedSeqs:
    if query.size < minsize:
      # Due to descending order, subsequent sizes will be smaller
      break

    var bestTargetIdx = -1
    var bestDiffs = int.high
    let querySizeTimes2 = query.size * 2
    let querySeq = query.seq

    for i in 0..<centroids.len:
      let targetSize = centroids[i].seqObj.size
      if targetSize < querySizeTimes2:
        # Centroids are added in descending size order, so remaining
        # targets will not satisfy the UNOISE skew threshold either.
        break

      let skew = float(targetSize) / float(query.size)

      let maxDiffFloat = (log2(skew) - 1.0) / alpha
      if maxDiffFloat < 0.0:
        break

      let maxDiff = int(maxDiffFloat)

      let targetSeq = centroids[i].seqObj.seq
      if abs(querySeq.len - targetSeq.len) > maxDiff:
        continue

      let diff = editDistanceLimitImpl(querySeq, targetSeq, maxDiff, row0, row1)
      if diff <= maxDiff:
        if diff < bestDiffs:
          bestDiffs = diff
          bestTargetIdx = i
          if bestDiffs <= 1:
            break # Good enough match

    if bestTargetIdx != -1:
      centroids[bestTargetIdx].totalSize += query.size
    else:
      centroids.add(Centroid(seqObj: query, totalSize: query.size))

  centroids.sort(proc(a, b: Centroid): int = cmp(b.totalSize, a.totalSize))
  return centroids
