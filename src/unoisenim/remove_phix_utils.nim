## PhiX174 contamination detection utilities.
##
## Provides k-mer containment scoring against the PhiX174 bacteriophage
## genome, compiled into the binary at build time via ``staticRead``.
## Intended for filtering PhiX spike-in reads from Illumina sequencing data.
##
## **Algorithm:** For a query read, compute what fraction of its 8-mers appear
## in the PhiX k-mer index.  For identity threshold ``p``, the expected
## fraction of 8-mers intact is approximately ``p^8`` (e.g. 0.97^8 ≈ 0.784).

import strutils, math

proc baseBitsPhix(c: char): int {.inline.} =
  case c:
    of 'A', 'a': return 0
    of 'C', 'c': return 1
    of 'G', 'g': return 2
    of 'T', 't', 'U', 'u': return 3
    else: return -1

proc reverseComplementStatic(s: string): string {.compileTime.} =
  var res = newString(s.len)
  for i in 0 ..< s.len:
    let c = s[s.len - 1 - i]
    case c:
      of 'A': res[i] = 'T'
      of 'T', 'U': res[i] = 'A'
      of 'C': res[i] = 'G'
      of 'G': res[i] = 'C'
      of 'a': res[i] = 't'
      of 't', 'u': res[i] = 'a'
      of 'c': res[i] = 'g'
      of 'g': res[i] = 'c'
      else: res[i] = c
  return res

proc buildPhixKmers(): array[65536, bool] {.compileTime.} =
  let fasta = staticRead("../../data/phi.fasta")
  var seq = ""
  for line in fasta.splitLines():
    let l = line.strip()
    if l.len > 0 and l[0] != '>':
      seq.add(l)

  # Forward strand
  if seq.len >= 8:
    var k: uint16 = 0
    var valid = 0
    for j in 0 ..< seq.len:
      let b = baseBitsPhix(seq[j])
      if b < 0:
        k = 0
        valid = 0
        continue
      k = (k shl 2) or uint16(b)
      if valid < 7:
        inc(valid)
        continue
      result[k] = true

  # Reverse complement strand
  let rcSeq = reverseComplementStatic(seq)
  if rcSeq.len >= 8:
    var k2: uint16 = 0
    var valid2 = 0
    for j in 0 ..< rcSeq.len:
      let b = baseBitsPhix(rcSeq[j])
      if b < 0:
        k2 = 0
        valid2 = 0
        continue
      k2 = (k2 shl 2) or uint16(b)
      if valid2 < 7:
        inc(valid2)
        continue
      result[k2] = true

proc getPhixSeqLen(): int {.compileTime.} =
  let fasta = staticRead("../../data/phi.fasta")
  var seqLen = 0
  for line in fasta.splitLines():
    let l = line.strip()
    if l.len > 0 and l[0] != '>':
      seqLen += l.len
  return seqLen

const phixKmers* = buildPhixKmers()
  ## Compile-time array mapping every 16-bit 8-mer hash to ``true`` if that
  ## 8-mer occurs in the PhiX174 genome (forward or reverse complement).

const phixSeqLen* = getPhixSeqLen()
  ## Length of the PhiX174 reference sequence in bases.

proc phixScore*(query: string): float =
  ## Returns the fraction of valid 8-mers in ``query`` found in the PhiX174
  ## k-mer index.
  ##
  ## A score near 1.0 indicates the read is likely PhiX contamination.
  ## Uses the same 2-bit base encoding as ``sintax_algo`` (A=0, C=1, G=2, T/U=3).
  ## Ambiguous bases (N etc.) reset the k-mer accumulator.
  ##
  ## Returns 0.0 for sequences shorter than 8 bases.
  if query.len < 8:
    return 0.0
  var total = 0
  var found = 0
  var k: uint16 = 0
  var valid = 0
  for j in 0 ..< query.len:
    let b = baseBitsPhix(query[j])
    if b < 0:
      k = 0
      valid = 0
      continue
    k = (k shl 2) or uint16(b)
    if valid < 7:
      inc(valid)
      continue
    inc(total)
    if phixKmers[k]:
      inc(found)
  if total == 0:
    return 0.0
  return float(found) / float(total)

proc isPhix*(query: string, minId: float = 0.97, minKmers: int = 8): bool =
  ## Returns ``true`` if ``query`` is identified as PhiX174 contamination.
  ##
  ## Classification is based on 8-mer containment: a read is called PhiX if
  ## the fraction of its k-mers present in ``phixKmers`` is at least
  ## ``minId^8`` (mirroring the USEARCH ``-id`` threshold convention for
  ## 8-mer pre-filtering).  Reads with fewer than ``minKmers`` valid k-mers
  ## are always kept (returned ``false``).
  ##
  ## **Parameters**
  ## * ``query``    — DNA/RNA query sequence
  ## * ``minId``    — identity threshold (default ``0.97``, matching usearch ``-closed_ref``)
  ## * ``minKmers`` — minimum valid 8-mers required to make a call (default ``8``)
  if query.len < 8:
    return false
  var total = 0
  var found = 0
  var k: uint16 = 0
  var valid = 0
  for j in 0 ..< query.len:
    let b = baseBitsPhix(query[j])
    if b < 0:
      k = 0
      valid = 0
      continue
    k = (k shl 2) or uint16(b)
    if valid < 7:
      inc(valid)
      continue
    inc(total)
    if phixKmers[k]:
      inc(found)
  if total < minKmers:
    return false
  return float(found) / float(total) >= pow(minId, 8.0)
