import math, algorithm

type
  UnoiseSeq* = object
    id*: string
    seq*: string
    size*: int

  Centroid* = object
    seqObj*: UnoiseSeq
    totalSize*: int

proc editDistanceLimit*(s1, s2: string, limit: int): int =
  let len1 = s1.len
  let len2 = s2.len
  if abs(len1 - len2) > limit: return limit + 1
  if len1 == 0: return if len2 <= limit: len2 else: limit + 1
  if len2 == 0: return if len1 <= limit: len1 else: limit + 1

  var v0 = newSeq[int](len2 + 1)
  var v1 = newSeq[int](len2 + 1)
  for i in 0..len2: v0[i] = i
  for i in 0..<len1:
    v1[0] = i + 1
    var minVal = v1[0]
    for j in 0..<len2:
      let cost = if s1[i] == s2[j]: 0 else: 1
      v1[j + 1] = min(v1[j] + 1, min(v0[j + 1] + 1, v0[j] + cost))
      if v1[j + 1] < minVal: minVal = v1[j + 1]

    if minVal > limit: return limit + 1
    swap(v0, v1)

  if v0[len2] <= limit: return v0[len2]
  return limit + 1

proc extractKmers*(s: string, kmers: var seq[uint16]) =
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
  # 1. Sort sequences by size, descending
  var sortedSeqs = sequences
  sortedSeqs.sort(proc(a, b: UnoiseSeq): int = cmp(b.size, a.size))

  var centroids = newSeqOfCap[Centroid](1024)

  for query in sortedSeqs:
    if query.size < minsize:
      # Due to descending order, subsequent sizes will be smaller
      break

    var bestTargetIdx = -1
    var bestDiffs = int.high

    for i in 0..<centroids.len:
      let targetSize = centroids[i].seqObj.size
      let skew = float(targetSize) / float(query.size)

      let maxDiffFloat = (log2(skew) - 1.0) / alpha
      if maxDiffFloat < 0.0:
        continue

      let maxDiff = int(maxDiffFloat)

      if abs(query.seq.len - centroids[i].seqObj.seq.len) > maxDiff:
        continue

      let diff = editDistanceLimit(query.seq, centroids[i].seqObj.seq, maxDiff)
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
