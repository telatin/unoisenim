const
  AlignBand = 16
  AlignInf = 32000'i16

var traceMat {.threadvar.}: seq[seq[uint8]]
var scoreRowA {.threadvar.}: seq[int16]
var scoreRowB {.threadvar.}: seq[int16]

import math
import std/threadpool
import unoise_algo

type
  AlignPath* = seq[char]

proc ensureAlignCapacity(lenQ, lenT: int) =
  if traceMat.len == 0:
    traceMat = newSeq[seq[uint8]](2000)
    for i in 0..<traceMat.len:
      traceMat[i] = newSeq[uint8](2000)
    scoreRowA = newSeq[int16](2000)
    scoreRowB = newSeq[int16](2000)

  if lenQ + 2 >= traceMat.len:
    let oldLen = traceMat.len
    traceMat.setLen(lenQ + 500)
    for i in oldLen..<traceMat.len:
      traceMat[i] = newSeq[uint8](scoreRowA.len)

  if lenT + 2 >= scoreRowA.len:
    let newLen = lenT + 500
    scoreRowA.setLen(newLen)
    scoreRowB.setLen(newLen)
    for i in 0..<traceMat.len:
      traceMat[i].setLen(newLen)

proc globalAlign*(q, t: string, path: var AlignPath): int =
  let lenQ = q.len
  let lenT = t.len
  ensureAlignCapacity(lenQ, lenT)

  # Same band used by original implementation.
  if abs(lenQ - lenT) > AlignBand:
    path.setLen(0)
    return int(AlignInf)

  var prev = scoreRowA
  var curr = scoreRowB

  for j in 0..lenT:
    if j <= AlignBand:
      prev[j] = int16(j)
      traceMat[0][j] = 2'u8
    else:
      prev[j] = AlignInf
      traceMat[0][j] = 255'u8

  traceMat[0][0] = 0'u8

  for i in 1..lenQ:
    let minJ = max(1, i - AlignBand)
    let maxJ = min(lenT, i + AlignBand)
    if minJ > maxJ:
      path.setLen(0)
      return int(AlignInf)

    curr[0] = if i <= AlignBand: int16(i) else: AlignInf
    traceMat[i][0] = 1'u8
    if minJ > 1:
      curr[minJ - 1] = AlignInf

    for j in minJ..maxJ:
      let cost = if q[i - 1] == t[j - 1]: 0'i16 else: 1'i16

      # Tie-breaking order mirrors original traceback logic:
      # diagonal > deletion > insertion.
      var best = prev[j - 1] + cost
      var dir = 0'u8

      let del = prev[j] + 1'i16
      if del < best:
        best = del
        dir = 1'u8

      let ins = curr[j - 1] + 1'i16
      if ins < best:
        best = ins
        dir = 2'u8

      curr[j] = best
      traceMat[i][j] = dir

    if maxJ < lenT:
      curr[maxJ + 1] = AlignInf

    swap(prev, curr)

  let score = prev[lenT]
  if score >= AlignInf:
    path.setLen(0)
    return int(score)

  path.setLen(0)
  var i = lenQ
  var j = lenT
  while i > 0 or j > 0:
    if i == 0:
      path.add('I')
      j -= 1
      continue
    if j == 0:
      path.add('D')
      i -= 1
      continue

    case traceMat[i][j]:
      of 0'u8:
        path.add('M')
        i -= 1
        j -= 1
      of 1'u8:
        path.add('D')
        i -= 1
      of 2'u8:
        path.add('I')
        j -= 1
      else:
        path.add('M')
        i -= 1
        j -= 1

  for k in 0 ..< path.len div 2:
    swap(path[k], path[path.len - 1 - k])

  return int(score)

proc getLeftRight*(q, t: string, path: AlignPath): tuple[diffs, posL0, posL1,
    posR0, posR1: int] =
  # Mimics USEARCH GetLeftRight
  var diffs = 0
  var posL0 = -1
  var posL1 = -1
  var posR0 = -1
  var posR1 = -1

  var qPos = 0
  var tPos = 0

  for c in path:
    if c == 'M':
      if q[qPos] != t[tPos]:
        diffs += 1
      if diffs == 0: posL0 = qPos
      elif diffs == 1: posL1 = qPos
      qPos += 1
      tPos += 1
    elif c == 'D':
      diffs += 1
      if diffs == 0: posL0 = qPos
      elif diffs == 1: posL1 = qPos
      qPos += 1
    elif c == 'I':
      # USEARCH says I increases diffs only if inside "InternalColRange"
      # We'll just increase diffs for all gaps for simplicity in this minimal impl
      diffs += 1
      if diffs == 0: posL0 = qPos
      elif diffs == 1: posL1 = qPos
      tPos += 1

  var diffsR = 0
  qPos = q.len
  tPos = t.len
  for i in countdown(path.len - 1, 0):
    let c = path[i]
    if c == 'M':
      qPos -= 1
      tPos -= 1
      if q[qPos] != t[tPos]: diffsR += 1
      if diffsR == 0: posR0 = qPos
      elif diffsR == 1: posR1 = qPos
    elif c == 'D':
      qPos -= 1
      diffsR += 1
      if diffsR == 0: posR0 = qPos
      elif diffsR == 1: posR1 = qPos
    elif c == 'I':
      tPos -= 1
      diffsR += 1
      if diffsR == 0: posR0 = qPos
      elif diffsR == 1: posR1 = qPos

  return (diffs, posL0, posL1, posR0, posR1)


proc uchime*(centroids: seq[Centroid], minAbSkew: float = 16.0,
    threads: int = 1): seq[bool] =
  type
    UchimeResult = object
      chimera: bool
      aligns: int

  proc classifyOne(i: int, c: ptr UncheckedArray[Centroid],
      parentFlags: ptr UncheckedArray[bool], minSkew: float): UchimeResult {.gcsafe.} =
    let qObj = c[i]
    let q = qObj.seqObj.seq
    let qSize = qObj.totalSize
    let minParentSize = int(ceil(float(qSize) * minSkew))

    var posBestL0 = -1
    var posBestL1 = -1
    var posBestR0 = q.len + 1
    var posBestR1 = q.len + 1

    var bestL0 = -1
    var bestL1 = -1
    var bestR0 = -1
    var bestR1 = -1

    var minDiffsQT = int.high
    var totalAligns = 0
    var path = newSeqOfCap[char](2048)

    for j in 0 ..< i:
      let tObj = c[j]
      if tObj.totalSize < minParentSize:
        break
      if parentFlags != nil and parentFlags[j]:
        continue

      let t = tObj.seqObj.seq
      inc totalAligns
      let alnDiffs = globalAlign(q, t, path)

      if alnDiffs < minDiffsQT:
        minDiffsQT = alnDiffs

      if alnDiffs == 0:
        return UchimeResult(chimera: false, aligns: totalAligns)

      let lr = getLeftRight(q, t, path)

      if lr.posL0 != -1 and lr.posL0 > posBestL0:
        posBestL0 = lr.posL0
        bestL0 = j
      if lr.posL1 != -1 and lr.posL1 > posBestL1:
        posBestL1 = lr.posL1
        bestL1 = j
      if lr.posR0 != -1 and lr.posR0 < posBestR0:
        posBestR0 = lr.posR0
        bestR0 = j
      if lr.posR1 != -1 and lr.posR1 < posBestR1:
        posBestR1 = lr.posR1
        bestR1 = j

    var chimera = false
    if posBestL0 > 2 and posBestR0 != q.len + 1 and posBestL0 + 1 >=
        posBestR0 and bestL0 != bestR0:
      chimera = true
    elif minDiffsQT > 4:
      if posBestL1 > 2 and posBestR0 != q.len + 1 and posBestL1 + 1 >=
          posBestR0 and bestL1 != bestR0:
        chimera = true
      elif posBestL0 > 2 and posBestR1 != q.len + 1 and posBestL0 + 1 >=
          posBestR1 and bestL0 != bestR1:
        chimera = true

    return UchimeResult(chimera: chimera, aligns: totalAligns)

  proc classifyRange(startIdx, count: int, c: ptr UncheckedArray[Centroid],
      parentFlags: ptr UncheckedArray[bool], minSkew: float,
      outRes: ptr UncheckedArray[UchimeResult]) {.gcsafe.} =
    for k in 0 ..< count:
      let i = startIdx + k
      outRes[i] = classifyOne(i, c, parentFlags, minSkew)

  var isChimera = newSeq[bool](centroids.len)
  if centroids.len == 0:
    return isChimera

  var totalAligns = 0

  let cPtr = cast[ptr UncheckedArray[Centroid]](unsafeAddr centroids[0])
  let flagsPtr = cast[ptr UncheckedArray[bool]](addr isChimera[0]) # sequential mode

  # threads semantics:
  # 1 => sequential (parent-filtered, closest to USEARCH de novo behavior)
  # 0 => threaded auto pool size (speed mode, independent-query evaluation)
  # >1 => threaded with fixed pool size (speed mode)
  let threaded = threads != 1

  if threaded and centroids.len >= 64:
    if threads > 1:
      let t = max(1, min(threads, MaxThreadPoolSize))
      setMaxPoolSize(t)
      setMinPoolSize(t)

    const chunkSize = 32
    var outRes = newSeq[UchimeResult](centroids.len)
    let outPtr = cast[ptr UncheckedArray[UchimeResult]](addr outRes[0])
    var start = 0
    while start < centroids.len:
      let thisStart = start
      let thisCount = min(chunkSize, centroids.len - thisStart)
      # Threaded mode intentionally does not use evolving chimera flags
      # to keep queries independent and parallelizable.
      spawn classifyRange(thisStart, thisCount, cPtr, nil, minAbSkew, outPtr)
      start += thisCount
    sync()

    for i in 0 ..< centroids.len:
      isChimera[i] = outRes[i].chimera
      totalAligns += outRes[i].aligns
  else:
    for i in 0 ..< centroids.len:
      let r = classifyOne(i, cPtr, flagsPtr, minAbSkew)
      isChimera[i] = r.chimera
      totalAligns += r.aligns

  echo "Total UCHIME global alignments: ", totalAligns
  return isChimera
