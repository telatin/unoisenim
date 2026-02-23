var scoreMat = newSeq[seq[int16]](2000)
for i in 0..<2000: scoreMat[i] = newSeq[int16](2000)

import unoise_algo

type
  AlignPath* = seq[char]

proc globalAlign*(q, t: string, path: var AlignPath): int =
  let lenQ = q.len
  let lenT = t.len

  if lenQ + 2 >= scoreMat.len:
    scoreMat.setLen(lenQ + 500)
    for i in 0..<scoreMat.len:
      if scoreMat[i].len == 0:
        scoreMat[i] = newSeq[int16](lenT + 500)

  for i in 0..lenQ + 1:
    if lenT + 2 >= scoreMat[i].len:
      scoreMat[i].setLen(lenT + 500)

  let band = 16
  let INF = 32000'i16

  for i in 0..lenQ:
    let minJ = max(0, i - band - 1)
    let maxJ = min(lenT, i + band + 1)
    for j in minJ..maxJ: scoreMat[i][j] = INF

  for i in 0..min(band, lenQ): scoreMat[i][0] = int16(i)
  for j in 0..min(band, lenT): scoreMat[0][j] = int16(j)

  for i in 1..lenQ:
    let minJ = max(1, i - band)
    let maxJ = min(lenT, i + band)
    for j in minJ..maxJ:
      let cost = if q[i-1] == t[j-1]: 0'i16 else: 1'i16

      # Ensure accesses are within valid banded initialized range
      let m = scoreMat[i-1][j-1] + cost
      let d = if j >= i - band + 1: scoreMat[i-1][j] + 1'i16 else: INF
      let ins = if j - 1 >= i - band - 1: scoreMat[i][j-1] + 1'i16 else: INF

      scoreMat[i][j] = min(m, min(d, ins))

  path.setLen(0)
  var i = lenQ
  var j = lenT
  while i > 0 or j > 0:
    if i > 0 and j > 0:
      let cost = if q[i-1] == t[j-1]: 0'i16 else: 1'i16
      if scoreMat[i][j] == scoreMat[i-1][j-1] + cost:
        path.add('M')
        i -= 1
        j -= 1
        continue
    if i > 0 and scoreMat[i][j] == scoreMat[i-1][j] + 1'i16:
      path.add('D')
      i -= 1
      continue
    if j > 0 and scoreMat[i][j] == scoreMat[i][j-1] + 1'i16:
      path.add('I')
      j -= 1
      continue

  for k in 0 ..< path.len div 2:
    swap(path[k], path[path.len - 1 - k])

  return int(scoreMat[lenQ][lenT])

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


proc uchime*(centroids: seq[Centroid], minAbSkew: float = 16.0): seq[bool] =
  var isChimera = newSeq[bool](centroids.len)

  var totalAligns = 0

  for i in 0 ..< centroids.len:
    let qObj = centroids[i]
    let q = qObj.seqObj.seq
    let qSize = qObj.totalSize

    var posBestL0 = -1
    var posBestL1 = -1
    var posBestR0 = q.len + 1
    var posBestR1 = q.len + 1

    var bestL0 = -1
    var bestL1 = -1
    var bestR0 = -1
    var bestR1 = -1

    var minDiffsQT = int.high

    var path = newSeqOfCap[char](2048)

    for j in 0 ..< i:
      let tObj = centroids[j]

      if float(tObj.totalSize) < float(qSize) * minAbSkew:
        break

      let t = tObj.seqObj.seq
      totalAligns += 1
      let alnDiffs = globalAlign(q, t, path)

      if alnDiffs < minDiffsQT:
        minDiffsQT = alnDiffs

      if alnDiffs == 0:
        break # Exact match, cannot be chimera

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

    if minDiffsQT == 0:
      isChimera[i] = false
      continue

    var chimera = false
    # 0 diff crossover
    if posBestL0 > 2 and posBestR0 != q.len + 1 and posBestL0 + 1 >=
        posBestR0 and bestL0 != bestR0:
      chimera = true
    # 1 diff crossovers
    elif minDiffsQT > 4:
      if posBestL1 > 2 and posBestR0 != q.len + 1 and posBestL1 + 1 >=
          posBestR0 and bestL1 != bestR0:
        chimera = true
      elif posBestL0 > 2 and posBestR1 != q.len + 1 and posBestL0 + 1 >=
          posBestR1 and bestL0 != bestR1:
        chimera = true

    isChimera[i] = chimera

  echo "Total UCHIME global alignments: ", totalAligns
  return isChimera
