import strutils, tables

type
  SintaxIndex* = object
    postings*: array[65536, seq[int32]]
    taxStrings*: seq[seq[string]]

  SintaxHit* = object
    rankNames*: seq[string]
    rankProbs*: seq[float]
    strand*: char

proc extractTaxRanks*(tax: string): seq[string] =
  # Assuming tax string like "d:Bacteria,p:Firmicutes..."
  for p in tax.split(','):
    result.add(p)

proc buildIndex*(seqs: seq[string], taxStrings: seq[string]): SintaxIndex =
  var idx = SintaxIndex()
  idx.taxStrings = newSeq[seq[string]](taxStrings.len)
  for i, ts in taxStrings:
    idx.taxStrings[i] = extractTaxRanks(ts)

  var lastSeen: array[65536, int32]
  for k in 0 ..< 65536: lastSeen[k] = -1

  for i, s in seqs:
    let i32 = int32(i)
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

        if lastSeen[k] != i32:
          idx.postings[k].add(i32)
          lastSeen[k] = i32

  return idx

# LCG Random as in USEARCH
var r_seed = 1'u32
proc nextRand(): uint32 =
  r_seed = 1664525'u32 * r_seed + 1013904223'u32
  return r_seed

proc rc*(s: string): string =
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

proc classifyOneDir*(query: string, idx: SintaxIndex, bootSubset: int = 32,
    bootIters: int = 100): tuple[topCount: int, hit: SintaxHit] =
  var queryWords = newSeq[uint16]()
  var seen: array[65536, bool]

  if query.len >= 8:
    var k: uint16 = 0
    for j in 0 ..< 7:
      k = (k shl 2)
      case query[j]:
        of 'C', 'c': k = k or 1
        of 'G', 'g': k = k or 2
        of 'T', 't', 'U', 'u': k = k or 3
        else: discard

    for j in 7 ..< query.len:
      k = (k shl 2)
      case query[j]:
        of 'C', 'c': k = k or 1
        of 'G', 'g': k = k or 2
        of 'T', 't', 'U', 'u': k = k or 3
        else: discard

      if not seen[k]:
        queryWords.add(k)
        seen[k] = true

  if queryWords.len < 8:
    return (0, SintaxHit())

  var taxStrToCount = initTable[int, int]()
  let dbSize = idx.taxStrings.len
  var u = newSeq[int32](dbSize)
  var modifiedTargets = newSeqOfCap[int](1024)

  var topTotalWordCount = 0

  for boot in 0 ..< bootIters:
    # SetUShuffle
    for t in modifiedTargets: u[t] = 0
    modifiedTargets.setLen(0)

    let M = min(bootSubset, queryWords.len)

    for k in 0 ..< M:
      let ri = nextRand() mod uint32(queryWords.len)
      let word = queryWords[ri]
      for targetIdx in idx.postings[word]:
        if u[targetIdx] == 0: modifiedTargets.add(targetIdx)
        u[targetIdx] += 1

    var topU = 0'i32
    var topTargets: seq[int] = @[]

    for targetIdx in modifiedTargets:
      let v = u[targetIdx]
      if v > topU:
        topU = v
        topTargets = @[targetIdx]
      elif v == topU and topU > 0:
        topTargets.add(targetIdx)

    if topTargets.len > 0:
      let topTargetIndex = topTargets[nextRand() mod uint32(topTargets.len)]
      if int(topU) > topTotalWordCount:
        topTotalWordCount = int(topU)

      let taxIdx = topTargetIndex
      taxStrToCount[taxIdx] = taxStrToCount.getOrDefault(taxIdx, 0) + 1

  # Find most frequent tax string index
  var maxTaxIdx = -1
  var maxCount = 0
  for k, v in taxStrToCount:
    if v > maxCount:
      maxCount = v
      maxTaxIdx = k

  if maxTaxIdx == -1:
    return (0, SintaxHit())

  let topTaxRanks = idx.taxStrings[maxTaxIdx]

  var hit = SintaxHit()
  hit.rankNames = topTaxRanks
  hit.rankProbs = newSeq[float](topTaxRanks.len)

  var prodP = 1.0
  for depth in 0 ..< topTaxRanks.len:
    let predName = topTaxRanks[depth]
    var predNameCount = 0
    for k, v in taxStrToCount:
      let tr = idx.taxStrings[k]
      if predName in tr:
        predNameCount += v

    var p = float(predNameCount) / float(bootIters)
    prodP *= p # SINTAX multiplies probabilities as depth increases
    p = prodP
    hit.rankProbs[depth] = p

  return (topTotalWordCount, hit)

proc sintax*(query: string, idx: SintaxIndex, bootSubset: int = 32,
    bootIters: int = 100): SintaxHit =
  let (countFwd, hitFwd) = classifyOneDir(query, idx, bootSubset, bootIters)
  let (countRev, hitRev) = classifyOneDir(rc(query), idx, bootSubset, bootIters)

  if countFwd >= countRev:
    var res = hitFwd
    res.strand = '+'
    return res
  else:
    var res = hitRev
    res.strand = '-'
    return res
