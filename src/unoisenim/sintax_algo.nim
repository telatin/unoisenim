import strutils, tables, algorithm

type
  SintaxIndex* = object
    postings*: array[65536, seq[int32]]
    taxStrings*: seq[seq[string]]
    seqToTaxIndex: seq[int32]
    uniqTaxStrings: seq[string]
    uniqTaxRanks: seq[seq[string]]

  SintaxHit* = object
    rankNames*: seq[string]
    rankProbs*: seq[float]
    strand*: char

  SintaxState* = object
    ws: SintaxWorkspace

  SintaxWorkspace = object
    queryWords: seq[uint16]
    seenWords: array[65536, int32]
    seenMark: int32
    u: seq[int32]
    modifiedTargets: seq[int]
    topTargets: seq[int]
    taxIdxToCount: Table[int32, int]
    tieRng: MwcRng
    randSeed: uint32

  MwcRng = object
    x: array[5, uint32]
    initialized: bool

proc extractTaxRanks*(tax: string): seq[string] =
  # Assuming tax string like "d:Bacteria,p:Firmicutes..."
  for p in tax.split(','):
    result.add(p)

proc baseBits(c: char): int =
  case c:
    of 'A', 'a': return 0
    of 'C', 'c': return 1
    of 'G', 'g': return 2
    of 'T', 't', 'U', 'u': return 3
    else: return -1

proc nameIsInTaxStr(taxStr: string, name: string): bool =
  let n = taxStr.find(name)
  if n < 0:
    return false
  let m = n + name.len
  if m >= taxStr.len:
    return true
  return taxStr[m] == ','

proc slcgReset(state: var uint32, seed: uint32) =
  state = seed
  for _ in 0 ..< 10:
    state = state * 214013'u32 + 2531011'u32

proc slcgNext(state: var uint32): uint32 =
  state = state * 214013'u32 + 2531011'u32
  return state

proc resetMwc(rng: var MwcRng, seed: uint32) =
  var s = seed
  slcgReset(s, seed)
  for i in 0 ..< 5:
    rng.x[i] = slcgNext(s)
  for _ in 0 ..< 100:
    let sum = 2111111111'u64 * uint64(rng.x[3]) +
        1492'u64 * uint64(rng.x[2]) +
        1776'u64 * uint64(rng.x[1]) +
        5115'u64 * uint64(rng.x[0]) +
        uint64(rng.x[4])
    rng.x[3] = rng.x[2]
    rng.x[2] = rng.x[1]
    rng.x[1] = rng.x[0]
    rng.x[4] = uint32(sum shr 32)
    rng.x[0] = uint32(sum)
  rng.initialized = true

proc nextMwc(rng: var MwcRng): uint32 =
  if not rng.initialized:
    resetMwc(rng, 1'u32)
  let sum = 2111111111'u64 * uint64(rng.x[3]) +
      1492'u64 * uint64(rng.x[2]) +
      1776'u64 * uint64(rng.x[1]) +
      5115'u64 * uint64(rng.x[0]) +
      uint64(rng.x[4])
  rng.x[3] = rng.x[2]
  rng.x[2] = rng.x[1]
  rng.x[1] = rng.x[0]
  rng.x[4] = uint32(sum shr 32)
  rng.x[0] = uint32(sum)
  return rng.x[0]

proc buildIndex*(seqs: seq[string], taxStrings: seq[string]): SintaxIndex =
  var idx = SintaxIndex()
  idx.taxStrings = newSeq[seq[string]](taxStrings.len)
  idx.seqToTaxIndex = newSeq[int32](taxStrings.len)
  var taxToIdx = initTable[string, int32]()
  for i, ts in taxStrings:
    idx.taxStrings[i] = extractTaxRanks(ts)
    if not taxToIdx.hasKey(ts):
      let ti = int32(idx.uniqTaxStrings.len)
      taxToIdx[ts] = ti
      idx.uniqTaxStrings.add(ts)
      idx.uniqTaxRanks.add(idx.taxStrings[i])
    idx.seqToTaxIndex[i] = taxToIdx[ts]

  var lastSeen: array[65536, int32]
  for k in 0 ..< 65536: lastSeen[k] = -1

  for i, s in seqs:
    let i32 = int32(i)
    if s.len >= 8:
      var k: uint16 = 0
      var valid = 0
      for j in 0 ..< s.len:
        let b = baseBits(s[j])
        if b < 0:
          k = 0
          valid = 0
          continue

        k = (k shl 2) or uint16(b)
        if valid < 7:
          inc(valid)
          continue

        if lastSeen[k] != i32:
          idx.postings[k].add(i32)
          lastSeen[k] = i32

  return idx

# LCG Random as in USEARCH
proc nextRand(seed: var uint32): uint32 =
  seed = 1664525'u32 * seed + 1013904223'u32
  return seed

proc initWorkspace(dbSize: int): SintaxWorkspace =
  result.queryWords = newSeqOfCap[uint16](256)
  result.seenMark = 1
  result.u = newSeq[int32](dbSize)
  result.modifiedTargets = newSeqOfCap[int](1024)
  result.topTargets = newSeqOfCap[int](dbSize)
  result.taxIdxToCount = initTable[int32, int]()
  result.taxIdxToCount.clear()
  result.randSeed = 1'u32
  resetMwc(result.tieRng, result.randSeed)

proc nextSeenMark(ws: var SintaxWorkspace): int32 =
  inc ws.seenMark
  if ws.seenMark <= 0:
    for i in 0 ..< 65536:
      ws.seenWords[i] = 0
    ws.seenMark = 1
  return ws.seenMark

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

proc classifyOneDirImpl(query: string, idx: SintaxIndex, ws: var SintaxWorkspace,
    bootSubset: int = 32, bootIters: int = 100): tuple[topCount: int, hit: SintaxHit] =
  ws.queryWords.setLen(0)
  let seenMark = nextSeenMark(ws)

  if query.len >= 8:
    var k: uint16 = 0
    var valid = 0
    for j in 0 ..< query.len:
      let b = baseBits(query[j])
      if b < 0:
        k = 0
        valid = 0
        continue

      k = (k shl 2) or uint16(b)
      if valid < 7:
        inc(valid)
        continue

      if ws.seenWords[k] != seenMark:
        ws.queryWords.add(k)
        ws.seenWords[k] = seenMark

  if ws.queryWords.len < 8:
    return (0, SintaxHit())

  ws.taxIdxToCount.clear()
  var sampleSeed = ws.randSeed

  var topTotalWordCount = 0

  for _ in 0 ..< bootIters:
    # SetUShuffle
    for t in ws.modifiedTargets:
      ws.u[t] = 0
    ws.modifiedTargets.setLen(0)

    let M = bootSubset

    for _ in 0 ..< M:
      let ri = int(nextRand(sampleSeed) mod uint32(ws.queryWords.len))
      let word = ws.queryWords[ri]
      for targetIdx in idx.postings[word]:
        if ws.u[targetIdx] == 0:
          ws.modifiedTargets.add(targetIdx)
        ws.u[targetIdx] += 1

    var topU = 0'i32
    ws.topTargets.setLen(0)

    if ws.modifiedTargets.len == 0:
      let topTargetIndex = int(nextMwc(ws.tieRng) mod uint32(ws.u.len))
      let taxIdx = idx.seqToTaxIndex[topTargetIndex]
      ws.taxIdxToCount[taxIdx] = ws.taxIdxToCount.getOrDefault(taxIdx, 0) + 1
    else:
      for targetIdx in ws.modifiedTargets:
        let v = ws.u[targetIdx]
        if v > topU:
          topU = v
          ws.topTargets.setLen(0)
          ws.topTargets.add(targetIdx)
        elif v == topU:
          ws.topTargets.add(targetIdx)

      # USEARCH scans targets in increasing index order before tie-break.
      if ws.topTargets.len > 1:
        ws.topTargets.sort(system.cmp[int])

      let topTargetIndex = ws.topTargets[int(nextMwc(ws.tieRng) mod uint32(ws.topTargets.len))]
      if int(topU) > topTotalWordCount:
        topTotalWordCount = int(topU)

      let taxIdx = idx.seqToTaxIndex[topTargetIndex]
      ws.taxIdxToCount[taxIdx] = ws.taxIdxToCount.getOrDefault(taxIdx, 0) + 1

  var taxCounts = newSeq[tuple[taxIdx: int32, count: int]]()
  taxCounts.setLen(ws.taxIdxToCount.len)
  var p = 0
  for k, v in ws.taxIdxToCount:
    taxCounts[p] = (k, v)
    inc p
  if taxCounts.len == 0:
    return (0, SintaxHit())

  taxCounts.sort(proc(a, b: tuple[taxIdx: int32, count: int]): int =
    if a.count != b.count:
      return cmp(b.count, a.count)
    return cmp(idx.uniqTaxStrings[a.taxIdx], idx.uniqTaxStrings[b.taxIdx]))

  let topTaxIdx = taxCounts[0].taxIdx
  let topCount = taxCounts[0].count
  let topTaxRanks = idx.uniqTaxRanks[topTaxIdx]

  var hit = SintaxHit()
  hit.rankNames = topTaxRanks
  hit.rankProbs = newSeq[float](topTaxRanks.len)

  var prodP = 1.0
  for depth in 0 ..< topTaxRanks.len:
    let predName = topTaxRanks[depth]
    var predNameCount = topCount
    for j in 1 ..< taxCounts.len:
      let taxStr = idx.uniqTaxStrings[taxCounts[j].taxIdx]
      if nameIsInTaxStr(taxStr, predName):
        predNameCount += taxCounts[j].count

    var p = float(predNameCount) / float(bootIters)
    prodP *= p # SINTAX multiplies probabilities as depth increases
    p = prodP
    hit.rankProbs[depth] = p

  return (topTotalWordCount, hit)

proc classifyOneDir*(query: string, idx: SintaxIndex, bootSubset: int = 32,
    bootIters: int = 100): tuple[topCount: int, hit: SintaxHit] =
  var ws = initWorkspace(idx.taxStrings.len)
  return classifyOneDirImpl(query, idx, ws, bootSubset, bootIters)

proc initSintaxState*(idx: SintaxIndex): SintaxState =
  result.ws = initWorkspace(idx.taxStrings.len)

proc sintaxWithState*(query: string, idx: SintaxIndex, state: var SintaxState,
    bootSubset: int = 32, bootIters: int = 100): SintaxHit =
  let (countFwd, hitFwd) = classifyOneDirImpl(query, idx, state.ws, bootSubset, bootIters)
  let (countRev, hitRev) = classifyOneDirImpl(rc(query), idx, state.ws, bootSubset, bootIters)

  if countFwd >= countRev:
    var res = hitFwd
    res.strand = '+'
    return res
  else:
    var res = hitRev
    res.strand = '-'
    return res

proc sintax*(query: string, idx: SintaxIndex, bootSubset: int = 32,
    bootIters: int = 100): SintaxHit =
  var state = initSintaxState(idx)
  return sintaxWithState(query, idx, state, bootSubset, bootIters)
