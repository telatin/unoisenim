## SINTAX rapid taxonomic classifier for amplicon sequences.
##
## Implements the SINTAX algorithm which uses 8-mer word matching and
## bootstrap resampling to classify sequences against a reference database
## without requiring a training step.
##
## **Reference:** Edgar RC (2016). SINTAX: a simple non-Bayesian taxonomy
## classifier for 16S and ITS sequences. *bioRxiv* doi:10.1101/074161

import strutils, tables, algorithm

type
  SintaxIndex* = object
    ## Prebuilt k-mer posting-list index for SINTAX classification.
    ##
    ## Constructed from a reference database by `buildIndex`.
    ## ``taxStrings`` holds per-sequence rank arrays; the posting lists
    ## allow O(1) lookup of which sequences contain a given 8-mer.
    postingStarts: array[65536, int32]
    postingLens: array[65536, int32]
    postingData: seq[int32]
    taxStrings*: seq[seq[string]]  ## Per-sequence parsed taxonomy ranks
    seqToTaxIndex: seq[int32]
    uniqTaxStrings: seq[string]
    uniqTaxRanks: seq[seq[string]]
    uniqTaxRankIds: seq[seq[int32]]

  SintaxHit* = object
    ## Result of a SINTAX taxonomic classification.
    rankNames*: seq[string]  ## Taxonomy rank strings (e.g. ``"d:Bacteria"``)
    rankProbs*: seq[float]   ## Bootstrap confidence for each rank (0.0–1.0)
    strand*: char            ## Matched strand: ``'+'`` (forward) or ``'-'`` (reverse complement)

  SintaxState* = object
    ## Reusable per-thread workspace for batch SINTAX classification.
    ##
    ## Initialise once with `initSintaxState` and pass to `sintaxWithState`
    ## for each query to avoid repeated heap allocation.
    ws: SintaxWorkspace

  SintaxWorkspace = object
    queryWords: seq[uint16]
    seenWords: array[65536, int32]
    seenMark: int32
    u: seq[int32]
    modifiedTargets: seq[int]
    topTargets: seq[int]
    taxVoteIdxs: seq[int32]
    taxVoteCounts: seq[int]
    tieRng: MwcRng
    randSeed: uint32

  MwcRng = object
    x: array[5, uint32]
    initialized: bool

proc postingLen*(idx: SintaxIndex, word: uint16): int =
  ## Returns the number of database sequences that contain 8-mer ``word``.
  return int(idx.postingLens[word])

proc postingFirst*(idx: SintaxIndex, word: uint16): int32 =
  ## Returns the index of the first database sequence containing ``word``,
  ## or ``-1`` if the posting list is empty.
  if idx.postingLens[word] == 0:
    return -1
  return idx.postingData[idx.postingStarts[word]]

proc extractTaxRanks*(tax: string): seq[string] =
  ## Parses a comma-delimited taxonomy string into a list of rank strings.
  ##
  ## Expects a string of the form ``"d:Bacteria,p:Firmicutes,c:Bacilli"``
  ## and returns each ``rank:name`` token as a separate element.
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
  ## Builds a k-mer posting-list index from reference sequences and taxonomy strings.
  ##
  ## Each sequence in ``seqs`` is paired with the corresponding entry in
  ## ``taxStrings``.  All overlapping 8-mers in each sequence are indexed
  ## into packed posting lists for fast word-match lookups during classification.
  ##
  ## **Parameters**
  ## * ``seqs``       — reference DNA sequences
  ## * ``taxStrings`` — taxonomy strings matching each sequence (e.g. ``"d:Bacteria,p:Firmicutes"``)
  ##
  ## Returns a `SintaxIndex` ready for use with `sintax` or `sintaxWithState`.
  var idx = SintaxIndex()
  idx.taxStrings = newSeq[seq[string]](taxStrings.len)
  idx.seqToTaxIndex = newSeq[int32](taxStrings.len)
  var taxToIdx = initTable[string, int32]()
  var rankToId = initTable[string, int32]()
  for i, ts in taxStrings:
    idx.taxStrings[i] = extractTaxRanks(ts)
    if not taxToIdx.hasKey(ts):
      let ti = int32(idx.uniqTaxStrings.len)
      taxToIdx[ts] = ti
      idx.uniqTaxStrings.add(ts)
      idx.uniqTaxRanks.add(idx.taxStrings[i])
      var rankIds = newSeq[int32](idx.taxStrings[i].len)
      for r, rankName in idx.taxStrings[i]:
        if not rankToId.hasKey(rankName):
          rankToId[rankName] = int32(rankToId.len)
        rankIds[r] = rankToId[rankName]
      idx.uniqTaxRankIds.add(rankIds)
    idx.seqToTaxIndex[i] = taxToIdx[ts]

  var counts: array[65536, int32]
  var lastSeen: array[65536, int32]
  for k in 0 ..< 65536:
    counts[k] = 0
    lastSeen[k] = -1

  # Pass 1: count postings per word.
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
          inc counts[k]
          lastSeen[k] = i32

  var total = 0
  for w in 0 ..< 65536:
    idx.postingStarts[w] = int32(total)
    idx.postingLens[w] = counts[w]
    total += int(counts[w])
  idx.postingData = newSeq[int32](total)

  var writePos: array[65536, int32]
  for w in 0 ..< 65536:
    writePos[w] = idx.postingStarts[w]
    lastSeen[w] = -1

  # Pass 2: fill packed postings.
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
          idx.postingData[writePos[k]] = i32
          inc writePos[k]
          lastSeen[k] = i32

  return idx

# LCG Random as in USEARCH
proc nextRand(seed: var uint32): uint32 =
  seed = 1664525'u32 * seed + 1013904223'u32
  return seed

proc kthSmallestInPlace(values: var seq[int], k: int): int =
  var left = 0
  var right = values.len - 1
  var kk = k
  while true:
    if left == right:
      return values[left]

    let pivot = values[(left + right) shr 1]
    var i = left
    var j = right
    while i <= j:
      while values[i] < pivot:
        inc i
      while values[j] > pivot:
        dec j
      if i <= j:
        if i != j:
          swap(values[i], values[j])
        inc i
        dec j

    if kk <= j:
      right = j
    elif kk >= i:
      left = i
    else:
      return values[kk]

proc initWorkspace(dbSize: int): SintaxWorkspace =
  result.queryWords = newSeqOfCap[uint16](256)
  result.seenMark = 1
  result.u = newSeq[int32](dbSize)
  result.modifiedTargets = newSeqOfCap[int](1024)
  result.topTargets = newSeqOfCap[int](dbSize)
  result.taxVoteIdxs = newSeqOfCap[int32](128)
  result.taxVoteCounts = newSeqOfCap[int](128)
  result.randSeed = 1'u32
  resetMwc(result.tieRng, result.randSeed)

proc nextSeenMark(ws: var SintaxWorkspace): int32 =
  inc ws.seenMark
  if ws.seenMark <= 0:
    for i in 0 ..< 65536:
      ws.seenWords[i] = 0
    ws.seenMark = 1
  return ws.seenMark

proc incTaxVote(ws: var SintaxWorkspace, taxIdx: int32) =
  for i in 0 ..< ws.taxVoteIdxs.len:
    if ws.taxVoteIdxs[i] == taxIdx:
      inc ws.taxVoteCounts[i]
      return
  ws.taxVoteIdxs.add(taxIdx)
  ws.taxVoteCounts.add(1)

proc rc*(s: string): string =
  ## Returns the reverse complement of DNA/RNA sequence ``s``.
  ##
  ## Both upper- and lower-case bases are handled; non-ACGTU characters are
  ## passed through unchanged.
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

proc fillQueryWords(query: string, ws: var SintaxWorkspace, reverseComp: bool) =
  ws.queryWords.setLen(0)
  let seenMark = nextSeenMark(ws)

  if query.len >= 8:
    var k: uint16 = 0
    var valid = 0
    for j in 0 ..< query.len:
      let srcIdx = (if reverseComp: query.len - 1 - j else: j)
      var b = baseBits(query[srcIdx])
      if reverseComp and b >= 0:
        # A<->T and C<->G.
        b = b xor 3
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

proc classifyOneDirImpl(query: string, idx: SintaxIndex, ws: var SintaxWorkspace,
    reverseComp: bool = false, bootSubset: int = 32,
    bootIters: int = 100): tuple[topCount: int, hit: SintaxHit] =
  fillQueryWords(query, ws, reverseComp)

  if ws.queryWords.len < 8:
    return (0, SintaxHit())

  ws.taxVoteIdxs.setLen(0)
  ws.taxVoteCounts.setLen(0)
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
      let start = idx.postingStarts[word]
      let n = idx.postingLens[word]
      var pos = start
      let endPos = start + n
      while pos < endPos:
        let t = idx.postingData[pos]
        let targetIdx = int(t)
        if ws.u[targetIdx] == 0:
          ws.modifiedTargets.add(targetIdx)
        ws.u[targetIdx] += 1
        inc pos

    var topU = 0'i32
    ws.topTargets.setLen(0)

    if ws.modifiedTargets.len == 0:
      let topTargetIndex = int(nextMwc(ws.tieRng) mod uint32(ws.u.len))
      let taxIdx = idx.seqToTaxIndex[topTargetIndex]
      incTaxVote(ws, taxIdx)
    else:
      for targetIdx in ws.modifiedTargets:
        let v = ws.u[targetIdx]
        if v > topU:
          topU = v
          ws.topTargets.setLen(0)
          ws.topTargets.add(targetIdx)
        elif v == topU:
          ws.topTargets.add(targetIdx)

      # USEARCH tie-breaks on index-ordered top targets.
      let r = int(nextMwc(ws.tieRng) mod uint32(ws.topTargets.len))
      let topTargetIndex = kthSmallestInPlace(ws.topTargets, r)
      if int(topU) > topTotalWordCount:
        topTotalWordCount = int(topU)

      let taxIdx = idx.seqToTaxIndex[topTargetIndex]
      incTaxVote(ws, taxIdx)

  var taxCounts = newSeq[tuple[taxIdx: int32, count: int]]()
  taxCounts.setLen(ws.taxVoteIdxs.len)
  for i in 0 ..< ws.taxVoteIdxs.len:
    taxCounts[i] = (ws.taxVoteIdxs[i], ws.taxVoteCounts[i])
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
  let topTaxRankIds = idx.uniqTaxRankIds[topTaxIdx]
  for depth in 0 ..< topTaxRanks.len:
    let predRankId = topTaxRankIds[depth]
    var predNameCount = topCount
    for j in 1 ..< taxCounts.len:
      let otherRankIds = idx.uniqTaxRankIds[taxCounts[j].taxIdx]
      for rid in otherRankIds:
        if rid == predRankId:
          predNameCount += taxCounts[j].count
          break

    var p = float(predNameCount) / float(bootIters)
    prodP *= p # SINTAX multiplies probabilities as depth increases
    p = prodP
    hit.rankProbs[depth] = p

  return (topTotalWordCount, hit)

proc classifyOneDir*(query: string, idx: SintaxIndex, bootSubset: int = 32,
    bootIters: int = 100): tuple[topCount: int, hit: SintaxHit] =
  ## Classifies ``query`` against ``idx`` on the forward strand only.
  ##
  ## Allocates a fresh workspace on each call; prefer `sintaxWithState` for
  ## batch processing to reuse allocations.
  ##
  ## **Parameters**
  ## * ``query``      — DNA query sequence
  ## * ``idx``        — prebuilt SINTAX index from `buildIndex`
  ## * ``bootSubset`` — number of words sampled per bootstrap iteration (default: ``32``)
  ## * ``bootIters``  — number of bootstrap iterations (default: ``100``)
  ##
  ## Returns a tuple of ``(topCount, hit)`` where ``topCount`` is the maximum
  ## word-match count seen and ``hit`` contains rank names and bootstrap confidences.
  var ws = initWorkspace(idx.taxStrings.len)
  return classifyOneDirImpl(query, idx, ws, false, bootSubset, bootIters)

proc initSintaxState*(idx: SintaxIndex): SintaxState =
  ## Allocates and returns a reusable `SintaxState` workspace for the given index.
  ##
  ## Pass the returned state to `sintaxWithState` for each query to avoid
  ## repeated heap allocations in batch classification loops.
  result.ws = initWorkspace(idx.taxStrings.len)

proc sintaxWithState*(query: string, idx: SintaxIndex, state: var SintaxState,
    bootSubset: int = 32, bootIters: int = 100): SintaxHit =
  ## Classifies ``query`` using a pre-allocated `SintaxState` workspace.
  ##
  ## Both strands are evaluated; the strand with the higher top word-match
  ## count wins.  Use this proc in batch loops for best performance.
  ##
  ## **Parameters**
  ## * ``query``      — DNA query sequence
  ## * ``idx``        — prebuilt SINTAX index from `buildIndex`
  ## * ``state``      — reusable workspace from `initSintaxState`
  ## * ``bootSubset`` — number of words sampled per bootstrap iteration (default: ``32``)
  ## * ``bootIters``  — number of bootstrap iterations (default: ``100``)
  let (countFwd, hitFwd) = classifyOneDirImpl(query, idx, state.ws, false,
      bootSubset, bootIters)
  let (countRev, hitRev) = classifyOneDirImpl(query, idx, state.ws, true,
      bootSubset, bootIters)

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
  ## Convenience single-query SINTAX classification.
  ##
  ## Allocates a fresh workspace, classifies ``query`` on both strands, and
  ## returns the best hit.  For classifying many queries, prefer `initSintaxState`
  ## + `sintaxWithState` to reuse workspace allocations.
  ##
  ## **Parameters**
  ## * ``query``      — DNA query sequence
  ## * ``idx``        — prebuilt SINTAX index from `buildIndex`
  ## * ``bootSubset`` — number of words sampled per bootstrap iteration (default: ``32``)
  ## * ``bootIters``  — number of bootstrap iterations (default: ``100``)
  var state = initSintaxState(idx)
  return sintaxWithState(query, idx, state, bootSubset, bootIters)
