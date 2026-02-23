import math, strutils, tables

const
  NbcKmerSize* = 8
  NbcWordSpace = 65536

type
  NbcNode = object
    name: string
    parent: int
    depth: int
    children: seq[int]
    seqCount: int
    wordCounts: Table[uint16, int]

  NbcIndex* = object
    nodes: seq[NbcNode]

  NbcHit* = object
    rankNames*: seq[string]
    rankProbs*: seq[float]
    strand*: char

proc parseTaxRanks*(tax: string): seq[string] =
  for p in tax.split(','):
    let v = p.strip()
    if v.len > 0:
      result.add(v)

proc reverseComplement*(s: string): string =
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

proc nextRand(seed: var uint32): uint32 =
  seed = 1664525'u32 * seed + 1013904223'u32
  return seed

proc baseBits(c: char): int =
  case c:
    of 'A', 'a': return 0
    of 'C', 'c': return 1
    of 'G', 'g': return 2
    of 'T', 't', 'U', 'u': return 3
    else: return -1

proc collectUniqueWords(s: string, seen: var array[NbcWordSpace, int32],
    mark: int32): seq[uint16] =
  if s.len < NbcKmerSize:
    return

  var kmer: uint16 = 0
  var valid = 0

  for c in s:
    let b = baseBits(c)
    if b < 0:
      kmer = 0
      valid = 0
      continue

    kmer = (kmer shl 2) or uint16(b)
    if valid < NbcKmerSize - 1:
      inc(valid)
      continue

    if seen[kmer] != mark:
      seen[kmer] = mark
      result.add(kmer)

proc buildNbcIndex*(seqs: seq[string], taxStrings: seq[string]): NbcIndex =
  var idx = NbcIndex()
  idx.nodes = @[
    NbcNode(name: "*root*", parent: -1, depth: 0, children: @[], seqCount: 0,
      wordCounts: initTable[uint16, int]())
  ]

  var nodeByParentTax = initTable[string, int]()
  var seen: array[NbcWordSpace, int32]
  var mark = 1'i32

  let n = min(seqs.len, taxStrings.len)
  for i in 0 ..< n:
    let ranks = parseTaxRanks(taxStrings[i])
    if ranks.len == 0:
      continue

    let words = collectUniqueWords(seqs[i], seen, mark)
    inc(mark)
    if mark == high(int32):
      for j in 0 ..< NbcWordSpace:
        seen[j] = 0
      mark = 1

    var path = newSeq[int]()
    var parent = 0
    for rank in ranks:
      let key = $parent & "\t" & rank
      var nodeId = nodeByParentTax.getOrDefault(key, -1)
      if nodeId < 0:
        nodeId = idx.nodes.len
        idx.nodes.add(
          NbcNode(name: rank, parent: parent, depth: idx.nodes[parent].depth + 1,
            children: @[], seqCount: 0, wordCounts: initTable[uint16, int]())
        )
        idx.nodes[parent].children.add(nodeId)
        nodeByParentTax[key] = nodeId
      parent = nodeId
      path.add(nodeId)

    for nodeId in path:
      inc idx.nodes[nodeId].seqCount
      for w in words:
        var c = idx.nodes[nodeId].wordCounts.mgetOrPut(w, 0)
        inc c

  return idx

proc choosePath(words: seq[uint16], idx: NbcIndex, randomTie: bool,
    seed: var uint32): tuple[path: seq[int], score: float] =
  var parent = 0
  var totalScore = 0.0

  while idx.nodes[parent].children.len > 0:
    let children = idx.nodes[parent].children
    var sumSeq = 0
    for child in children:
      sumSeq += idx.nodes[child].seqCount
    if sumSeq == 0:
      break

    var bestScore = NegInf
    var bestChildren = newSeq[int]()

    for child in children:
      let node = idx.nodes[child]
      let prior = ln((float(node.seqCount) + 1.0) /
          (float(sumSeq) + float(children.len)))
      let denom = float(node.seqCount) + 2.0
      var score = prior
      for w in words:
        let c = node.wordCounts.getOrDefault(w, 0)
        let p = (float(c) + 1.0) / denom
        score += ln(p)

      if score > bestScore:
        bestScore = score
        bestChildren.setLen(0)
        bestChildren.add(child)
      elif score == bestScore:
        bestChildren.add(child)

    if bestChildren.len == 0:
      break

    let chosen =
      if randomTie and bestChildren.len > 1:
        bestChildren[nextRand(seed) mod uint32(bestChildren.len)]
      else:
        bestChildren[0]

    result.path.add(chosen)
    totalScore += bestScore
    parent = chosen

  result.score = totalScore

proc sampleWords(words: seq[uint16], sampleSize: int, seed: var uint32): seq[uint16] =
  if words.len == 0 or sampleSize <= 0:
    return
  result = newSeq[uint16](sampleSize)
  for i in 0 ..< sampleSize:
    result[i] = words[nextRand(seed) mod uint32(words.len)]

proc classifyOneDir(query: string, idx: NbcIndex, bootIters: int,
    minWords: int): tuple[path: seq[int], probs: seq[float], score: float] =
  var seen: array[NbcWordSpace, int32]
  let words = collectUniqueWords(query, seen, 1)
  if words.len == 0:
    return (@[], @[], NegInf)

  var seed = 1'u32
  let det = choosePath(words, idx, randomTie = false, seed)
  if det.path.len == 0:
    return (@[], @[], det.score)

  let wordsByFrac = words.len div 8
  var sampleSize = max(minWords, wordsByFrac)
  if sampleSize < 1:
    sampleSize = 1
  if sampleSize > words.len:
    sampleSize = words.len

  var agree = newSeq[int](det.path.len)
  if bootIters > 0:
    for _ in 0 ..< bootIters:
      let sampled = sampleWords(words, sampleSize, seed)
      let bt = choosePath(sampled, idx, randomTie = true, seed)
      var allPrevMatch = true
      for d in 0 ..< det.path.len:
        if allPrevMatch and d < bt.path.len and bt.path[d] == det.path[d]:
          inc agree[d]
        else:
          allPrevMatch = false

  var probs = newSeq[float](det.path.len)
  if bootIters > 0:
    for d in 0 ..< det.path.len:
      probs[d] = float(agree[d]) / float(bootIters)

  return (det.path, probs, det.score)

proc nbc*(query: string, idx: NbcIndex, bootIters: int = 100,
    minWords: int = 5): NbcHit =
  let fwd = classifyOneDir(query, idx, bootIters, minWords)
  let rev = classifyOneDir(reverseComplement(query), idx, bootIters, minWords)

  var chosenPath: seq[int]
  var chosenProbs: seq[float]
  if rev.path.len > fwd.path.len or (rev.path.len == fwd.path.len and rev.score >
      fwd.score):
    chosenPath = rev.path
    chosenProbs = rev.probs
    result.strand = '-'
  else:
    chosenPath = fwd.path
    chosenProbs = fwd.probs
    result.strand = '+'

  for nodeId in chosenPath:
    result.rankNames.add(idx.nodes[nodeId].name)
  result.rankProbs = chosenProbs
