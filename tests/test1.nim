import unittest
import unoisenim/utils
import unoisenim/unoise_algo
import unoisenim/uchime2_algo
import unoisenim/sintax_algo
import unoisenim/nbc_algo
import unoisenim/remove_phix_utils

suite "utils":
  test "parseSize extracts abundance from FASTA label":
    check parseSize("Zotu1;size=42;") == 42
    check parseSize("Zotu2;foo=bar;size=9;") == 9
    check parseSize("Zotu3;size=bad;") == 0
    check parseSize("Zotu4;") == 0

suite "unoise":
  test "editDistanceLimit respects threshold early stop":
    check editDistanceLimit("ACGT", "ACGT", 0) == 0
    check editDistanceLimit("ACGT", "ACGA", 0) == 1
    check editDistanceLimit("ACGT", "ACGA", 1) == 1
    check editDistanceLimit("ACGT", "ACGTA", 1) == 1
    check editDistanceLimit("AAAA", "TTTT", 2) == 3

  test "unoise clusters by abundance skew and minsize":
    let seqs = @[
      UnoiseSeq(id: "a;size=80;", seq: "ACGTACGTAA", size: 80),
      UnoiseSeq(id: "b;size=10;", seq: "ACGTACGTAT", size: 10),
      UnoiseSeq(id: "c;size=9;", seq: "TGCATGCATT", size: 9),
      UnoiseSeq(id: "d;size=7;", seq: "TGCATGCATA", size: 7)
    ]

    let centroids = unoise(seqs, alpha = 2.0, minsize = 8)
    check centroids.len == 2
    check centroids[0].seqObj.id == "a;size=80;"
    check centroids[0].totalSize == 90
    check centroids[1].seqObj.id == "c;size=9;"
    check centroids[1].totalSize == 9

suite "uchime":
  proc mkCentroid(id, seq: string, size: int): Centroid =
    Centroid(seqObj: UnoiseSeq(id: id, seq: seq, size: size), totalSize: size)

  test "empty input returns empty flags":
    check uchime(@[], minAbSkew = 16.0, threads = 1).len == 0
    check uchime(@[], minAbSkew = 16.0, threads = 0).len == 0
    check uchime(@[], minAbSkew = 16.0, threads = 4).len == 0

  test "exact low-abundance match is not called chimera":
    let centroids = @[
      mkCentroid("parent", "ACGTACGTAA", 200),
      mkCentroid("query", "ACGTACGTAA", 10)
    ]
    let flags = uchime(centroids, minAbSkew = 16.0)
    check flags.len == 2
    check flags[0] == false
    check flags[1] == false

  test "threaded options fall back to sequential for small input":
    let centroids = @[
      mkCentroid("parentA", "ACGTACGTAA", 300),
      mkCentroid("parentB", "ACGTACGTAT", 40),
      mkCentroid("query", "ACGTACGTAA", 10)
    ]
    let seqFlags = uchime(centroids, minAbSkew = 16.0, threads = 1)
    let autoFlags = uchime(centroids, minAbSkew = 16.0, threads = 0)
    let fixedFlags = uchime(centroids, minAbSkew = 16.0, threads = 8)
    check autoFlags == seqFlags
    check fixedFlags == seqFlags

  test "threaded auto and fixed modes are equivalent on larger inputs":
    const n = 96
    var centroids = newSeq[Centroid](n)
    for i in 0..<n:
      var seq = "ACGTACGTACGTACGT"
      seq[^1] = char(ord('A') + (i mod 4))
      centroids[i] = mkCentroid("c" & $i, seq, n - i + 100)

    let autoFlags = uchime(centroids, minAbSkew = 16.0, threads = 0)
    let fixedFlags = uchime(centroids, minAbSkew = 16.0, threads = 4)
    check autoFlags.len == n
    check fixedFlags.len == n
    check autoFlags == fixedFlags

  test "chimera constructed from two parents is detected":
    # Parent A covers positions 0-19, Parent B covers 20-39.
    # The chimera splices them at position 20: the left half matches Parent A
    # (posL0=19) and the right half matches Parent B (posR0=20), satisfying
    # posBestL0+1 >= posBestR0 with two different parents.
    # minParentSize = ceil(10 * 16.0) = 160; both parents qualify.
    let parentASeq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" # 40 A's
    let parentBSeq = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT" # 40 T's
    let chimeraSeq = "AAAAAAAAAAAAAAAAAAAA" & "TTTTTTTTTTTTTTTTTTTT" # 20A+20T
    let centroids = @[
      mkCentroid("parentA", parentASeq, 1000),
      mkCentroid("parentB", parentBSeq, 800),
      mkCentroid("chimera", chimeraSeq, 10),
    ]
    let flags = uchime(centroids, minAbSkew = 16.0, threads = 1)
    check flags.len == 3
    check flags[0] == false  # parentA: no higher-abundance parents to check against
    check flags[1] == false  # parentB: posL0 never set vs parentA (first char mismatches)
    check flags[2] == true   # chimera: breakpoint detected between parentA and parentB

  test "threaded mode is deterministic across repeated runs":
    const n = 96
    var centroids = newSeq[Centroid](n)
    for i in 0..<n:
      var seq = "TGCATGCATGCATGCA"
      seq[^2] = char(ord('A') + (i mod 4))
      centroids[i] = mkCentroid("d" & $i, seq, n - i + 100)

    let run1 = uchime(centroids, minAbSkew = 16.0, threads = 0)
    let run2 = uchime(centroids, minAbSkew = 16.0, threads = 0)
    check run1 == run2

suite "sintax":
  test "extractTaxRanks splits comma-delimited taxonomy":
    let ranks = extractTaxRanks("d:Bacteria,p:Firmicutes,g:Testus")
    check ranks == @["d:Bacteria", "p:Firmicutes", "g:Testus"]

  test "rc computes reverse complement preserving case":
    check rc("ACGTNacgtu") == "aacgtNACGT"

  test "buildIndex stores unique kmer postings per target":
    let idx = buildIndex(@["AAAAAAAAAA"], @["d:Bacteria"])
    check idx.postingLen(0'u16) == 1
    check idx.postingFirst(0'u16) == 0

  test "sintax classifies a known query with full confidence":
    let dbSeq = "ACGTCAGTGCATGACCTGTAAG"
    let idx = buildIndex(@[dbSeq], @["d:Bacteria,p:Firmicutes,g:Testus"])
    let hit = sintax(dbSeq, idx, bootSubset = 16, bootIters = 50)

    check hit.strand == '+'
    check hit.rankNames == @["d:Bacteria", "p:Firmicutes", "g:Testus"]
    check hit.rankProbs.len == 3
    for p in hit.rankProbs:
      check p >= 0.99

  test "sintax handles short queries as unclassified":
    let dbSeq = "ACGTCAGTGCATGACCTGTAAG"
    let idx = buildIndex(@[dbSeq], @["d:Bacteria,p:Firmicutes,g:Testus"])
    let hit = sintax("ACGTACG", idx)
    check hit.rankNames.len == 0
    check hit.rankProbs.len == 0
    check hit.strand == '+'

  test "sintax uses reverse-complement strand when better":
    let dbSeq = "ACGTCAGTGCATGACCTGTAAG"
    let idx = buildIndex(@[dbSeq], @["d:Bacteria,p:Firmicutes,g:Testus"])
    let hit = sintax(rc(dbSeq), idx, bootSubset = 16, bootIters = 50)
    check hit.strand == '-'
    check hit.rankNames == @["d:Bacteria", "p:Firmicutes", "g:Testus"]

  test "sintax ignores words spanning ambiguous letters":
    let dbSeq = "ACGTCAGTGCATGACCTGTAAGACGT"
    let idx = buildIndex(@[dbSeq], @["d:Bacteria,p:Firmicutes,g:Testus"])
    let badQuery = "ACGNACGNACGNACGNACGNACGN"
    let hit = sintax(badQuery, idx)
    check hit.rankNames.len == 0
    check hit.rankProbs.len == 0

  test "sintax groups duplicated taxonomy strings":
    let s1 = "ACGTCAGTGCATGACCTGTAAGACGT"
    let s2 = "ACGTCAGTGCATGACCTGTAAGACGA"
    let s3 = "ACGTCAGTGCATGACCTGTAAGACGC"
    let idx = buildIndex(
      @[s1, s2, s3],
      @[
        "d:Bacteria,p:Firmicutes,g:Alpha",
        "d:Bacteria,p:Firmicutes,g:Alpha",
        "d:Bacteria,p:Actinobacteria,g:Beta"
      ]
    )
    let hit = sintax(s1, idx, bootSubset = 32, bootIters = 100)
    check hit.rankNames.len == 3
    check hit.rankNames[2] == "g:Alpha"

suite "nbc":
  test "parseTaxRanks splits comma-delimited taxonomy":
    let ranks = parseTaxRanks("d:Bacteria,p:Firmicutes,g:Testus")
    check ranks == @["d:Bacteria", "p:Firmicutes", "g:Testus"]

  test "nbc classifies a known query with high confidence":
    let seqA = "ACGTCAGTGCATGACCTGTAAG"
    let seqB = "TTGGAACCGTTTCCGGAACTTG"
    let idx = buildNbcIndex(
      @[seqA, seqB],
      @["d:Bacteria,p:Firmicutes,g:Alpha", "d:Bacteria,p:Actinobacteria,g:Beta"]
    )

    let hit = nbc(seqA, idx, bootIters = 80, minWords = 5)
    check hit.strand == '+'
    check hit.rankNames == @["d:Bacteria", "p:Firmicutes", "g:Alpha"]
    check hit.rankProbs.len == 3
    check hit.rankProbs[1] >= 0.9
    check hit.rankProbs[2] >= 0.9

  test "nbc handles short queries as unclassified":
    let idx = buildNbcIndex(@["ACGTCAGTGCATGACCTGTAAG"],
      @["d:Bacteria,p:Firmicutes,g:Alpha"])
    let hit = nbc("ACGTACG", idx)
    check hit.rankNames.len == 0
    check hit.rankProbs.len == 0
    check hit.strand == '+'

  test "nbc uses reverse-complement strand when better":
    let seqA = "ACGTCAGTGCATGACCTGTAAG"
    let idx = buildNbcIndex(@[seqA], @["d:Bacteria,p:Firmicutes,g:Alpha"])
    let hit = nbc(reverseComplement(seqA), idx, bootIters = 80, minWords = 5)
    check hit.strand == '-'
    check hit.rankNames == @["d:Bacteria", "p:Firmicutes", "g:Alpha"]

suite "remove_phix":
  # First 140 bases of the PhiX174 genome (NC_001422.1)
  const phixSnip =
    "GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTT" &
    "GATAAAGCAGGAATTACTACTGCTTGTTTACGAATTAAATCGAAGTGGACTGCTGGCGGAAAATGAGAAA"

  # E. coli 16S rRNA gene opening region — no overlap with PhiX
  const rnd16s =
    "AGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCATGCTTAACACATGCAAGTCGAACGGTAACAGGA" &
    "AGAAGCTTGCTTCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGG"

  test "phixScore returns high score for PhiX-derived sequence":
    check phixScore(phixSnip) >= 0.7

  test "phixScore returns low score for random 16S":
    # ~8% overlap expected for random reads; real PhiX threshold is ~0.78 (0.97^8)
    check phixScore(rnd16s) < 0.3

  test "isPhix returns true for PhiX sequence at default threshold":
    check isPhix(phixSnip) == true

  test "isPhix returns false for 16S sequence":
    check isPhix(rnd16s) == false

  test "isPhix respects minKmers guard for short sequences":
    check isPhix("ACGT", minId = 0.97, minKmers = 8) == false

  test "phixSeqLen matches expected PhiX174 genome length":
    check phixSeqLen == 5386
