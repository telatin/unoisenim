import unittest
import unoisenim/utils
import unoisenim/unoise_algo
import unoisenim/uchime2_algo
import unoisenim/sintax_algo

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
  test "exact low-abundance match is not called chimera":
    let centroids = @[
      Centroid(seqObj: UnoiseSeq(id: "parent", seq: "ACGTACGTAA", size: 200),
        totalSize: 200),
      Centroid(seqObj: UnoiseSeq(id: "query", seq: "ACGTACGTAA", size: 10),
        totalSize: 10)
    ]
    let flags = uchime(centroids, minAbSkew = 16.0)
    check flags.len == 2
    check flags[0] == false
    check flags[1] == false

suite "sintax":
  test "extractTaxRanks splits comma-delimited taxonomy":
    let ranks = extractTaxRanks("d:Bacteria,p:Firmicutes,g:Testus")
    check ranks == @["d:Bacteria", "p:Firmicutes", "g:Testus"]

  test "rc computes reverse complement preserving case":
    check rc("ACGTNacgtu") == "aacgtNACGT"

  test "buildIndex stores unique kmer postings per target":
    let idx = buildIndex(@["AAAAAAAAAA"], @["d:Bacteria"])
    check idx.postings[0].len == 1
    check idx.postings[0][0] == 0

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
