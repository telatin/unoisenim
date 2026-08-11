import argparse, readfx, strutils, os, times, malebolgia
import unoisenim/sintax_algo

const BatchSize = 4096
const ChunkSize = 8

proc dada2HeaderToTax(header: string): string =
  const rankPrefixes = ["d:", "p:", "c:", "o:", "f:", "g:", "s:"]
  var parts: seq[string]
  for p in header.split(';'):
    let trimmed = p.strip()
    if trimmed.len > 0:
      parts.add(trimmed)
  if parts.len == 0:
    return ""
  var taxedParts: seq[string]
  for i, part in parts:
    if i < rankPrefixes.len:
      taxedParts.add(rankPrefixes[i] & part)
    else:
      taxedParts.add("r" & $(i + 1) & ":" & part)
  return taxedParts.join(",")

proc taxDepth(tax: string): int =
  for p in tax.split(','):
    if p.strip().len > 0:
      inc result

proc reportInitialTaxonomy(n, noTax: int, sample: string, dada2Format: bool) =
  let withTax = n - noTax
  if withTax > 0:
    echo "Taxonomy recognised in ", withTax, "/", n, " initial records."
    if sample != "":
      echo "  Example: ", sample
  else:
    stderr.writeLine("Warning: no taxonomy recognised in the first ", n,
        " records. Check the database format",
        (if not dada2Format: " or try --dada2-format" else: ""), ".")

proc formatHit(hit: SintaxHit, cutoff: float,
    outTax, passedTax: var string): bool =
  ## Fills ``outTax``/``passedTax``; returns true when the domain-level
  ## prediction passes the confidence cutoff.
  if hit.rankNames.len > 0:
    var taxes = newSeqOfCap[string](hit.rankNames.len)
    for i in 0 ..< hit.rankNames.len:
      taxes.add(hit.rankNames[i] & "(" & formatFloat(hit.rankProbs[i],
          format = ffDecimal, precision = 4) & ")")
    outTax = taxes.join(",")

    var passed = newSeqOfCap[string](hit.rankNames.len)
    for i in 0 ..< hit.rankNames.len:
      if hit.rankProbs[i] >= cutoff:
        passed.add(hit.rankNames[i])
      else:
        break
    passedTax = passed.join(",")
    return hit.rankProbs[0] >= cutoff
  else:
    outTax = "*"
    passedTax = "*"
    return false

# Module-level worker: classify a contiguous slice of a query batch.
# Must be {.gcsafe.} and take no closures so malebolgia can spawn it safely.
# Each task owns its SintaxState; the index is shared read-only.
proc classifySintaxRange(seqsPtr: ptr UncheckedArray[string],
    outTaxPtr, passedTaxPtr: ptr UncheckedArray[string],
    strandPtr: ptr UncheckedArray[char],
    domPtr: ptr UncheckedArray[bool],
    start, count: int, idxPtr: ptr SintaxIndex, cutoff: float) {.gcsafe.} =
  var state = initSintaxState(idxPtr[])
  for i in start ..< start + count:
    let hit = sintaxWithState(seqsPtr[i], idxPtr[], state)
    strandPtr[i] = hit.strand
    domPtr[i] = formatHit(hit, cutoff, outTaxPtr[i], passedTaxPtr[i])

# Classify n queries; threads semantics match uchime/remove_phix:
#   1  => sequential
#   0  => auto (malebolgia scheduler, all chunks in one awaitAll)
#   >1 => fixed max concurrent chunk tasks per awaitAll round
proc classifyBatch(qseqs: var seq[string],
    outTaxArr, passedTaxArr: var seq[string],
    strandArr: var seq[char], domArr: var seq[bool],
    n: int, idxPtr: ptr SintaxIndex, cutoff: float,
    threads: int, state: var SintaxState) =
  if threads == 1 or n < ChunkSize:
    for i in 0 ..< n:
      let hit = sintaxWithState(qseqs[i], idxPtr[], state)
      strandArr[i] = hit.strand
      domArr[i] = formatHit(hit, cutoff, outTaxArr[i], passedTaxArr[i])
    return

  let sp = cast[ptr UncheckedArray[string]](addr qseqs[0])
  let op = cast[ptr UncheckedArray[string]](addr outTaxArr[0])
  let pp = cast[ptr UncheckedArray[string]](addr passedTaxArr[0])
  let stp = cast[ptr UncheckedArray[char]](addr strandArr[0])
  let dp = cast[ptr UncheckedArray[bool]](addr domArr[0])
  var start = 0

  if threads > 1:
    while start < n:
      var m = createMaster()
      m.awaitAll:
        var launched = 0
        while start < n and launched < threads:
          let c = min(ChunkSize, n - start)
          m.spawn classifySintaxRange(sp, op, pp, stp, dp, start, c,
              idxPtr, cutoff)
          start += c
          inc launched
  else: # threads == 0: auto-scheduler
    var m = createMaster()
    m.awaitAll:
      while start < n:
        let c = min(ChunkSize, n - start)
        m.spawn classifySintaxRange(sp, op, pp, stp, dp, start, c,
            idxPtr, cutoff)
        start += c

proc main() =
  let p = newParser("sintax"):
    help("USEARCH sintax algorithm implementation in Nim")
    option("-i", "--input", help = "Input query FASTA file")
    option("-d", "--database", help = "Reference database FASTA file with taxonomy annotations")
    option("-t", "--tabbedout", help = "Output tabbed file")
    option("-c", "--cutoff", default = some("0.8"),
        help = "Confidence cutoff (default 0.8)")
    option("--threads", default = some("0"),
        help = "Threads: 1=sequential, 0=auto, N=fixed max concurrent tasks (default: 0)")
    flag("--dada2-format",
        help = "Database uses DADA2 taxonomy format (header is semicolon-delimited taxonomy)")

  var opts: typeof(p.parse())
  try:
    opts = p.parse()
  except ShortCircuit as e:
    if e.flag == "argparse_help":
      echo p.help
      quit(0)
  except UsageError:
    stderr.writeLine getCurrentExceptionMsg()
    echo p.help
    quit(1)

  if opts.input == "" or opts.database == "":
    echo "Error: Need both --input and --database"
    quit(1)

  let cutoff = parseFloat(opts.cutoff)
  let dada2Format = opts.dada2_format
  let threads = parseInt(opts.threads)

  # Load DB
  var dbSeqs = newSeq[string]()
  var taxStrings = newSeq[string]()
  var noTaxCount = 0
  var totalDbSeqs = 0
  var sampleTax = ""
  var fDb = readfx.xopen[GzFile](opts.database, mode = fmRead)
  var record: FQRecord
  while readFastx(fDb, record):
    dbSeqs.add(record.sequence)
    inc totalDbSeqs
    var tax = ""
    if dada2Format:
      tax = dada2HeaderToTax(record.name)
    else:
      for p in record.name.split(';'):
        if p.startsWith("tax="):
          tax = p[4..^1]
          break
    if tax == "":
      inc noTaxCount
    elif sampleTax == "":
      sampleTax = tax
    taxStrings.add(tax)
    if totalDbSeqs == 10:
      reportInitialTaxonomy(totalDbSeqs, noTaxCount, sampleTax, dada2Format)

  if totalDbSeqs > 0 and totalDbSeqs < 10:
    reportInitialTaxonomy(totalDbSeqs, noTaxCount, sampleTax, dada2Format)

  if noTaxCount > 0:
    stderr.writeLine("Warning: ", noTaxCount, "/", totalDbSeqs,
        " database sequences had no taxonomy detected. ",
        "Ensure database is correctly formatted or use --dada2-format.")

  var speciesLevel = 0
  var genusLevel = 0
  var belowGenus = 0
  for ts in taxStrings:
    let depth = taxDepth(ts)
    if depth >= 7:
      inc speciesLevel
    elif depth == 6:
      inc genusLevel
    elif depth > 0:
      inc belowGenus
  echo "Taxonomy summary: ", speciesLevel, " sequences at species level, ",
      genusLevel, " at genus level, ", belowGenus, " below genus, ",
      noTaxCount, " without taxonomy."

  echo "Building index from ", dbSeqs.len, " database sequences..."
  let t0 = epochTime()
  var idx = buildIndex(dbSeqs, taxStrings)
  let t1 = epochTime()
  echo "Index built in ", formatFloat(t1 - t0, format = ffDecimal,
      precision = 2), " seconds."

  # Raw DB data is fully absorbed by the index: release it before
  # classifying to keep the memory footprint small.
  dbSeqs = newSeq[string]()
  taxStrings = newSeq[string]()

  var state = initSintaxState(idx)
  let idxPtr = unsafeAddr idx

  var tabF: File
  let writeTab = opts.tabbedout != ""
  if writeTab:
    tabF = open(opts.tabbedout, fmWrite)

  var fQ = readfx.xopen[GzFile](opts.input, mode = fmRead)
  if threads == 1:
    echo "Classifying queries..."
  elif threads == 0:
    echo "Classifying queries (threaded: auto)..."
  else:
    echo "Classifying queries (threaded: max ", threads, " concurrent tasks)..."
  let t2 = epochTime()

  var names = newSeq[string](BatchSize)
  var qseqs = newSeq[string](BatchSize)
  var outTaxArr = newSeq[string](BatchSize)
  var passedTaxArr = newSeq[string](BatchSize)
  var strandArr = newSeq[char](BatchSize)
  var domArr = newSeq[bool](BatchSize)

  var queryCount = 0
  var domainClassified = 0
  var lastReportAt = 0
  var lastReportTime = t2
  var atEOF = false

  while not atEOF:
    var batchLen = 0
    while batchLen < BatchSize:
      if readFastx(fQ, record):
        names[batchLen] = record.name
        qseqs[batchLen] = record.sequence
        inc batchLen
      else:
        atEOF = true
        break

    if batchLen == 0:
      break

    classifyBatch(qseqs, outTaxArr, passedTaxArr, strandArr, domArr,
        batchLen, idxPtr, cutoff, threads, state)

    for i in 0 ..< batchLen:
      if writeTab:
        tabF.writeLine(names[i] & "\t" & outTaxArr[i] & "\t" &
            strandArr[i] & "\t" & passedTaxArr[i])
      if domArr[i]:
        inc domainClassified
    queryCount += batchLen

    let now = epochTime()
    if queryCount - lastReportAt >= 25000 or now - lastReportTime >= 5.0:
      echo queryCount, " queries processed (", domainClassified,
          " classified at least at domain level)..."
      lastReportAt = queryCount
      lastReportTime = now

  let t3 = epochTime()
  if writeTab:
    tabF.close()
  echo "Classification took ", formatFloat(t3 - t2, format = ffDecimal,
      precision = 2), " seconds."
  echo queryCount, " queries processed, ", domainClassified,
      " classified at least at domain level."
  echo "Done."

when isMainModule:
  main()
