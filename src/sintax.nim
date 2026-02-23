import argparse, readfx, strutils, os, times
import std/cpuinfo
import malebolgia
import unoisenim/sintax_algo

proc main() =
  let p = newParser("sintax"):
    help("USEARCH sintax algorithm implementation in Nim")
    option("-i", "--input", help = "Input query FASTA file")
    option("-d", "--database", help = "Reference database FASTA file with taxonomy annotations")
    option("-t", "--tabbedout", help = "Output tabbed file")
    option("-c", "--cutoff", default = some("0.8"),
        help = "Confidence cutoff (default 0.8)")
    option("-w", "--workers", default = some("1"),
        help = "Number of worker threads (default 1)")

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
  var nThreads = parseInt(opts.workers)
  if nThreads <= 0:
    nThreads = countProcessors()

  # Load DB
  var dbSeqs = newSeq[string]()
  var taxStrings = newSeq[string]()
  var fDb = readfx.xopen[GzFile](opts.database, mode = fmRead)
  var record: FQRecord
  while readFastx(fDb, record):
    dbSeqs.add(record.sequence)
    # Extract tax=... from label
    var tax = ""
    for p in record.name.split(';'):
      if p.startsWith("tax="):
        tax = p[4..^1]
        break
    taxStrings.add(tax)

  echo "Building index from ", dbSeqs.len, " database sequences..."
  let t0 = cpuTime()
  let idx = buildIndex(dbSeqs, taxStrings)
  let t1 = cpuTime()
  echo "Index built in ", (t1 - t0), " seconds."

  var tabF: File
  let writeTab = opts.tabbedout != ""
  if writeTab:
    tabF = open(opts.tabbedout, fmWrite)

  var fQ = readfx.xopen[GzFile](opts.input, mode = fmRead)
  echo "Classifying queries using ", nThreads, " threads..."
  let t2 = cpuTime()

  type OutQuery = tuple[name, outTax, strand, passedTax: string]
  var batch = newSeq[string]()
  var names = newSeq[string]()
  let batchSize = 1000 * nThreads

  proc processSeq(q: string, qName: string, idxPtr: ptr SintaxIndex): OutQuery =
    let hit = sintax(q, idxPtr[])
    var outTax = ""
    var passedTax = ""
    if hit.rankNames.len > 0:
      var taxes = newSeq[string]()
      for i in 0 ..< hit.rankNames.len:
        taxes.add(hit.rankNames[i] & "(" & formatFloat(hit.rankProbs[i],
            format = ffDecimal, precision = 4) & ")")
      outTax = taxes.join(",")

      var passed = newSeq[string]()
      for i in 0 ..< hit.rankNames.len:
        if hit.rankProbs[i] >= cutoff:
          passed.add(hit.rankNames[i])
        else:
          break
      passedTax = passed.join(",")
    else:
      outTax = "*"
      passedTax = "*"

    return (qName, outTax, $hit.strand, passedTax)

  var m = createMaster()

  while true:
    var record: FQRecord
    let hasMore = readFastx(fQ, record)
    if hasMore:
      batch.add(record.sequence)
      names.add(record.name)

    if batch.len >= batchSize or (not hasMore and batch.len > 0):
      # Process this batch concurrently
      var results = newSeq[OutQuery](batch.len)
      m.awaitAll:
        for i in 0 ..< batch.len:
          m.spawn processSeq(batch[i], names[i], addr idx) -> results[i]

      if writeTab:
        for res in results:
          tabF.writeLine(res.name & "\t" & res.outTax & "\t" & res.strand &
              "\t" & res.passedTax)

      batch.setLen(0)
      names.setLen(0)

    if not hasMore:
      break

  let t3 = cpuTime()
  if writeTab:
    tabF.close()
  echo "Classification took ", (t3 - t2), " seconds."
  echo "Done."

when isMainModule:
  main()
