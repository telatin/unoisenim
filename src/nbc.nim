import argparse, readfx, strutils, times
import unoisenim/nbc_algo

proc extractTaxFromLabel(label: string): string =
  for p in label.split(';'):
    if p.startsWith("tax="):
      return p[4..^1]
  return ""

proc main() =
  let p = newParser("nbc"):
    help("RDP-style Naive Bayesian taxonomy classifier in Nim")
    option("-i", "--input", help = "Input query FASTA file")
    option("-d", "--database", help = "Reference database FASTA file with taxonomy annotations")
    option("-t", "--tabbedout", help = "Output tabbed file")
    option("-c", "--cutoff", default = some("0.8"),
        help = "Confidence cutoff (default 0.8)")
    option("--boot-iters", default = some("100"),
        help = "Bootstrap iterations (default 100)")
    option("--min-words", default = some("5"),
        help = "Minimum words sampled per bootstrap (default 5)")

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
  let bootIters = parseInt(opts.boot_iters)
  let minWords = parseInt(opts.min_words)

  var dbSeqs = newSeq[string]()
  var taxStrings = newSeq[string]()
  var fDb = readfx.xopen[GzFile](opts.database, mode = fmRead)
  var record: FQRecord
  while readFastx(fDb, record):
    dbSeqs.add(record.sequence)
    taxStrings.add(extractTaxFromLabel(record.name))

  echo "Building NBC index from ", dbSeqs.len, " database sequences..."
  let t0 = cpuTime()
  let idx = buildNbcIndex(dbSeqs, taxStrings)
  let t1 = cpuTime()
  echo "Index built in ", (t1 - t0), " seconds."

  var tabF: File
  let writeTab = opts.tabbedout != ""
  if writeTab:
    tabF = open(opts.tabbedout, fmWrite)

  var fQ = readfx.xopen[GzFile](opts.input, mode = fmRead)
  echo "Classifying queries..."
  let t2 = cpuTime()

  while readFastx(fQ, record):
    let hit = nbc(record.sequence, idx, bootIters = bootIters, minWords = minWords)
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

    if writeTab:
      tabF.writeLine(record.name & "\t" & outTax & "\t" & hit.strand & "\t" & passedTax)

  let t3 = cpuTime()
  if writeTab:
    tabF.close()
  echo "Classification took ", (t3 - t2), " seconds."
  echo "Done."

when isMainModule:
  main()
