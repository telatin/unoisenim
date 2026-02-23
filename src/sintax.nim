import argparse, readfx, strutils, os, times
import unoisenim/sintax_algo

proc main() =
  let p = newParser("sintax"):
    help("USEARCH sintax algorithm implementation in Nim")
    option("-i", "--input", help = "Input query FASTA file")
    option("-d", "--database", help = "Reference database FASTA file with taxonomy annotations")
    option("-t", "--tabbedout", help = "Output tabbed file")
    option("-c", "--cutoff", default = some("0.8"),
        help = "Confidence cutoff (default 0.8)")

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
  echo "Classifying queries..."
  let t2 = cpuTime()

  while readFastx(fQ, record):
    let hit = sintax(record.sequence, idx)
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
