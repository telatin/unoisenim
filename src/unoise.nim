import argparse, readfx, strutils, os, times
import unoisenim/utils
import unoisenim/unoise_algo
import unoisenim/uchime2_algo

proc main() =
  let p = newParser("unoise"):
    help("USEARCH unoise algorithm implementation in Nim")
    arg("input", help = "Input FASTA file with size=X; annotations")
    option("-z", "--zotus", help = "Output ZOTUs FASTA file")
    option("-a", "--alpha", default = some("2.0"),
        help = "Alpha parameter for skew calculation (default 2.0)")
    option("-m", "--minsize", default = some("8"),
        help = "Minimum abundance for a sequence to be processed (default 8)")
    option("--uchime-threads", default = some("0"),
        help = "UCHIME threads: 0=auto threaded, 1=sequential, >1 fixed threadpool size")

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

  let inputFile = opts.input
  let alpha = parseFloat(opts.alpha)
  let minsize = parseInt(opts.minsize)
  let uchimeThreads = parseInt(opts.uchime_threads)

  var seqs = newSeq[UnoiseSeq]()
  var f = readfx.xopen[GzFile](inputFile, mode = fmRead)
  var record: FQRecord

  echo "Reading file: ", inputFile
  while readFastx(f, record):
    let sz = parseSize(record.name)
    seqs.add(UnoiseSeq(id: record.name, seq: record.sequence, size: sz))

  echo "Read ", seqs.len, " sequences."
  echo "Running UNOISE..."
  let t0 = cpuTime()
  let centroids = unoise(seqs, alpha, minsize)
  let t1 = cpuTime()
  echo "UNOISE produced ", centroids.len, " amplicons before chimera filtering."
  echo "UNOISE algorithm took ", (t1 - t0), " seconds."

  let t2 = cpuTime()
  if uchimeThreads == 1:
    echo "UCHIME mode: sequential"
  elif uchimeThreads == 0:
    echo "UCHIME mode: threaded (auto pool)"
  else:
    echo "UCHIME mode: threaded (pool=", uchimeThreads, ")"
  let chimeras = uchime(centroids, 16.0, uchimeThreads)
  let t3 = cpuTime()
  echo "UCHIME chimera filtering took ", (t3 - t2), " seconds."
  var zotuCount = 0
  for i, c in centroids:
    if not chimeras[i]:
      zotuCount += 1

  echo "Found ", centroids.len - zotuCount, " chimeras."
  echo "Final: ", zotuCount, " ZOTUs."

  if opts.zotus != "":
    var outF = open(opts.zotus, fmWrite)
    var zIndex = 1
    for i, c in centroids:
      if not chimeras[i]:
        let outLabel = "Zotu" & $zIndex
        outF.writeLine(">" & outLabel)
        outF.writeLine(c.seqObj.seq)
        zIndex += 1
    outF.close()
    echo "ZOTUs written to ", opts.zotus

when isMainModule:
  main()
