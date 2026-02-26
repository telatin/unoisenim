import argparse, readfx, strutils, os, times, algorithm
import unoisenim/utils
import unoisenim/unoise_algo
import unoisenim/uchime2_algo

proc main() =
  let p = newParser("uchime"):
    help("UCHIME2 chimera detection on a dereplicated FASTA file")
    option("-i", "--input", help = "Input dereplicated FASTA file with ;size=N annotations")
    option("-o", "--output", help = "Output FASTA file (non-chimeric sequences, original headers)")
    option("-s", "--summary", help = "Output TSV with chimera status for every input sequence")
    option("--min-skew", default = some("16.0"),
        help = "Minimum abundance ratio between query and parent (default 16.0)")
    option("--threads", default = some("1"),
        help = "Threading: 1=sequential, 0=auto threaded, >1=fixed max concurrent chunk tasks")

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

  if opts.input == "":
    stderr.writeLine "Error: --input is required"
    echo p.help
    quit(1)

  let minSkew = parseFloat(opts.min_skew)
  let threads = parseInt(opts.threads)

  # Read input FASTA
  var seqs = newSeq[UnoiseSeq]()
  var f = readfx.xopen[GzFile](opts.input, mode = fmRead)
  var record: FQRecord

  echo "Reading: ", opts.input
  while readFastx(f, record):
    let sz = parseSize(record.name)
    seqs.add(UnoiseSeq(id: record.name, seq: record.sequence, size: sz))

  echo "Read ", seqs.len, " sequences."

  # Build centroids sorted descending by abundance.
  # totalSize == size since there is no prior UNOISE merging.
  var centroids = newSeq[Centroid](seqs.len)
  for i, s in seqs:
    centroids[i] = Centroid(seqObj: s, totalSize: s.size)
  centroids.sort(proc(a, b: Centroid): int = cmp(b.totalSize, a.totalSize))

  if threads == 1:
    echo "UCHIME mode: sequential"
  elif threads == 0:
    echo "UCHIME mode: threaded (auto scheduler)"
  else:
    echo "UCHIME mode: threaded (max concurrent chunk tasks=", threads, ")"

  let t0 = cpuTime()
  let chimeras = uchime(centroids, minSkew, threads)
  let t1 = cpuTime()
  echo "UCHIME filtering took ", (t1 - t0), " seconds."

  var chimeraCount = 0
  for c in chimeras:
    if c: inc chimeraCount
  echo "Chimeras: ", chimeraCount, " / ", centroids.len, " sequences."

  # Write non-chimeric FASTA (preserves original headers)
  if opts.output != "":
    var outF = open(opts.output, fmWrite)
    for i, c in centroids:
      if not chimeras[i]:
        outF.writeLine(">" & c.seqObj.id)
        outF.writeLine(c.seqObj.seq)
    outF.close()
    echo "Non-chimeric sequences written to: ", opts.output

  # Write summary TSV: one row per input sequence, abundance-sorted
  if opts.summary != "":
    var sumF = open(opts.summary, fmWrite)
    sumF.writeLine("id\tsize\tstatus")
    for i, c in centroids:
      let status = if chimeras[i]: "chimera" else: "ok"
      sumF.writeLine(c.seqObj.id & "\t" & $c.totalSize & "\t" & status)
    sumF.close()
    echo "Summary written to: ", opts.summary

when isMainModule:
  main()
