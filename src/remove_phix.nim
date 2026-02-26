import argparse, readfx, strutils, malebolgia
import unoisenim/remove_phix_utils

const BatchSize = 1024
const ChunkSize = 64

# Module-level worker: classify a contiguous slice of a sequence array.
# Must be {.gcsafe.} and take no closures so malebolgia can spawn it safely.
proc classifyPhixRange(seqs: ptr UncheckedArray[string],
                       flags: ptr UncheckedArray[bool],
                       start, count: int,
                       minId: float, minKmers: int) {.gcsafe.} =
  for i in start ..< start + count:
    flags[i] = isPhix(seqs[i], minId, minKmers)

# Classify n entries of seqs[] into flags[].
# threads semantics match uchime:
#   1  => sequential
#   0  => auto (malebolgia scheduler, all chunks in one awaitAll)
#   >1 => fixed max concurrent tasks per awaitAll round
proc classifyBatch(seqs: var seq[string], flags: var seq[bool],
                   n, threads: int, minId: float, minKmers: int) =
  if threads == 1 or n < ChunkSize:
    for i in 0 ..< n:
      flags[i] = isPhix(seqs[i], minId, minKmers)
    return

  let sp = cast[ptr UncheckedArray[string]](addr seqs[0])
  let fp = cast[ptr UncheckedArray[bool]](addr flags[0])
  var start = 0

  if threads > 1:
    while start < n:
      var m = createMaster()
      m.awaitAll:
        var launched = 0
        while start < n and launched < threads:
          let s = start
          let c = min(ChunkSize, n - s)
          m.spawn classifyPhixRange(sp, fp, s, c, minId, minKmers)
          start += c
          inc launched
  else: # threads == 0: auto-scheduler
    var m = createMaster()
    m.awaitAll:
      while start < n:
        let s = start
        let c = min(ChunkSize, n - s)
        m.spawn classifyPhixRange(sp, fp, s, c, minId, minKmers)
        start += c

proc writeFastq(f: File, rec: FQRecord) =
  let nameLine = if rec.comment.len > 0: rec.name & " " & rec.comment else: rec.name
  f.writeLine("@" & nameLine)
  f.writeLine(rec.sequence)
  f.writeLine("+")
  f.writeLine(rec.quality)

proc main() =
  let p = newParser("remove_phix"):
    help("Remove PhiX174 contamination from Illumina FASTQ reads")
    option("-1", "--read1", help = "Input R1 (or single-end) FASTQ [.gz supported]")
    option("-2", "--read2", help = "Input R2 FASTQ (paired-end)")
    option("-o", "--out1", help = "Output R1 FASTQ")
    option("-O", "--out2", help = "Output R2 FASTQ (required when --read2 is given)")
    option("--min-id", default = some("0.97"),
        help = "Identity threshold for PhiX removal (default: 0.97)")
    option("--min-kmers", default = some("8"),
        help = "Min valid 8-mers required to classify a read (default: 8)")
    option("--paired-mode", default = some("strict"),
        help = "strict=remove pair if EITHER is PhiX; lenient=remove only if BOTH (default: strict)")
    option("--threads", default = some("1"),
        help = "Threads: 1=sequential, 0=auto, N=fixed concurrency (default: 1)")
    option("--report", help = "Write TSV summary (reads_in, reads_removed, pct) to FILE")
    flag("-v", "--verbose", help = "Print per-read scores to stderr")

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

  if opts.read1 == "" or opts.out1 == "":
    stderr.writeLine "Error: --read1 and --out1 are required"
    echo p.help
    quit(1)

  if opts.read2 != "" and opts.out2 == "":
    stderr.writeLine "Error: --out2 is required when --read2 is given"
    echo p.help
    quit(1)

  let minId = parseFloat(opts.min_id)
  let minKmers = parseInt(opts.min_kmers)
  let strictMode = opts.paired_mode != "lenient"
  let threads = parseInt(opts.threads)

  var readsIn = 0
  var readsRemoved = 0

  # Reusable batch buffers
  var seqs1 = newSeq[string](BatchSize)
  var seqs2 = newSeq[string](BatchSize)
  var flags1 = newSeq[bool](BatchSize)
  var flags2 = newSeq[bool](BatchSize)

  if opts.read2 == "":
    # Single-end mode
    var batch = newSeq[FQRecord](BatchSize)
    var fOut = open(opts.out1, fmWrite)
    var fR1 = readfx.xopen[GzFile](opts.read1, mode = fmRead)
    var rec: FQRecord
    var atEOF = false

    while not atEOF:
      var batchLen = 0
      while batchLen < BatchSize:
        if readFastx(fR1, rec):
          batch[batchLen] = rec
          seqs1[batchLen] = rec.sequence
          inc batchLen
          inc readsIn
        else:
          atEOF = true
          break

      if batchLen == 0:
        break

      classifyBatch(seqs1, flags1, batchLen, threads, minId, minKmers)

      for i in 0 ..< batchLen:
        if opts.verbose:
          stderr.writeLine batch[i].name & "\t" &
              formatFloat(phixScore(batch[i].sequence), ffDecimal, 4) & "\t" &
              (if flags1[i]: "phix" else: "keep")
        if flags1[i]:
          inc readsRemoved
        else:
          writeFastq(fOut, batch[i])

    fOut.close()
  else:
    # Paired-end mode
    var batch1 = newSeq[FQRecord](BatchSize)
    var batch2 = newSeq[FQRecord](BatchSize)
    var fOut1 = open(opts.out1, fmWrite)
    var fOut2 = open(opts.out2, fmWrite)
    var fR1 = readfx.xopen[GzFile](opts.read1, mode = fmRead)
    var fR2 = readfx.xopen[GzFile](opts.read2, mode = fmRead)
    var rec1, rec2: FQRecord
    var atEOF = false

    while not atEOF:
      var batchLen = 0
      while batchLen < BatchSize:
        if readFastx(fR1, rec1) and readFastx(fR2, rec2):
          batch1[batchLen] = rec1
          batch2[batchLen] = rec2
          seqs1[batchLen] = rec1.sequence
          seqs2[batchLen] = rec2.sequence
          inc batchLen
          inc readsIn
        else:
          atEOF = true
          break

      if batchLen == 0:
        break

      classifyBatch(seqs1, flags1, batchLen, threads, minId, minKmers)
      classifyBatch(seqs2, flags2, batchLen, threads, minId, minKmers)

      for i in 0 ..< batchLen:
        let remove = if strictMode: flags1[i] or flags2[i] else: flags1[i] and flags2[i]
        if opts.verbose:
          stderr.writeLine batch1[i].name & "\t" &
              formatFloat(phixScore(batch1[i].sequence), ffDecimal, 4) & "/" &
              formatFloat(phixScore(batch2[i].sequence), ffDecimal, 4) & "\t" &
              (if remove: "phix" else: "keep")
        if remove:
          inc readsRemoved
        else:
          writeFastq(fOut1, batch1[i])
          writeFastq(fOut2, batch2[i])

    fOut1.close()
    fOut2.close()

  let pct = if readsIn > 0: 100.0 * float(readsRemoved) / float(readsIn) else: 0.0
  echo "reads_in=", readsIn, " reads_removed=", readsRemoved,
      " pct=", formatFloat(pct, ffDecimal, 2), "%"

  if opts.report != "":
    var repF = open(opts.report, fmWrite)
    repF.writeLine("reads_in\treads_removed\tpct_removed")
    repF.writeLine($readsIn & "\t" & $readsRemoved & "\t" &
        formatFloat(pct, ffDecimal, 4))
    repF.close()

when isMainModule:
  main()
