# Package

version       = "0.2.0"
author        = "Andrea Telatin"
description   = "Nim library for USEARCH algorithms unoise and sintax."
license       = "GPL-3.0"
srcDir        = "src"
binDir        = "bin"
bin           = @["unoise", "sintax", "nbc", "uchime", "remove_phix"]


# Dependencies

requires "nim >= 2.0.0"
requires "readfx >= 0.3.1"
requires "argparse >= 4.0.0"
requires "malebolgia >= 1.3.0"

task docs, "Generate HTML documentation":
  exec "nim doc --project --index:on --outdir:docs src/unoisenim.nim"
