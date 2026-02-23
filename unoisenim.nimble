# Package

version       = "0.1.0"
author        = "Andrea Telatin"
description   = "Nim library for USEARCH algorithms unoise and sintax."
license       = "MIT"
srcDir        = "src"
binDir        = "bin"
bin           = @["unoise", "sintax"]


# Dependencies

requires "nim >= 2.0.0"
requires "readfx >= 0.1.0"
requires "argparse >= 4.0.0"
requires "malebolgia >= 1.3.0"
