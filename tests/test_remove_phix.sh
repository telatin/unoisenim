#!/usr/bin/env bash
# test_remove_phix.sh — integration test for bin/remove_phix
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
BIN="$REPO_ROOT/bin/remove_phix"
DATA="$REPO_ROOT/data/sim-reads"

PASS=0
FAIL=0

pass() { echo "  [PASS] $1"; ((PASS++)) || true; }
fail() { echo "  [FAIL] $1"; ((FAIL++)) || true; }

# ── 1. Check required inputs ──────────────────────────────────────────────────
echo "[check] Input files present"
for f in "$DATA/phix_R1.fastq.gz" "$DATA/phix_R2.fastq.gz" \
          "$DATA/eco_R1.fastq.gz"  "$DATA/eco_R2.fastq.gz"; do
  if [[ -f "$f" ]]; then
    pass "$f"
  else
    fail "$f not found"
  fi
done

if [[ ! -x "$BIN" ]]; then
  echo "ERROR: binary not found at $BIN — run 'nimble build' first"
  exit 1
fi

# ── 2. Temporary working directory ───────────────────────────────────────────
TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT
echo "[tmpdir] $TMPDIR"

# ── 3. Run PhiX removal ───────────────────────────────────────────────────────
echo
echo "[run] PhiX-spiked paired-end reads (expect 100% removal)"
PHIX_OUT=$("$BIN" \
  -1 "$DATA/phix_R1.fastq.gz" \
  -2 "$DATA/phix_R2.fastq.gz" \
  -o "$TMPDIR/phix_nophix_R1.fastq" \
  -O "$TMPDIR/phix_nophix_R2.fastq" \
  --report "$TMPDIR/phix_report.tsv" 2>/dev/null)
echo "  $PHIX_OUT"

echo
echo "[run] E. coli paired-end reads (expect 0% removal)"
ECO_OUT=$("$BIN" \
  -1 "$DATA/eco_R1.fastq.gz" \
  -2 "$DATA/eco_R2.fastq.gz" \
  -o "$TMPDIR/eco_nophix_R1.fastq" \
  -O "$TMPDIR/eco_nophix_R2.fastq" \
  --report "$TMPDIR/eco_report.tsv" 2>/dev/null)
echo "  $ECO_OUT"

# ── 4. Check results ──────────────────────────────────────────────────────────
echo
echo "[check] Results"

# Helper: extract a field from the reads_in=N reads_removed=N pct=X% line
reads_in()      { echo "$1" | grep -oP 'reads_in=\K[0-9]+'; }
reads_removed() { echo "$1" | grep -oP 'reads_removed=\K[0-9]+'; }

PHIX_IN=$(reads_in      "$PHIX_OUT")
PHIX_RM=$(reads_removed "$PHIX_OUT")
ECO_IN=$(reads_in       "$ECO_OUT")
ECO_RM=$(reads_removed  "$ECO_OUT")

# PhiX sample: all reads should be removed
if [[ "$PHIX_RM" -eq "$PHIX_IN" && "$PHIX_IN" -gt 0 ]]; then
  pass "PhiX sample: $PHIX_RM/$PHIX_IN reads removed (100%)"
else
  fail "PhiX sample: expected 100% removal, got $PHIX_RM/$PHIX_IN"
fi

# E. coli sample: no reads should be removed
if [[ "$ECO_RM" -eq 0 && "$ECO_IN" -gt 0 ]]; then
  pass "E. coli sample: 0/$ECO_IN reads removed (0%)"
else
  fail "E. coli sample: expected 0% removal, got $ECO_RM/$ECO_IN"
fi

# Output FASTQ for PhiX run should be empty (all discarded)
PHIX_LINES=$(wc -l < "$TMPDIR/phix_nophix_R1.fastq")
if [[ "$PHIX_LINES" -eq 0 ]]; then
  pass "PhiX R1 output is empty"
else
  fail "PhiX R1 output unexpectedly has $PHIX_LINES lines"
fi

# Output FASTQ for E. coli run should have records (4 lines per read)
ECO_LINES=$(wc -l < "$TMPDIR/eco_nophix_R1.fastq")
ECO_EXPECTED=$(( ECO_IN * 4 ))
if [[ "$ECO_LINES" -eq "$ECO_EXPECTED" ]]; then
  pass "E. coli R1 output has $ECO_LINES lines ($ECO_IN reads × 4)"
else
  fail "E. coli R1 output: expected $ECO_EXPECTED lines, got $ECO_LINES"
fi

# Report TSV should have a header + one data line
PHIX_REPORT_LINES=$(wc -l < "$TMPDIR/phix_report.tsv")
if [[ "$PHIX_REPORT_LINES" -eq 2 ]]; then
  pass "PhiX report has 2 lines (header + data)"
else
  fail "PhiX report: expected 2 lines, got $PHIX_REPORT_LINES"
fi

# ── Summary ───────────────────────────────────────────────────────────────────
echo
echo "Results: $PASS passed, $FAIL failed"
[[ "$FAIL" -eq 0 ]]
