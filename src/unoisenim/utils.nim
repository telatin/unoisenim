import strutils

proc parseSize*(label: string): int =
  let parts = label.split(";")
  for p in parts:
    if p.startsWith("size="):
      try:
        return parseInt(p[5..^1])
      except ValueError:
        return 0
  return 0
