# Parse samtools flagstat into a coarse summary:
# mapper  mm  total  mapped  pct_mapped

from pathlib import Path
import re
import pandas as pd

def parse_flagstat(path: str):
    total = mapped = None
    with open(path, "r") as fh:
        for line in fh:
            m = re.search(r"^(\d+)\s+\+\s+\d+\s+in total", line)
            if m: total = int(m.group(1))
            m = re.search(r"^(\d+)\s+\+\s+\d+\s+mapped", line)
            if m: mapped = int(m.group(1))
    if total is None or mapped is None:
        raise ValueError(f"Could not parse totals from {path}")
    return total, mapped

rows = []
for path in snakemake.input:
    p = Path(path)
    # Expect: map/<mapper>/mmXX/flagstat.txt
    mapper = p.parts[1]
    mm = int(re.search(r"mm(\d+)", path).group(1))
    total, mapped = parse_flagstat(path)
    pct = (100.0 * mapped / total) if total else 0.0
    rows.append((mapper, mm, total, mapped, pct))

df = pd.DataFrame(rows, columns=["mapper", "mm", "total", "mapped", "pct_mapped"])
df.sort_values(["mapper", "mm"], inplace=True)
df.to_csv(snakemake.output[0], sep="\t", index=False, float_format="%.4f")
