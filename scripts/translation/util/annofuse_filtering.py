# util/annofuse_filtering.py
# Shared helpers to gate fusion calls by *.annoFuse_filter.tsv

import csv

# ---------- Canonicalization (handles "/" vs "," and token order) ----------

def _canon_side(side: str) -> str:
    """
    Normalize one side of FusionName:
    - uppercase
    - treat "/" and "," as equivalent
    - strip spaces
    - deduplicate tokens
    - sort tokens for stable, order-insensitive comparison
    """
    side = (side or "").upper().replace("/", ",")
    toks = [t.strip() for t in side.split(",") if t.strip()]
    return ",".join(sorted(set(toks)))
    
def canonical_fusionname(fn: str) -> str:
    """
    Canonicalize FusionName "LEFT--RIGHT" so formatting differences
    ("/" vs ",", token order) don't affect equality.
    """
    s = (fn or "").strip().upper()
    if "--" not in s:
        return _canon_side(s)
    left, right = s.split("--", 1)
    return f"{_canon_side(left)}--{_canon_side(right)}"
  
# ---------- Load AnnoFuse keep list ----------

def load_annofuse(path: str) -> set[str]:
    """
    Load FusionName values from *.annoFuse_filter.tsv and canonicalize them.
    Requires a 'FusionName' column.
    """
    keep = set()
    with open(path, "r", newline="") as f:
        rdr = csv.DictReader(f, delimiter="\t")
        if "FusionName" not in rdr.fieldnames:
            raise ValueError(f"AnnoFuse file {path} must contain 'FusionName'")
        for row in rdr:
            v = (row.get("FusionName") or "").strip()
            if v:
                keep.add(canonical_fusionname(v))
    print(f"Loaded {len(keep)} canonical FusionName entries from {path}")
    return keep

# ---------- Filters (return header + kept rows; print summary) ----------

def filter_arriba_by_annofuse(rows, annofuse_keep_set):
    """Keep Arriba rows whose ordered FusionName (gene1--gene2) appears in AnnoFuse keep set"""
    if not rows or not annofuse_keep_set:
        return rows
    header = rows[0]
    i1 = header.index("gene1")
    i2 = header.index("gene2")
    kept, removed = [header], 0
    for r in rows[1:]:
        name = canonical_fusionname(f"{r[i1]}--{r[i2]}")
        if name not in annofuse_keep_set:
            removed += 1
            continue
        kept.append(r)
    print(f"AnnoFuse gate: removed {removed}, kept {len(kept)-1}")
    return kept

def filter_starfusion_by_annofuse(rows, annofuse_keep_set):
    """Keep STAR-Fusion rows whose FusionName is present in the AnnoFuse keep set"""
    if not rows or not annofuse_keep_set:
        return rows
    header = rows[0]
    iF = header.index("FusionName")
    kept, removed = [header], 0
    for r in rows[1:]:
        name = canonical_fusionname(r[iF])
        if name not in annofuse_keep_set:
            removed += 1
            continue
        kept.append(r)
    print(f"AnnoFuse gate: removed {removed}, kept {len(kept)-1}")
    return kept
