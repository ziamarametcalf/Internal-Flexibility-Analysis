import csv, glob, os, sys

# Match both no-extension and .txt variants:
PATTERNS = ["2kdq_model*-torsions", "2kdq_model*-torsions.*"]
ANGLES = ["alpha","beta","gamma","delta","epsilon","zeta","chi"]
OUT_FILE = "2kdq_master_torsions.csv"

def parse_torsions(path):
    rows = []
    header = False
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if not header:
                if s.lower().startswith("nt "):  # start of main torsion block
                    header = True
                continue
            # stop at next section
            if s.startswith("*") or "Virtual eta/theta" in s:
                break
            parts = s.split()
            nt = parts[0]
            vals = []
            for tok in parts[1:1+7]:
                try:
                    vals.append(float(tok))
                except:
                    vals.append(None)
            rows.append((nt, *vals))
    return rows

# find files using both patterns
files = []
for pat in PATTERNS:
    files.extend(glob.glob(pat))
# dedupe while preserving order
seen = set(); files = [f for f in files if not (f in seen or seen.add(f))]

print(f"[INFO] cwd = {os.getcwd()}")
print(f"[INFO] found {len(files)} torsion files")
if not files:
    print("❌ No files matched. Try these commands:")
    print("   dir /b 2kdq_model*-torsions*")
    sys.exit(1)

with open(OUT_FILE, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["structure","index_in_chain"] + ANGLES)
    for file in sorted(files):
        # keep the base name (no extension) as structure label
        struct = os.path.splitext(file)[0]
        data = parse_torsions(file)
        for i, row in enumerate(data, start=1):
            _, *angles = row
            w.writerow([struct, i] + angles)

print(f"✅ Master torsion table written to {OUT_FILE}")
