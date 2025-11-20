import os, sys, re
src = sys.argv[1]  # e.g., 1anr.pdb
base = os.path.splitext(os.path.basename(src))[0]
with open(src, encoding="utf-8", errors="ignore") as f:
    text = f.read()
blocks = re.split(r'(?m)^MODEL\b.*?\n', text)[1:]
if not blocks:
    print("No MODEL blocks found. Is this an NMR ensemble PDB?")
    raise SystemExit(1)
blocks = [re.split(r'(?m)^ENDMDL\b.*?\n', b)[0] for b in blocks]
outdir = "models"
os.makedirs(outdir, exist_ok=True)
for i, b in enumerate(blocks, 1):
    path = os.path.join(outdir, f"{base}_model{i:02d}.pdb")
    with open(path, "w", encoding="utf-8") as w:
        w.write(b if b.endswith("\n") else b+"\n")
    print("wrote", path)
print(f"Done. Wrote {len(blocks)} model files to ./{outdir}")
