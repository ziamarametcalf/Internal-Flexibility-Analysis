python
import os
from pymol import cmd

# ---- CONFIG ----
ref_prefix = "ref_"   # expects ref_1 ... ref_6
num_refs   = 6
conf_count = 20
# RNA backbone atoms
backbone = "name P+O5'+C5'+C4'+C3'+O3'+O1P+O2P+OP1+OP2"
core_sel = f"resi 19-43 and ({backbone})"
# ---------------

# Safety cleanup: remove any old selections that conflict
for n in cmd.get_names("selections"):
    if n.startswith("conf_") or n.startswith("ref_"):
        cmd.delete(n)

# Check conformers
missing = []
for i in range(1, conf_count+1):
    if cmd.count_atoms(f"conf_{i}") == 0:
        missing.append(i)
if len(missing) == conf_count:
    raise SystemExit("[ERROR] No conf_# objects found. Load conf_1..conf_20 first.")
elif missing:
    print(f"[WARN] Missing confs (skipped): {missing}")

# Determine output file path
out_path = os.path.join(os.getcwd(), "RMSD_results_all.txt")
out = open(out_path, "w", encoding="utf-8")

out.write("=== RMSD RESULTS ===\n")
out.write(f"Core selection: {core_sel}\n\n")

for r in range(1, num_refs+1):
    ref = f"{ref_prefix}{r}"
    if cmd.count_atoms(ref) == 0:
        print(f"[WARN] {ref} not found; skipping.")
        continue

    ref_states = cmd.count_states(ref)
    if ref_states < 1:
        print(f"[WARN] {ref} has 0 states; skipping.")
        continue

    out.write(f"\n=== {ref} ===\n")

    for rs in range(1, ref_states+1):
        # Skip if no atoms in this state
        if cmd.count_atoms(f"{ref} and ({core_sel})", state=rs) == 0:
            continue

        out.write(f"\nState {rs}:\n")
        out.write("Conformer\tRMSD (Å)\n")
        out.write("-" * 30 + "\n")

        best_rmsd = float("inf")
        best_conf = None

        for i in range(1, conf_count+1):
            obj = f"conf_{i}"
            if cmd.count_atoms(obj) == 0 or cmd.count_atoms(f"{obj} and ({core_sel})") == 0:
                continue

            mobile_sel = f"{obj} and ({core_sel})"
            target_sel = f"{ref} and ({core_sel})"

            # conf_i (state 1) vs ref (state rs)
            rms = cmd.super(mobile_sel, target_sel,
                            cycles=5, quiet=1,
                            mobile_state=1, target_state=rs,
                            transform=0)[0]

            out.write(f"{obj}\t{rms:.3f}\n")

            if rms < best_rmsd:
                best_rmsd, best_conf = rms, obj

        if best_conf is not None:
            out.write(f"\nBest match: {best_conf} (RMSD {best_rmsd:.3f} Å)\n")

out.close()

print("\n✅ All RMSD results written to:")
print(out_path)
print("You can open this file directly in Notepad or any text editor.")
python end
