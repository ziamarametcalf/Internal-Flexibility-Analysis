python
import os
from pymol import cmd

# ---- CONFIG ----
ref_prefix = "ref_"   # expects ref_1 ... ref_N
num_refs   = 6
# RNA backbone atoms (OP1/OP2 and O1P/O2P tolerant)
backbone = "name P+O5'+C5'+C4'+C3'+O3'+O1P+O2P+OP1+OP2"
core_sel = f"resi 19-43 and ({backbone})"
# ---------------

# Helper: count core atoms in a specific state
def core_atoms(obj, st):
    return cmd.count_atoms(f"{obj} and ({core_sel})", state=st)

# Prepare output
out_path = os.path.join(os.getcwd(), "RMSD_ref_to_ref.txt")
f = open(out_path, "w", encoding="utf-8")
f.write("=== Pairwise RMSD between references (best state-to-state over core 19–43 backbone) ===\n")
f.write(f"Core: {core_sel}\n\n")
f.write("Format: refA vs refB  ->  best(A_state, B_state)  RMSD\n\n")

# Gather refs and their state counts
refs = []
for i in range(1, num_refs+1):
    rname = f"{ref_prefix}{i}"
    if cmd.count_atoms(rname) == 0:
        continue
    S = cmd.count_states(rname)
    if S < 1:
        continue
    # Skip refs that have zero atoms in core across all states
    if cmd.count_atoms(f"{rname} and ({core_sel})") == 0:
        continue
    refs.append((rname, S))

if len(refs) < 2:
    print("[WARN] Fewer than two valid references with core atoms; nothing to compare.")
    f.write("[WARN] Fewer than two valid references with core atoms; nothing to compare.\n")
    f.close()
    print("Wrote:", out_path)
else:
    # Compute pairwise best RMSD (min over all state pairs)
    for i in range(len(refs)):
        A, Sa = refs[i]
        for j in range(i+1, len(refs)):
            B, Sb = refs[j]

            best = (float("inf"), None, None)  # (rmsd, a_state, b_state)

            for a in range(1, Sa+1):
                if core_atoms(A, a) == 0:
                    continue
                for b in range(1, Sb+1):
                    if core_atoms(B, b) == 0:
                        continue

                    rms = cmd.super(f"{A} and ({core_sel})", f"{B} and ({core_sel})",
                                    cycles=5, quiet=1,
                                    mobile_state=a, target_state=b,
                                    transform=0)[0]
                    if rms < best[0]:
                        best = (rms, a, b)

            if best[1] is None:
                line = f"{A} vs {B}  ->  no overlapping atoms in core; skipped\n"
            else:
                line = f"{A} vs {B}  ->  best(A_state={best[1]}, B_state={best[2]})  RMSD={best[0]:.3f} Å\n"

            f.write(line)
            print(line, end="")

    f.close()
    print("\n✅ Pairwise RMSDs written to:")
    print(out_path)
python end
