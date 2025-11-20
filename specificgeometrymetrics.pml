python
import os, math
from pymol import cmd

# --- Make sure the apo (1anr) is called ref_0 ---
if "1anr" in cmd.get_names("objects"):
    cmd.set_name("1anr", "ref_0")
    print('Renamed "1anr" → ref_0')

# --- Label your references (update if known) ---
CLASS = {
    'ref_0': 'apo',        # apo baseline
    'ref_1': 'unknown',
    'ref_2': 'unknown',
    'ref_3': 'unknown',
    'ref_4': 'unknown',
    'ref_5': 'unknown',
    'ref_6': 'unknown',
}

# --- Define residue regions and atom groups ---
backbone = "name P+O5'+C5'+C4'+C3'+O3'+O1P+O2P+OP1+OP2"
bulge_loop = "resi 22-25+39-43"
lower_stem_P = "resi 19-22 and name P"
upper_stem_P = "resi 40-43 and name P"
groove_pairs = [(20,42),(21,41),(22,40)]

# ---------- helper functions ----------
def chi_atom_names(resn):
    resn=(resn or '').upper()
    return ("O4'","C1'","N1","C2") if resn in ("C","U","CYT","URA","CMP","UMP") else ("O4'","C1'","N9","C4")

def dihedral_for_resi(obj, state, resi):
    model = cmd.get_model(f"model {obj} and resi {resi} and name C1'", state=state)
    if not model.atom:
        return None
    resn = model.atom[0].resn
    a1,a2,a3,a4 = chi_atom_names(resn)
    sels=[f"model {obj} and resi {resi} and name {a}" for a in (a1,a2,a3,a4)]
    for s in sels:
        if cmd.count_atoms(s, state=state)==0: 
            return None
    return cmd.get_dihedral(sels[0],sels[1],sels[2],sels[3], state=state)

def avg_P_centroid(obj, st, sel):
    m=cmd.get_model(f"model {obj} and ({sel})", state=st)
    pts=[a.coord for a in m.atom if a.name=='P']
    if not pts: return None
    x=sum(p[0] for p in pts)/len(pts)
    y=sum(p[1] for p in pts)/len(pts)
    z=sum(p[2] for p in pts)/len(pts)
    return (x,y,z)

def angle_between(v1,v2):
    if not v1 or not v2: return None
    ax,ay,az=v1; bx,by,bz=v2
    a=(ax*ax+ay*ay+az*az)**0.5; b=(bx*bx+by*by+bz*bz)**0.5
    if a==0 or b==0: return None
    cos=(ax*bx+ay*by+az*bz)/(a*b)
    cos=max(-1,min(1,cos))
    return math.degrees(math.acos(cos))

def atom_coord(obj, st, resi, name):
    m=cmd.get_model(f"model {obj} and resi {resi} and name {name}", state=st)
    if not m.atom: return None
    a=m.atom[0]; return (a.coord[0],a.coord[1],a.coord[2])

def dist(a,b):
    if not a or not b: return None
    return ((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2)**0.5

def mean_or_none(vals):
    vals=[v for v in vals if v is not None]
    return sum(vals)/len(vals) if vals else None

# ---------- choose apo baseline ----------
apo_ref="ref_0"
apo_state=None
if cmd.count_atoms(apo_ref):
    for s in range(1,cmd.count_states(apo_ref)+1):
        if cmd.count_atoms(f"model {apo_ref} and ({bulge_loop})",state=s):
            apo_state=s
            break
if apo_state is None:
    print("[WARN] Could not find bulge+loop atoms in any state of ref_0; M1 RMSD will be NA")

# ---------- which objects to include ----------
objects=[o for o in cmd.get_names("objects") if o.startswith("ref_") and cmd.count_atoms(o)]

# ================== M1–M4 metrics ==================
metrics_path=os.path.join(os.getcwd(),"TAR_metrics_by_state.tsv")
with open(metrics_path,"w",encoding="utf-8") as f:
    f.write("ref\tclass\tstate\tM1_bulgeLoop_RMSD(A)\tM2_bend_angle(deg)\tM3_majorGroove(A)\tM4_chi22\tM4_chi23\tM4_chi24\tM4_chi25\tM4_chi39\tM4_chi40\tM4_chi41\tM4_chi42\tM4_chi43\n")
    for ref in objects:
        rclass=CLASS.get(ref,"unknown")
        S=cmd.count_states(ref)
        for st in range(1,S+1):
            if cmd.count_atoms(f"model {ref} and ({bulge_loop})",state=st)==0:
                continue

            # --- M1 bulge+loop RMSD vs apo baseline ---
            rms=None
            if apo_state:
                if ref == apo_ref:
                    rms = 0.0  # skip self-alignment for apo
                else:
                    sel_mobile = f"model {ref} and ({bulge_loop})"
                    sel_target = f"model {apo_ref} and ({bulge_loop})"
                    mob_n = cmd.count_atoms(sel_mobile, state=st)
                    tgt_n = cmd.count_atoms(sel_target, state=apo_state)
                    if mob_n > 0 and tgt_n > 0:
                        try:
                            rms = cmd.super(
                                sel_mobile, sel_target,
                                mobile_state=st, target_state=apo_state,
                                quiet=1, transform=0, cycles=5
                            )[0]
                        except Exception as e:
                            print(f"[WARN] super() failed for {ref} state {st} vs {apo_ref} state {apo_state}: {e}")
                            rms = None
                    else:
                        print(f"[WARN] Empty alignment selection: {ref} state {st} (n={mob_n}) vs {apo_ref} state {apo_state} (n={tgt_n})")

            # --- M2 inter-stem bend angle ---
            lowC=avg_P_centroid(ref,st,lower_stem_P)
            upC=avg_P_centroid(ref,st,upper_stem_P)
            bend=None
            if lowC and upC:
                mid=[(lowC[i]+upC[i])/2 for i in range(3)]
                vL=[lowC[i]-mid[i] for i in range(3)]
                vU=[upC[i]-mid[i] for i in range(3)]
                bend=angle_between(vL,vU)

            # --- M3 major-groove width ---
            groove=mean_or_none([dist(atom_coord(ref,st,a,'P'),atom_coord(ref,st,b,'P')) for a,b in groove_pairs])

            # --- M4 chi angles ---
            chi_resis=[22,23,24,25,39,40,41,42,43]
            chis={r:dihedral_for_resi(ref,st,r) for r in chi_resis}

            row=[ref,rclass,str(st),
                 f"{rms:.3f}" if rms is not None else "NA",
                 f"{bend:.2f}" if bend is not None else "NA",
                 f"{groove:.2f}" if groove is not None else "NA"
                ]+[f"{chis[r]:.1f}" if chis[r] is not None else "NA" for r in chi_resis]
            f.write("\t".join(row)+"\n")

# ================== base-stacking geometry ==================
PAIRS=[(23,27),(24,40),(25,39)]

def base_centroid(obj,st,resi):
    m=cmd.get_model(f"model {obj} and resi {resi} and name C1'+C2+C4+C5+C6+C8+N1+N3+N7+N9+O2+O6",state=st)
    if not m.atom: return None
    x,y,z=[sum(a.coord[i] for a in m.atom)/len(m.atom) for i in range(3)]
    return (x,y,z)

def base_normal(obj,st,resi):
    m=cmd.get_model(f"model {obj} and resi {resi} and name C1'+C2+C4+C5+C6+C8+N1+N3+N7+N9",state=st)
    if len(m.atom)<3: return None
    p1,p2,p3=[a.coord for a in m.atom[:3]]
    v1=[p2[i]-p1[i] for i in range(3)]
    v2=[p3[i]-p1[i] for i in range(3)]
    n=[v1[1]*v2[2]-v1[2]*v2[1],v1[2]*v2[0]-v1[0]*v2[2],v1[0]*v2[1]-v1[1]*v2[0]]
    ln=(n[0]**2+n[1]**2+n[2]**2)**0.5
    if ln==0: return None
    return [c/ln for c in n]

def plane_angle(n1,n2):
    if not n1 or not n2: return None
    dot=sum(n1[i]*n2[i] for i in range(3))
    dot=max(-1.0,min(1.0,abs(dot)))
    return math.degrees(math.acos(dot))

stack_path=os.path.join(os.getcwd(),"TAR_base_stacking.tsv")
with open(stack_path,"w",encoding="utf-8") as g:
    g.write("ref\tstate\tpair\tplane_angle_deg\tcentroid_dist_A\n")
    for ref in objects:
        for st in range(1,cmd.count_states(ref)+1):
            for a,b in PAIRS:
                ca,cb=base_centroid(ref,st,a),base_centroid(ref,st,b)
                na,nb=base_normal(ref,st,a),base_normal(ref,st,b)
                ang=plane_angle(na,nb) if na and nb else None
                dd=dist(ca,cb)
                if (ang is not None) or (dd is not None):
                    g.write(f"{ref}\t{st}\t{a}-{b}\t{(f'{ang:.1f}' if ang is not None else 'NA')}\t{(f'{dd:.2f}' if dd is not None else 'NA')}\n")

# ================== ion proximity ==================
IONS=["MG","K","NA"]; RADIUS=4.0
ion_path=os.path.join(os.getcwd(),"TAR_ion_counts.tsv")
with open(ion_path,"w",encoding="utf-8") as h:
    h.write("ref\tstate\tMG_within4\tK_within4\tNA_within4\n")
    for ref in objects:
        for st in range(1,cmd.count_states(ref)+1):
            row=[ref,str(st)]
            for ion in IONS:
                n=cmd.count_atoms(f"model {ref} and resn {ion} within {RADIUS} of (model {ref} and ({bulge_loop}))",state=st)
                row.append(str(n))
            h.write("\t".join(row)+"\n")

print("\n✅ Finished. Files written to:")
print(" ", os.path.join(os.getcwd(),"TAR_metrics_by_state.tsv"))
print(" ", os.path.join(os.getcwd(),"TAR_base_stacking.tsv"))
print(" ", os.path.join(os.getcwd(),"TAR_ion_counts.tsv"))
python end
