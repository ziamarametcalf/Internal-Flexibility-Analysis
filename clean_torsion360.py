# clean_torsion360.py
import re, glob, pathlib
import pandas as pd
from io import StringIO

# Run this script from the folder that contains *.torsion360.txt
ROOT = pathlib.Path(".")

HEADER_PAT = re.compile(r"(?m)^\s*nt\s+(?:id\s+res\s+)?alpha\s+beta\s+gamma\s+delta\s+epsilon\s+zeta", re.I)

def extract_main_table_text(text: str) -> str:
    m = HEADER_PAT.search(text)
    if not m:
        raise ValueError("Couldn't find torsion table header (line starting with 'nt ... alpha').")
    start = m.start()

    # End of table: a line of ***** or the next titled section (robust fallback to *****)
    star_rel = re.search(r"(?m)^\*{5,}\s*$", text[start:])
    stop = start + (star_rel.start() if star_rel else len(text))
    return text[start:stop].strip()

def parse_table_to_df(tbl: str, file_code: str) -> pd.DataFrame:
    # collapse 2+ spaces to delimiter
    df = pd.read_csv(StringIO(tbl), sep=r"\s{2,}", engine="python")
    df.columns = [c.strip().lower() for c in df.columns]

    # normalize column names we care about
    # possible sets: nt, id, res, alpha...chi, phase-angle
    want_order = ["nt","id","res","alpha","beta","gamma","delta","epsilon","zeta","chi","phase-angle"]
    have = [c for c in want_order if c in df.columns]
    df = df[have]

    # split nt like "A.C19" into chain/resn/resi
    def parse_nt(s):
        s = str(s)
        m = re.search(r"(?P<chain>[^.\s]+)\.(?P<resn>[ACGU])(?P<resi>\d+)", s)
        if m:
            return m.group("chain"), m.group("resn"), int(m.group("resi"))
        m2 = re.search(r"(?P<resn>[ACGU])(?P<resi>\d+)", s)
        if m2:
            return None, m2.group("resn"), int(m2.group("resi"))
        return None, None, None

    chain,resn,resi = zip(*df["nt"].map(parse_nt))
    df.insert(0, "chain", chain)
    df.insert(1, "resn",  resn)
    df.insert(2, "resi",  resi)
    df.drop(columns=["nt"], inplace=True)

    # clean numeric columns (strip annotations like "(anti)" or "(C3'-endo)")
    num_cols = [c for c in df.columns if c in {"alpha","beta","gamma","delta","epsilon","zeta","chi","phase-angle"}]
    for c in num_cols:
        df[c] = (
            df[c].astype(str)
                 .str.extract(r"(-?\d+(?:\.\d+)?)")[0]
                 .astype(float)
        )

    # provenance
    df.insert(0, "model", file_code)
    return df

def main():
    files = sorted(glob.glob(str(ROOT / "*.torsion360.txt")))
    if not files:
        print("No '*.torsion360.txt' files found. Run this in the folder with those files.")
        return

    all_df = []
    for f in files:
        path = pathlib.Path(f)
        text = path.read_text(encoding="utf-8", errors="ignore")
        try:
            table = extract_main_table_text(text)
            # model code from filename like 1anr_model01.torsion360.txt → 1anr_model01
            code = path.name.replace(".torsion360.txt","")
            df = parse_table_to_df(table, code)
            out_csv = path.with_name(f"{code}.torsion360_clean.csv")
            df.to_csv(out_csv, index=False)
            print(f"✓ {code}: {len(df)} rows -> {out_csv.name}")
            all_df.append(df)
        except Exception as e:
            print(f"× Problem with {path.name}: {e}")

    if all_df:
        big = pd.concat(all_df, ignore_index=True)
        big.to_csv(ROOT / "ALLMODELS_torsion360_clean.csv", index=False)
        print(f"\nCombined -> ALLMODELS_torsion360_clean.csv ({len(big)} total rows)")

if __name__ == "__main__":
    main()
