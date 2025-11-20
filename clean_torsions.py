# clean_torsions.py
import re, glob, pathlib
import pandas as pd

# ---- change this to your torsions folder OR leave as "." and run from the folder
ROOT = pathlib.Path(r".")


def extract_main_table_text(text: str) -> str:
    """Return the main-chain torsion table text including header line starting with 'nt'."""
    # start at the header line that contains 'nt' and 'alpha'
    start_match = re.search(r"(?m)^\s*nt\s+alpha\s+beta\s+gamma\s+delta\s+epsilon\s+zeta", text)
    if not start_match:
        raise ValueError("Couldn't locate main table header (line starting with 'nt ... alpha').")
    start = start_match.start()

    # stop at the next section (Virtual eta/theta) or a line of asterisks
    stop_match = re.search(r"\*{10,}\s*\n\s*Virtual eta/theta", text)
    if not stop_match:
        # fallback: stop at first line of 10+ asterisks after start
        stop_match = re.search(r"(?s)\*{10,}.*?$", text[start:])
        stop = start + (stop_match.start() if stop_match else len(text))
    else:
        stop = stop_match.start()

    return text[start:stop].strip()


def parse_table_to_df(tbl: str, file_code: str) -> pd.DataFrame:
    """Parse DSSR main-chain table (whitespace aligned) to a DataFrame with clean numeric columns."""
    # compress multiple spaces to two+ spaces so pandas can split; keep columns with >=2 spaces
    df = pd.read_csv(pd.io.common.StringIO(tbl), sep=r"\s{2,}", engine="python")
    df.columns = [c.strip().lower() for c in df.columns]

    # keep only expected columns if present
    want = ["nt", "alpha", "beta", "gamma", "delta", "epsilon", "zeta", "chi", "phase-angle"]
    have = [c for c in want if c in df.columns]
    df = df[have]

    # split nt like "A.C19" (chain=A, resn=C, resi=19)
    def parse_nt(s):
        # handles things like "A.C19", "A.U23", "B.G5" etc.
        m = re.search(r"(?P<chain>[^\.]+)\.(?P<resn>[ACGU])(?P<resi>\d+)", str(s))
        if m:
            return m.group("chain"), m.group("resn"), int(m.group("resi"))
        # fallback: try just "C19"
        m2 = re.search(r"(?P<resn>[ACGU])(?P<resi>\d+)", str(s))
        if m2:
            return None, m2.group("resn"), int(m2.group("resi"))
        return None, None, None

    parsed = df["nt"].apply(parse_nt).apply(pd.Series)
    parsed.columns = ["chain", "resn", "resi"]
    df = pd.concat([parsed, df.drop(columns=["nt"])], axis=1)

    # remove text like "(anti)" "(C3'-endo)" from numeric columns
    num_cols = [c for c in df.columns if c in {"alpha","beta","gamma","delta","epsilon","zeta","chi","phase-angle"}]
    for c in num_cols:
        df[c] = (
            df[c]
            .astype(str)
            .str.extract(r"(-?\d+(?:\.\d+)?)")  # first number
            .astype(float)
        )

    # add pdb code column for provenance
    df.insert(0, "pdb", file_code.lower())
    return df


def clean_one(path: pathlib.Path) -> pd.DataFrame:
    text = path.read_text(encoding="utf-8", errors="ignore")
    table = extract_main_table_text(text)
    code = path.name.split("-torsions")[0]  # e.g., "1anr"
    df = parse_table_to_df(table, code)
    # write per-file csv next to the source
    out = path.with_name(f"{code}-torsions_clean.csv")
    df.to_csv(out, index=False)
    print(f"✓ {code}: {len(df)} rows -> {out.name}")
    return df


def main():
    files = sorted(glob.glob(str(ROOT / "*-torsions.txt")))
    if not files:
        print("No '*-torsions.txt' files found. Adjust ROOT or run from the torsions folder.")
        return

    all_df = []
    for f in files:
        try:
            all_df.append(clean_one(pathlib.Path(f)))
        except Exception as e:
            print(f"× Problem with {f}: {e}")

    if all_df:
        big = pd.concat(all_df, ignore_index=True)
        big.to_csv(ROOT / "all_torsions_clean.csv", index=False)
        print(f"\nCombined table -> all_torsions_clean.csv ({len(big)} total rows)")


if __name__ == "__main__":
    main()
