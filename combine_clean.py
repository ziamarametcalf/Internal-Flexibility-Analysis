import pandas as pd
import glob
import os

# Folder where your cleaned files are located
folder = os.path.dirname(__file__)  # same folder as this script

# Find all cleaned torsion CSVs in this folder
files = glob.glob(os.path.join(folder, '*.torsion360_clean.csv'))

if not files:
    raise FileNotFoundError(f"No .torsion360_clean.csv files found in {folder}")

# Combine all CSVs into one DataFrame
dfs = []
for f in files:
    df = pd.read_csv(f)
    df.insert(0, "source_file", os.path.basename(f))  # keep filename info
    dfs.append(df)

# ✅ This is the correct line:
combined = pd.concat(dfs, ignore_index=True)

# Save combined master file in the same folder
out_path = os.path.join(folder, 'master_torsion360.csv')
combined.to_csv(out_path, index=False)

print(f"✅ Combined {len(files)} files into:\n{out_path}")
