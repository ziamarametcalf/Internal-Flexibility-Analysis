# HIV-1 TAR RNA Torsion-Angle and Conformational Analysis

This project analyzes conformational flexibility in HIV-1 TAR RNA across multiple PDB ensembles (1ANR, 6MCE, 2KDQ, 9DE5). Torsion angles (α, β, γ, δ, ε, ζ, χ) are extracted using DSSR, aligned using PyMOL, and visualized using Python (rose plots, RMSD heatmaps, and torsion distributions).

## Project Goals
- Compare torsion-angle flexibility across bulge/loop regions (U25, C24, A35, etc.)
- Identify structural patterns relevant for peptide binding (Tat, L22, MV2003)
- Build analyzable CSV datasets for future peptide design work

## Repository Structure
main=codes
pdb=files and trimmed ones
results=results and data

## 🛠️ Tools Used
- **DSSR / x3DNA** – Backbone torsion extraction
- **PyMOL** – Structure alignment & visualization
- **Python (NumPy, Pandas, Matplotlib)** – Data analysis & plots
- **Excel** – Pivot tables for angle statistics

## ▶️ How to Run
1. Install dependencies  
   ```bash
   pip install numpy pandas matplotlib

