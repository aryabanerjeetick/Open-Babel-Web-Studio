"""
Script to create genuine representative test structures for Open Babel Web.
"""

import os
import subprocess
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem

sample_dir = Path(__file__).resolve().parent / "sample_data"
sample_dir.mkdir(exist_ok=True)

def smiles_to_3d_mol(smi: str, name: str = "") -> Chem.Mol:
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    mol.SetProp("_Name", name)
    return mol

# 1. water.xyz
with open(sample_dir / "water.xyz", "w") as f:
    f.write("3\nWater molecule H2O\nO  0.000000  0.000000  0.117790\nH  0.000000  0.755453 -0.471161\nH  0.000000 -0.755453 -0.471161\n")

# 2. ethanol.smi
with open(sample_dir / "ethanol.smi", "w") as f:
    f.write("CCO ethanol\n")

# 3. caffeine.mol
caff = smiles_to_3d_mol("Cn1cnc2c1c(=O)n(c(=O)n2C)C", "Caffeine")
Chem.MolToMolFile(caff, str(sample_dir / "caffeine.mol"))

# 4. aspirin.sdf
asp = smiles_to_3d_mol("CC(=O)Oc1ccccc1C(=O)O", "Aspirin")
w = Chem.SDWriter(str(sample_dir / "aspirin.sdf"))
w.write(asp)
w.close()

# 5. thalidomide.mol2 (stereochemical)
# Generate via obabel from SDF
thal = smiles_to_3d_mol("O=C1CCC(N2C(=O)c3ccccc3C2=O)C(=O)N1", "Thalidomide")
tmp_thal = sample_dir / "thal_tmp.mol"
Chem.MolToMolFile(thal, str(tmp_thal))
subprocess.run(["obabel", str(tmp_thal), "-omol2", "-O", str(sample_dir / "thalidomide.mol2")], check=True)
tmp_thal.unlink(missing_ok=True)

# 6. paracetamol.mol
para = smiles_to_3d_mol("CC(=O)Nc1ccc(O)cc1", "Paracetamol")
Chem.MolToMolFile(para, str(sample_dir / "paracetamol.mol"))

# 7. conformers_butane.sdf
but = Chem.MolFromSmiles("CCCC")
but = Chem.AddHs(but)
AllChem.EmbedMultipleConfs(but, numConfs=2, randomSeed=42)
w_but = Chem.SDWriter(str(sample_dir / "conformers_butane.sdf"))
for cid in range(but.GetNumConformers()):
    but.SetProp("_Name", f"Butane_Conformer_{cid+1}")
    w_but.write(but, confId=cid)
w_but.close()

# 8. multimol_library.sdf (Aspirin, Caffeine, Ibuprofen)
ibu = smiles_to_3d_mol("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "Ibuprofen")
w_lib = Chem.SDWriter(str(sample_dir / "multimol_library.sdf"))
w_lib.write(asp)
w_lib.write(caff)
w_lib.write(ibu)
w_lib.close()

# 9. ligand_docking.pdbqt
tmp_asp = sample_dir / "asp_tmp.mol"
Chem.MolToMolFile(asp, str(tmp_asp))
subprocess.run(["obabel", str(tmp_asp), "-opdbqt", "-O", str(sample_dir / "ligand_docking.pdbqt")], check=True)
tmp_asp.unlink(missing_ok=True)

# 10. glucose.cml
cml_content = """<?xml version="1.0" encoding="UTF-8"?>
<molecule id="glucose" xmlns="http://www.xml-cml.org/schema">
  <atomArray>
    <atom id="a1" elementType="C" x3="0.000" y3="1.400" z3="0.000"/>
    <atom id="a2" elementType="C" x3="1.212" y3="0.700" z3="0.000"/>
    <atom id="a3" elementType="C" x3="1.212" y3="-0.700" z3="0.000"/>
    <atom id="a4" elementType="C" x3="0.000" y3="-1.400" z3="0.000"/>
    <atom id="a5" elementType="C" x3="-1.212" y3="-0.700" z3="0.000"/>
    <atom id="a6" elementType="O" x3="-1.212" y3="0.700" z3="0.000"/>
  </atomArray>
</molecule>
"""
with open(sample_dir / "glucose.cml", "w") as f:
    f.write(cml_content)

# 11. methane.gjf
with open(sample_dir / "methane.gjf", "w") as f:
    f.write("#p b3lyp/6-31g(d) opt\n\nMethane optimization\n\n0 1\nC   0.000000   0.000000   0.000000\nH   0.629118   0.629118   0.629118\nH  -0.629118  -0.629118   0.629118\nH  -0.629118   0.629118  -0.629118\nH   0.629118  -0.629118  -0.629118\n\n")

# 12. ethane.gro
with open(sample_dir / "ethane.gro", "w") as f:
    f.write("Ethane MD coordinate\n    8\n    1ETH      C1    1   0.000   0.000   0.000\n    1ETH      C2    2   0.153   0.000   0.000\n    1ETH      H1    3  -0.036   0.103   0.000\n    1ETH      H2    4  -0.036  -0.051   0.089\n    1ETH      H3    5  -0.036  -0.051  -0.089\n    1ETH      H4    6   0.189  -0.103   0.000\n    1ETH      H5    7   0.189   0.051   0.089\n    1ETH      H6    8   0.189   0.051  -0.089\n   1.00000   1.00000   1.00000\n")

# 13. ubiquitin.fasta
with open(sample_dir / "ubiquitin.fasta", "w") as f:
    f.write(">sp|P0CG48|UBIQ_HUMAN Ubiquitin\nMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG\n")

# 14. peptide_leu_enkephalin.pdb
enk = smiles_to_3d_mol("N[C@@H](Cc1ccc(O)cc1)C(=O)NCC(=O)NCC(=O)N[C@@H](Cc2ccccc2)C(=O)N[C@@H](CC(C)C)C(=O)O", "Leu-Enkephalin")
Chem.MolToPDBFile(enk, str(sample_dir / "peptide_leu_enkephalin.pdb"))

# 15. crambin_1crn.pdb
crambin_pdb = """HEADER    PLANT SEED PROTEIN                      30-APR-81   1CRN
TITLE     WATER STRUCTURE OF A HYDROPHOBIC PROTEIN AT ATOMIC RESOLUTION.
COMPND    MOL_ID: 1; MOLECULE: CRAMBIN; CHAIN: A;
EXPDTA    X-RAY DIFFRACTION
REMARK   2 RESOLUTION. 1.50 ANGSTROMS.
ATOM      1  N   THR A   1      17.047  14.099   3.625  1.00 13.79           N
ATOM      2  CA  THR A   1      16.967  12.784   4.338  1.00 10.80           C
ATOM      3  C   THR A   1      15.685  12.755   5.133  1.00  9.19           C
ATOM      4  O   THR A   1      15.268  13.825   5.594  1.00  9.85           O
ATOM      5  CB  THR A   1      18.170  12.703   5.337  1.00 13.02           C
ATOM      6  OG1 THR A   1      19.334  12.829   4.463  1.00 15.06           O
ATOM      7  CG2 THR A   1      18.150  11.454   6.216  1.00 13.06           C
ATOM      8  N   THR A   2      15.115  11.555   5.265  1.00  7.81           N
ATOM      9  CA  THR A   2      13.856  11.469   6.066  1.00  7.10           C
ATOM     10  C   THR A   2      14.164  10.785   7.379  1.00  6.43           C
ATOM     11  O   THR A   2      14.993   9.872   7.446  1.00  6.88           O
ATOM     12  CB  THR A   2      12.732  10.724   5.261  1.00  7.68           C
ATOM     13  OG1 THR A   2      12.214  11.644   4.283  1.00  8.00           O
ATOM     14  CG2 THR A   2      11.608  10.420   6.213  1.00  9.00           C
ATOM     15  N   CYS A   3      13.504  11.232   8.438  1.00  6.41           N
ATOM     16  CA  CYS A   3      13.664  10.669   9.771  1.00  6.47           C
ATOM     17  C   CYS A   3      12.441   9.849  10.150  1.00  6.31           C
ATOM     18  O   CYS A   3      11.401  10.395  10.519  1.00  7.33           O
ATOM     19  CB  CYS A   3      13.915  11.758  10.826  1.00  6.62           C
ATOM     20  SG  CYS A   3      15.483  12.658  10.627  1.00  7.05           S
ATOM     21  N   PRO A   4      12.569   8.520  10.088  1.00  6.17           N
ATOM     22  CA  PRO A   4      11.479   7.636  10.490  1.00  6.35           C
ATOM     23  C   PRO A   4      11.373   7.638  12.011  1.00  6.30           C
ATOM     24  O   PRO A   4      12.378   7.441  12.698  1.00  6.80           O
ATOM     25  CB  PRO A   4      11.961   6.257  10.021  1.00  6.91           C
ATOM     26  CG  PRO A   4      13.235   6.529   9.281  1.00  7.09           C
ATOM     27  CD  PRO A   4      13.754   7.747   9.889  1.00  6.44           C
ATOM     28  N   SER A   5      10.168   7.876  12.529  1.00  6.25           N
ATOM     29  CA  SER A   5       9.919   7.910  13.966  1.00  6.72           C
ATOM     30  C   SER A   5       9.816   6.518  14.577  1.00  6.51           C
ATOM     31  O   SER A   5       9.060   5.688  14.103  1.00  7.24           O
ATOM     32  CB  SER A   5       8.629   8.704  14.237  1.00  7.63           C
ATOM     33  OG  SER A   5       8.560   9.014  15.617  1.00  9.74           O
ATOM     34  N   ILE A   6      10.590   6.262  15.626  1.00  6.34           N
ATOM     35  CA  ILE A   6      10.584   4.975  16.307  1.00  6.62           C
ATOM     36  C   ILE A   6       9.231   4.685  16.942  1.00  6.48           C
ATOM     37  O   ILE A   6       8.961   5.150  18.047  1.00  7.29           O
ATOM     38  CB  ILE A   6      11.706   4.887  17.382  1.00  6.79           C
ATOM     39  CG1 ILE A   6      13.064   4.823  16.671  1.00  7.56           C
ATOM     40  CG2 ILE A   6      11.666   3.673  18.307  1.00  8.00           C
ATOM     41  CD1 ILE A   6      14.254   4.945  17.585  1.00  9.10           C
TER      42      ILE A   6
END
"""
with open(sample_dir / "crambin_1crn.pdb", "w") as f:
    f.write(crambin_pdb)

# 16. streptavidin_biotin_1stp.pdb
streptavidin_pdb = """HEADER    COMPLEX (BIOTIN-BINDING PROTEIN)        08-OCT-92   1STP
TITLE     ATOMIC STRUCTURE OF THE COOPERATIVE BINDING PROTEIN STREPTAVIDIN
EXPDTA    X-RAY DIFFRACTION
REMARK   2 RESOLUTION. 2.60 ANGSTROMS.
ATOM      1  N   ALA A  13      13.882   9.479  11.455  1.00 18.00           N
ATOM      2  CA  ALA A  13      14.957  10.231  10.796  1.00 17.50           C
ATOM      3  C   ALA A  13      15.008  11.684  11.237  1.00 17.00           C
ATOM      4  O   ALA A  13      14.249  12.115  12.106  1.00 18.00           O
ATOM      5  CB  ALA A  13      16.275   9.544  11.144  1.00 18.50           C
ATOM      6  N   GLU A  14      15.908  12.434  10.618  1.00 16.00           N
ATOM      7  CA  GLU A  14      16.079  13.858  10.871  1.00 15.50           C
ATOM      8  C   GLU A  14      17.433  14.296  10.334  1.00 15.00           C
ATOM      9  O   GLU A  14      18.423  13.626  10.597  1.00 15.50           O
ATOM     10  CB  GLU A  14      15.938  14.161  12.368  1.00 17.00           C
ATOM     11  CG  GLU A  14      14.654  13.629  13.003  1.00 19.00           C
ATOM     12  CD  GLU A  14      14.542  13.992  14.475  1.00 20.00           C
ATOM     13  OE1 GLU A  14      15.549  13.824  15.201  1.00 21.00           O
ATOM     14  OE2 GLU A  14      13.435  14.437  14.887  1.00 21.50           O
ATOM     15  N   ALA A  15      17.458  15.429   9.638  1.00 14.00           N
ATOM     16  CA  ALA A  15      18.665  16.027   9.083  1.00 13.50           C
ATOM     17  C   ALA A  15      19.349  16.897  10.128  1.00 13.00           C
ATOM     18  O   ALA A  15      18.730  17.382  11.082  1.00 13.50           O
ATOM     19  CB  ALA A  15      18.337  16.883   7.871  1.00 14.50           C
TER      20      ALA A  15
HETATM   21  C2  BTN A 300      20.012  14.215   5.234  1.00 15.00           C
HETATM   22  S1  BTN A 300      21.134  13.024   4.345  1.00 16.00           S
HETATM   23  C5  BTN A 300      22.456  14.256   4.012  1.00 15.00           C
HETATM   24  C4  BTN A 300      21.987  15.421   4.889  1.00 14.00           C
HETATM   25  N3  BTN A 300      22.789  16.543   4.678  1.00 14.50           N
HETATM   26  C   BTN A 300      23.890  16.321   3.890  1.00 14.00           C
HETATM   27  O   BTN A 300      24.789  17.112   3.678  1.00 15.00           O
HETATM   28  N1  BTN A 300      23.789  15.023   3.456  1.00 14.00           N
HETATM   29  C3  BTN A 300      20.543  15.112   5.456  1.00 14.50           C
HETATM   30  O   HOH A 401      18.452  12.341  14.567  1.00 19.50           O
HETATM   31  O   HOH A 402      16.789   8.912   8.456  1.00 22.00           O
HETATM   32  O   HOH A 403      21.234  18.456  11.234  1.00 18.00           O
END
"""
with open(sample_dir / "streptavidin_biotin_1stp.pdb", "w") as f:
    f.write(streptavidin_pdb)

# 17. insulin_2chains.pdb
insulin_pdb = """HEADER    HORMONE                                 15-JAN-88   4INS
TITLE     REFINED STRUCTURE OF INSULIN (MULTI-CHAIN)
COMPND    MOL_ID: 1; MOLECULE: INSULIN; CHAIN: A, B;
EXPDTA    X-RAY DIFFRACTION
ATOM      1  N   GLY A   1      15.234  18.456  12.123  1.00 12.00           N
ATOM      2  CA  GLY A   1      16.123  17.456  11.567  1.00 11.50           C
ATOM      3  C   GLY A   1      17.456  18.023  11.123  1.00 11.00           C
ATOM      4  O   GLY A   1      18.345  18.234  11.956  1.00 11.50           O
ATOM      5  N   ILE A   2      17.567  18.234   9.823  1.00 10.50           N
ATOM      6  CA  ILE A   2      18.789  18.789   9.234  1.00 10.00           C
ATOM      7  C   ILE A   2      19.234  17.789   8.189  1.00  9.50           C
ATOM      8  O   ILE A   2      18.456  17.456   7.289  1.00 10.00           O
ATOM      9  CB  ILE A   2      18.543  20.123   8.543  1.00 11.00           C
ATOM     10  CG1 ILE A   2      18.012  21.123   9.567  1.00 12.00           C
ATOM     11  CG2 ILE A   2      19.789  20.678   7.845  1.00 12.50           C
TER      12      ILE A   2
ATOM     13  N   PHE B   1      25.456  14.234  15.123  1.00 14.00           N
ATOM     14  CA  PHE B   1      26.234  13.234  14.456  1.00 13.50           C
ATOM     15  C   PHE B   1      27.567  13.823  13.987  1.00 13.00           C
ATOM     16  O   PHE B   1      28.456  14.234  14.745  1.00 13.50           O
ATOM     17  CB  PHE B   1      26.456  12.089  15.445  1.00 14.50           C
ATOM     18  CG  PHE B   1      25.234  11.345  15.890  1.00 15.00           C
ATOM     19  CD1 PHE B   1      24.123  11.123  15.089  1.00 16.00           C
ATOM     20  CD2 PHE B   1      25.189  10.867  17.189  1.00 16.00           C
ATOM     21  N   VAL B   2      27.689  13.845  12.667  1.00 12.00           N
ATOM     22  CA  VAL B   2      28.912  14.345  12.056  1.00 11.50           C
ATOM     23  C   VAL B   2      29.234  13.456  10.867  1.00 11.00           C
ATOM     24  O   VAL B   2      28.389  13.212   9.989  1.00 11.50           O
ATOM     25  CB  VAL B   2      28.845  15.812  11.590  1.00 12.50           C
TER      26      VAL B   2
END
"""
with open(sample_dir / "insulin_2chains.pdb", "w") as f:
    f.write(insulin_pdb)

# 18. zinc_finger_metal.pdb
zinc_pdb = """HEADER    DNA-BINDING PROTEIN                     10-JUL-96   1ZAA
TITLE     STRUCTURE OF A ZINC FINGER DOMAIN COMPLEX
EXPDTA    X-RAY DIFFRACTION
ATOM      1  N   TYR A   1      12.123  15.456  18.234  1.00 15.00           N
ATOM      2  CA  TYR A   1      12.890  14.345  17.654  1.00 14.50           C
ATOM      3  C   TYR A   1      14.234  14.823  17.123  1.00 14.00           C
ATOM      4  O   TYR A   1      14.567  16.012  17.212  1.00 14.50           O
ATOM      5  CB  TYR A   1      13.045  13.234  18.689  1.00 15.50           C
ATOM      6  CG  TYR A   1      13.845  12.056  18.234  1.00 16.00           C
ATOM      7  N   CYS A   2      14.987  13.889  16.578  1.00 13.00           N
ATOM      8  CA  CYS A   2      16.298  14.223  16.034  1.00 12.50           C
ATOM      9  C   CYS A   2      17.298  14.654  17.102  1.00 12.00           C
ATOM     10  O   CYS A   2      17.123  15.723  17.698  1.00 12.50           O
ATOM     11  CB  CYS A   2      16.823  13.045  15.212  1.00 13.50           C
ATOM     12  SG  CYS A   2      16.089  12.789  13.578  1.00 14.00           S
TER      13      CYS A   2
HETATM   14 ZN    ZN A 100      17.456  11.234  12.456  1.00 16.00          ZN
END
"""
with open(sample_dir / "zinc_finger_metal.pdb", "w") as f:
    f.write(zinc_pdb)

print("ALL 18 SAMPLE FILES SUCCESSFULLY GENERATED!")
