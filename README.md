# Open Babel Web Studio

A complete, production-grade web application bringing the power of Open Babel and RDKit to modern browsers with full 2D/3D molecular visualization, format conversion, validation auditing, and property calculations.

---

## 1. Architecture Overview

```
                      +---------------------------------------+
                      |       Modern Web Client (UI)          |
                      |  - Vanilla ES6 JavaScript Modules     |
                      |  - 3Dmol.js WebGL (Interactive 3D)    |
                      |  - SVG / Canvas Rendering (2D)        |
                      |  - Atom/Bond Interactive Sync Table   |
                      +-------------------+-------------------+
                                          |
                                HTTP / REST API
                                          |
                      +-------------------v-------------------+
                      |      FastAPI Backend Engine           |
                      |          (Python 3.12)                |
                      +---------+-------------------+---------+
                                |                   |
                 +--------------v-----+       +-----v---------------+
                 |  Open Babel 3.1.1  |       |   RDKit 2026.03.6   |
                 | - CLI Subprocess   |       | - 2D Depiction (SVG)|
                 | - 18+ File Formats |       | - ETKDG v3 3D Embed |
                 | - Clean Ring Gen   |       | - MMFF94 ForceField |
                 | - Format Auditing  |       | - Property Descript |
                 +--------------------+       +---------------------+
```

- **Frontend**: Clean, responsive, dark/light themed UI built with standard web technologies. Integrates local offline-capable `3Dmol.js` (WebGL) for both small ligands (Ball & Stick, Stick, Space-Filling, Wireframe) and biological macromolecules (Cartoon, Ribbon, Surfaces, Chains, B-factor coloring). Two-way atom picking synchronizes the 3D viewport with an atom inspection table.
- **Backend**: Built with **FastAPI** (`backend/main.py`), offering endpoints for format registry, parsing, conversions, comparative validation, structural manipulation, 2D/3D coordinate generation, batch conversion with ZIP packaging, and pre-packaged sample files.
- **Conversion Core**: Authentic Open Babel 3.1.1 CLI engine executed via safe subprocesses with argument arrays (preventing shell injection), timeout guards, isolated temporary directories, and strict error handling.
- **Cheminformatics Enhancements**: Integrated RDKit engine for robust 2D SVG/PNG depiction, stereochemical preservation, exact monoisotopic mass and Lipinski descriptor calculations, and ETKDG 3D coordinate generation fallback (mitigating Open Babel 3.1.1 Windows ring geometry bugs).

---

## 2. Supported Formats Capability Matrix

| Format | Extension | Read | Write | 2D Support | 3D Support | Ligand Support | Macromolecule Support | Description / Notes |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---|
| **SMILES** | `.smi`, `.smiles` | Yes | Yes | Yes | Via Gen3D | Yes | Limited | Simplified Molecular Input Line Entry System |
| **Canonical SMILES** | `.can` | Yes | Yes | Yes | Via Gen3D | Yes | Limited | Canonical representation for molecule indexing |
| **SDF / MDL Molfile** | `.sdf`, `.mol`, `.sd` | Yes | Yes | Yes | Yes | Yes | Limited | V2000 / V3000 format with explicit bonds & stereo |
| **Tripos MOL2** | `.mol2` | Yes | Yes | Yes | Yes | Yes | Yes | Full atom types, SYBYL bonds, partial charges |
| **Protein Data Bank** | `.pdb`, `.ent` | Yes | Yes | Yes | Yes | Yes | Yes | Macromolecular structures, chains, residues, B-factors |
| **AutoDock PDBQT** | `.pdbqt` | Yes | Yes | Yes | Yes | Yes | Limited | PDB format with AutoDock atom types & Gasteiger charges |
| **XYZ Coordinate** | `.xyz` | Yes | Yes | Via Gen2D | Yes | Yes | Limited | Simple Cartesian coordinates (no explicit bonds) |
| **Crystallographic CIF** | `.cif` | Yes | Yes | Via Gen2D | Yes | Yes | Yes | Crystallographic Information Framework |
| **Macromolecular CIF** | `.mcif`, `.mmcif` | Yes | Yes | Via Gen2D | Yes | Yes | Yes | PDBx/mmCIF standard for high-res macromolecules |
| **Chemical Markup Language**| `.cml` | Yes | Yes | Yes | Yes | Yes | Yes | XML-based chemical structure specification |
| **ChemDraw CDX** | `.cdx` | Yes | No | Yes | Via Gen3D | Yes | No | ChemDraw binary format (read-only input) |
| **MOPAC Format** | `.mop`, `.mopout` | Yes | Yes | Via Gen2D | Yes | Yes | No | Semi-empirical quantum chemistry input/output |
| **Gaussian Input** | `.gjf`, `.com` | Yes | Yes | Via Gen2D | Yes | Yes | No | Gaussian Z-matrix / Cartesian input format |
| **GROMACS Structure** | `.gro` | Yes | Yes | Via Gen2D | Yes | Yes | Yes | GROMACS molecular dynamics coordinates (nm converted) |
| **FASTA Sequence** | `.fasta`, `.fa` | Yes | Yes | No | Via Model | Yes (peptides) | Yes | Amino acid / nucleic acid sequence format |
| **InChI** | `.inchi` | Yes | Yes | Yes | Via Gen3D | Yes | No | IUPAC International Chemical Identifier |
| **InChIKey** | `.inchikey` | Yes | No | No | No | Yes | No | 27-character hashed InChI key (lookup identifier) |

---

## 3. Installation & Getting Started

### Prerequisites
- Python 3.10+ (tested on Python 3.12)
- Open Babel 3.1.1 installed locally (e.g. `D:\OpenBabel-3.1.1` or on system `PATH`)

### Virtual Environment Setup
```powershell
# Clone or navigate to the repository
cd "c:\Users\Anneswa Das\Downloads\openbabel mockup"

# Create virtual environment if not already present
python -m venv .venv

# Activate virtual environment
.\.venv\Scripts\activate

# Install required dependencies
pip install fastapi uvicorn pydantic rdkit pytest pytest-asyncio httpx
```

### Launching the Application
Run the FastAPI development server:
```powershell
.\.venv\Scripts\python.exe -m uvicorn backend.main:app --host 127.0.0.1 --port 8000 --reload
```

Once running:
- **Web Interface**: Open `http://localhost:8000` in your web browser.
- **Interactive API Documentation (Swagger)**: Open `http://localhost:8000/docs`.
- **Health Check Endpoint**: Open `http://localhost:8000/api/health`.

---

## 4. Automated Testing Suite

The project includes an exhaustive, multi-tier test suite of **315 automated tests** executing authentic Open Babel conversions, validations, edge-case evaluations, and API integration checks.

### Running the Tests
```powershell
# Run the complete test suite
.\.venv\Scripts\python.exe -m pytest -v

# Run only the 200+ multi-format conversion matrix tests
.\.venv\Scripts\python.exe -m pytest tests/test_matrix_200.py -v

# Run with concise summary
.\.venv\Scripts\python.exe -m pytest -q
```

### Test Suite Structure
| Test Module | Test Count | Scope |
|:---|:---:|:---|
| `tests/test_formats.py` | 21 | Heuristic detection, extension lookup, format capability matrix |
| `tests/test_parsing.py` | 24 | PDB, PDBQT, MOL2, SDF, XYZ, GJF, GRO, macromolecule detection |
| `tests/test_properties.py`| 16 | MW, exact mass, formula, charge, Lipinski rules, TPSA, LogP |
| `tests/test_conversion.py`| 28 | Authentic Open Babel conversions, stereochemistry, conformers |
| `tests/test_roundtrip.py` | 12 | PDB->MOL2->PDB, SDF->MOL2->SDF, SMI->SDF->SMI consistency |
| `tests/test_validation.py`| 9 | Comparative auditing, atom conservation, charge loss warnings |
| `tests/test_manipulation.py`| 10 | Add/remove H, 3D coordinate generation, water stripping, MMFF94 |
| `tests/test_depiction.py` | 8 | 2D vector SVG generation, raster PNG, stereochemical rendering |
| `tests/test_batch.py` | 5 | Multi-file conversions, isolated error handling, ZIP archive |
| `tests/test_errors.py` | 11 | Malformed structures, empty input, timeouts, corrupted syntax |
| `tests/test_api.py` | 19 | FastAPI endpoints, CORS, file upload, error responses |
| `tests/test_matrix_200.py`| 152 | Systematic cross-format matrix across 18 chemical classes |
| **Total** | **315** | **100% Passed** |

---

## 5. Honest Format Preservation & Limitation Analysis

The application features a built-in **Conversion Validation Auditor** (`backend/validator.py`) that transparently compares input and output structures and alerts users when information is inevitably discarded due to format specifications:

1. **PDB / MOL2 -> SMILES**:
   - *Limitation*: Ordinary SMILES syntax encodes only chemical connectivity and configuration. 3D Cartesian coordinates (`X, Y, Z`), macromolecular chain identifiers, residue numbers, B-factors, and crystal unit cell parameters are **not preserved**.
2. **SDF / MOL2 -> PDB**:
   - *Limitation*: The standard PDB format does not explicitly record bond orders or formal charges for standard residues. Downstream parsers must infer connectivity from interatomic distances unless `CONECT` records are explicitly supplied.
3. **MOL / SDF / PDB -> XYZ**:
   - *Limitation*: XYZ files contain solely element symbols and coordinates. All chemical bond orders, formal charges, aromaticity flags, and stereochemical descriptors are omitted.
4. **AutoDock PDBQT -> SDF / PDB**:
   - *Limitation*: AutoDock atom types (e.g., `A` for aromatic carbon, `OA` for hydrogen bond acceptor oxygen) and Gasteiger partial charges are mapped back to standard elements; dock scoring metadata is not retained in generic PDB.
5. **Small Molecule Ring Systems & Open Babel 3.1.1 on Windows**:
   - *Issue*: Open Babel 3.1.1 for Windows contains a known crash in `rigid-fragments.txt` when running `--gen3d` on certain aromatic or aliphatic rings.
   - *Mitigation*: Our converter automatically employs RDKit's ETKDG v3 conformer generation with MMFF94 energy minimization as a reliable 3D coordinator, ensuring robust execution.

---

## 6. Sample Workflows & Step-by-Step Instructions

### Workflow 1: Docking Ligand (PDBQT -> MOL2 / SDF / SMILES)
1. Open the application at `http://localhost:8000`.
2. Click **Load Sample Molecule** in the top bar and select **Docking Ligand (PDBQT)** (or upload `sample_data/ligand_docking.pdbqt`).
3. Notice:
   - Format is auto-detected as **PDBQT**.
   - The 3D viewer displays the ligand with Ball & Stick representation and AutoDock partial charges.
   - The **Structural Information** panel displays formula `C10H14N2` and exact molecular weight.
   - The **Atom & Bond Inspection Table** lists each atom with its Cartesian coordinates and charge.
4. Select **MOL2** in the **Convert To** dropdown and click **Convert Structure**.
5. The **Validation Audit** confirms 100% heavy atom conservation (12 heavy atoms preserved).
6. Click **Download Converted File** to obtain `ligand.mol2`.

### Workflow 2: Macromolecule / Protein Analysis (Crambin 1CRN / Streptavidin 1STP)
1. In the sample loader, select **Crambin (PDB 1CRN)** or **Streptavidin + Biotin (PDB 1STP)**.
2. The application automatically classifies the file as a **Macromolecule**:
   - Renders in **Cartoon** representation colored by **Secondary Structure** or **Chain**.
   - B-factor statistics (min, max, mean) are displayed.
   - Chains (`A`, `B`, etc.) and bound ligands (`BTN` in 1STP) are identified with dedicated visibility toggles.
3. Click an atom in the 3D viewer to select and highlight the corresponding residue in the Atom Table.
4. Use the **Remove Water Molecules** manipulation tool to strip crystallization waters and clean the structure.

### Workflow 3: Multi-Conformer / Multi-Molecule Library
1. Upload `sample_data/multimol_library.sdf` (containing multiple compounds).
2. The viewer detects multiple records and provides a molecule stepper.
3. Convert the entire library to multi-line SMILES or multi-record MOL2.

### Workflow 4: Multi-File Batch Conversion
1. Scroll to the **Batch File Conversion** section.
2. Drag and drop multiple files (e.g. `caffeine.mol`, `aspirin.sdf`, `water.xyz`).
3. Select target format (e.g., `pdb` or `mol2`).
4. Click **Start Batch Conversion**.
5. Each file is processed in isolation; once complete, click **Download Converted Package (.zip)** to obtain all converted structures.

---

## 7. License & Acknowledgements
- Powered by the **Open Babel** chemical toolbox project (http://openbabel.org).
- Powered by **RDKit** cheminformatics toolkit (https://www.rdkit.org).
- 3D rendering powered by **3Dmol.js** (https://3dmol.csb.pitt.edu).
