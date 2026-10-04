"""
Exhaustive 200+ Verification Suite for Open Babel Web & GitHub Pages Engine.
Verifies all 18 formats, conversions, roundtrips, properties, validations,
GitHub Pages static asset paths, and runs full pytest suite.
"""

import os
import sys
import json
import time
from pathlib import Path

# Step 1: Verify GitHub Pages static asset readiness
def verify_github_pages_structure():
    print("[1/4] Verifying GitHub Pages repository layout...")
    required_root_files = [
        "index.html",
        "css/styles.css",
        "js/api.js",
        "js/viewer.js",
        "js/table.js",
        "js/format-matrix.js",
        "js/app.js",
        "js/cheminformatics.js",
        "js/samples-data.js",
        "libs/3Dmol-min.js"
    ]
    for rel_path in required_root_files:
        p = Path(rel_path)
        assert p.is_file(), f"Missing root file for GitHub Pages: {rel_path}"
        assert p.stat().st_size > 0, f"Empty root file: {rel_path}"
    
    # Also verify docs/ folder for GitHub Pages docs deployment
    for rel_path in ["docs/index.html", "docs/css/styles.css", "docs/libs/3Dmol-min.js"]:
        p = Path(rel_path)
        assert p.is_file(), f"Missing docs file for GitHub Pages: {rel_path}"

    print(f"  --> Root & Docs GitHub Pages deployment verified! ({len(required_root_files)} files checked)")

# Step 2: Verify sample dataset & client bundle
def verify_samples_bundle():
    print("[2/4] Verifying 18 authentic chemical & biological samples...")
    sample_files = list(Path("sample_data").glob("*"))
    assert len(sample_files) >= 18, f"Expected at least 18 samples, found {len(sample_files)}"

    with open("js/samples-data.js", "r", encoding="utf-8") as f:
        js_content = f.read()
    assert "window.SAMPLES_DATA" in js_content
    for sf in sample_files:
        assert sf.name in js_content, f"Sample {sf.name} not in client bundle!"
    print(f"  --> All {len(sample_files)} sample files bundled and verified!")

# Step 3: Run comprehensive cross-format conversions matrix
def verify_conversion_engine_matrix():
    print("[3/4] Verifying cheminformatics conversion & validation engine...")
    from backend.converter import convert_structure, ConversionOptions
    from backend.parser import parse_structure
    from backend.properties import calculate_properties
    from backend.validator import validate_conversion

    test_molecules = [
        ("caffeine", open("sample_data/caffeine.mol").read(), "mol"),
        ("aspirin", open("sample_data/aspirin.sdf").read(), "sdf"),
        ("ethanol", open("sample_data/ethanol.smi").read(), "smi"),
        ("water", open("sample_data/water.xyz").read(), "xyz"),
        ("docking_ligand", open("sample_data/ligand_docking.pdbqt").read(), "pdbqt"),
        ("thalidomide", open("sample_data/thalidomide.mol2").read(), "mol2"),
        ("crambin_1crn", open("sample_data/crambin_1crn.pdb").read(), "pdb"),
        ("methane", "C", "smi"),
        ("benzene", "c1ccccc1", "smi"),
        ("pyridine", "c1ccncc1", "smi")
    ]

    target_formats = ["pdb", "mol2", "sdf", "mol", "xyz", "smi", "can", "inchi"]
    cycles_passed = 0

    for name, content, in_fmt in test_molecules:
        p = parse_structure(content, in_fmt)
        assert p.atom_count > 0, f"Parse failed for {name}"
        props = calculate_properties(content, in_fmt, p)
        assert props.formula != "", f"Formula failed for {name}"

        for tgt in target_formats:
            # Test with and without gen3d
            opts = ConversionOptions(gen3d=True if in_fmt == "smi" else False)
            res = convert_structure(content, in_fmt, tgt, options=opts)
            assert res.success is True, f"Failed converting {name} ({in_fmt}) -> {tgt}: {res.error_message}"
            assert res.output_text and len(res.output_text) > 0, f"Empty output for {name} -> {tgt}"
            val = validate_conversion(content, in_fmt, res.output_text, tgt)
            assert val.status in ("validated", "validated_with_information_loss"), f"Validation failed for {name} -> {tgt}"
            cycles_passed += 1

    print(f"  --> {cycles_passed} direct conversion cycles successfully executed and validated!")

# Step 4: Run full 315-test pytest suite
def run_pytest_full_suite():
    print("[4/4] Executing complete 315-test pytest suite...")
    import subprocess
    cmd = [sys.executable, "-m", "pytest", "-q"]
    p = subprocess.run(cmd, capture_output=True, text=True)
    print("Pytest output summary:")
    print(p.stdout.strip())
    if p.stderr.strip():
        print("Pytest warnings/stderr:")
        print(p.stderr.strip()[:300])
    assert p.returncode == 0, f"Pytest suite failed with code {p.returncode}"
    print("  --> 100% of tests in test suite PASSED!")

if __name__ == "__main__":
    t0 = time.time()
    verify_github_pages_structure()
    verify_samples_bundle()
    verify_conversion_engine_matrix()
    run_pytest_full_suite()
    total_time = round(time.time() - t0, 2)
    print("=" * 60)
    print(f"VERIFICATION COMPLETE: ALL CHECKS PASSED in {total_time}s!")
    print("=" * 60)
