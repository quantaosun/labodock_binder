#!/usr/bin/env python3
"""
Covalent Docking Workflow using GNINA

Converted from covalent_dock_gnina_colab.ipynb to a local Python script.

Primary input is a **protein-ligand complex PDB file**. The script auto-extracts
the receptor and native ligand, then uses the ligand's position to set the grid
box center and size — no manual coordinate entry needed.

Workflow:
  1. Parse complex PDB  — extract receptor + native ligand, compute grid box
  2. Environment         — verify RDKit and GNINA are available
  3. Ligand prep         — SMILES → 3D conformer, identify warhead via SMARTS
  4. Covalent docking    — run GNINA with covalent constraints
  5. Analysis            — parse SDF output, show docking affinities

Usage:
    conda activate rdkit
    python covalent_dock_gnina.py \
        --complex protein_ligand_complex.pdb \
        --native_ligand XEY \
        --reactive_res "A:145:SG" \
        --smiles "C1(=O)[C@H]...CCN1"

GNINA must be installed and on PATH, or use --gnina_bin /path/to/gnina.
RDKit must be installed in the active Python environment.
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

# ── Optional: silence RDKit warnings ──────────────────────────────────────
try:
    from rdkit import RDLogger
    RDLogger.logger().setLevel(RDLogger.ERROR)
except Exception:
    pass


# ── PDB parsing constants ─────────────────────────────────────────────────

# Standard 20 amino acid 3-letter codes (plus common protonation variants)
STANDARD_RESIDUES: Set[str] = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "CYX", "CYM",
    "GLN", "GLU", "GLY", "HIS", "HID", "HIE", "HIP",
    "ILE", "LEU", "LYS", "LYN", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "TYM", "VAL",
}

# Common non-ligand HETATM entries to ignore
NON_LIGAND_RESIDUES: Set[str] = {
    "HOH", "WAT",        # water
    "NA",  "K",   "CL",  "MG",  "CA",  "ZN",  "FE",  "MN",
    "CU",  "CO",  "NI",  "CD",  "HG",  "PT",  "AU",  "AG",
    "BR",  "I",   "F",            # ions
    "SO4", "PO4", "GOL", "EDO", "ACT", "DMS", "MPD", "PEG",
    "BME", "TRS", "CIT", "EPE", "FMT", "IPA", "PGE", "TAM",
    "NO3", "NH4", "IMD", "ACM", "BCT", "BMA", "BOG", "BNG",
    "C8E", "DDM", "DMF", "DTT", "EGL", "LDA", "LMT", "MES",
    "MRD", "MYR", "N8E", "OCT", "OGL", "OLC", "PLC", "SDS",
    "SPM", "TMA", "TRT", "URE",
}


# ============================================================================
# Step 1 — Parse complex PDB
# ============================================================================

def _calc_center_and_size(coords: List[Tuple[float, float, float]],
                           padding: float = 8.0
                           ) -> Tuple[Tuple[float, float, float],
                                      Tuple[float, float, float]]:
    """
    Compute grid-box center (centroid) and size from a list of (x, y, z)
    coordinates.  `padding` is added to each side (so total extra = 2×padding).
    """
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]

    center = (
        (min(xs) + max(xs)) / 2.0,
        (min(ys) + max(ys)) / 2.0,
        (min(zs) + max(zs)) / 2.0,
    )
    size = (
        max(1.0, max(xs) - min(xs) + 2 * padding),
        max(1.0, max(ys) - min(ys) + 2 * padding),
        max(1.0, max(zs) - min(zs) + 2 * padding),
    )
    return center, size


def parse_complex_pdb(
    complex_pdb: str,
    output_dir: str = ".",
    box_padding: float = 8.0,
    protein_chain: Optional[str] = None,
    native_ligand: Optional[str] = None,
) -> Tuple[str, Tuple[float, float, float], Tuple[float, float, float]]:
    """
    Parse a protein-ligand complex PDB file.

    Separates protein atoms from non-standard (ligand) atoms, writes
    ``receptor.pdb`` and ``native_ligand.pdb`` into *output_dir*.

    Parameters
    ----------
    native_ligand : str or None
        3-letter residue name of the native ligand (e.g. ``"XEY"``, ``"GTP"``).
        When provided, only atoms with this residue name are used for grid-box
        computation.  When ``None``, all non-standard residues are used.

    Returns
    -------
    receptor_path : str
        Path to the written receptor PDB file.
    grid_center : (float, float, float)
        Centroid of the detected ligand atoms.
    grid_size : (float, float, float)
        Bounding-box extent of the ligand plus *box_padding* on each side.
    """
    from rdkit import Chem

    # ── Read raw lines ─────────────────────────────────────────────────
    with open(complex_pdb, "r") as fh:
        lines = fh.readlines()

    rec_atoms: List[str] = []
    all_lig_atoms: List[str] = []        # all candidate ligand atoms
    all_lig_coords: List[Tuple[float, float, float]] = []

    het_resnames_seen: Set[str] = set()

    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            if line.startswith(("TER", "END", "MODEL", "ENDMDL")):
                rec_atoms.append(line)
            continue

        resname = line[17:20].strip()
        chain = line[21:22].strip()
        atom_type = line[:6].strip()

        # If a specific protein chain is requested, route other chains
        if protein_chain and chain != protein_chain:
            if atom_type == "HETATM" and resname not in NON_LIGAND_RESIDUES:
                all_lig_atoms.append(line)
                try:
                    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    all_lig_coords.append((x, y, z))
                except ValueError:
                    pass
            elif atom_type == "ATOM" and resname not in STANDARD_RESIDUES:
                all_lig_atoms.append(line)
                try:
                    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    all_lig_coords.append((x, y, z))
                except ValueError:
                    pass
            continue

        if atom_type == "HETATM":
            if resname in NON_LIGAND_RESIDUES:
                continue
            het_resnames_seen.add(resname)
            all_lig_atoms.append(line)
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                all_lig_coords.append((x, y, z))
            except ValueError:
                pass
        elif atom_type == "ATOM":
            if resname in STANDARD_RESIDUES:
                rec_atoms.append(line)
            else:
                all_lig_atoms.append(line)
                try:
                    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    all_lig_coords.append((x, y, z))
                except ValueError:
                    pass

    # ── Filter to requested native ligand (if specified) ───────────────
    if native_ligand:
        native_upper = native_ligand.upper()
        lig_atoms = [ln for ln in all_lig_atoms
                     if ln[17:20].strip().upper() == native_upper]
        lig_coords_raw = [
            (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
            for ln in lig_atoms
        ]
        lig_coords = []
        for c in lig_coords_raw:
            try:
                lig_coords.append(c)
            except ValueError:
                pass
        if not lig_atoms:
            all_seen = sorted({ln[17:20].strip() for ln in all_lig_atoms})
            print(
                f"⚠️  Native ligand '{native_ligand}' not found in PDB.\n"
                f"   HETATM residues present: {all_seen or 'none'}\n"
                f"   Falling back to all non-standard residues for grid box."
            )
            lig_atoms = all_lig_atoms
            lig_coords = all_lig_coords
        else:
            print(f"🔬 Filtered to native ligand '{native_ligand}' "
                  f"({len(lig_atoms)} atoms)")
    else:
        lig_atoms = all_lig_atoms
        lig_coords = all_lig_coords
        if len(het_resnames_seen) > 1:
            print(
                f"💡 Multiple non-standard residues found: {sorted(het_resnames_seen)}.\n"
                f"   Using ALL of them for grid box. To target a specific ligand,"
                f" pass --native_ligand <3-letter-code>."
            )

    # ── Sanity checks ──────────────────────────────────────────────────
    if not rec_atoms:
        sys.exit(
            f"❌ No protein atoms found in '{complex_pdb}'.\n"
            f"   Standard residues checked: {sorted(STANDARD_RESIDUES)}\n"
            f"   If your protein uses unusual residue names, prepare the receptor"
            f"   separately and pass it via --receptor instead of --complex."
        )

    if not lig_coords:
        print(
            f"⚠️  No ligand coordinates could be extracted from '{complex_pdb}'.\n"
            f"   HETATM residues seen: {het_resnames_seen or 'none'}\n"
            f"   Non-ligand filter: {sorted(NON_LIGAND_RESIDUES)}\n"
            f"   Use --center and --size manually for the grid box."
        )

    # ── Write receptor PDB ─────────────────────────────────────────────
    os.makedirs(output_dir, exist_ok=True)
    rec_path = os.path.join(output_dir, "receptor.pdb")
    with open(rec_path, "w") as fh:
        fh.writelines(rec_atoms)
    print(f"✅ Receptor ({len(rec_atoms)} atoms) → {rec_path}")

    # ── Write native ligand PDB (for reference / inspection) ───────────
    lig_path = os.path.join(output_dir, "native_ligand.pdb")
    with open(lig_path, "w") as fh:
        fh.writelines(lig_atoms)
    displayed = native_ligand.upper() if native_ligand else \
                ", ".join(sorted(het_resnames_seen)) if het_resnames_seen else "none"
    print(f"✅ Native ligand ({len(lig_atoms)} atoms, residue(s): {displayed})"
          f" → {lig_path}")

    # ── Compute grid box ───────────────────────────────────────────────
    if lig_coords:
        center, size = _calc_center_and_size(lig_coords, padding=box_padding)
    else:
        center, size = (0.0, 0.0, 0.0), (20.0, 20.0, 20.0)

    print(f"📐 Grid center: ({center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f})")
    print(f"📐 Grid size:   ({size[0]:.1f}, {size[1]:.1f}, {size[2]:.1f}) "
          f"(ligand extent + {box_padding} Å padding)")

    return rec_path, center, size


# ============================================================================
# Step 2 — Environment Setup
# ============================================================================

def check_rdkit() -> str:
    """Verify RDKit is importable and return its version."""
    try:
        from rdkit import rdBase
        ver = rdBase.rdkitVersion
        print(f"✅ RDKit {ver} is available")
        return ver
    except ImportError:
        sys.exit(
            "❌ RDKit is not installed. Install it with:\n"
            "   conda activate rdkit && conda install -c conda-forge rdkit\n"
            "   or  conda activate rdkit && pip install rdkit-pypi"
        )


GNINA_URL = (
    "https://github.com/gnina/gnina/releases/download/v1.1/gnina"
)

RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"


def download_file(url: str, dest: str, desc: str = "") -> None:
    """Download a file with a progress bar via urllib."""
    import urllib.request

    dest_path = Path(dest)
    if dest_path.exists():
        print(f"📦 {desc or dest} already exists, skipping download.")
        return

    print(f"⬇️  Downloading {desc or url} → {dest} ...")
    try:
        from tqdm import tqdm as _tqdm  # optional progress bar
        import urllib.request

        class _TqdmUpTo(_tqdm):
            def update_to(self, b=1, bsize=1, tsize=None):
                if tsize is not None:
                    self.total = tsize
                self.update(b * bsize - self.n)

        with _TqdmUpTo(unit="B", unit_scale=True, unit_divisor=1024,
                        miniters=1, desc=dest_path.name) as t:
            urllib.request.urlretrieve(url, dest, reporthook=t.update_to)
    except ImportError:
        # tqdm not available — simple download
        urllib.request.urlretrieve(url, dest)
        print(f"✅ Downloaded to {dest}")

    if not dest_path.exists() or dest_path.stat().st_size == 0:
        sys.exit(f"❌ Download failed: {url}")


def download_pdb(pdb_id: str, dest_dir: str = ".") -> str:
    """Download a PDB file from RCSB by its 4-character ID."""
    pdb_id = pdb_id.upper().strip()
    dest = os.path.join(dest_dir, f"{pdb_id}.pdb")
    url = RCSB_URL.format(pdb_id=pdb_id)
    download_file(url, dest, desc=f"PDB {pdb_id}")
    return dest


def check_gnina(gnina_bin: str, download: bool = False) -> str:
    """Verify the GNINA binary exists; optionally download it."""
    path = Path(gnina_bin)
    if not path.exists():
        if download:
            download_file(GNINA_URL, gnina_bin, desc="GNINA binary")
        else:
            sys.exit(
                f"❌ GNINA not found at '{gnina_bin}'.\n"
                f"   Run with --download_gnina to fetch it automatically, or\n"
                f"   download from {GNINA_URL}\n"
                f"   or pass --gnina_bin /path/to/gnina"
            )
    if not os.access(path, os.X_OK):
        path.chmod(0o755)
        print(f"🔧 Made '{gnina_bin}' executable")
    print(f"✅ GNINA binary found at '{gnina_bin}'")
    return str(path.resolve())


# ============================================================================
# Step 3 — Ligand Preparation
# ============================================================================

def prepare_ligand(
    smiles: str,
    output_sdf: str,
    warhead_smarts: str,
) -> Tuple[Optional["Chem.Mol"], List[int]]:
    """
    Generate a 3D conformer from SMILES, write to SDF, and identify the
    warhead atom indices via SMARTS matching.

    Returns (mol, warhead_atom_indices).  ``mol`` is ``None`` on failure.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = None
    warhead_indices: List[int] = []

    # ── Generate 3D conformer from SMILES ──────────────────────────────
    if smiles:
        print(f"⚙️  Generating 3D conformer from SMILES: {smiles}")
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print("❌ Invalid SMILES string. Please check your input.")
            return None, []

        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)

        writer = Chem.SDWriter(output_sdf)
        writer.write(mol)
        writer.close()
        print(f"✅ Saved 3D ligand to {output_sdf}")
    elif os.path.exists(output_sdf):
        print(f"📥 Loading existing ligand from {output_sdf}")
        mol = Chem.MolFromMolFile(output_sdf, removeHs=False)
        if mol is None:
            print(f"❌ Could not parse {output_sdf}")
            return None, []
    else:
        print("❌ Please provide a SMILES string (--smiles) or an existing SDF file.")
        return None, []

    # ── Identify warhead atom(s) ───────────────────────────────────────
    if mol is not None:
        pattern = Chem.MolFromSmarts(warhead_smarts)
        if pattern is None:
            print(f"❌ Invalid warhead SMARTS pattern: {warhead_smarts}")
            return mol, []

        matches = mol.GetSubstructMatches(pattern)
        if matches:
            warhead_indices = list(matches[0])
            print(f"✅ Warhead found at atom indices: {warhead_indices}")

            try:
                from rdkit.Chem import Draw
                mol_2d = Chem.Mol(mol)
                AllChem.Compute2DCoords(mol_2d)
                img = Draw.MolToImage(mol_2d, highlightAtoms=warhead_indices)
                png_path = Path(output_sdf).with_suffix(".png")
                img.save(str(png_path))
                print(f"📸 2D structure with warhead highlight saved to {png_path}")
            except Exception as e:
                print(f"⚠️  Could not generate 2D image: {e}")
        else:
            print(
                f"⚠️  Warhead pattern '{warhead_smarts}' not found in the ligand.\n"
                f"   GNINA will not apply covalent constraints automatically."
            )

    return mol, warhead_indices


# ============================================================================
# Step 4 — Covalent Docking
# ============================================================================

def run_docking(
    gnina_bin: str,
    receptor_file: str,
    ligand_file: str,
    reactive_residue: str,
    warhead_smarts: str,
    center_x: float,
    center_y: float,
    center_z: float,
    size_x: float,
    size_y: float,
    size_z: float,
    output_file: str = "covalent_poses.sdf",
    exhaustiveness: int = 8,
    cnn_scoring: str = "none",
    autobox_ligand: Optional[str] = None,
    extra_flags: Optional[List[str]] = None,
) -> int:
    """Run GNINA covalent docking.  Returns the exit code."""
    cmd = [
        gnina_bin,
        "--receptor", receptor_file,
        "--ligand", ligand_file,
        "--out", output_file,
        "--center_x", str(center_x),
        "--center_y", str(center_y),
        "--center_z", str(center_z),
        "--size_x", str(size_x),
        "--size_y", str(size_y),
        "--size_z", str(size_z),
        "--covalent_rec_atom", reactive_residue,
        "--covalent_lig_atom_pattern", warhead_smarts,
        "--exhaustiveness", str(exhaustiveness),
        "--cnn_scoring", cnn_scoring,
    ]

    if autobox_ligand:
        cmd += ["--autobox_ligand", autobox_ligand]

    if extra_flags:
        cmd += extra_flags

    display_cmd = " ".join(
        f"'{a}'" if (" " in a or "(" in a or "=" in a) else a
        for a in cmd
    )
    print(f"🚀 Running: {display_cmd}")
    print("-" * 72)

    result = subprocess.run(cmd, capture_output=False)
    return result.returncode


# ============================================================================
# Step 5 — Analysis
# ============================================================================

def analyze_results(sdf_file: str) -> "pd.DataFrame":
    """Parse GNINA output SDF and return a DataFrame of minimized affinities."""
    import pandas as pd

    affinities: List[float] = []
    with open(sdf_file, "r") as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if "<minimizedAffinity>" in line:
                affinities.append(float(lines[i + 1].strip()))

    if not affinities:
        print(f"⚠️  No '<minimizedAffinity>' tags found in {sdf_file}")
        return pd.DataFrame()

    df = pd.DataFrame(affinities, columns=["Minimized Affinity (kcal/mol)"])
    df.index.name = "Pose"
    df.index = df.index + 1
    print(df.to_string())

    best = df["Minimized Affinity (kcal/mol)"].min()
    print(f"\n🏆 Best affinity: {best:.2f} kcal/mol "
          f"(pose {df['Minimized Affinity (kcal/mol)'].idxmin()})")

    # Save bar chart
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 4))
        df.plot(kind="bar", ax=ax, legend=False, color="steelblue")
        ax.set_xlabel("Pose")
        ax.set_ylabel("Minimized Affinity (kcal/mol)")
        ax.set_title("GNINA Covalent Docking — Affinities")
        fig.tight_layout()
        chart_path = Path(sdf_file).with_suffix(".png")
        fig.savefig(str(chart_path), dpi=150)
        print(f"📊 Affinity chart saved to {chart_path}")
    except Exception as e:
        print(f"⚠️  Could not save chart: {e}")

    return df


# ============================================================================
# CLI
# ============================================================================

def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Covalent Docking Workflow using GNINA",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Auto-extract receptor + grid box from a complex PDB:
  python covalent_dock_gnina.py \\
      --complex 6lze_complex.pdb \\
      --native_ligand XEY \\
      --reactive_res "A:145:SG" \\
      --smiles "C1(=O)[C@H](...)CCN1"

  # Download PDB directly by ID (Colab-friendly, no upload needed):
  python covalent_dock_gnina.py \\
      --pdb_id 5XV7 \\
      --native_ligand EMH \\
      --reactive_res "A:227:OH" \\
      --smiles "O=S(...)(F)=O" \\
      --warhead_smarts "S(=O)(=O)F" \\
      --download_gnina

  # Manual receptor + grid box (skip complex parsing):
  python covalent_dock_gnina.py \\
      --receptor rec.pdb \\
      --reactive_res "A:145:SG" \\
      --smiles "C1(=O)[C@H](...)CCN1" \\
      --center 30.9 -12.0 7.0 --size 20 20 20
        """,
    )

    # ── Input source ───────────────────────────────────────────────────
    p.add_argument(
        "--complex", default=None,
        help="Path to protein-ligand complex PDB file. "
             "The script will extract the receptor and native ligand, "
             "and auto-compute the grid box from the ligand position.",
    )
    p.add_argument(
        "--pdb_id", default=None,
        help="RCSB PDB ID (e.g., '5XV7'). The script will download "
             "{ID}.pdb from RCSB and treat it as --complex. "
             "Convenient for Colab — no manual upload needed.",
    )
    p.add_argument(
        "-r", "--receptor", default=None,
        help="Path to receptor PDB file. Use this if you already have a "
             "prepared receptor (skips --complex parsing).",
    )

    # ── Required for covalent docking ──────────────────────────────────
    p.add_argument(
        "-c", "--reactive_res", required=True,
        help="Reactive residue in GNINA format: Chain:ResNum:AtomName "
             "(e.g., A:145:SG for Cys145-SG).",
    )
    p.add_argument(
        "-s", "--smiles", default=None,
        help="SMILES string of the covalent inhibitor to dock. "
             "If omitted, --ligand_sdf must point to an existing SDF file.",
    )

    # ── Ligand options ─────────────────────────────────────────────────
    p.add_argument(
        "--ligand_sdf", default="ligand_prepared.sdf",
        help="Path to write (or read) the prepared ligand SDF.  "
             "Default: ligand_prepared.sdf",
    )
    p.add_argument(
        "-w", "--warhead_smarts", default="C(=O)N",
        help="SMARTS pattern to identify the covalent warhead atom(s) on "
             "the ligand.  Default: C(=O)N",
    )

    # ── Output ─────────────────────────────────────────────────────────
    p.add_argument(
        "-o", "--output", default="covalent_poses.sdf",
        help="Output SDF for docking poses.  Default: covalent_poses.sdf",
    )
    p.add_argument(
        "--out_dir", default=".",
        help="Directory for intermediate files (receptor.pdb, native_ligand.pdb).  "
             "Default: current directory",
    )

    # ── Grid box (optional overrides; auto-computed from --complex) ────
    p.add_argument(
        "--center", nargs=3, type=float, default=None,
        metavar=("X", "Y", "Z"),
        help="Grid box center (x y z).  Overrides auto-computation from --complex.",
    )
    p.add_argument(
        "--size", nargs=3, type=float, default=None,
        metavar=("X", "Y", "Z"),
        help="Grid box size (x y z).  Overrides auto-computation from --complex.",
    )
    p.add_argument(
        "-n", "--native_ligand", default=None,
        help="3-letter residue name of the native ligand in the complex PDB "
             "(e.g., 'XEY', 'GTP'). When provided, only this residue is used "
             "for grid-box computation. When omitted, all non-standard residues "
             "are used.",
    )
    p.add_argument(
        "--box_padding", type=float, default=8.0,
        help="Padding added to each side of the ligand bounding box for the "
             "grid size (Å).  Only used with --complex.  Default: 8.0",
    )
    p.add_argument(
        "--protein_chain", default=None,
        help="Limit receptor extraction to a specific chain (e.g., 'A').  "
             "Useful for multi-chain complexes.",
    )

    # ── Docking parameters ─────────────────────────────────────────────
    p.add_argument(
        "--exhaustiveness", type=int, default=8,
        help="Docking exhaustiveness.  Default: 8",
    )
    p.add_argument(
        "--cnn_scoring", default="none",
        choices=["none", "rescore", "refinement", "all"],
        help="CNN scoring mode.  Default: none",
    )
    p.add_argument(
        "--autobox_ligand", default=None,
        help="If set, use GNINA's --autobox_ligand to center the box on "
             "this ligand file (instead of explicit coordinates).",
    )

    # ── GNINA binary ───────────────────────────────────────────────────
    p.add_argument(
        "--gnina_bin", default="./gnina",
        help="Path to the GNINA executable.  Default: ./gnina",
    )
    p.add_argument(
        "--download_gnina", action="store_true",
        help="Auto-download the GNINA binary if not found. "
             "Useful on Colab where nothing is pre-installed.",
    )
    p.add_argument(
        "--gpu", action="store_true",
        help="Pass --gpu flag to GNINA (use GPU acceleration if available).",
    )
    p.add_argument(
        "--extra", nargs="*", default=[],
        help="Extra flags to pass directly to GNINA "
             "(e.g., --extra --scoring vinardo).",
    )

    return p.parse_args(argv)


# ============================================================================
# Main
# ============================================================================

def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)

    # ── Validate input source ──────────────────────────────────────────
    if args.pdb_id:
        args.complex = download_pdb(args.pdb_id, dest_dir=args.out_dir)

    if not args.complex and not args.receptor:
        sys.exit(
            "❌ Either --complex (protein-ligand complex PDB), --pdb_id, or "
            "--receptor (prepared receptor PDB) is required."
        )
    if args.complex and args.receptor:
        print("⚠️  Both --complex and --receptor provided; using --receptor directly.")

    # ── 1 — Environment ────────────────────────────────────────────────
    print("=" * 72)
    print("Step 1 — Environment Setup")
    print("=" * 72)
    check_rdkit()
    gnina = check_gnina(args.gnina_bin, download=args.download_gnina)
    print(f"🎯 Reactive residue: {args.reactive_res}")

    # ── 2 — Receptor & grid box ────────────────────────────────────────
    print()
    print("=" * 72)
    print("Step 2 — Receptor & Grid Box")
    print("=" * 72)

    center = args.center
    size = args.size

    if args.receptor:
        # Use pre-prepared receptor
        receptor_path = args.receptor
        if not Path(receptor_path).exists():
            sys.exit(f"❌ Receptor file not found: {receptor_path}")
        print(f"✅ Receptor: {receptor_path}")
        if not center or not size:
            sys.exit(
                "❌ --center and --size are required when using --receptor "
                "(they cannot be auto-computed without --complex)."
            )
    else:
        # Parse complex PDB → extract receptor + compute grid box
        receptor_path, auto_center, auto_size = parse_complex_pdb(
            complex_pdb=args.complex,
            output_dir=args.out_dir,
            box_padding=args.box_padding,
            protein_chain=args.protein_chain,
            native_ligand=args.native_ligand,
        )
        # Manual overrides take precedence
        if center is None:
            center = auto_center
        if size is None:
            size = auto_size

    print(f"📍 Final grid center: ({center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f})")
    print(f"📍 Final grid size:   ({size[0]:.1f}, {size[1]:.1f}, {size[2]:.1f})")

    # ── 3 — Ligand ─────────────────────────────────────────────────────
    print()
    print("=" * 72)
    print("Step 3 — Ligand Preparation")
    print("=" * 72)
    mol, warhead_idx = prepare_ligand(
        smiles=args.smiles,
        output_sdf=args.ligand_sdf,
        warhead_smarts=args.warhead_smarts,
    )
    if mol is None:
        return 1

    # ── 4 — Docking ────────────────────────────────────────────────────
    print()
    print("=" * 72)
    print("Step 4 — Covalent Docking")
    print("=" * 72)

    extra_flags = list(args.extra)
    if args.gpu:
        extra_flags.append("--gpu")

    rc = run_docking(
        gnina_bin=gnina,
        receptor_file=receptor_path,
        ligand_file=args.ligand_sdf,
        reactive_residue=args.reactive_res,
        warhead_smarts=args.warhead_smarts,
        center_x=center[0],
        center_y=center[1],
        center_z=center[2],
        size_x=size[0],
        size_y=size[1],
        size_z=size[2],
        output_file=args.output,
        exhaustiveness=args.exhaustiveness,
        cnn_scoring=args.cnn_scoring,
        autobox_ligand=args.autobox_ligand,
        extra_flags=extra_flags if extra_flags else None,
    )

    if rc != 0:
        print(f"❌ GNINA exited with non-zero status: {rc}")
        return rc

    # ── 5 — Analysis ───────────────────────────────────────────────────
    print()
    print("=" * 72)
    print("Step 5 — Analysis")
    print("=" * 72)
    if not os.path.exists(args.output):
        print(f"❌ Output file '{args.output}' was not created. Docking may have failed.")
        return 1

    analyze_results(args.output)
    print("\n✅ Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
