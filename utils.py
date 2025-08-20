import os
import shutil
import time
from pathlib import Path

from dockstring.utils import convert_pdbqt_to_pdb, get_bin_dir, get_vina_filename


VINA_WRAPPER_SCRIPT = """#!/usr/bin/env python3
import sys
import subprocess

def main():
    args = sys.argv[1:]
    log_path = None
    filtered_args = []
    
    # Parse arguments to extract --log flag
    i = 0
    while i < len(args):
        if args[i] == '--log' and i + 1 < len(args):
            log_path = args[i + 1]
            i += 2
        else:
            filtered_args.append(args[i])
            i += 1
    
    # Run vina with filtered arguments
    proc = subprocess.Popen(
        ['/usr/local/bin/vina'] + filtered_args, 
        stdout=subprocess.PIPE, 
        stderr=subprocess.STDOUT
    )
    output, _ = proc.communicate()
    
    # Write log file if requested
    if log_path:
        try:
            with open(log_path, 'wb') as f:
                f.write(output)
        except Exception:
            pass  # Silently fail on log write errors
    
    # Output to stdout
    sys.stdout.buffer.write(output)
    sys.exit(proc.returncode)

if __name__ == '__main__':
    main()
"""


def setup_vina():
    pkg_vina = get_bin_dir() / get_vina_filename()

    if not Path("/usr/local/bin/vina").exists():
        raise RuntimeError("System vina not found at /usr/local/bin/vina")

    # Create wrapper script
    pkg_vina.parent.mkdir(parents=True, exist_ok=True)
    pkg_vina.write_text(VINA_WRAPPER_SCRIPT)
    os.chmod(pkg_vina, 0o755)


# Viewer HTML is generated on the fly in the endpoint using 3Dmol.js


# --- HELPERS --- #
def save_docking_results(target_name: str, inchi_key: str) -> dict:
    """Save the docked ligand PDB file to persistent storage"""

    # Ensure output directory exists
    output_dir = f"/data/docking_results/{target_name}"
    os.makedirs(output_dir, exist_ok=True)

    # Wait for docking output file
    source_file = "/data/docked_ligand.pdb"
    for _ in range(50):
        if os.path.exists(source_file):
            break
        time.sleep(0.1)
    else:
        raise FileNotFoundError("Docking output not found")

    # Move to permanent location (use copy+remove; rename not supported on this FS)
    destination = f"{output_dir}/{inchi_key}.pdb"
    shutil.copyfile(source_file, destination)
    try:
        os.remove(source_file)
    except Exception:
        pass

    return {"ligand": destination}


def create_complex_pdb(target_name: str, inchi_key: str, ligand_pdb: str) -> str:
    """Combine protein and ligand into a single PDB file"""

    # Convert protein from PDBQT to PDB format
    protein_pdbqt = f"/data/targets/{target_name}_target.pdbqt"
    temp_protein_pdb = f"/tmp/{target_name}_target.pdb"
    convert_pdbqt_to_pdb(protein_pdbqt, temp_protein_pdb)

    # Read protein and ligand structures
    with open(temp_protein_pdb, "r") as f:
        protein_lines = f.readlines()
    with open(ligand_pdb, "r") as f:
        ligand_lines = f.readlines()

    # Extract best-scoring ligand pose (first MODEL if multiple)
    ligand_atoms = _extract_best_pose(ligand_lines)

    # Combine structures with proper atom numbering
    output_path = f"/data/docking_results/{target_name}/{inchi_key}_complex.pdb"
    _write_complex_pdb(protein_lines, ligand_atoms, output_path)

    return output_path


def _extract_best_pose(ligand_lines: list[str]) -> list[str]:
    """Extract atoms from the first MODEL block (best Vina pose)."""

    atoms: list[str] = []
    in_model = False

    for line in ligand_lines:
        if line.startswith("MODEL") and not in_model:
            in_model = True
            continue
        if in_model:
            if line.startswith("ENDMDL"):
                break
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atoms.append(line)

    if atoms:
        return atoms

    # Fallback: files without MODEL blocks – collect all atoms
    return [
        line
        for line in ligand_lines
        if line.startswith("ATOM") or line.startswith("HETATM")
    ]


def _write_complex_pdb(
    protein_lines: list[str], ligand_atoms: list[str], output_path: str
):
    """Write combined protein-ligand complex with proper formatting"""

    with open(output_path, "w") as out:
        # Write protein atoms
        atom_counter = 0
        for line in protein_lines:
            if line.startswith(("ATOM", "HETATM")):
                atom_counter += 1
                # Update atom serial number
                line = line[:6] + f"{atom_counter:5d}" + line[11:]
            if not line.startswith("END"):
                out.write(line)

        # Write ligand atoms with chain ID "L"
        for line in ligand_atoms:
            atom_counter += 1
            # Update atom serial and chain ID
            line = line[:6] + f"{atom_counter:5d}" + line[11:21] + "L" + line[22:]
            out.write(line)

        out.write("TER\nEND\n")


def generate_visualizations(
    target_name: str, inchi_key: str, ligand_pdb: str, complex_pdb: str
) -> list[str]:
    """Return viewer URLs for ligand and complex without pre-rendering HTML files"""

    base_url = "https://onepot-ai--awesome-docking-view-structure.modal.run"

    return [
        f"{base_url}?structure_type=ligand&target={target_name}&molecule_id={inchi_key}",
        f"{base_url}?structure_type=complex&target={target_name}&molecule_id={inchi_key}",
    ]
