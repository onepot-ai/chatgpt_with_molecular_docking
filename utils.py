import os
import time
from pathlib import Path

import nglview as nv

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


# --- PDB visualization (returns standalone HTML) ---
def plot_pdb_to_html(
    pdb_path: str,
    background: str = "#FDFDFD",
    out_path: str = "view.html",
    width: str = "100%",
    height: str = "900px",
) -> str:
    v = nv.show_structure_file(pdb_path, ext="pdb")
    v.background = background
    v.clear_representations()
    v.add_representation("cartoon", selection="polymer", color="sstruc")
    v.add_representation("licorice", selection="not polymer")
    v.center()
    try:
        v._set_size(width, height)
    except Exception:
        try:
            v.layout.width = width
            v.layout.height = height
        except Exception:
            pass
    try:
        nv.write_html(out_path, [v])
    except Exception:
        from nglview import write_html

        write_html(out_path, [v])


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

    # Move to permanent location
    destination = f"{output_dir}/{inchi_key}.pdb"
    os.rename(source_file, destination)

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
    """Extract the best-scoring pose from multi-model PDB"""

    atom_lines = []
    in_first_model = True

    for line in ligand_lines:
        if line.startswith("MODEL") and line.strip() != "MODEL        1":
            in_first_model = False
        if in_first_model and (line.startswith("ATOM") or line.startswith("HETATM")):
            atom_lines.append(line)
        if line.startswith("ENDMDL"):
            break

    return (
        atom_lines
        if atom_lines
        else [l for l in ligand_lines if l.startswith("ATOM") or l.startswith("HETATM")]
    )


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
) -> dict:
    """Generate HTML visualizations using NGL viewer"""

    base_url = "https://onepot-ai--molecular-docking-demo-view-structure.modal.run"

    # Generate visualization files
    plot_pdb_to_html(
        ligand_pdb, output_path=f"/data/docking_results/{target_name}/{inchi_key}.html"
    )

    plot_pdb_to_html(
        complex_pdb,
        output_path=f"/data/docking_results/{target_name}/{inchi_key}_complex.html",
    )

    # Force filesystem sync
    os.sync()

    return {
        "ligand": f"{base_url}?structure_type=ligand&target={target_name}&molecule_id={inchi_key}",
        "complex": f"{base_url}?structure_type=complex&target={target_name}&molecule_id={inchi_key}",
    }


def read_file_with_retry(file_path: str, max_attempts: int = 50) -> Optional[str]:
    """Read file with retries to handle volume sync delays"""

    for attempt in range(max_attempts):
        if os.path.exists(file_path):
            try:
                # Force metadata sync
                os.stat(file_path)
                time.sleep(0.05)

                with open(file_path, "r", encoding="utf-8") as f:
                    content = f.read()

                # Verify content is valid
                if content and len(content) > 100:
                    return content
            except Exception:
                pass

        # Exponential backoff
        delay = min(0.1 * (1.5 ** (attempt // 10)), 2.0)
        time.sleep(delay)

    return None
