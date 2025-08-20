import os

import modal
from dotenv import load_dotenv
from fastapi.responses import HTMLResponse
from pydantic import BaseModel, Field
from rdkit import Chem

# ==================== Modal Setup ====================

load_dotenv()
image = (
    modal.Image.debian_slim(python_version="3.10")
    .apt_install("wget", "openbabel", "libboost-all-dev")
    .run_commands(
        "wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64 -O /usr/local/bin/vina",
        "chmod +x /usr/local/bin/vina",
    )
    .pip_install("rdkit", "dockstring", "fastapi", "pydantic", "dotenv")
    .add_local_file("utils.py", "/root/utils.py")
)

app = modal.App("awesome-docking", image=image)
volume = modal.Volume.from_name("awesome-docking-volume", create_if_missing=True)

s3_bucket_name = os.getenv("BUCKET_NAME")
s3_access_credentials = modal.Secret.from_name("aws-s3-credentials")

with image.imports():
    from dockstring import load_target
    from utils import (
        create_complex_pdb,
        generate_visualizations,
        save_docking_results,
        setup_vina,
    )


# ==================== Data Models ====================


class DockingRequest(BaseModel):
    smiles: str = Field(..., description="SMILES string of the molecule to dock")
    target_name: str = Field("DRD2", description="Name of the protein target")


class DockingResult(BaseModel):
    smiles: str
    target_name: str
    score: float = Field(..., description="AutoDock Vina docking score (kcal/mol)")
    visualization_urls: list[str] = Field(..., description="URLs for 3D visualizations")


# ==================== Core Docking Function ====================


@app.function(
    timeout=300,
    volumes={
        "/data": modal.CloudBucketMount(s3_bucket_name, secret=s3_access_credentials)
    },
    min_containers=1,
    cpu=4,
)
@modal.fastapi_endpoint(method="POST")
def dock_molecule(request: DockingRequest) -> DockingResult:
    """
    Main docking endpoint that performs molecular docking and generates visualizations.

    Args:
        request: DockingRequest containing SMILES and target name

    Returns:
        DockingResult with docking score and visualization URLs
    """

    setup_vina()

    target = load_target(
        request.target_name,
        working_dir="/data",
        targets_dir="/data/targets",
    )

    docking_score, auxiliary_data = target.dock(request.smiles, num_cpus=4)

    if docking_score is None:
        raise ValueError(f"Docking failed for {request.smiles}")

    inchi_key = Chem.MolToInchiKey(auxiliary_data["ligand"])

    result_paths = save_docking_results(
        target_name=request.target_name,
        inchi_key=inchi_key,
    )

    combined_pdb_path = create_complex_pdb(
        target_name=request.target_name,
        inchi_key=inchi_key,
        ligand_pdb=result_paths["ligand"],
    )

    visualizations = generate_visualizations(
        target_name=request.target_name,
        inchi_key=inchi_key,
        ligand_pdb=result_paths["ligand"],
        complex_pdb=combined_pdb_path,
    )

    return DockingResult(
        smiles=request.smiles,
        target_name=request.target_name,
        score=float(docking_score),
        visualization_urls=visualizations,
    )


# ==================== Visualization Endpoint ====================


@app.function(
    timeout=60,
    volumes={
        "/data": modal.CloudBucketMount(s3_bucket_name, secret=s3_access_credentials)
    },
)
@modal.fastapi_endpoint(method="GET")
def view_structure(structure_type: str, target: str, molecule_id: str):
    """
    Serve interactive 3D visualizations of docked structures.

    Args:
        structure_type: "ligand" or "complex"
        target: Protein target name
        molecule_id: InChI key of the docked molecule
    """

    if structure_type not in ("ligand", "complex"):
        return HTMLResponse("Invalid structure_type", status_code=400)

    pdb_filename = (
        f"{molecule_id}_complex.pdb"
        if structure_type == "complex"
        else f"{molecule_id}.pdb"
    )
    pdb_path = f"/data/docking_results/{target}/{pdb_filename}"

    # pdb_text = read_file_with_retry(pdb_path, max_attempts=50)
    with open(pdb_path, "r", encoding="utf-8") as f:
        pdb_text = f.read()

    if not pdb_text:
        return HTMLResponse(
            content=f"<h2>Visualization not found</h2><p>File: {pdb_filename}</p>",
            status_code=404,
        )

    # Self-contained 3Dmol.js viewer embedding the PDB text directly
    style_script = (
        "viewer.setStyle({chain:'L'},{stick:{radius:0.2}});\n"
        "viewer.setStyle({not:{chain:'L'}},{cartoon:{color:'spectrum'}});\n"
        if structure_type == "complex"
        else "viewer.setStyle({}, {stick:{}});\n"
    )

    html = f"""
<!DOCTYPE html>
<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n<title>Docking Viewer</title>\n<script src=\"https://unpkg.com/3dmol@2.0.4/build/3Dmol-min.js\"></script>\n<style>html,body{{margin:0;padding:0;height:100%;}}#viewer{{width:100vw;height:100vh;}}</style>\n</head>\n<body>\n<div id=\"viewer\"></div>\n<script id=\"pdb\" type=\"text/plain\">{pdb_text}</script>\n<script>\n(function(){{\n  const pdb = document.getElementById('pdb').textContent;\n  const container = document.getElementById('viewer');\n  if (!window.$3Dmol || typeof $3Dmol.createViewer !== 'function') {{\n    container.innerHTML = '<div style=\\'padding:16px;font-family:sans-serif;\\'>Failed to load 3Dmol.js</div>';\n    return;\n  }}\n  const viewer = $3Dmol.createViewer(container, {{ backgroundColor: 'white' }});\n  viewer.addModel(pdb, 'pdb');\n  {style_script}viewer.zoomTo();\n  viewer.render();\n  window.addEventListener('resize', () => viewer.resize());\n}})();\n</script>\n</body>\n</html>\n"""

    return HTMLResponse(
        content=html,
        headers={
            "Content-Type": "text/html",
            "Cache-Control": "public, max-age=3600",
        },
    )
