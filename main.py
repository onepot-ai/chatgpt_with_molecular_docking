import modal
from fastapi.responses import HTMLResponse
from pydantic import BaseModel, Field
from rdkit import Chem

# ==================== Modal Setup ====================


image = (
    modal.Image.debian_slim(python_version="3.10")
    .apt_install("wget", "openbabel", "libboost-all-dev")
    .run_commands(
        "wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64 -O /usr/local/bin/vina",
        "chmod +x /usr/local/bin/vina",
    )
    .pip_install("rdkit", "dockstring", "nglview", "fastapi", "pydantic")
    .add_local_file("utils.py", "/root/utils.py")
)

app = modal.App("molecular-docking-demo", image=image)
volume = modal.Volume.from_name("docking-results-volume", create_if_missing=True)

with image.imports():
    from dockstring import load_target
    from utils import (
        create_complex_pdb,
        generate_visualizations,
        read_file_with_retry,
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
    docking_score: float = Field(
        ..., description="AutoDock Vina docking score (kcal/mol)"
    )
    visualization_urls: dict = Field(..., description="URLs for 3D visualizations")


# ==================== Core Docking Function ====================


@app.function(
    timeout=300,
    volumes={"/data": volume},
    min_containers=1,
    cpu=4,
)
@modal.web_endpoint(method="POST")
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
        docking_score=float(docking_score),
        visualization_urls=visualizations,
    )


# ==================== Visualization Endpoint ====================


@app.function(
    timeout=60,
    volumes={"/data": volume},
    min_containers=2,
    allow_concurrent_inputs=100,
)
@modal.web_endpoint(method="GET")
def view_structure(structure_type: str, target: str, molecule_id: str):
    """
    Serve interactive 3D visualizations of docked structures.

    Args:
        structure_type: "ligand" or "complex"
        target: Protein target name
        molecule_id: InChI key of the docked molecule
    """

    filename = (
        f"{molecule_id}_complex.html"
        if structure_type == "complex"
        else f"{molecule_id}.html"
    )
    file_path = f"/data/docking_results/{target}/{filename}"

    html_content = read_file_with_retry(file_path, max_attempts=50)
    if html_content:
        return HTMLResponse(
            content=html_content,
            headers={
                "Content-Type": "text/html",
                "Cache-Control": "public, max-age=3600",
            },
        )
    else:
        return HTMLResponse(
            content=f"<h2>Visualization not found</h2><p>File: {filename}</p>",
            status_code=404,
        )
