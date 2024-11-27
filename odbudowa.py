import argparse
from Bio import PDB
import os
import numpy as np

def load_template(template_dir, residue_name):
    pdb_parser = PDB.PDBParser(QUIET=True)
    template_path = os.path.join(template_dir, f"{residue_name}.pdb")
    if not os.path.exists(template_path):
        raise FileNotFoundError(f"Template file for {residue_name} not found: {template_path}")
    structure = pdb_parser.get_structure(residue_name, template_path)
    return structure

def iteratively_superimpose(fixed_atoms, template_atoms):
    max_iterations = 5
    rmsd_threshold = 0.5

    for iteration in range(max_iterations):
        if len(fixed_atoms) != len(template_atoms):
            raise ValueError("Mismatched atom counts between fixed and template atoms during iteration.")

        fixed_coords = np.array([atom.get_coord() for atom in fixed_atoms])
        template_coords = np.array([atom.get_coord() for atom in template_atoms])

        sup = PDB.Superimposer()
        sup.set_atoms(fixed_atoms, template_atoms)
        sup.apply(template_atoms)

        distances = np.linalg.norm(fixed_coords - np.array([atom.get_coord() for atom in template_atoms]), axis=1)
        rmsd = np.sqrt(np.mean(distances**2))
        print(f"Iteration {iteration + 1}: RMSD = {rmsd:.4f}")

        if rmsd <= rmsd_threshold or len(fixed_atoms) <= 3:
            break

    return sup.rotran, rmsd

def rebuild_structure(input_path, output_path, template_dir):
    pdb_parser = PDB.PDBParser(QUIET=True)
    io = PDB.PDBIO()

    structure = pdb_parser.get_structure("cg_structure", input_path)
    output_structure = PDB.Structure.Structure("rebuilt_structure")

    processed_residues = 0
    skipped_residues = 0
    rmsd_values = []

    for model in structure:
        output_model = PDB.Model.Model(model.id)
        for chain in model:
            output_chain = PDB.Chain.Chain(chain.id)
            existing_ids = set()
            for residue in chain:
                print(f"Residue {residue.resname} at {residue.id} contains atoms: {[atom.get_name() for atom in residue.get_atoms()]}")

                if residue.resname in ["A", "C", "G", "U"]:
                    template_structure = load_template(template_dir, residue.resname)
                    template_residue = next(template_structure.get_residues())


                    if residue.resname in ["A", "G"]:  #purines
                        coarse_atoms = ["N9", "C2", "C6"]
                    elif residue.resname in ["C", "U"]:  #pyrimidines
                        coarse_atoms = ["N1", "C2", "C4"]
                    else:
                        print(f"Skipping unsupported residue: {residue.resname}")
                        skipped_residues += 1
                        continue

                    coarse_atoms += ["P", "C4'", "O5'", "C3'", "O3'"]

                    fixed_atoms = []
                    for atom in residue.get_atoms():
                        if atom.get_name() in coarse_atoms:
                            fixed_atoms.append(atom)

                    template_atoms = []
                    for atom in template_residue.get_atoms():
                        if atom.get_name() in coarse_atoms:
                            template_atoms.append(atom)

                    if len(fixed_atoms) != len(template_atoms):
                        print(f"Warning: Residue {residue.resname} at {residue.id} has mismatched atom counts.")
                        print(f"Fixed atoms ({len(fixed_atoms)}): {[atom.get_name() for atom in fixed_atoms]}")
                        print(f"Template atoms ({len(template_atoms)}): {[atom.get_name() for atom in template_atoms]}")
                        matched_atoms = min(len(fixed_atoms), len(template_atoms))
                        fixed_atoms = fixed_atoms[:matched_atoms]
                        template_atoms = template_atoms[:matched_atoms]

                    if len(fixed_atoms) < 3:
                        print(f"Skipping residue {residue.resname} at {residue.id} due to insufficient atoms for alignment.")
                        skipped_residues += 1
                        continue

                    rotran, rmsd = iteratively_superimpose(fixed_atoms, template_atoms)
                    rmsd_values.append(rmsd)

                    new_residue = PDB.Residue.Residue(residue.id, residue.resname, residue.segid)
                    for atom in template_residue.get_atoms():
                        new_atom = atom.copy()
                        new_atom.transform(rotran[0], rotran[1])
                        new_residue.add(new_atom)

                    for atom in residue.get_atoms():
                        if atom.get_name() not in coarse_atoms:
                            new_residue.add(atom.copy())

                    output_chain.add(new_residue)
                    print(f"Residue {residue.resname} at {residue.id} successfully rebuilt with {len(new_residue)} atoms and RMSD {rmsd:.4f}.")
                    processed_residues += 1
                else:
                    print(f"Unsupported residue: {residue.resname}")
                    skipped_residues += 1

            output_model.add(output_chain)
        output_structure.add(output_model)

    io.set_structure(output_structure)
    io.save(output_path)

    if rmsd_values:
        min_rmsd = min(rmsd_values)
        max_rmsd = max(rmsd_values)
        avg_rmsd = sum(rmsd_values) / len(rmsd_values)
        print(f"RMSD statistics: Min = {min_rmsd:.4f}, Max = {max_rmsd:.4f}, Avg = {avg_rmsd:.4f}")
    print("Rebuilding complete.")
    print(f"Processed residues: {processed_residues}")
    print(f"Skipped residues: {skipped_residues}")

def main():
    parser = argparse.ArgumentParser(description="Rebuild full-atom structure from coarse-grained model using templates.")
    parser.add_argument("--input", required=True, help="Input coarse-grained PDB file.")
    parser.add_argument("--output", required=True, help="Output rebuilt PDB file.")
    parser.add_argument("--templates", required=True, help="Directory containing template PDB files.")

    args = parser.parse_args()

    rebuild_structure(args.input, args.output, args.templates)

if __name__ == "__main__":
    main()
