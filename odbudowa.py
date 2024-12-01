import argparse
from Bio.PDB import PDBParser, Superimposer, PDBIO, Residue, Model, Chain
from Bio.PDB.Structure import Structure
import os
from copy import deepcopy
from typing import Dict, List

def duplicate_residue(original: Residue.Residue, new_id, new_segid):
    cloned_residue = Residue.Residue(new_id, original.resname, new_segid)
    for atom in original:
        cloned_residue.add(deepcopy(atom))
    return cloned_residue

def classify_residue_type(res: Residue.Residue):
    if res.get_resname() in ['A', 'G']:
        return 'purine'
    elif res.get_resname() in ['C', 'T', 'U']:
        return 'pyrimidine'
    return None

def parse_structure(file_path: str, structure_id: str = 'cg_model'):
    pdb_parser = PDBParser(QUIET=True)
    return pdb_parser.get_structure(structure_id, file_path)

def read_templates(template_paths: Dict[str, List[str]]) -> Dict[str, List[Structure]]:
    pdb_parser = PDBParser(QUIET=True)
    loaded_templates = {}
    for residue_type, files in template_paths.items():
        loaded_templates[residue_type] = [pdb_parser.get_structure(f"template_{residue_type}_{idx}", file) for idx, file
                                          in enumerate(files)]
    return loaded_templates

def fetch_initial_residue(structure: Structure):
    for model in structure:
        for chain in model:
            for residue in chain:
                return residue

def common_atom_names(res1: Residue.Residue, res2: Residue.Residue):
    atoms1 = {atom.name for atom in res1}
    atoms2 = {atom.name for atom in res2}
    return atoms1 & atoms2

def construct_element(res: Residue.Residue, template_data: Dict[str, List[Structure]], is_backbone: bool = False):
    template_res = fetch_initial_residue(template_data['backbone'][0]) if is_backbone else fetch_initial_residue(
        template_data[res.resname][0])

    new_res = duplicate_residue(template_res, res.id, res.segid)
    shared_atoms = common_atom_names(res, new_res)

    fixed_atoms = [atom for atom in res if atom.name in shared_atoms]
    movable_atoms = [atom for atom in new_res if atom.name in shared_atoms]

    superimposer = Superimposer()
    superimposer.set_atoms(fixed_atoms, movable_atoms)
    superimposer.apply(new_res.get_atoms())

    return new_res

def regenerate_residue(res: Residue.Residue, template_data: Dict[str, List[Structure]]):
    base = construct_element(res, template_data)
    backbone = construct_element(res, template_data, is_backbone=True)

    rebuilt_residue = duplicate_residue(base, res.id, res.segid)
    existing_atoms = {atom.name for atom in rebuilt_residue}

    for atom in backbone.get_atoms():
        if atom.name not in existing_atoms:
            rebuilt_residue.add(atom)

    return rebuilt_residue

def regenerate_structure(input_structure: Structure, templates: Dict[str, List[Structure]],
                         output_structure_id: str = 'rebuilt_model') -> Structure:
    rebuilt_structure = Structure(output_structure_id)

    for model in input_structure:
        rebuilt_model = Model.Model(model.id)
        rebuilt_structure.add(rebuilt_model)

        for chain in model:
            rebuilt_chain = Chain.Chain(chain.id)
            rebuilt_model.add(rebuilt_chain)

            for res in chain:
                if res.resname in templates:
                    rebuilt_residue = regenerate_residue(res, templates)
                    rebuilt_chain.add(rebuilt_residue)
                else:
                    print(f"Residue {res.resname} not supported; copying without changes.")
                    rebuilt_chain.add(duplicate_residue(res, res.id, res.segid))

    return rebuilt_structure

def main():
    arg_parser = argparse.ArgumentParser(
        description="Reconstruct full-atom structure from a coarse-grained model using provided templates.")
    arg_parser.add_argument("--input", required=True, help="Path to the input coarse-grained PDB file.")
    arg_parser.add_argument("--output", required=True, help="Path for saving the reconstructed PDB file.")
    arg_parser.add_argument("--templates", required=True, help="Directory containing PDB template files.")

    args = arg_parser.parse_args()

    template_files = {
        'A': [os.path.join(args.templates, 'A.pdb')],
        'C': [os.path.join(args.templates, 'C.pdb')],
        'G': [os.path.join(args.templates, 'G.pdb')],
        'U': [os.path.join(args.templates, 'U.pdb')],
        'backbone': [os.path.join(args.templates, 'backbone.pdb')]
    }

    templates = read_templates(template_files)
    coarse_grained_structure = parse_structure(args.input)
    rebuilt_structure = regenerate_structure(coarse_grained_structure, templates, args.output)

    io = PDBIO()
    io.set_structure(rebuilt_structure)
    io.save(args.output)


if __name__ == "__main__":
    main()
