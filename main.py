import argparse
from Bio import PDB

def filter_atoms(structure):
    cg_atoms = set()

    purine_atoms = {'N9', 'C2', 'C6'}
    pyrimidine_atoms = {'N1', 'C2', 'C4'}
    backbone_atoms = {'P', "C4'"}

    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname()

                if resname in {'G', 'A', 'C', 'U'}:
                    for atom in residue:
                        if atom.get_name() in backbone_atoms:
                            cg_atoms.add(atom)

                    if resname in {'A', 'G'}:
                        for atom in residue:
                            if atom.get_name() in purine_atoms:
                                cg_atoms.add(atom)
                    elif resname in {'C', 'U'}:
                        for atom in residue:
                            if atom.get_name() in pyrimidine_atoms:
                                cg_atoms.add(atom)

    return cg_atoms

def save_structure(structure, atoms, output_filename):
    io = PDB.PDBIO()
    io.set_structure(structure)

    class Select(PDB.Select):
        def accept_atom(self, atom):
            return atom in atoms

    io.save(output_filename, select=Select())

def main():
    parser = argparse.ArgumentParser(
        description="Konwersja pełnoatomowej struktury RNA do reprezentacji gruboziarnistej.")
    parser.add_argument('--input', required=True, help="Plik wejściowy z pełnoatomową strukturą PDB")
    parser.add_argument('--output', required=True, help="Plik wyjściowy z gruboziarnistą strukturą PDB")

    args = parser.parse_args()

    pdb_parser = PDB.PDBParser(QUIET=True)
    structure = pdb_parser.get_structure('RNA', args.input)

    cg_atoms = filter_atoms(structure)

    save_structure(structure, cg_atoms, args.output)
    print(f"Zapisano gruboziarnistą strukturę do pliku: {args.output}")

if __name__ == "__main__":
    main()
