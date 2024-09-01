import openmm.app as app
from openmm import *
from openmm.unit import *
import numpy as np
import random
class membrane_generator_master:
    def __init__(self, MainMainWindow):
        self.lipid_files = None
        self.lipid_ratios = None
        self.output_file = None
        self.n_lipids_x = None
        self.n_lipids_y = None
        self.spacing = None
        self.box_z = None
        self.final_box_z = None 
    def generate_bilipid_membrane(self, lipid_files, lipid_ratios, output_file, n_lipids_x, n_lipids_y, spacing, box_z, final_box_z):
        self.lipid_files = lipid_files
        self.lipid_ratios = lipid_ratios
        self.output_file = output_file
        self.n_lipids_x = n_lipids_x
        self.n_lipids_y = n_lipids_y
        self.spacing = spacing
        self.box_z = box_z
        self.final_box_z = final_box_z

        try:
            lipid_pdbs = [app.PDBFile(file) for file in self.lipid_files]
            lipid_topologies = [pdb.topology for pdb in lipid_pdbs]
            lipid_positions_list = [pdb.positions for pdb in lipid_pdbs]
            total_ratio = sum(self.lipid_ratios)
            self.lipid_ratios = [ratio / total_ratio for ratio in self.lipid_ratios]
            combined_positions = []
            combined_topology = app.Topology()
            all_atoms = []
            atom_counter = 0
       
            def add_lipid_to_combined_topology(lipid_choice, translation, is_lower_layer=False):
                nonlocal atom_counter
                lipid_positions = lipid_positions_list[lipid_choice]
                lipid_atoms = list(lipid_topologies[lipid_choice].atoms())
                lipid_bonds = list(lipid_topologies[lipid_choice].bonds())

                # Translate and add positions
                for pos in lipid_positions:
                    new_pos = pos + translation
                    if is_lower_layer:
                        new_pos = Vec3(new_pos[0], new_pos[1], -new_pos[2])  # Flip the z-coordinate to invert the lipid
                        new_pos += Vec3(0, 0, -box_z) * nanometers  # Adjust the position to be in the lower layer, away from the center
                    combined_positions.append(new_pos)

                # Add atoms and bonds to combined topology
                chain = combined_topology.addChain()
                residue = combined_topology.addResidue(lipid_atoms[0].residue.name, chain)
                new_atoms = []
                for atom in lipid_atoms:
                    new_atom = combined_topology.addAtom(atom.name, atom.element, residue)
                    new_atoms.append(new_atom)
                    all_atoms.append(new_atom)

                print(f"Added {len(lipid_atoms)} atoms. Total atoms: {len(all_atoms)}")

                # Add bonds
                for bond in lipid_bonds:
                    atom1 = new_atoms[bond.atom1.index]
                    atom2 = new_atoms[bond.atom2.index]
                    combined_topology.addBond(atom1, atom2)

                atom_counter += len(lipid_atoms)

            # Create the upper layer
            for i in range(int(self.n_lipids_x)):
                for j in range(int(self.n_lipids_y)):
                    print(lipid_ratios)
                    if not lipid_ratios:
                        raise ValueError("The lipid_ratios list is empty")

                    lipid_choice = random.choices(range(len(self.lipid_ratios)), weights=self.lipid_ratios, k=1)[0]
                    translation = Vec3(i * self.spacing, j * self.spacing, self.box_z / 2) * nanometers  # Position the upper layer half the distance
                    add_lipid_to_combined_topology(lipid_choice, translation)

            # Create the lower layer
            for i in range(int(self.n_lipids_x)):
                for j in range(int(self.n_lipids_y)):
                    lipid_choice = random.choices(range(len(self.lipid_ratios)), weights=self.lipid_ratios, k=1)[0]
                    translation = Vec3(i * self.spacing, j * self.spacing, -self.box_z / 2) * nanometers  # Position the lower layer half the distance
                    add_lipid_to_combined_topology(lipid_choice, translation, is_lower_layer=True)

            final_positions = combined_positions

            # Verify the number of positions matches the number of atoms
            print(f"Number of final positions: {len(final_positions)}")
            print(f"Number of all atoms: {len(all_atoms)}")
            if len(final_positions) != len(all_atoms):
                raise ValueError("The number of positions must match the number of atoms")

            # Create a new system with all the lipids
            system = System()
            for pos in final_positions:
                system.addParticle(39.9 * daltons)  # Example mass for a typical lipid atom

            # Define the box dimensions
            system.setDefaultPeriodicBoxVectors(Vec3(self.n_lipids_x * self.spacing, 0, 0) * nanometers,
                                                Vec3(0, self.n_lipids_y * self.spacing, 0) * nanometers,
                                                Vec3(0, 0, self.final_box_z) * nanometers)

            integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
            context = Context(system, integrator)
            context.setPositions(final_positions)

            # Save the new bilipid membrane to a PDB file
            with open(self.output_file, 'w') as f:
                app.PDBFile.writeFile(combined_topology, context.getState(getPositions=True).getPositions(), f)
        except Exception as e:
            print(f"Error in generate_bilipid_membrane: {e}")
