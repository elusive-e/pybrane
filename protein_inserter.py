import openmm.app as app
from openmm import *
from openmm.unit import *
import numpy as np
class protein_inserter_master:
    def __init__(self, MainMainWindow):
        self.output_file = None
        self.n_lipids_x = None
        self.n_lipids_y = None
        self.spacing = None
        self.box_z = None
        self.final_box_z = None
        self.protein_files = None
        self.rotation_matrices = None
        self.membrane_size = None
        self.membrane_file = None
        self.final_box_z = None
        self.output_file = None
        self.protein_spacing = None
    def load_pdb(file):
        return app.PDBFile(file)

    def center_protein(protein_positions):
        positions_array = np.array([pos.value_in_unit(nanometers) for pos in protein_positions])
        centroid = np.mean(positions_array, axis=0)
        centered_positions = positions_array - centroid
        centered_positions = [Vec3(*pos) * nanometers for pos in centered_positions]
        return centered_positions

    def rotate_protein(protein_positions, rotation_matrix):
        positions_array = np.array([pos.value_in_unit(nanometers) for pos in protein_positions])
        rotated_positions = np.dot(positions_array, rotation_matrix.T)
        rotated_positions = [Vec3(*pos) * nanometers for pos in rotated_positions]
        return rotated_positions

    def translate_protein(protein_positions, translation_vector):
        translated_positions = [Vec3(pos[0] + translation_vector[0], pos[1] + translation_vector[1], pos[2] + translation_vector[2]) for pos in protein_positions]
        return translated_positions

    def combine_structures(membrane_topology, membrane_positions, protein_topologies, protein_positions_list):
        combined_topology = app.Topology()
        combined_positions = []

        membrane_atoms = list(membrane_topology.atoms())
        for atom, pos in zip(membrane_atoms, membrane_positions):
            chain = combined_topology.addChain()
            residue = combined_topology.addResidue(atom.residue.name, chain)
            combined_topology.addAtom(atom.name, atom.element, residue)
            combined_positions.append(pos)

        for protein_topology, protein_positions in zip(protein_topologies, protein_positions_list):
            protein_atoms = list(protein_topology.atoms())
            for atom, pos in zip(protein_atoms, protein_positions):
                chain = combined_topology.addChain()
                residue = combined_topology.addResidue(atom.residue.name, chain)
                combined_topology.addAtom(atom.name, atom.element, residue)
                combined_positions.append(pos)

        return combined_topology, combined_positions

    def orient_proteins_in_membrane(protein_files, membrane_file, output_file, rotation_matrices, membrane_size, spacing):
        self.output_file = output_file
        self.n_lipids_x = n_lipids_x
        self.n_lipids_y = n_lipids_y
        self.spacing = spacing
        self.box_z = box_z
        self.final_box_z = final_box_z
        self.protein_files = []
        self.protein_spacing = self.pro_spin_ratio.value()
        self.rotation_matrices = [np.eye(3) for _ in self.protein_files]
        self.membrane_size = (self.n_lipids_x, self.n_lipids_y, self.box_z)
        self.membrane_file = 'membrane_maker_output.pdb'
        self.final_box_z = self.main_window.final_z_spin.value()
        self.protein_choice = self.main_window.radioButton.isChecked()
        self.output_file = 'membrane_maker_output.pdb'
        membrane_pdb = load_pdb(membrane_file)

        protein_pdbs = [load_pdb(file) for file in self.protein_files]
        protein_positions_list = [pdb.positions for pdb in protein_pdbs]

        # Center and rotate each protein
        centered_rotated_protein_positions_list = []
        for protein_positions, rotation_matrix in zip(protein_positions_list, self.rotation_matrices):
            centered_positions = center_protein(protein_positions)
            rotated_positions = rotate_protein(centered_positions, rotation_matrix)
            centered_rotated_protein_positions_list.append(rotated_positions)

        # Calculate translation vectors to space proteins evenly in the membrane
        num_proteins = len(protein_files)
        membrane_center = Vec3(self.membrane_size[0] / 2, self.membrane_size[1] / 2, self.membrane_size[2] / 2) * nanometers
        total_spacing = self.protein_spacing * (num_proteins - 1) * nanometers
        start_x = membrane_center[0] - (total_spacing / 2)

        translated_protein_positions_list = []
        for i, rotated_protein_positions in enumerate(centered_rotated_protein_positions_list):
            z_positions = np.array([pos[2].value_in_unit(nanometers) for pos in rotated_protein_positions])
            protein_center_z = (z_positions.max() + z_positions.min()) / 2 * nanometers

            translation_vector = Vec3((start_x + i * spacing * nanometers).value_in_unit(nanometers), 
                                      membrane_center[1].value_in_unit(nanometers), 
                                      (membrane_center[2] - protein_center_z).value_in_unit(nanometers)) * nanometers
            translated_positions = translate_protein(rotated_protein_positions, translation_vector)
            translated_protein_positions_list.append(translated_positions)
            #print(f"Protein {i+1} translation vector: {translation_vector}")  # Debug print

        # Combine the membrane and protein structures
        combined_topology, combined_positions = combine_structures(
            membrane_pdb.topology, membrane_pdb.positions,
            [pdb.topology for pdb in protein_pdbs], translated_protein_positions_list
        )

        # Create a new system with all the lipids and protein
        system = System()
        for pos in combined_positions:
            system.addParticle(39.9 * daltons)  # Example mass for a typical atom

        # Define the box dimensions
        system.setDefaultPeriodicBoxVectors(Vec3(self.membrane_size[0], 0, 0) * nanometers,
                                            Vec3(0, self.membrane_size[1], 0) * nanometers,
                                            Vec3(0, 0, self.membrane_size[2]) * nanometers)

        integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
        context = Context(system, integrator)
        context.setPositions(combined_positions)

        # Save the new combined structure to a PDB file
        with open(self.output_file, 'w') as f:
            app.PDBFile.writeFile(combined_topology, context.getState(getPositions=True).getPositions(), f)



