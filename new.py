import openmm.app as app
from openmm import *
from openmm.unit import *
import numpy as np
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
        self.protein_files = None
        self.rotation_matrices = None
        self.membrane_size = None
        self.membrane_file = None
        self.final_box_z = None
        self.output_file = None
        self.protein_spacing = None
    def start_generator(self, lipid_files, lipid_ratios, output_file, n_lipids_x, n_lipids_y, spacing, box_z, final_box_z):
        if protein_spacing or rotation_matrices == None:
            
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
            #load pdbs
            lipid_pdbs = [app.PDBFile(file) for file in self.lipid_files]            
            #grab pdb data
            lipid_topology = pdb.lipid_chosen
            lipid_positions = pdb.lipid_chosen
            #data on which lipid i want
            array_size = int(self.n_lipids_x * self.n_lipids_y)
            lipid_chosen = random.choices(lipid_pdbs, weights=self.lipid_ratios, k=array_size)
            #fill a list with the appropritarte positions
            grid_positions = [
    Vec3(i * self.spacing, j * self.spacing, self.box_z / 2) * nanometers
    for i in range(int(self.n_lipids_x))
    for j in range(int(self.n_lipids_y))
]
            #"spawn in" the lipid i want
            for lipid_pdb, grid_position in zip(lipid_chosen, grid_positions):
                translated_positions = [pos + grid_position for pos in lipid_pdb.positions]
                
                if combined_topology is None:
                    combined_topology = lipid_pdb.topology
                else:
                    combined_topology = combined_topology.merge(lipid_pdb.topology)
                
                combined_positions.extend(translated_positions)
            #i gotta figure out a way to add the lower leaf too bc rn i have upper
            with open('combined_membrane.pdb', 'w') as f:
                PDBFile.writeFile(combined_topology, combined_positions, file=f)
                
                
            if varible == not None:
                #protein_spacing = len(proteinpdbs)
                protein_pdbs = [app.PDBFile(file) for file in self.protein_files]     
                #take oriented rpoteins
            #calculate approptiare plaec and delete the corresponding lipds
                self.membrane_size = (self.n_lipids_x, self.n_lipids_y, self.box_z)
                num_proteins = len(protein_files)
                membrane_center = Vec3(self.membrane_size[0] / 2, self.membrane_size[1] / 2, self.membrane_size[2] / 2) * nanometers
                total_spacing = self.protein_spacing * (num_proteins - 1) * nanometers
                start_x = membrane_center[0] - (total_spacing / 2)
                #translate the protein by membrane_division + and - it 
            #insert it
                for i, rotated_protein_positions in enumerate(centered_rotated_protein_positions_list):
                    z_positions = np.array([pos[2].value_in_unit(nanometers) for pos in rotated_protein_positions])
                    protein_center_z = (z_positions.max() + z_positions.min()) / 2 * nanometers

                    translation_vector = Vec3((start_x + i * spacing * nanometers).value_in_unit(nanometers), 
                                              membrane_center[1].value_in_unit(nanometers), 
                                              (membrane_center[2] - protein_center_z).value_in_unit(nanometers)) * nanometers
                    translated_protein_positions_list.append(translated_positions)
                for protein_pdb, translated_position in zip(protein_pdbs, translated_positions):
                    combined_topology = combined_topology.merge(lipid_pdb.topology)
                    combined_positions.extend(translated_positions)
                    #minimize energy
                    system.setDefaultPeriodicBoxVectors(Vec3(self.membrane_size[0], 0, 0) * nanometers,
                                                    Vec3(0, self.membrane_size[1], 0) * nanometers,
                                                    Vec3(0, 0, self.membrane_size[2]) * nanometers)

                integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
                context = Context(system, integrator)
                context.setPositions(combined_positions)
                with open('membrane_with_protein', 'w') as f:
                    PDBFile.writeFile(combined_topology, combined_positions, file=f)
                
           
    
    

            #delete membrane liipds
                #find the number atom
#                         so probabily make a copy of the membrane file,
#                         delete the mathcing clash text,
#                         if its not found then whatever,
#                         resave the membrane file,
#                         add the protein in the same exact spot as the ones in the og file
        except Exception as e:
            print(f"Error in generate_bilipid_membrane: {e}")

