from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QApplication, QWidget, QPushButton, QLineEdit, QVBoxLayout, QDialog, QDialogButtonBox, QMessageBox
from PyQt5 import QtOpenGL
from PyQt5 import QtWidgets
from OpenGL import GLU
#from openmm.app import app
from openmm import *
from openmm.unit import *
from openmm.app import *
from openmm import *
from pdbfixer import PDBFixer
from openmm.unit import *
from sys import stdout
import pandas as pd
import MDAnalysis
import numpy as np
import random
import sys
import re
import subprocess
import pdbfixer
import requests
#import yaml
from opengl_generic import genericOpenGLWidget
from protein_inserter import protein_inserter_master
from membrane_maker import membrane_generator_master
from pybraneui12 import Ui_MainWindow
import webbrowser
import sys
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import QDialog, QApplication, QPushButton, QVBoxLayout
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import random
from rdkit.Chem import rdDepictor as rdd
from rdkit.Chem.Draw import rdMolDraw2D as draw2D
from rdkit import Chem
from rdkit.Chem.Draw import DrawingOptions
from rdkit.Chem import Draw, AllChem
import wikipedia 
from chemspipy import ChemSpider
import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
from michaelai import michaelai
matplotlib.use('Qt5Agg')
import mdtraj as md



class MainMainWindow(QtWidgets.QMainWindow, Ui_MainWindow):
    def __init__(self):
        super(MainMainWindow, self).__init__()
        self.setupUi(self)
        self.mmmm = self.molecule_maker_tab(self)
        self.mf = self.MenuFunctionsss(self)
        self.saf = self.search_and_files(self)
        self.mem_sim = self.MemMaker(self)
        self.ana = self.analysis_tab(self)
        self.simsim = self.simulation_tab(self)
        self.mem_load_but1.clicked.connect(self.mem_sim.load_file_mem_maker)
        self.start_mem_make_but.clicked.connect(self.mem_sim.start_membrane_maker)
        self.pushButton.clicked.connect(self.activate_michael)
        self.pushButton_16.clicked.connect(self.execute_command)
        self.pro_load_but1.clicked.connect(self.mem_sim.load_file_protein)
        
        self.pushButton_14.clicked.connect(self.mmmm.start_bonds)
        self.pushButton_13.clicked.connect(self.mmmm.start_bonds)
        self.pushButton_12.clicked.connect(self.mmmm.start_atoms)
        self.pushButton_11.clicked.connect(self.mmmm.smiles_import)
        self.pushButton_9.clicked.connect(self.mmmm.load_pdbmmmm)
        self.pushButton_10.clicked.connect(self.mmmm.save_file)

        self.pushButton_17.clicked.connect(self.simulation_tab.load_sim_ff)
        self.pushButton_18.clicked.connect(self.simulation_tab.load_sim_pro)
        self.pushButton_19.clicked.connect(self.simulation_tab.load_sim_mem)
        self.pushButton_20.clicked.connect(self.simulation_tab.start_simulation)
        self.pushButton_21.clicked.connect(self.simulation_tab.load_sim_thing)
        self.pushButton_22.clicked.connect(self.simulation_tab.load_sim_sol)
        
        self.pushButton_23.clicked.connect(self.saf.search_lipid_maps)
        self.pushButton_24.clicked.connect(self.saf.search_pdb)
        self.pushButton_25.clicked.connect(self.saf.search_wiki)
        self.pushButton_26.clicked.connect(self.saf.search_chemspi)
        
        self.pushButton_27.clicked.connect(self.mmmm.switch_display)
        self.opengl_with_buttons.Button1.clicked.connect(self.ana.open_viewer)
        self.opengl_with_buttons.Button2.clicked.connect(self.ana.update_graphs)
        self.opengl_with_buttons.Button3.clicked.connect(self.ana.settings_graphs)
        self.actionNew.triggered.connect(self.mf.openNew)
        self.actionOpen.triggered.connect(self.mf.openOpen)
        self.actionSave.triggered.connect(self.mf.openSave)
        self.actionSave_as.triggered.connect(self.mf.openSaveAs)
        self.actionQuit.triggered.connect(self.mf.openQuit)
        self.actionGraphs.triggered.connect(self.mf.adjustGraphs)
        self.actionViewer.triggered.connect(self.mf.adjustViewer)
        self.actionInsert_CS_key.triggered.connect(self.mf.acceptKey)
        self.actionSource_Code.triggered.connect(self.mf.openCode)
        self.actionAbout_PyBRANE.triggered.connect(self.mf.openDes)
        self.actionMaximize.triggered.connect(self.mf.makeBig)
        self.actionMinimize.triggered.connect(self.mf.makeSmall)
        self.actionOpen_Mol_Viewer.triggered.connect(self.mf.molView)
        self.actionOpen_Mol_Editor.triggered.connect(self.mf.molEdit)

        DrawingOptions.bondLineWidth=1.8
        DrawingOptions.atomLabelFontSize=14
        DrawingOptions.includeAtomNumbers=True
        
        twod_view = True
    class MenuFunctionsss:
        def __init__(self, MainMainWindow):
            self.MainMainWindow = MainMainWindow
        def openSaveAs(self):
            tabNumber = window.tabWidget.currentIndex()
            if tabNumber == 0:
                window.mmmm.save_fileAS()
            if tabNumber == 1:
                window.textBrowser_13.append("NOTICE: There is no current way to save the settings. The membrane will save automatically once it is completed. Please keep monitoring the text field for output updates.")
            if tabNumber == 2:
                window.textBrowser_13.append("NOTICE: There is no current way to save the settings. The simulation will save automatically once it is completed. Please keep monitoring the text field for output updates.")
            if tabNumber == 3:
                window.ana.save_fileAS()
        def openSave(self):
            tabNumber = window.tabWidget.currentIndex()
            if tabNumber == 0:
                window.mmmm.save_file()
            if tabNumber == 1:
                window.textBrowser_13.append("NOTICE: There is no current way to save the settings. The membrane will save automatically once it is completed. Please keep monitoring the text field for output updates.")
            if tabNumber == 2:
                window.textBrowser_13.append("NOTICE: There is no current way to save the settings. The simulation will save automatically once it is completed. Please keep monitoring the text field for output updates.")
            if tabNumber == 3:
                window.ana.save_file()
        def openOpen(self):
            pass
            #see tab index, update the appropriate thing to be opened (maybe open optional screen to ask where they want it opened?)
        def openNew(self):
            tabNumber = window.tabWidget.currentIndex()
            if tabNumber == 0:
                pass
                #clear the file that's made (blank pdb with a hydrogen prob)
            if tabNumber == 1:
                pass 
                # stop the simulation if theres one running and clear all the settings
            if tabNumber == 2:
                pass
                #clear all the settings and stop the simulaiton idk if theres a way to do that maybe like clear all the subprocesses
            if tabNumber == 3:
                pass
                #clear the file and make a blank thingie
            pass
        def openQuit(self):
            tabNumber = window.tabWidget.currentIndex()
            if tabNumber == 1 or tabNumber == 2:
                window.textBrowser_13.append("WARNING: Please verify all proccesses have finished and swtich to the analysis or molecule editor tab to close window.")
            else:
                window.close()
        def openCode(self):
            webbrowser.open('https://github.com/elusive-e/pybrane/')
        def openDes(self):
            webbrowser.open('https://github.com/elusive-e/pybrane/wiki')
        def makeSmall(self):
            window.showMinimized()
        def makeBig(self):
            window.showMaximized()
        def acceptKey(self):
            # TO DO: open the settings window
            cswindow.show()
        def adjustGraphs(self):
            #launch the graphs setting window, continue as normal
            # TO DO: open the settings window
            pass
        def adjustViewer(self):
            #launch the viewer settngs with style changes and stuff
            pass
        def molView(self):
            #window.tabWidget.setCurrentIndex(1)
            viewerwindow.show()
            #add actual names 
            pass
        def molEdit(self):
            #window.tabWidget.setCurrentIndex(1)
            pass
    class search_and_files:
        def __init__(self, MainMainWindow):
            self.MainMainWindow = MainMainWindow
            self.pdb_request = MainMainWindow.lineEdit_3.text().strip()
            self.cs = None
            
        def search_lipid_maps(self):
            try:
                urls = []
                url = f"https://www.lipidmaps.org/rest/compound/pubchem_cid/{self.pdb_request}"
                urls.append(url)
                url = f"https://www.lipidmaps.org/rest/compound/formula/{self.pdb_request}/smiles"
                urls.append(url)
                url = f"https://www.lipidmaps.org/rest/compound/formula/{self.pdb_request}/pubchem_cid"
                urls.append(url)
                url = f"https://www.lipidmaps.org/rest/compound/formula/{self.pdb_request}/all"
                urls.append(url)
                url = f"https://www.lipidmaps.org/rest/compound/formula/{self.pdb_request}/molfile"
                urls.append(url)
                url = f"https://www.lipidmaps.org/rest/compound/smiles/{self.pdb_request}/all"
                urls.append(url)
                url = f"https://www.lipidmaps.org/rest/protein/uniprot_id/{self.pdb_request}/all"
                urls.append(url)
                url = f"https://www.lipidmaps.org/rest/protein/protein_name/{self.pdb_request}/all"
                urls.append(url)
                for url in urls:
                    response = requests.get(url)
                    data = response.json()
                    with open(f"{entry_id}.json", 'w') as file:
                        file.write(response.text)
                window.textBrowser_18.append(f"Your lipid maps request yieled results saved in the file {entry_id}.json. Please be advised not all the results may be what you're looking for, as the query was ran through smiles, formula, and other ids, which will not match perfectly.")
            except:
                window.textBrowser_18.append(f"ERROR: Your lipid maps search resulted in 0 results or encountered an error.")

        def search_chemspi(self):
            data = []
            pdb_request = window.lineEdit_3.text().strip()
            if self.cs == None:
                self.cs = ChemSpider('3G9Vc3ehSXmA8dggOLpbNkxMYbYwQsWe')
            r = self.cs.search(self.pdb_request)
            try:
                try:
                    for result in self.cs.search(self.pdb_request):
                        with open(f"{entry_id}.json", 'w') as file:
                            file.write(result.text)
                    data.append(r.message)
                except:
                    window.textBrowser_18.append("ERROR: Your query did not return a result")
                else:
                    chemspi = cs.get_compound(pdb_request)
                    data.append(chemspi.molecular_formula)
                    data.append(chemspi.molecular_weight)
                    data.append(chemspi.smiles)
                    
                    data.append(c.search_by_mass)
                    data.append(c.search_by_formula)
                    for i in range(len(data)):
                        with open(f"{entry_id}.json", 'w') as file:
                                file.write(data.text)
                    window.textBrowser_18.append(f"Your Chem Spider request yieled results saved in the file {entry_id}.json. Please be advised not all the results may be what you're looking for, as the query was ran through smiles, formula, and other ids, which will not match perfectly. If it is the result your are lookign for, please save the file with a unique name so it is not overwritten if this query is searched again.")
               
            except:
                window.textBrowser_18.append(f"ERROR: Your Chem Spider request encountered an error or found 0 results. Please try again.")

        def search_pdb(self):
            try:
                self.pdb_request = window.lineEdit_3.text().strip()
                url = f"https://files.rcsb.org/download/{self.pdb_request}.cif"

                url2 = f"https://data.rcsb.org/rest/v1/core/entry/{self.pdb_request}"

                response = requests.get(url2)

                data = response.json()
                data.keys()
                with open(f"{self.pdb_request}.cif", 'w') as file:
                    file.write(response.text)
                window.textBrowser_18.append(f"Protein Data Bank Search resulted in a succesfu; retrival of info and structure. Structure file saved as {self.pdb_request}.cif")
            except:
                window.textBrowser_18.append("ERROR: Your Protein Data BAnk Search resulted in 0 results or the search has encountered an error.")
        def search_wiki(self):
            self.pdb_request = window.lineEdit_3.text().strip()
            pdb_request = self.pdb_request
            try:
                output = wikipedia.summary(pdb_request)
            except:
                pdb_request = wikipedia.suggest(pdb_request)
                output = wikipedia.summary(pdb_request)
            window.textBrowser_18.append(f"--> Wikipedia result for {pdb_request}:")
            window.textBrowser_18.append(output)
        def update_cskey(self):
            self.cs = cswindow.inputField.text.strip()
    class analysis_tab:
        def __init__(self, MainMainWindow):
            self.MainMainWindow = MainMainWindow
            self.atom_positions=[]
        def zoom_in_max(self):
            viewerwindow.openGLWidget.zoomin()
        def auto_rotate(self):
            viewerwindow.openGLWidget.animatefun(self.atom_positions)
        def open_viewer(self):
            viewerwindow.show()
        def zoom_out_max(self):
            viewerwindow.openGLWidget.zoomout()
        def change_style(self):
#           for bond in mol.GetBonds():
#     			start_atom = bond.GetBeginAtomIdx()
#    			end_atom = bond.GetEndAtomIdx()
#    			bond_type = bond.GetBondTypeAsDouble()
            pass
        def load_file(self):
            file_path, _ = QtWidgets.QFileDialog.getOpenFileName(window, "Open File", "", "All Files (*)")
            if file_path:
                with open(file_path, 'r') as file:
                    content = file.read()
                self.mol_input = Chem.MolFromPDBBlock(content)
                if self.mol_input.GetNumConformers() == 0:
                    self.mol_input = Chem.Mol(self.mol_input)
                    self.mol_input = AllChem.AddHs(self.mol_input, addCoords=True)
                    AllChem.EmbedMolecule(self.mol_input, AllChem.ETKDG())
                    AllChem.MMFFOptimizeMolecule(self.mol_input)
                mol_input = Chem.AddHs(self.mol_input)
                AllChem.EmbedMolecule(self.mol_input)
                conf = self.mol_input.GetConformer()
                atom_positions = []
                for atom in self.mol_input.GetAtoms():
                    pos = conf.GetAtomPosition(atom.GetIdx())
                    atom_positions.append({
                            'element': atom.GetSymbol(),
                            'radius': Chem.GetPeriodicTable().GetRvdw(atom.GetAtomicNum()),
                            'position': np.array([pos.x, pos.y, pos.z])
                        })
                print(atom_positions)
                self.atom_positions = atom_positions
                viewerwindow.openGLWidget.setFocus()
                viewerwindow.openGLWidget.raise_()
                viewerwindow.openGLWidget.set_coordinates(atom_positions) 
                
            else:
                window.textBrowser_13.append("WARNING: The file chosen is not the correct format or has not been picked.")
        def save_file(self):
            self.mol_input = Chem.Mol(self.mol_input)
            content = Chem.rdmolfiles.MolToPDBBlock(self.mol_input)
            with open(f"molecule_analysis_export.pdb", 'w') as file:
                        file.write(content)
        def minimize_viewer(self):
            viewerwindow.hide()
        def settings_graphs(self):
            #add settings logic
            pass
        def expand_viewer(self):
            viewerwindow.showMaximized()
        def save_fileAS(self):
            mol_input = self.mol_input
            mol_input = Chem.MolToPDBBlock(self.mol_input)
            file_name, _ = QFileDialog.getSaveFileName(self,"Save File","","All Files(*);;Text Files(*.txt)",options = options)
            if file_name:
                f = open(file_name, 'w')
                f.write(mol_input)
                self.setWindowTitle(str(os.path.basename(file_name)))
                f.close()
            else:
                window.Text_browser13.append("WARNING: No file name provided, file HAS NOT BEEN SAVED. Please try again.")
        def update_graphs(self, output):            
            def update_radiusOfGyration(self):
                data = md.compute_rg(traj)
                parameter = 0
                window.genericGraph_4.update(data, parameter)
           ####
            def update_centerOfGeometry(self):
                data = md.compute_rg(traj)
                parameter = 1
                window.genericGraph_5.update(data, parameter)
            def update_totalMass(self):
                data = md.compute_rg(traj)
                parameter = 2
                window.genericGraph_3.update(data, parameter)
            def update_centerOfMass(self):
                data = md.compute_rg(traj)
                parameter = 3
                window.genericGraph_1.update(data, parameter)
            ####
            def update_totalCharge(self):
                data = md.compute_rg(traj)
                parameter = 4
                window.genericGraph_5.update(data, parameter)
            def update_rmsf(self):
                data = md.compute_rg(traj)
                parameter = 5
                window.genericGraph_5.update(data, parameter)
            def update_gyration_tensor(self):
                data = md.compute_gyration_tensor(traj)
                parameter = 6
                window.genericGraph_5.update(data, parameter)
            def update_rmsf(self):
                data = md.asphericity(traj)
                parameter = 7
                window.genericGraph_5.update(data, parameter)
            def update_gyration_tensor(self):
                data = md.acylindricity(traj)
                parameter = 8
                window.genericGraph_5.update(data, parameter)
            def update_rmsf(self):
                data = md.principal_moments(traj)
                parameter = 9
                window.genericGraph_5.update(data, parameter)
            def update_gyration_tensor(self):
                data = md.relative_shape_antisotropy(traj)
                parameter = 10
                window.genericGraph_5.update(data, parameter)
            def update_rmsf(self):
                data = md.principal_moments(traj)
                parameter = 11
                window.genericGraph_5.update(data, parameter)
            try:
                traj = md.load('trajectory.dcd',top='output.pdb')
                # add logic that gets the stuff from settings
                if user_graph_choice == 1:
                    update_radiusOfGyration(self)
                elif user_graph_choice == 2:
                    update_centerOfGeometry(self)
                elif user_graph_choice == 3:
                    update_totalMass(self)
                elif user_graph_choice == 4:
                    update_centerOfMass(self)
                elif user_graph_choice == 5:
                    update_totalCharge(self)
                
            except Exception as e:
                window.textBrowser_13.append(f"ERROR: {e}")
            
    class simulation_tab:
        def __init__(self, MainMainWindow):
            self.MainMainWindow = MainMainWindow
            #super(MainMainWindow, self).__init__()
            self.forcefield = None
            self.pdb = None
            self.mem = None
            self.macro_other = []
            
        def start_simulation(self):
            forcefield = window.simsim.forcefield
            pdb = window.simsim.pdb
            
            window.textBrowser_13.append(f"NOTICE: The simulation was started. Inputs: {pdb, forcefield}")
            try: 
                modeller = Modeller(pdb.topology, pdb.positions)
                mem = window.simsim.mem
                if mem == None:
                    window.textBrowser_13.append("NOTICE: No mmebrane was selected, creating a default.")
                    modeller.addMembrane(forcefield, lipidType='POPC', minimumPadding=1*nanometer)
                else:
                    modeller.add(mem.topology, mem.positions)
                    modeller.addSolvent(forcefield, model='tip3p', padding=1.0*nanometers)
                print(f"Number of atoms in topology: {pdb.topology.atoms()}")
                print(f"Number of positions: {pdb.positions}")
                window.textBrowser_13.append("NOTICE: The modeller was created.")

                system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)

                window.textBrowser_13.append("NOTICE: The system was created.")
                print(f"Number of atoms in modeller: {modeller.topology.getNumAtoms()}")
                print(f"Number of particles in system: {system.getNumParticles()}")

                integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
                window.textBrowser_13.append("NOTICE: The intergator was created.")
                simulation = Simulation(modeller.topology, system, integrator)
                window.textBrowser_13.append("NOTICE: The simulation was created.")
                simulation.context.setPositions(modeller.positions)
                window.textBrowser_13.append("NOTICE: The position were set.")
                simulation.minimizeEnergy()
                window.textBrowser_13.append("NOTICE: Energy was minimized.")
                simulation.context.setVelocitiesToTemperature(300*kelvin)
                window.textBrowser_13.append("NOTICE: The temperature was set.")
                
                simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))
                simulation.reporters.append(DCDReporter('trajectory.dcd', 1000))
                
                window.textBrowser_13.append("NOTICE: The simulation is running...")

                simulation.step(10000)
                positions = simulation.context.getState(getPositions=True).getPositions()
                with open('output.pdb', 'w') as output_file:
                    PDBFile.writeFile(simulation.topology, positions, output_file)
                window.textBrowser_13.append("NOTICE: The simulation is finished.")
            except Exception as e:
                window.textBrowser_13.append(f"ERROR: {e}")
        def load_sim_ff(self):
            forcefield = ForceField('charmm36.xml','charmm36/water.xml')            #forcefield.loadFile('lipid17.xml')
           # forcefield.loadFile('tip3p-pme-b.xml')
            window.simsim.forcefield = forcefield
            window.textBrowser_13.append("NOTICE: Default forcefield was loaded, you may now pick custom force field or cancel the operation to use default.")
            file_path, _ = QtWidgets.QFileDialog.getOpenFileName(window, "Open File", "", "All Files (*)")
            if file_path:
                with open(file_path, 'r') as file:
                    content = file.read()
                forcefield = ForceField(file_path)
                window.simsim.forcefield = forcefield
                window.textBrowser_13.append("Forcefield file was succesfully loaded.")
            else:
                window.simsim.forcefield = forcefield
                window.textBrowser_13.append("WARNING: The file chosen is not the correct format or has not been picked, the default forcefield is being invoked.")
            
        def load_sim_pro(self):
            file_path, _ = QtWidgets.QFileDialog.getOpenFileName(window, "Open File", "", "All Files (*)")
            if file_path:
                with open(file_path, 'r') as file:
                    content = file.read()
                fixer = PDBFixer(file_path)
                fixer.findMissingResidues()
                fixer.findMissingAtoms()
                fixer.addMissingAtoms()
                fixer.addMissingHydrogens()
                with open('fixed_structure.pdb', 'w') as f:
                    PDBFile.writeFile(fixer.topology, fixer.positions, f)
                pdb = PDBFile('fixed_structure.pdb')
                window.simsim.pdb = pdb
                window.textBrowser_13.append("PDB file was succesfully loaded.")
            else:
                window.textBrowser_13.append("WARNING: The file chosen is not the correct format or has not been picked.")
    
        def load_sim_mem(self):
            file_path, _ = QtWidgets.QFileDialog.getOpenFileName(window, "Open File", "", "All Files (*)")
            if file_path:
                with open(file_path, 'r') as file:
                    content = file.read()
                fixer = PDBFixer(file_path)
                window.textBrowser_13.append("fixer was activated :)")
    
                fixer.findMissingResidues()
                window.textBrowser_13.append("we found nemo (missing residues) :)")
                fixer.findMissingAtoms()
                window.textBrowser_13.append("we found dory (missing atoms) :)")
                fixer.addMissingAtoms()
                window.textBrowser_13.append("oh helm yeah we added the missing atoms :)")
                #fixer.addMissingHydrogens()
                window.textBrowser_13.append("we added the hydrogens! chat :)")
                with open('fixed_mem_structure.pdb', 'w') as f:
                    PDBFile.writeFile(fixer.topology, fixer.positions, f)
                mem = PDBFile('fixed_mem_structure.pdb')
                window.simsim.mem = mem
                window.textBrowser_13.append("Membrane file was succesfully loaded.")
            else:
                window.textBrowser_13.append("WARNING: The file chosen is not the correct format or has not been picked.")
        
        def load_sim_thing(self):
            file_path, _ = QtWidgets.QFileDialog.getOpenFileName(window, "Open File", "", "All Files (*)")
            if file_path:
                with open(file_path, 'r') as file:
                    content = file.read()
                macro = PDBFile(file_path)
                macro_other.append(macro)
                window.simsim.macro_other = macro_other
                window.textBrowser_13.append("Macromolecule file was succesfully loaded.")
            else:
                window.textBrowser_13.append("WARNING: The file chosen is not the correct format or has not been picked.")
                     
        def load_sim_sol(self):
            file_path, _ = QtWidgets.QFileDialog.getOpenFileName(window, "Open File", "", "All Files (*)")
            if file_path:
                with open(file_path, 'r') as file:
                    content = file.read()
                sol = PDBFile(file_path)
                window.simsim.sol = sol
                window.textBrowser_13.append("Solvent file was succesfully loaded.")
            
            else:
                window.textBrowser_13.append("WARNING: The file chosen is not the correct format or has not been picked.")
    
    
    class molecule_maker_tab:
        def __init__(self, MainMainWindow):
            self.MainMainWindow = MainMainWindow
            self.twod_view = False
            self.mol_input = None
            DrawingOptions.bondLineWidth=1.2
            DrawingOptions.atomLabelFontSize=10
            DrawingOptions.includeAtomNumbers=True
        def save_file(self):
            self.mol_input = Chem.Mol(self.mol_input)
            content = Chem.rdmolfiles.MolToPDBBlock(self.mol_input)
            with open(f"molecule_maker_export.pdb", 'w') as file:
                        file.write(content) 
        def save_fileAS(self):
            mol_input = self.mol_input
            mol_input = Chem.MolToPDBBlock(self.mol_input)
            file_name, _ = QFileDialog.getSaveFileName(self,"Save File","","All Files(*);;Text Files(*.txt)",options = options)
            if file_name:
                f = open(file_name, 'w')
                f.write(mol_input)
                self.setWindowTitle(str(os.path.basename(file_name)))
                f.close()
            else:
                window.Text_browser13.append("WARNING: No file name provided, file HAS NOT BEEN SAVED. Please try again.")
        def switch_display(self, twod_view):
            self.twod_view = not self.twod_view
            
            mol_input = self.mol_input
            print(self.mol_input)
            self.update_display(mol_input)
            
                
        def update_display(self, mol_input):
            self.mol_input = mol_input
            #self.mol_input = self.mol_input
            DrawingOptions.bondLineWidth=1.2   
            DrawingOptions.atomLabelFontSize=10
            DrawingOptions.includeAtomNumbers=True
            if self.twod_view == True:
                print(self.mol_input)
                window.tabWidget.setCurrentIndex(1)
                backup_mol = self.mol_input
                print(backup_mol)
                if self.mol_input.GetNumConformers() == 0:
                    AllChem.Compute2DCoords(self.mol_input)
                print(backup_mol)
                drawer = Draw.rdMolDraw2D.MolDraw2DSVG(300, 300)
                Draw.rdMolDraw2D.PrepareAndDrawMolecule(drawer, self.mol_input)
                drawer.FinishDrawing()
                svg = drawer.GetDrawingText()
                print(self.mol_input)
                with open('initial_mol.svg', 'w') as f:
                    f.write(svg)
                print(backup_mol)
                window.graphicsView.load('initial_mol.svg')
                self.mol_input = backup_mol
            if self.twod_view == False:
                print(self.mol_input)
                try:
                    window.tabWidget.setCurrentIndex(0)
                    if self.mol_input.GetNumConformers() == 0:
                        self.mol_input = Chem.Mol(self.mol_input)
                        self.mol_input = AllChem.AddHs(self.mol_input, addCoords=True)
                        AllChem.EmbedMolecule(self.mol_input, AllChem.ETKDG())
                        AllChem.MMFFOptimizeMolecule(self.mol_input)
                    mol_input = Chem.AddHs(self.mol_input)
                    print(self.mol_input)
                    AllChem.EmbedMolecule(self.mol_input)
                    print(self.mol_input)
                    conf = self.mol_input.GetConformer()
                    atom_positions = []
                    for atom in self.mol_input.GetAtoms():
                        pos = conf.GetAtomPosition(atom.GetIdx())
                        atom_positions.append({
                            'element': atom.GetSymbol(),
                            'radius': Chem.GetPeriodicTable().GetRvdw(atom.GetAtomicNum()),
                            'position': np.array([pos.x, pos.y, pos.z])
                        })
                    print(atom_positions)
                    self.atom_positions = atom_positions
                    window.openGLWidget.setFocus()
                    window.openGLWidget.raise_()
                    window.openGLWidget.set_coordinates(atom_positions)
                except:
                    window.textBrowser_18.append("ERROR: Your molecule could not be visualized. Please amke sure there are no errors in the file.")
        def load_pdbmmmm(self, mol_input):
            file_path, _ = QtWidgets.QFileDialog.getOpenFileName(window, "Open File", "", "All Files (*)")
            if file_path:
                with open(file_path, 'r') as file:
                    content = file.read()
                mol_input = Chem.MolFromPDBBlock(content)
                window.molecule_maker_tab.update_display(self, mol_input)
            
            else:
                window.textBrowser_13.append("WARNING: The file chosen is not the correct format or has not been picked.")
    
        def smiles_import(self):
            smilesdialog.show()
            
        def smiles_read(self):
            smiles_input = smilesdialog.inputField.text().strip()
            mol_input = Chem.MolFromSmiles(smiles_input)
            window.mmmm.update_display(mol_input)
        def start_atoms(self):
            atomdialog.show()
#         def add_atom(self):
#             rw_user_mol = Chem.RWMol(self.mol_input)
#             user_atom_number = self.spinBox_5.value()
# 
#             new_atom = rw_user_mol.AddAtom(Chem.Atom(user_atom_number))
#             self.update_display(self, self.mol_input)
#             rw_user_mol = Chem.RWMol(self.mol_input)
#             user_atom_index = self.spinBox_6.value()
#             rw_user_mol.RemoveAtom(user_atom_index)
# 
#             mol_input = rw_user_mol.GetMol()
#             self.update_display(self, mol_input)
            
        def start_bonds(self):
            bonddialog.show()
            
        def manage_bonds(self):
            atom_1 = bonddialog.spinBox.value()
            atom_2 = bonddialog.spinBox_2.value()
            addbond = bonddialog.radioButton.isChecked()
            removebond = bonddialog.radioButton_2.isChecked()
            bond_choice = bonddialog.comboBox.choice()
            rw_user_mol = Chem.RWMol(self.mol_input)
            try:
                if addbond and removebond == True:
                    window.textBox_18.append("WARNING: You must pick one bond action, not both.")
                if addbond == True:
                    rw_user_mol.AddBond(atom_1,atom_2,order=Chem.rdchem.BondType.bond_choice)
                if removebond == True:
                    rw_user_mol.RemoveBond(atom_1,atom_2,order=Chem.rdchem.BondType.bond_choice)
                else:
                    window.textBox_18.append("WARNING: You must pick one bond action, not none")
                self.mol_input = rw_user_mol
                self.update_display()
            except:
                winodw.textBrowser_18.append("ERROR: The bond picked could not be updated.")
        def manage_atoms(self):
            atom_1 = atomdialog.spinBox.value()
            addatom = atomdialog.radioButton.isChecked()
            removeatom = atomdialog.radioButton_2.isChecked()
            mol_input = window.mmmm.mol_input
            rw_user_mol = Chem.RWMol(mol_input)
            user_atom_index = atomdialog.spinBox.value()
            try:
                if addatom and removeatom:
                    window.textBrowser_18.append("WARNING: You must pick one atom action, not both.")
                if addatom:
                    new_atom = Chem.Atom(user_atom_index)  # 6 is the atomic number for Carbon (C)
                    rw_user_mol.AddAtom(user_atom_index)
                    window.textBrowser_18.append(f"Atom added at index: {rw_user_mol.GetNumAtoms() - 1}")

                if removeatom:
                    rw_user_mol.RemoveAtom(user_atom_index)
                    window.textBrowser_18.append(f"Atom removed at index: {user_atom_index}")
                else:
                    window.textBrowser_18.append("WARNING: You must pick one atom action, not none")
                self.mol_input = rw_user_mol
                window.mmmm.update_display(rw_user_mol)
            except Exception as e:
                window.textBrowser_18.append(f"Error: {str(e)}")
                
                
    class MemMaker:
        def __init__(self, MainMainWindow):
            self.mm = membrane_generator_master(self)
            self.pp = protein_inserter_master(self)
            self.lipid_files = []
            self.lipid_ratios = []
            self.n_lipids_x = None
            self.n_lipids_y = None
            self.box_z = None
            self.spacing = None
            self.membrane_size = None
            self.membrane_file = 'membrane_maker_output.pdb'
            self.final_box_z = None
            self.protein_choice = None
            self.output_file = 'membrane_maker_output.pdb'
            self.protein_spacing = None
            self.protein_files = []  
            self.rotation_matrices = []
            
        def start_membrane_maker(self):
            self.n_lipids_x = window.x_spin_box.value()
            self.n_lipids_y = window.y_spin_box.value()
            self.box_z = window.z_spin_box.value()
            self.spacing = window.spacing_spin.value()
            self.membrane_size = (self.n_lipids_x, self.n_lipids_y, self.box_z)
            self.membrane_file = 'membrane_maker_output.pdb'
            self.final_box_z = window.final_z_spin.value()
            self.protein_choice = window.radioButton.isChecked()
            self.output_file = 'membrane_maker_output.pdb'
            self.protein_files = []  
            self.rotation_matrices = []
            self.protein_spacing = window.pro_spin_ratio.value()
            window.textBrowser_7.append("Variables for Membrane Making were set! Beginning process..")
            window.textBrowser_7.append("Lipid files:", self.lipid_files)
            window.textBrowser_7.append("Lipid ratios:", self.lipid_ratios)
            window.textBrowser_7.append("Output file:", self.output_file)
            window.textBrowser_7.append("Number of lipids in X:", self.n_lipids_x)
            window.textBrowser_7.append("Number of lipids in Y:", self.n_lipids_y)
            window.textBrowser_7.append("Spacing:", self.spacing)
            window.textBrowser_7.append("Initial box Z dimension:", self.box_z)
            window.textBrowser_7.append("Final box Z dimension:", self.final_box_z)
            window.textBrowser_7.append("ratio", self.lipid_ratios)
            if not self.protein_choice:
                try: 
                    self.mm.generate_bilipid_membrane(
                    self.lipid_files, 
                    self.lipid_ratios, 
                    self.output_file, 
                    self.n_lipids_x, 
                    self.n_lipids_y, 
                    self.spacing, 
                    self.box_z, 
                    self.final_box_z
                )
                    
                    window.textBrowser_7.append("Membrane Maker Has Finished!")
                except:
                    window.textBrowser_7.append("ERROR: The Membrane Maker encountered an error. Please make sure all inputs are correct and try again.")
            else:
                try:
                    self.mm.generate_bilipid_membrane(lipid_files, lipid_ratios, output_file, n_lipids_x, n_lipids_y, spacing, box_z, final_box_z)
                    self.pp.orient_proteins_in_membrane(protein_files, membrane_file, output_file, rotation_matrices, membrane_size, protein_spacing)
                except:
                    window.textBroswer_7.append("ERROR: The Membrane Maker and Protein inserters ran into an error.")
        def load_file_mem_maker(self):
            lipid_ratio_input = window.memratiospin.value()
            if lipid_ratio_input == 0 or lipid_ratio_input >= 1.0:
                window.textBroswer_7.append("WARNING: The ratio must be between 0-1. Your file will not be imported.")
            else:    
                try: 
                    self.lipid_ratios.append(float(lipid_ratio_input))
                        
                    lipid_input, _ = QtWidgets.QFileDialog.getOpenFileName(window, "Open File", "", "All Files (*)")
                    if lipid_input:
                        with open(lipid_input, 'r') as file:
                            content = file.read()
                    if content not in self.lipid_files:
                        self.lipid_files.append(str(lipid_input))
                        window.textBrowser_7.append("Succefully imported path: "+ lipid_input+" with the following ratio: " + str(lipid_ratio_input))
                except:
                    window.textBrowser_7.append("ERROR: the lipid was not succesfully imported.")
            
        def load_file_protein(self):
            if self.protein_choice:
                try:
                    file_path, _ = QtWidgets.QFileDialog.getOpenFileName(window, "Open File", "", "All Files (*)")
                    if file_path:
                        self.protein_files.append(file_path)
                except:
                    window.textBrowser_7.append("ERROR: Your protein file was not loaded. Please make sure the file is valid.")


    def execute_command(self):
        command = window.lineEdit_2.text().strip()
        if command:
            window.textBrowser_13.append(f'> {command}')
            window.lineEdit_2.clear()
            
            try:
                result = subprocess.run(command, shell=True, capture_output=True, text=True, timeout=10)
                output = result.stdout if result.stdout else result.stderr
            except subprocess.TimeoutExpired:
                output = 'Command timed out.'
            except Exception as e:
                output = f'Error: {e}'
            
            self.textBrowser_13.append(output)
            
    def activate_michael(self):
        exit_conditions = (":q","quit","exit")
        PROMPT = self.lineEdit.text().strip()
        window.textBrowser_3.append(f'You: {PROMPT}')
        window.lineEdit.clear()
        if query in exit_conditions:
            window.textBrowser_3.append("Michael: Goodbye! I hope I was able to help you today.")
        if query == '':
            window.textBrowser_3.append("Michael: Your response cannot be blank.")
        else:
            window.textBrowser_3.append(f"Michael: {michaelai.comp(PROMPT, MaxToken=3000, outputs=3)}")
             
class CSInputWindow(QtWidgets.QDialog):
    def __init__(self):
        super().__init__()
        self.setupUi()

    def setupUi(self):
        self.setWindowTitle('Chemistry Spider Input Window')
        self.setGeometry(150, 150, 300, 200)

        self.layout = QVBoxLayout()

        self.inputField = QLineEdit(self)
        self.layout.addWidget(self.inputField)

        self.buttonBox = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, self)
        self.buttonBox.accepted.connect(self.acceptInput)
        self.buttonBox.rejected.connect(self.reject)
        self.layout.addWidget(self.buttonBox)

        self.setLayout(self.layout)

    def acceptInput(self):
        window.saf.update_cskey()
        self.close()
class SMILESInputWindow(QtWidgets.QDialog):
    def __init__(self):
        super().__init__()
        self.setupUi()

    def setupUi(self):
        self.setWindowTitle('SMILES Input Window')
        self.setGeometry(150, 150, 300, 200)

        self.layout = QVBoxLayout()

        self.inputField = QLineEdit(self)
        self.layout.addWidget(self.inputField)

        self.buttonBox = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, self)
        self.buttonBox.accepted.connect(self.acceptInput)
        self.buttonBox.rejected.connect(self.reject)
        self.layout.addWidget(self.buttonBox)

        self.setLayout(self.layout)

    def acceptInput(self):
        smiles_input = self.inputField.text()
        self.mol_input = Chem.MolFromSmiles(smiles_input)
        window.molecule_maker_tab.smiles_read(self)
        self.close()

class BondManager(QtWidgets.QDialog):
    def __init__(self):
        super().__init__()
        self.setupUi()
    def setupUi(self):
        self.setWindowTitle('SMILES Input Window')
        self.setGeometry(150, 150, 300, 200)
        self.gridLayout = QtWidgets.QGridLayout(self)
        self.gridLayout.setObjectName("gridLayout")
        self.label = QtWidgets.QLabel(self)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)
        self.label_2 = QtWidgets.QLabel(self)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 0, 1, 1, 1)
        self.spinBox = QtWidgets.QSpinBox(self)
        self.spinBox.setObjectName("spinBox")
        self.gridLayout.addWidget(self.spinBox, 1, 0, 1, 1)
        self.spinBox_2 = QtWidgets.QSpinBox(self)
        self.spinBox_2.setObjectName("spinBox_2")
        self.gridLayout.addWidget(self.spinBox_2, 1, 1, 1, 1)
        self.label_3 = QtWidgets.QLabel(self)
        self.label_3.setObjectName("label_3")
        self.gridLayout.addWidget(self.label_3, 2, 0, 1, 1)
        self.radioButton = QtWidgets.QRadioButton(self)
        self.radioButton.setObjectName("radioButton")
        self.gridLayout.addWidget(self.radioButton, 2, 1, 2, 1)
        self.comboBox = QtWidgets.QComboBox(self)
        self.comboBox.setObjectName("comboBox")
        self.gridLayout.addWidget(self.comboBox, 3, 0, 1, 1)
        self.radioButton_2 = QtWidgets.QRadioButton(self)
        self.radioButton_2.setObjectName("radioButton_2")
        self.gridLayout.addWidget(self.radioButton_2, 4, 1, 1, 1)
        self.buttonBox = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, self)
        self.buttonBox.setObjectName("buttonBox")
        self.gridLayout.addWidget(self.buttonBox, 5, 0, 1, 1)
        self.buttonBox.accepted.connect(self.acceptInput)
        self.buttonBox.rejected.connect(self.reject)
        self.retranslateUi()


    def retranslateUi(self):
        _translate = QtCore.QCoreApplication.translate
        self.setWindowTitle(_translate("Dialog1", "Bond Manager"))
        self.label.setText(_translate("Dialog1", "Atom #1 Id Number"))
        self.label_2.setText(_translate("Dialog1", "Atom #2 Id Number"))
        self.label_3.setText(_translate("Dialog1", "Bond type"))
        self.radioButton.setText(_translate("Dialog1", "Add Bond"))
        self.radioButton_2.setText(_translate("Dialog1", "Remove Bond"))

    def acceptInput(self):
        window.molecule_maker_tab.manage_bonds(self)
        self.close()
class Viewer(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.setupUi()
    def setupUi(self):
        self.resize(620, 559)
        self.gridLayout = QtWidgets.QGridLayout(self)
        self.gridLayout.setObjectName("gridLayout")
        self.openGLWidget = genericOpenGLWidget()
        self.openGLWidget.setObjectName("openGLWidget")
        self.gridLayout.addWidget(self.openGLWidget, 0, 0, 1, 4)
        self.pushButton = QtWidgets.QPushButton(self)
        self.pushButton.setObjectName("pushButton")
        self.gridLayout.addWidget(self.pushButton, 1, 0, 1, 1)
        self.pushButton_3 = QtWidgets.QPushButton(self)
        self.pushButton_3.setObjectName("pushButton_3")
        self.gridLayout.addWidget(self.pushButton_3, 1, 1, 1, 1)
        self.pushButton_6 = QtWidgets.QPushButton(self)
        self.pushButton_6.setObjectName("pushButton_6")
        self.gridLayout.addWidget(self.pushButton_6, 1, 2, 1, 1)
        self.pushButton_7 = QtWidgets.QPushButton(self)
        self.pushButton_7.setObjectName("pushButton_7")
        self.gridLayout.addWidget(self.pushButton_7, 1, 3, 1, 1)
        self.pushButton_2 = QtWidgets.QPushButton(self)
        self.pushButton_2.setObjectName("pushButton_2")
        self.gridLayout.addWidget(self.pushButton_2, 2, 0, 1, 1)
        self.pushButton_4 = QtWidgets.QPushButton(self)
        self.pushButton_4.setObjectName("pushButton_4")
        self.gridLayout.addWidget(self.pushButton_4, 2, 1, 1, 1)
        self.pushButton_8 = QtWidgets.QPushButton(self)
        self.pushButton_8.setObjectName("pushButton_8")
        self.gridLayout.addWidget(self.pushButton_8, 2, 2, 1, 1)
        self.pushButton_5 = QtWidgets.QPushButton(self)
        self.pushButton_5.setObjectName("pushButton_5")
        self.gridLayout.addWidget(self.pushButton_5, 2, 3, 1, 1)

        self.retranslateUi()
        self.pushButton_2.clicked.connect(window.ana.zoom_in_max)
        self.pushButton.clicked.connect(window.ana.zoom_out_max)
        self.pushButton_5.clicked.connect(window.ana.minimize_viewer)
        self.pushButton_6.clicked.connect(window.ana.auto_rotate)
        self.pushButton_7.clicked.connect(window.ana.expand_viewer)
        self.pushButton_8.clicked.connect(window.ana.change_style)
        self.pushButton_3.clicked.connect(window.ana.load_file)
        self.pushButton_4.clicked.connect(window.ana.save_file)
    def retranslateUi(self):
        _translate = QtCore.QCoreApplication.translate
        self.setWindowTitle(_translate("Form", "Analytival viewer"))
        self.pushButton.setText(_translate("Form", "Max Zoom Out"))
        self.pushButton_3.setText(_translate("Form", "Load File"))
        self.pushButton_6.setText(_translate("Form", "Auto Rotate"))
        self.pushButton_7.setText(_translate("Form", "Full Screen"))
        self.pushButton_2.setText(_translate("Form", "Max Zoom In"))
        self.pushButton_4.setText(_translate("Form", "Save File"))
        self.pushButton_8.setText(_translate("Form", "Change view style"))
        self.pushButton_5.setText(_translate("Form", "Minimize"))
        
class AtomManager(QtWidgets.QDialog):
    def __init__(self):
        super().__init__()
        self.setupUi()
    def setupUi(self):
        self.setWindowTitle('Atom Manager')
        self.setGeometry(150, 150, 300, 200)
        self.gridLayout = QtWidgets.QGridLayout(self)
        self.gridLayout.setObjectName("gridLayout")
        self.label = QtWidgets.QLabel(self)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)
        self.spinBox = QtWidgets.QSpinBox(self)
        self.spinBox.setObjectName("spinBox")
        self.gridLayout.addWidget(self.spinBox, 1, 0, 1, 1)
        self.radioButton = QtWidgets.QRadioButton(self)
        self.radioButton.setObjectName("radioButton")
        self.gridLayout.addWidget(self.radioButton, 1, 1, 1, 1)
        self.label_2 = QtWidgets.QLabel(self)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 2, 0, 1, 1)
        self.spinBox_2 = QtWidgets.QSpinBox(self)
        self.spinBox_2.setObjectName("spinBox_2")
        self.gridLayout.addWidget(self.spinBox_2, 3, 0, 1, 1)
        self.radioButton_2 = QtWidgets.QRadioButton(self)
        self.radioButton_2.setObjectName("radioButton_2")
        self.gridLayout.addWidget(self.radioButton_2, 3, 1, 1, 1)
        self.gridLayout.addWidget(self.radioButton_2, 4, 1, 1, 1)
        self.buttonBox = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, self)
        self.buttonBox.setObjectName("buttonBox")
        self.gridLayout.addWidget(self.buttonBox, 5, 0, 1, 1)
        self.buttonBox.accepted.connect(self.acceptInput)
        self.buttonBox.rejected.connect(self.reject)

        self.retranslateUi()


    def retranslateUi(self):
        _translate = QtCore.QCoreApplication.translate
        self.label.setText(_translate("Dialog", "Atom ID"))
        self.radioButton.setText(_translate("Dialog", "Add"))
        self.label_2.setText(_translate("Dialog", "Atomic Number"))
        self.radioButton_2.setText(_translate("Dialog", "Remove"))
    def acceptInput(self):
        window.molecule_maker_tab.manage_atoms(self)
        self.close()
        
class settings_ui(QtWidgets.QDialog):
    def __init__(self):
        super().__init__()
        self.setupUi()
     def setupUi(self):
        settings_ui.setObjectName("settings_ui")
        settings_ui.resize(320, 638)
        self.verticalLayout = QtWidgets.QVBoxLayout(settings_ui)
        self.verticalLayout.setObjectName("verticalLayout")
        self.input_label = QtWidgets.QLabel(settings_ui)
        self.input_label.setObjectName("input_label")
        self.verticalLayout.addWidget(self.input_label)
        self.input_combo = QtWidgets.QComboBox(settings_ui)
        self.input_combo.setObjectName("input_combo")
        self.verticalLayout.addWidget(self.input_combo)
        self.output_label = QtWidgets.QLabel(settings_ui)
        self.output_label.setObjectName("output_label")
        self.verticalLayout.addWidget(self.output_label)
        self.output_combo = QtWidgets.QComboBox(settings_ui)
        self.output_combo.setObjectName("output_combo")
        self.verticalLayout.addWidget(self.output_combo)
        self.output_label_3 = QtWidgets.QLabel(settings_ui)
        self.output_label_3.setObjectName("output_label_3")
        self.verticalLayout.addWidget(self.output_label_3)
        self.lineEdit = QtWidgets.QLineEdit(settings_ui)
        self.lineEdit.setObjectName("lineEdit")
        self.verticalLayout.addWidget(self.lineEdit)
        self.output_label_4 = QtWidgets.QLabel(settings_ui)
        self.output_label_4.setObjectName("output_label_4")
        self.verticalLayout.addWidget(self.output_label_4)
        self.lineEdit_2 = QtWidgets.QLineEdit(settings_ui)
        self.lineEdit_2.setObjectName("lineEdit_2")
        self.verticalLayout.addWidget(self.lineEdit_2)
        self.output_label_2 = QtWidgets.QLabel(settings_ui)
        self.output_label_2.setObjectName("output_label_2")
        self.verticalLayout.addWidget(self.output_label_2)
        self.checkBox = QtWidgets.QCheckBox(settings_ui)
        self.checkBox.setObjectName("checkBox")
        self.verticalLayout.addWidget(self.checkBox)
        self.checkBox_2 = QtWidgets.QCheckBox(settings_ui)
        self.checkBox_2.setObjectName("checkBox_2")
        self.verticalLayout.addWidget(self.checkBox_2)
        self.checkBox_5 = QtWidgets.QCheckBox(settings_ui)
        self.checkBox_5.setObjectName("checkBox_5")
        self.verticalLayout.addWidget(self.checkBox_5)
        self.checkBox_6 = QtWidgets.QCheckBox(settings_ui)
        self.checkBox_6.setObjectName("checkBox_6")
        self.verticalLayout.addWidget(self.checkBox_6)
        self.checkBox_4 = QtWidgets.QCheckBox(settings_ui)
        self.checkBox_4.setObjectName("checkBox_4")
        self.verticalLayout.addWidget(self.checkBox_4)
        self.checkBox_3 = QtWidgets.QCheckBox(settings_ui)
        self.checkBox_3.setObjectName("checkBox_3")
        self.verticalLayout.addWidget(self.checkBox_3)
        self.checkBox_10 = QtWidgets.QCheckBox(settings_ui)
        self.checkBox_10.setObjectName("checkBox_10")
        self.verticalLayout.addWidget(self.checkBox_10)
        self.checkBox_11 = QtWidgets.QCheckBox(settings_ui)
        self.checkBox_11.setObjectName("checkBox_11")
        self.verticalLayout.addWidget(self.checkBox_11)
        self.checkBox_12 = QtWidgets.QCheckBox(settings_ui)
        self.checkBox_12.setObjectName("checkBox_12")
        self.verticalLayout.addWidget(self.checkBox_12)
        self.checkBox_7 = QtWidgets.QCheckBox(settings_ui)
        self.checkBox_7.setObjectName("checkBox_7")
        self.verticalLayout.addWidget(self.checkBox_7)
        self.checkBox_8 = QtWidgets.QCheckBox(settings_ui)
        self.checkBox_8.setObjectName("checkBox_8")
        self.verticalLayout.addWidget(self.checkBox_8)
        self.checkBox_9 = QtWidgets.QCheckBox(settings_ui)
        self.checkBox_9.setObjectName("checkBox_9")
        self.verticalLayout.addWidget(self.checkBox_9)
        self.buttonBox = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, self)
        self.buttonBox.setObjectName("buttonBox")
        self.gridLayout.addWidget(self.buttonBox, 5, 0, 1, 1)
        self.buttonBox.accepted.connect(self.acceptInput)
        self.buttonBox.rejected.connect(self.reject)
        

        self.retranslateUi(settings_ui)
        
        QtCore.QMetaObject.connectSlotsByName(settings_ui)

    def retranslateUi(self, settings_ui):
        _translate = QtCore.QCoreApplication.translate
        settings_ui.setWindowTitle(_translate("settings_ui", "PyBRANE Settings"))
        self.input_label.setText(_translate("settings_ui", "Default Input File Type"))
        self.output_label.setText(_translate("settings_ui", "Default Output File Type"))
        self.output_label_3.setText(_translate("settings_ui", "ChemSpider API Key"))
        self.output_label_4.setText(_translate("settings_ui", "OpenAI API Key"))
        self.output_label_2.setText(_translate("settings_ui", "Trajectory Graphs"))
        self.checkBox.setText(_translate("settings_ui", "CheckBox"))
        self.checkBox_2.setText(_translate("settings_ui", "CheckBox"))
        self.checkBox_5.setText(_translate("settings_ui", "CheckBox"))
        self.checkBox_6.setText(_translate("settings_ui", "CheckBox"))
        self.checkBox_4.setText(_translate("settings_ui", "CheckBox"))
        self.checkBox_3.setText(_translate("settings_ui", "CheckBox"))
        self.checkBox_10.setText(_translate("settings_ui", "CheckBox"))
        self.checkBox_11.setText(_translate("settings_ui", "CheckBox"))
        self.checkBox_12.setText(_translate("settings_ui", "CheckBox"))
        self.checkBox_7.setText(_translate("settings_ui", "CheckBox"))
        self.checkBox_8.setText(_translate("settings_ui", "CheckBox"))
        self.checkBox_9.setText(_translate("settings_ui", "CheckBox"))

        
if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = MainMainWindow()
    smilesdialog = SMILESInputWindow()
    bonddialog = BondManager()
    atomdialog = AtomManager()
    cswindow = CSInputWindow()
    cswindow.hide()
    viewerwindow = Viewer()
    settings_ui = Settings()
    Settings.hide()
    viewerwindow.hide()
    smilesdialog.hide()
    bonddialog.hide()
    atomdialog.hide()
    window.show()


    sys.exit(app.exec())
