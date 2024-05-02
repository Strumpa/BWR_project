## Python3 script for parsing MCNP input file
# Author : R. Guasch
# Idea : convert MCNP geom to Serpent2
# Try to keep it general : idea could be to use functions in centralized package for geometry handling
# Goals :
# Step (1) : Convert MCNP input file for ATRIUM10 with control blade to Serpent2 format
# Step (2) : define a general Geometry object that could be written in Serpent2 or Dragon5 native geom
# Step (3) : allow for surfacic geometry ?

import numpy as np
import math
from GeometricTools import geom_PIN
import MCNP_cards as MCNP
import Dragon5_geometry as D5
import Serpent2_cards as S2

# Class DMLG_Interface()

class DMLG_Interface:
    def __init__(self, input_deck, type, mode):
        """
        loading geometry object to be parsed according to :
        input file (str) : path to file to be parsed
        type (str) : "MCNP", "Native", "Serpent2"
        mode (str) : 'input' or 'output' or 'check volumes'
        """
        self.type = type
        self.input_deck = input_deck
        self.mode = mode
        if self.type == "MCNP" and self.mode == 'input' :
            self.parseMCNP_input()
            self.createMCNPcard_objects()
        elif self.type == "Dragon" and self.mode == 'input' :
            self.parseDragon_input()
        elif self.type == "Dragon" and self.mode == 'output' :
            self.parseDragon_output()
            self.createDragon5_geometry()
        elif self.type == "Serpent2" and self.mode == 'input' :
            self.parseSerpent2_input()
        elif self.type == "Serpent2" and self.mode == 'check_volumes' :
            self.parseSepent2_check_volumes()
            self.createSerpent2_Material_volumes()
        elif self.type == "Serpent2" and self.mode == 'output' : 
            self.parseSerpent2_output()
            self.createSerpent2_output_cards()
        else:
            print("Error: invalid Geometry type to be parsed")
        
    # MCNP related functions    
    def parseMCNP_input(self):
        """
        parsing the MCNP input file
        Its main structure contains a defintion of "Cell Cards", "Surface Cards" and "Material Cards" 
        Each section has a list associated in order to store the info about each card.
        The lists are :
        for Cell Cards : group of cells with associated data
        for Surface Cards : group of surfaces with asociated data
        for Material Cards : material name with associated data

        Classes defining the proper data structure for each card/ group of cards are called in this function.
        """
        print(f"Parsing {self.type} {self.mode} file, from input file : {self.input_deck}")
        # Define flags to track when to start/stop parsing each section
        parsing_cell_cards = False
        parsing_surface_cards = False
        parsing_material_cards = False

        # Initialize dictionaries to store information from each section
        cell_cards = []
        surface_cards = []
        material_cards = []
        cell_titles = []
        surface_titles = []
        material_titles = []

        # Initialize variables to store current card title and contents
        current_title = None

        # Open the input file
        with open(self.input_deck, 'r') as file:
            for line in file:
                line = line.strip()  # Remove leading/trailing whitespace

                # Check if the line indicates the start of a section
                if "Cell  Cards" in line:
                    parsing_cell_cards = True
                    parsing_surface_cards = False
                    parsing_material_cards = False
                elif "Surface Cards" in line:
                    parsing_cell_cards = False
                    parsing_surface_cards = True
                    parsing_material_cards = False
                elif "Material  Cards" in line:
                    parsing_cell_cards = False
                    parsing_surface_cards = False
                    parsing_material_cards = True
                elif "Tally Cards" in line:
                    parsing_cell_cards = False
                    parsing_surface_cards = False
                    parsing_material_cards = False
                
                elif parsing_cell_cards or parsing_surface_cards:
                    if line.startswith("c"):
                        if line.lstrip("c").strip():
                            current_title=line.lstrip("c").strip()

                # Parse and store information from each section
                elif parsing_cell_cards or parsing_surface_cards and line:
                    # Store previous card's contents (if any) in the corresponding dictionary
                    if parsing_cell_cards:
                        cell_cards.append(line.replace("  ", " ").replace("   "," "))
                        cell_titles.append(current_title)
                    elif parsing_surface_cards:
                        surface_cards.append(line.replace("  ", " ").replace("   "," "))
                        surface_titles.append(current_title)
                elif parsing_material_cards and line:
                    if line.startswith("m"):
                        line=line.replace("  ", " ").replace("   "," ")
                        current_title=line.split(" ")[0]
                        material_titles.append(current_title)
                    if not line.startswith("c"):
                        line=line.replace("  ", " ").replace("   "," ")
                        material_cards.append(line)
        file.close()              

 

        self.Cell_Cards = self.clean_list(cell_cards)
        self.Surface_Cards = self.clean_list(surface_cards)
        self.cell_titles = self.clean_list(cell_titles)
        self.surface_titles = self.clean_list(surface_titles)
        
        self.Material_Cards = self.combine_list(material_cards) 
        
    def clean_list(self,list_):
        cleaned_list=[]
        for entry in list_:
            if entry != '' or entry !="\n":
                cleaned_list.append(entry.strip())
        return cleaned_list
    def combine_list(self, list_):
        combined_list=[]
        current_list=[]
        for sub_list in list_:
            sub_list = sub_list.split(" ")
            if sub_list[0][0]!="m" and current_list:
                current_list.extend(sub_list)
            else:
                if current_list:
                    combined_list.append(current_list)
                current_list=sub_list
        if current_list:
            combined_list.append(current_list)
        return combined_list

    def createMCNPcard_objects(self):
        for i in range(len(self.Cell_Cards)):
            self.Cell_Cards[i] = MCNP.MCNP_Cell_Card(1, self.cell_titles[i], self.Cell_Cards[i])
        for i in range(len(self.Surface_Cards)):
            self.Surface_Cards[i] = MCNP.MCNP_Surface_Card(self.surface_titles[i], self.Surface_Cards[i],printlvl=False)
        for i in range(len(self.Material_Cards)):
            self.Material_Cards[i] = MCNP.MCNP_Material_Card(self.Material_Cards[i], printlvl=True)
    
    def getMCNP_card_data(self, print_cells, print_surfaces, print_materials):
        """
        print_cells/surfaces/materials = boolean to specify print level
        """
        if print_cells:
            print(f"There are a total of {len(self.Cell_Cards)} cell cards")
            for cellcard in self.Cell_Cards:
                print(f"Cell card attrtibutes are : group of cells name = {cellcard.cells_group}, \n user defined cell numbers = {cellcard.cell_number} \n material numbers : {cellcard.material_number} \n and material densities : {cellcard.material_density} ")
                print(f"neutron importance is : {cellcard.neutron_importance}")
        if print_surfaces:
            print(f"There are a total of {len(self.Surface_Cards)} surface cards")
            for surface_card in self.Surface_Cards:
                print(f"Surface card attributes are : group of surfaces name : = {surface_card.surfaces_group}")
                print(f"The corresonding surface equation is : {surface_card.surfaceEquation}")
        if print_materials:
            print(f"There are a total of {len(self.Material_Cards)} Material cards")
        
        return

    #Dragon5 related functions
    def parseDragon_input(self):
        print("Native geom type not suported yet")
        return

    def parseDragon_output(self):
        """
        parse Dragon output file for tracking outputs
        retrieve the analytical volumes calculated by the tracking module
        """
        volume_data=[]
        save_volume_values = False
        module =""
        print(f"Parsing {self.type} {self.mode} file, from input file : {self.input_deck}")
        with open(self.input_deck, 'r') as file:
            lines = file.readlines()
            for line in lines:
            # Check if the line contains volume data
                if "SALT" in self.input_deck:
                    module="SALT"
                    if line.strip().startswith("VOLUME"):
                        # Split the line to extract volume values
                        volume_values = line.split()[1:]  # Exclude the "VOLUME" keyword
                        # Convert volume values to floats and append to volume_data list
                        volume_data.extend([float(value.replace('D', 'E')) for value in volume_values])
                        # Set the flag to True to start saving volume values
                        save_volume_values = True
                    elif line.strip().startswith("--------------------") and save_volume_values:
                        # Stop saving volume values when encountering the "--------------------" delimiter
                        save_volume_values = False
                    elif save_volume_values:
                        # Save volume values from lines below the delimiter
                        volume_values = line.split()
                        volume_data.extend([float(value.replace('D', 'E')) for value in volume_values])


        self.geometric_data = self.clean_vol_data(volume_data)
        self.D5module = module
        return
    
    def createDragon5_geometry(self):
        self.Dragon5_geom = D5.Dragon5_geom(self.mode, self.D5module, self.geometric_data)
        return
    
    def clean_vol_data(self, voldata):
        cleaned_data = []
        for vol in voldata:
            if vol != 0.0 :
                cleaned_data.append(vol)
        return cleaned_data
    
    # Serpent2 related functions
    def parseSerpent2_input(self):
        print("Serpent2 geom type not supported yet")
        return
    def parseSepent2_check_volumes(self):
        """
        parsing Serpent2 .mvol file containing material volumes MC evaluation obtained with -checkvolumes option
        """
        Material_Volumes = {}
        with open(self.input_deck, 'r') as file:
            for line in file:
                line = line.strip() # Remove leading/trailing whitespace
                line = line.replace("         ", " ")
                line = line.replace("        ", " ")
                line = line.replace("    ", " ")
                line = line.replace("   ", " ")
                line = line.replace("  ", " ")
                print(line)  
                if "% (0." in line :
                    #print(line.split(" "))
                    if float(line.split(" ")[2]) != 0.00000E+00:
                        mat_name=line.split(" ")[0]
                        mc_volume=line.split(" ")[2]
                        Material_Volumes[mat_name] = mc_volume
        self.Material_Volumes_data = Material_Volumes
    def createSerpent2_Material_volumes(self):
        self.Material_Volumes = []
        for material_vol in self.Material_Volumes_data.keys():
            self.Material_Volumes.append(S2.S2_Material_Vol(material_vol,self.Material_Volumes_data[material_vol],nbDim=3))
    def parseSerpent2_output(self):
        """
        parsing Serpent2 output file
        """
        lines=[]
        materials=[]
        material_data_={}
        current_material=None
        parsing_new_mat = False
        parsing_same_mat = False
        print(f"Parsing {self.type} {self.mode} file, from input file : {self.input_deck}")
        with open(self.input_deck, 'r') as file:
            for line in file:
                line = line.strip()  # Remove leading/trailing whitespace
                if "Material" in line and ":" in line:
                    current_material=line.split(" ")[1].strip('"').strip(":").strip('"')
                    materials.append(current_material)
                    parsing_new_mat=True
                elif parsing_new_mat:
                    material_data_[current_material] = line
                    parsing_new_mat=False
                    parsing_same_mat=True
                elif parsing_same_mat:
                    if "density" in line :
                        data = material_data_[current_material] + line.split(' ')[3]+" "
                    elif "Volume" in line or "Mass" in line :
                        data = material_data_[current_material] + line.split(' ')[2]+" "
                    elif ".06c" in line or ".05c" in line or ".03c" in line:
                        data = material_data_[current_material] + line.split('  ')[0]+" "+line.split('  ')[4]+" "
                    else:
                        data = material_data_[current_material]
                    material_data_.pop(current_material)
                    material_data_[current_material]=data
                """
                if current_material in material_data_.keys():
                    data=material_data_[current_material]+line
                    material_data_.pop(current_material)
                    material_data_[current_material]=data
                elif line:
                    material_data_[current_material]=line
                """
        file.close()
        #print(materials)
        #print(material_data_)
        self.Serpent2_output = material_data_

        return
    def createSerpent2_output_cards(self):
        """
        calling Serpent2 cards class to create associated serpent2 objects
        """
        self.Serpent2_mat_cards=[]
        for mat in self.Serpent2_output.keys():
            self.Serpent2_mat_cards.append(S2.S2_mat_card(mat, self.Serpent2_output[mat], "output"))
    def createS2_mat_properties(self, name):
        """
        Create material properties object from Serpent2 cards. Useful functions for checking materials compositions and volumes are implemented in the S2_geom class
        """
        self.S2_mat_properties = S2.S2_mat_properties(name, self.Serpent2_mat_cards)
        return
    def createS2_geom(self, name, height=1):
        if height == 0:
            ndim = 2
        else:
            ndim = 3
        self.S2_geom = S2.S2_geom(name, self.Material_Volumes, ndim, height)
        



