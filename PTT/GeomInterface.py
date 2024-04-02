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

# Class GeometryInterface()

class GeometryInterface:
    def __init__(self, input_deck, type):
        """
        loading geometry object to be parsed according to :
        input file (str) : path to file to be parsed
        type (str) : "MCNP", "Native", "Serpent2"
        """
        self.type = type
        self.input_deck = input_deck
        if type == "MCNP":
            self.parseMCNP_deck()
            self.cleanupMCNPcards()
        elif type == "Native_Dragon":
            self.parseNative_Geom()
        elif type == "Serpent2":
            self.parseSerpent2_Geom()
        else:
            print("Error: invalid Geometry type to be parsed")
        
        
    def parseMCNP_deck(self):
        """
        parsing the MCNP input file (deck)
        Its main structure contains a defintion of "Cell Cards", "Surface Cards" and "Material Cards" 
        Each section has a dictionnary associated in order to store the info about each card.
        The dictionnary's values are :
        for Cell Cards : 
        for Surface Cards :
        for Material Cards : the material number/identifier
        """
        print(f"Parsing geometry of type : {self.type}, from input file : {self.input_deck}")
        # Define flags to track when to start/stop parsing each section
        parsing_cell_cards = False
        parsing_surface_cards = False
        parsing_material_cards = False

        # Initialize dictionaries to store information from each section
        cell_cards = {}
        surface_cards = {}
        material_cards = {}

        # Initialize variables to store current card title and contents
        current_title = None
        current_contents = []

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

                # Parse and store information from each section
                elif parsing_cell_cards or parsing_surface_cards or parsing_material_cards:
                    if line:  # Skip empty lines
                        # Check if the line indicates a new card title
                        if line.startswith("c"):
                            # Store previous card's contents (if any) in the corresponding dictionary
                            if current_title and current_contents:
                                if parsing_cell_cards:
                                    cell_cards[current_title] = current_contents.copy()
                                elif parsing_surface_cards:
                                    surface_cards[current_title] = current_contents.copy()
                                elif parsing_material_cards:
                                    material_cards[current_title] = current_contents.copy()

                            # Update current card title and reset contents list
                            current_title = line.split("    ")[-1]
                            current_contents = []
                        else:
                            # Add line to current card's contents
                            current_contents.append(line.replace("  ", " ").replace("   ", " ").replace("    ", " ").replace("     ", " ").replace("      "," "))

            # Store contents of the last card in each section
            if current_title and current_contents:
                if parsing_cell_cards:
                    cell_cards[current_title] = current_contents.copy()
                elif parsing_surface_cards:
                    surface_cards[current_title] = current_contents.copy()
                elif parsing_material_cards:
                    material_cards[current_title] = current_contents.copy()

        self.Cell_Cards = cell_cards
        self.Surface_Cards = surface_cards
        self.Material_Cards = material_cards

    def cleanupMCNPcards(self):
        for cell in self.Cell_Cards.keys():
            card=""
            for line in self.Cell_Cards[cell]:
                card+=line+" "
            self.Cell_Cards[cell] = card







    def parseNative_Geom(self):
        print("Native geom type not suported yet")
        return

    def parseSerpent2_Geom(self):
        print("Serpent2 geom type not supported yet")
        return


input_f = "MCNP_AT10_sanitized.inp"
AT10_CRTL_MCNP = GeometryInterface(input_f, type="MCNP")
print(AT10_CRTL_MCNP.Cell_Cards)