# Serpent2 card class structures in Python3
#  
# Instances of this class are initialized in the DMLG_Interface class when parsing Serpent2 cards
# Author: R. Guasch
# Part of the DMLG code Package.
#
# Contribution to BWR geometry handling for the Version5 environment.
# 
# Purpose : retrieve geometric info from parsed data and create Serpent2 cards object
# 

class S2_mat_card:
    def __init__(self, name, data, io_mode):
        """
        Definition of Serpent2 output geometry object.
        io_mode = output read or input read
        data = dict with key = material name 
        values (Atom density, Mass density, Volume, Mass, + composition pairs (iso, a. dens))
        """
        self.io_mode = io_mode
        #print(data)
        self.mat_name = name
        self.atom_dens = float(data.split(" ")[0])
        self.mass_dens = float(data.split(" ")[1])
        self.volume =  float(data.split(" ")[2])
        self.mass = float(data.split(" ")[3])
        # now retrieve compositions !
        print(f"Processing material name : {self.mat_name} with volume {self.volume}")



class S2_mat_properties:
    def __init__(self, name, S2_materials, nbDim):
        """
        name = name of the geometry, user defined
        S2_materials = list of S2_mat_card objects contained in the geometry
        """
        
        self.Geo_name = name
        self.S2materials = S2_materials
        self.nbDim=nbDim
        if self.nbDim == 2:
            self.units = "cm2"
        elif self.nbDim ==3 :
            self.units = "cm3"

        print(f"The Serpent2 geometry being processed is {self.nbDim}-Dimensional, the 'volume' units associated are {self.units}")
        self.ComputeTotalFuelMass()
        self.ComputeTotalFuelVolume()
        self.ComputeTotalFuelMassDens()
        

    def getMaterialsandVolumes(self):
        material_vol_dict = {}
        for material in self.S2materials:
            material_vol_dict[material.mat_name]=material.volume
        return material_vol_dict
    
    
    def getMaterialsandMass(self):
        material_Mass_dict = []
        for material in self.S2materials:
            material_Mass_dict[material.mat_name]=material.mass
        return material_Mass_dict
    
    def getMaterialsandMassDens(self):
        material_MassDens_dict = []
        for material in self.S2materials:
            material_MassDens_dict[material.mat_name]=material.mass_dens 
        return material_MassDens_dict


    # Compute and get total fuel volume for a given geometry
    def ComputeTotalFuelVolume(self):
        fuel_vol = 0
        for material in self.S2materials:
            if "UOx" in material.mat_name or "fuel" in material.mat_name or "FUEL" in material.mat_name:
                fuel_vol+=material.volume
        self.TotalFuelVol = fuel_vol
        return
    
    def getFuelVolume(self):
        return self.TotalFuelVol
    

    # Compute and get total fuel mass for the given geometry
    def ComputeTotalFuelMass(self):
        fuel_mass = 0
        for material in self.S2materials:
            if "UOx" in material.mat_name or "fuel" in material.mat_name or "FUEL" in material.mat_name:
                fuel_mass+=material.mass
        self.TotalFuelMass=fuel_mass
        return
    
    def getTotalFuelMass(self):
        return self.TotalFuelMass

    # Compute and get total fuel mass density for the given geometry
    def ComputeTotalFuelMassDens(self):
        self.TotalFuelMassDens = self.TotalFuelMass/self.TotalFuelVol
        return 
    
    def getTotalFuelMassDens(self):
        return self.TotalFuelMassDens
    
class S2_Material_Vol:
    def __init__(self, name, data, nbDim):
        self.mode = "MVol"
        self.mat_vol_name = name
        self.nbDim=nbDim
        if self.nbDim == 2:
            self.units = "cm2"
        elif self.nbDim ==3 :
            self.units = "cm3"
        self.volume = data
        print(f"{self.mat_vol_name} with volume {self.volume} {self.units} processed")
    
class S2_geom:
    def __init__(self, name, data, nbDim, height):
        """
        Minimal implementation of the S2 geometry :
        at this stage only material volumes are taken to be included in data
        """
        self.geo_name = name
        self.nbDim = nbDim
        self.height = height
        material_volumes = {}
        #print(data)
        #for mat in data.keys():
        for elem in data:
            if elem.mode == "MVol":
                material_volumes[elem.mat_vol_name] = float(elem.volume)
        self.material_volumes = dict(sorted(material_volumes.items(), key=lambda item: item[1]))
        #print(len(self.material_volumes))
    
    def getOrderedMaterialVols(self):
        return self.material_volumes


        