# Dragon5 geometry Python3 class
#  
# Instances of this class are initialized in the DMLGInterface class when parsing data of type Dragon
# Author: R. Guasch
# Part of the DMLG code Package.
#
# Contribution to BWR geometry handling for the Version5 environment.
# 
# Purpose : check geometric properties from Dragon5 definitions
# 

class Dragon5_geom:
    """
    Definition of Dragon5 geometry object.
    io_mode = output read or input read
    Dragon_module = module used to define input or output of Dragon5.
    """

    def __init__(self, io_mode, Dragon_module, data):
        self.io_mode = io_mode
        self.module = Dragon_module
        print(data)
        self.volumes = {}
        for i in range(len(data)):
            self.volumes[f"region {i+1}"] = data[i]

        print(self.volumes)
    
    def getVolumesAndRegions(self):
        return self.volumes
    def ComputeOrderedVolumesandRegions(self):
        self.ordered_volumes = dict(sorted(self.volumes.items(), key=lambda item: item[1]))
        print(len(self.ordered_volumes))
    
    def getOrderedVolumesandRegions(self):
        return self.ordered_volumes
