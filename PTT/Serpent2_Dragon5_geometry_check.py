import DMLGInterface as DMLG
import numpy as np
"""
input_f = "MCNP_AT10_sanitized.inp"
AT10_CRTL_MCNP = DMLG_Interface(input_f, type="MCNP",mode="input")
#print(AT10_CRTL_MCNP.Cell_Cards)
AT10_CRTL_MCNP.getMCNP_card_data(print_cells=False, print_surfaces=False, print_materials=True)
#print(AT10_CRTL_MCNP.Cell_Cards)
#print(AT10_CRTL_MCNP.Cell_Cards['water box centered at ( 8.267, 6.973)'].material_densities[2])
"""
cell_serpent_output = "Serpent2/data_AT10_24UOX_try/AT10_24_test_mc.out"
AT10_24UOX_cell_serp = DMLG.DMLG_Interface(cell_serpent_output, type="Serpent2",mode="output")
cell_dragon_output = "../Version5_ev3232/Dragon/Linux_x86_64/SALT_TSPC.result"
AT10_24UOX_cell_drag = DMLG.DMLG_Interface(cell_dragon_output, type="Dragon",mode="output")

material_assocoation_dict = {"UOx_A":1, "UOx_B":2, "UOx_C":3, "UOx_D":4, "gap":5, "clad":6, "H2O":7}
#print(AT10_24UOX_cell_serp.Serpent2_cards)
AT10_24UOX_cell_serp.createSerpent2_geometry("AT10_24UOX_cell",height=3.8E2)
print(f"Total fuel mass is : {AT10_24UOX_cell_serp.Serpent2_geom.getTotalFuelMass():.3f} g")
print(f"Total fuel mass dens is : {AT10_24UOX_cell_serp.Serpent2_geom.getTotalFuelMassDens():.3f} g/cm3")

def check_SerpentvsDragon_vols(Serpent2_case, Dragon5_case, material_assocoation_dict):
    """
    Serpent2_case : DMLG object corresponding to Serpent2 case
    Dragon5_case : DMLG object corresponding to equivalent Dragon5 case
    material_association_dict : dictionary associating Serpent2 materials to region numbers from Dragon5. 
    """
    D5_regions_volumes = Dragon5_case.Dragon5_geom.getVolumesAndRegions()
    S2_regions_volumes = Serpent2_case.Serpent2_geom.getMaterialsandVolumes()
    print(D5_regions_volumes)
    print(S2_regions_volumes)
    errors=[]
    for mat in material_assocoation_dict.keys():
        region=material_assocoation_dict[mat]
        if mat == "H2O":
            errors.append(D5_regions_volumes[f"region {region}"]-8.47125E-01)
        else:
            errors.append(D5_regions_volumes[f"region {region}"]-S2_regions_volumes[mat])

    sum=0
    for err in errors:
        sum+=err**2
    rms=np.sqrt(sum/len(errors))
    print(f"errors on volumes are {errors}, with an RMS error of {rms}")

    return

check_SerpentvsDragon_vols(AT10_24UOX_cell_serp, AT10_24UOX_cell_drag, material_assocoation_dict)

