## Setup a multiphysics case calling the starterDD library.

import os
import numpy as np
from starterDD.DDModel.DonjonModel import CoreModel
from starterDD.DDModel.DonjonCalculationScheme import DonjonCalculationScheme, ModelInitialisation, NeutronicsSolve
from starterDD.GeometryAnalysis.cartesian_geometry_analysis import CartesianGeometricAnalyser
from starterDD.InterfaceToDD.donjon_module_calls import DonjonTHM1DProcedure, INIT, NEUTRONICS

def get_solver_data_and_core_model(geometry_yaml, case_name, number_of_axial_meshes):
    GE14_DATA = "input_yaml_configs"
    yaml_path = os.path.join(GE14_DATA , geometry_yaml)
        
    if not os.path.exists(yaml_path):
        raise FileNotFoundError(f"Input geometry YAML file {yaml_path} not found.")
    
    # 1. Split folder and file name
    path_to_configs = os.path.dirname(yaml_path)
    core_desc_file = os.path.basename(yaml_path)
    
    # 2. Initiate core model
    core_model = CoreModel(
        name=case_name, 
        path_to_yaml_configs=path_to_configs, 
        core_description_yaml=core_desc_file
    )
    
    # 3. Créer les modèles d'assemblage en mémoire
    core_model.createAssemblyModels()
    
    # 4. On passe l'objet `core_model` à l'analyseur
    analyser = CartesianGeometricAnalyser(core_model=core_model, core_i=1, core_j=1)
    # -------------------------

    solver_data = analyser.run_THM_analysis(nz=number_of_axial_meshes, include_water_rods=False) 
    z_boundaries = np.linspace(0, solver_data["fuel_data"]["max_rod_length"], number_of_axial_meshes + 1)
    z_values = (z_boundaries[:-1] + z_boundaries[1:]) / 2  # Midpoints of control volumes
    solver_data["z_mesh"] = z_values
    return solver_data, core_model

def initialise_starterDD_case(geometry_yaml, case_name, number_of_axial_meshes, reactor_power):
    """ 
    Replace what would be yaml stored data in future starterDD implementations : 
    specify the calculation options for a mutliphysics case
    Inputs : 
    core_model : Donjon CoreModel created by starterDD, holding the 2D / 3D information necessary to initialise the model
    number_of_axial_meshes (int) : number of axial divisions to be used to mesh the geometry and define unique materials on each plane
    reactor_power (float) : total thermal power to be released in the multi-physics case, in W.
    """

    LCM_object_names = {
                        "Geometry" : "Geom",
                        "Fuel_Map" : "Fmap",
                        "Material_Indexation": "Matex",
                        "Tracking" : "Track",
                        "Thermal_hydraulics_data": "THData",
                        }

    data_for_pyTHM, core_model = get_solver_data_and_core_model(geometry_yaml, case_name, number_of_axial_meshes)

    if case_name == "GE14_bundle_DOM":
        DonjonModelInitialiser = ModelInitialisation(core_model=core_model, 
                                            radial_homogenization_strategy="by_assembly",
                                            boundary_conditions_dict={"x": "reflective",
                                                                    "y": "reflective",
                                                                    "z": "void"},
                                            number_of_axial_materials_per_slice=[number_of_axial_meshes],
                                            number_of_energy_groups=2,
                                            interpolation_variables=[("TFuel", "local", "fuel_temperature"), 
                                                                    ("TCool", "local", "coolant_temperature"),
                                                                    ("DCool", "local", "coolant_density")])
        
        DonjonModelInitialiser.resolve_mesh()
        DonjonModelInitialiser.set_slice_to_compodir_correspondance(xpos=0, ypos=0, z_bounds=[0.0, 347.1291], compo_name="CpoDOM", dir_name="EDI2G_DOM")
        LCM_object_names["MultiCompos"] = ["CompoDOM"]

    elif case_name == "GE14_bundle_DOM_VAN":
        DonjonModelInitialiser = ModelInitialisation(core_model=core_model, 
                                                    radial_homogenization_strategy="by_assembly",
                                                    boundary_conditions_dict={"x": "reflective",
                                                                            "y": "reflective",
                                                                            "z": "void"},
                                                    number_of_axial_materials_per_slice=[int(number_of_axial_meshes*250.0/400.0), int(number_of_axial_meshes*150.0/400.0)],
                                                    number_of_energy_groups=2,
                                                    interpolation_variables=[("TFuel", "local", "fuel_temperature"), 
                                                                            ("TCool", "local", "coolant_temperature"),
                                                                            ("DCool", "local", "coolant_density")])
            
        DonjonModelInitialiser.resolve_mesh()
        DonjonModelInitialiser.set_slice_to_compodir_correspondance(xpos=0, ypos=0, z_bounds=[0.0, 250.0], compo_name="CpoDOM", dir_name="EDI2G_DOM")
        DonjonModelInitialiser.set_slice_to_compodir_correspondance(xpos=0, ypos=0, z_bounds=[250.0, 400.0], compo_name="CpoVAN", dir_name="EDI2G_VAN")
        LCM_object_names["MultiCompos"] = ["CompoDOM", "CompoVAN"]

    DonjonModelInitialiser.set_reactor_power(reactor_power)
    DonjonModelInitialiser.set_initial_axial_power_form("sine")
    DonjonModelInitialiser.set_initial_parameters({"TFuel": 900, "TCool": 600.0, "DCool": 0.73669})
    DonjonModelInitialiser.set_fuel_mass(total_fuel_mass = 0.5) # "BWR Progression Problems", "for the GE-14 bundle with 500g of heavy metal ..."


    return data_for_pyTHM, LCM_object_names, DonjonModelInitialiser

def generate_c2m_procedures(donjon_model_initialiser, procedure_names):
    """
    Generate cle2000 procedures with starterDD
    """

    calc_scheme = DonjonCalculationScheme(name="GE14_multiphysics_scheme")
    calc_scheme.add_initialisation_step(donjon_model_initialiser)
    neutronics_step = NeutronicsSolve(donjon_model_initialiser, "diffusion", "linear", ["MCFD", 1])
    nz = len(donjon_model_initialiser.meshz) - 1
    neutronics_step.set_interpolation_parameter("TFuel", [900]*nz, ("THData", "TFuelList"))
    neutronics_step.set_interpolation_parameter("TCool", [600]*nz, ("THData", "TCoolList"))
    neutronics_step.set_interpolation_parameter("DCool", [0.73669]*nz, ("THData", "DCoolList"))
    neutronics_step.resolve_cpo_dir_to_mix()
    calc_scheme.add_neutronics_step(neutronics_step)

    initalisation_step = INIT(calc_scheme, proc_name = procedure_names["initialisation"])
    initalisation_step.write_to_c2m(".", procedure_names["initialisation"])

    solve_neutronics_step = NEUTRONICS(calc_scheme, proc_name = procedure_names["neutronics_solve"])
    solve_neutronics_step.write_to_c2m(".",  procedure_names["neutronics_solve"])