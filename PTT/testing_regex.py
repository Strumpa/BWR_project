import re

def parse_material(text):
    # Define patterns to extract relevant information
    material_pattern = r'Material "(.*?)"'
    property_patterns = {
        'Atom density': r'Atom density ([\d.E-]+) 1/barn\*cm',
        'Mass density': r'Mass density ([\d.E+]+) g/cm3',
        'Volume': r'Volume ([\d.E-]+) cm3',
        'Mass': r'Mass ([\d.E+]+) g'
    }

    # Initialize variables to store parsed data
    materials = {}
    current_material = None

    # Iterate through each line in the text
    for line in text.split('\n'):
        # Check if the line matches the material pattern
        match_material = re.match(material_pattern, line)
        if match_material:
            material_name = match_material.group(1)
            materials[material_name] = {}
            current_material = materials[material_name]
            continue

        # Check if the line contains a property
        for prop, pattern in property_patterns.items():
            match_property = re.match(pattern, line)
            if match_property and current_material is not None:
                current_material[prop] = match_property.group(1)
                break

    return materials

# Example text to parse
text = """
Material "UOx_A":

 - Material is burnable
 - Material is included in majorant
 - Material is included in geometry
 - Atom density 7.00058E-02 1/barn*cm
 - Mass density 1.04612E+01 g/cm3
 - Volume 3.08964E-01 cm3
 - Mass 3.23212E+00 g
 - Photon emission rate 0.00000E+00 1/s
 - Neutron emission rate 0.00000E+00 1/s
 - 1582 nuclides in composition
 - No nuclides associated with S(a,b) data
 - 2.73 Mb of memory allocated for data
"""

# Parse the text and print the results
materials_info = parse_material(text)
for material, properties in materials_info.items():
    print(f"Material: {material}")
    for prop, value in properties.items():
        print(f"{prop}: {value}")
    print()
