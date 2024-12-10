from itertools import product
from periodictable import elements
from math import isclose
import models.base as mod
from mendeleev import element
from tabulate import tabulate

############################################################################################################################
############################################## PARAMETERS TO BE MODIFIED ###################################################
############################################################################################################################

# Define ranges for each element (min, max)
atomRanges = {
    "C": (0, 50),  # Carbon atoms can range from 0 to 50
    "H": (0, 100),  # Hydrogen atoms can range from 0 to 100
    "N": (0, 20),   # Nitrogen atoms can range from 0 to 20
    "O": (0, 30),   # Oxygen atoms can range from 0 to 30
    "S": (0, 5),   # Sulfur atoms can range from 0 to 5
    "F": (0, 0),  # Fluorine atoms can range from 0 to
    "Br": (0, 0), # Bromine atoms can range from 0 to
    "Cl": (0, 0), # Chlorine atoms can range from 0 to
    "I": (0, 0)   # Iodine atoms can range from 0 to
}

# Range of acceptable unsaturation indices
unsaturationRange = (0, 100) 

############################################## PARAMETERS TO START ##########################################################

# Extract atomic masses dynamically from the `periodictable` package
# atomicMasses = {elem: getattr(elements, elem).mass for elem in ["C", "H", "N", "O", "S"]}

# Monoisotopic masses of each elements from this website: https://www.unimod.org/masses.html
monoisotopicMasses = {
    "H": 1.007825035,   # Hydrogen-1
    "C": 12.0000000,  # Carbon-12
    "N": 14.003074,  # Nitrogen-14
    "O": 15.99491463,  # Oxygen-16
    "S": 31.9720707,  # Sulfur-32
    "F": 18.99840322, # Fluorine
    "Br": 78.9183361, # Bromine
    "Cl": 34.96885272, # Chlorine
    "I": 126.904473 # Iodine
}

# Ionization cases with corresponding mass adjustments
ionizationCases = {
    "+": 0,                  # Neutral
    "H+": 1.007276,          # Single protonation
    "(H+)*2": 2 * 1.007276,  # Double protonation
    "(H+)*3": 3 * 1.007276,  # Triple protonation
    "-": -1.007276,          # Single deprotonation
    "(-)*2": 2*-1.007276,    # Double deprotonation
}

############################################################################################################################
########################################## RESULTS DISPLAY + COMMAND #######################################################
############################################################################################################################

print("\nFind the Molecular Formula of your molecule thanks to the m/z mass!\n")

# Ask user for the target mass and accuracy
targetMass = float(input("Enter the target molecular mass (with a point) [g/mol]: "))
accuracyppm = float(input("\nEnter the accuracy of the measure [ppm]: "))

# Display ionization options
print("\nSelect an ionization state:")
for idx, (name, value) in enumerate(ionizationCases.items(), 1):
    print(f"{idx}. {name} (adjusts mass by {value:+.6f} Da)")

# User selects ionization
ionizationChoice = int(input("Enter the number corresponding to your choice: "))
ionizationName = list(ionizationCases.keys())[ionizationChoice - 1]
ionizationValue = ionizationCases[ionizationName]

# Find molecular formulas
formulas = mod.find_molecular_formulas(
    targetMass, accuracyppm, atomRanges, monoisotopicMasses, ionizationValue, unsaturationRange
)

# Sort formulas by deviation
formulas.sort(key=lambda x: abs(x["deviation"]))

# Display results
if formulas:
    print(f"\nResults for ionization: {ionizationName}")
    print(f"Possible molecular formulas for mass {targetMass} (precision: {accuracyppm} ppm):\n")
    
   # Prepare data for the table
    table_data = []
    for idx, result in enumerate(formulas, 1):
        formula_str = ''.join(f"{element}{count}" for element, count in result['formula'].items() if count > 0)
        table_data.append([
            idx,
            formula_str,
            f"{result['mass']:.4f}",
            f"{result['m/z']:.4f}",
            f"{result['deviation']:.3f}",
            f"{result['unsaturation']:.1f}"
        ])

    # Define table headers
    headers = ["#", "Molecular Formula", "Monoisotopic Mass [g/mol]", "m/z [g/mol]", "Deviation [ppm]", "Unsaturation [-]"]

    # Generate the table as plain text
    table = tabulate(table_data, headers=headers, tablefmt="plain", disable_numparse=True).splitlines()

    # Add a line of dashes matching column widths under the header
    header_line = "-" * len(table[0]) 
    table.insert(1, header_line)

    # Print the table
    print("\n".join(table))
    
else:
    print("\nNo formulas found within the specified tolerance.")