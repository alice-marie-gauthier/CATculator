from itertools import product
from periodictable import elements
from math import isclose

def calculate_mass(formula, monoisotopiqueMasses):
    """
    Calculate the molecular mass for a given formula.
    :param formula: Dictionary of elements and their counts (e.g., {"C": 6, "H": 12, "O": 6, "N":14, "S":32}).
    :param monoisotopiqueMasses: Dictionary of element symbols and their monoisotopique masses.
    :return: Monoisotopic mass of the molecule.
    """
    return sum(monoisotopiqueMasses[elements] * count for elements, count in formula.items())

def calculate_unsaturation(formula):
    """
    Calculate the unsaturation index (double bond equivalents, DBE) for a given formula.
    Formula: DBE = C - (H + X)/2 + N/2 + 1
    Where X = F + Cl + Br + I (halogens)
    :param formula: Dictionary of elements and their counts.
    :return: Unsaturation index (float).
    """
    C = formula.get("C", 0)
    H = formula.get("H", 0)
    N = formula.get("N", 0)
    F = formula.get("F", 0)
    Cl = formula.get("Cl", 0)
    Br = formula.get("Br", 0)
    I = formula.get("I", 0)

    # Total halogens
    X = F + Cl + Br + I

    # Calculate DBE
    return C - (H + X) / 2 + N / 2 + 1

def is_valid_formula(formula):
    """
    Validate a molecular formula based on valency rules.
    Example: Ensure that the number of bonds does not exceed the capacity of the atoms involved.
    :param formula: Dictionary of elements and their counts.
    :return: True if the formula is valid, False otherwise.
    """
    # Define maximum valencies for the elements (C, H, N, O, S)
    valency = {"C": 4, "H": 1, "N": 3, "O": 2, "S": 2, "F": 1, "Cl": 7, "Br": 7, "I": 7}

    # Calculate the total number of bonds in the formula
    totalBonds = (
        4 * formula.get("C", 0) +
        1 * formula.get("H", 0) +
        3 * formula.get("N", 0) +
        2 * formula.get("O", 0) +
        2 * formula.get("S", 0) +
        1 * formula.get("F", 0) +
        7 * formula.get("Cl", 0) +
        7 * formula.get("Br", 0) +
        7 * formula.get("I", 0)
    )

    # The total bonds in a molecule must match its connectivity; adjust as needed.
    maxPossibleBonds = sum(valency[elements] * count for elements, count in formula.items())

    # Check if total bonds do not exceed the maximum possible
    return totalBonds <= maxPossibleBonds

def generate_formulas(atomRanges, monoisotopiqueMasses, minMass, maxMass, maxTotalAtoms=150):
    """
    Generate all molecular formulas within a specified mass range.
    :param atomRanges: Dictionary with element symbols as keys and (min, max) tuples as values.
    :param monoisotopiqueMasses: Dictionary of element symbols and their monoisotopique masses.
    :param minMass: Minimum mass for the formulas.
    :param maxMass: Maximum mass for the formulas.
    :param maxTotalAtoms: Maximum total number of atoms allowed in the formula.
    :return: Generator of valid molecular formulas.
    """
    # Create ranges for all elements
    ranges = [range(atomRanges[element][0], atomRanges[element][1] + 1) for element in atomRanges.keys()]

    for combo in product(*ranges):
        formula = {element: count for element, count in zip(atomRanges.keys(), combo)}
        totalAtoms = sum(formula.values())

        # Skip formulas with zero atoms or exceeding the max atom limit
        if totalAtoms == 0 or totalAtoms > maxTotalAtoms:
            continue

        # Calculate the mass of the formula
        mass = calculate_mass(formula, monoisotopiqueMasses)

        # Check if the mass is within the specified range
        if minMass <= mass <= maxMass:
            yield formula

def find_molecular_formulas(targetMass, ppmTolerance, atomRanges, monoisotopiqueMasses, ionization, unsaturationRange):
    """
    Find possible molecular formulas within a mass tolerance defined by ppm.
    :param targetMass: Target molecular mass (float).
    :param ppmTolerance: Tolerance in parts per million (ppm) (float).
    :param atomRanges: Dictionary with element symbols as keys and (min, max) tuples as values.
    :param monoisotopiqueMasses: Dictionary of element symbols and their monoisotopique masses.
    :param ionization: Mass adjustment for ionization (e.g., +1.007 for protonation).
    :param unsaturationRange: Tuple specifying the range of acceptable unsaturation indices.
    :return: List of possible molecular formulas with their masses, deviations, and unsaturation values.
    """
    # Calculate tolerance range
    tolerance = ppmTolerance / 1e6
    exactMass = targetMass - ionization
    minMass = exactMass - tolerance * exactMass
    maxMass = exactMass + tolerance * exactMass
    print('\nMass range from the expected value:','[',minMass,',',maxMass,']')

    possible_formulas = []

    # Generate formulas and check within the tolerance range
    for formula in generate_formulas(atomRanges, monoisotopiqueMasses, minMass, maxMass):
        mass = calculate_mass(formula, monoisotopiqueMasses)
        unsaturation = calculate_unsaturation(formula)

        if unsaturationRange[0] <= unsaturation <= unsaturationRange[1]:
            deviationppm = ((exactMass - mass) / mass) * 1e6
            mz = mass + ionization
            possible_formulas.append({
                "formula": formula,
                "mass": mass,
                "m/z": mz,
                "deviation": deviationppm,
                "unsaturation": unsaturation,
            })

    # Sort formulas by the smallest deviation
    possible_formulas.sort(key=lambda x: abs(x["deviation"]))
    return possible_formulas