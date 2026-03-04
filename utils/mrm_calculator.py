"""
MRM (Multiple Reaction Monitoring) Mass Spectrometry Calculator.
Calculates precursor ions, common adducts, and predicted fragment ions.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Fragments
import re


# Common adduct masses
ADDUCTS_POSITIVE = {
    "[M+H]+": 1.007276,
    "[M+Na]+": 22.989218,
    "[M+K]+": 38.963158,
    "[M+NH4]+": 18.034164,
    "[M+2H]2+": 1.007276,  # divide by 2
    "[M+Li]+": 7.016003,
    "[M+CH3OH+H]+": 33.033491,
}

ADDUCTS_NEGATIVE = {
    "[M-H]-": -1.007276,
    "[M+Cl]-": 34.969402,
    "[M+HCOO]-": 44.998201,
    "[M+CH3COO]-": 59.013851,
    "[M-2H]2-": -1.007276,  # divide by 2
    "[M+Na-2H]-": 20.974666,
}

# Common neutral losses and their masses
NEUTRAL_LOSSES = {
    "H2O": 18.010565,
    "CO": 27.994915,
    "CO2": 43.989829,
    "NH3": 17.026549,
    "HCN": 27.010899,
    "CH3OH": 32.026215,
    "HCOOH (Formic acid)": 46.005480,
    "CH3COOH (Acetic acid)": 60.021130,
    "H3PO4 (Phosphoric acid)": 97.976895,
    "SO3": 79.956815,
    "HF": 20.006229,
    "HCl": 35.976678,
    "HBr": 79.926160,
    "C6H10O5 (Glucose)": 162.052824,
    "C6H8O6 (Glucuronic acid)": 176.032088,
}


def calculate_mrm_data(mol):
    """Calculate comprehensive MRM mass spectrometry data."""
    if mol is None:
        return None

    data = {}

    # Exact monoisotopic mass
    exact_mass = Descriptors.ExactMolWt(mol)
    data["Monoisotopic Mass"] = round(exact_mass, 6)
    data["Average Mass"] = round(Descriptors.MolWt(mol), 4)
    data["Molecular Formula"] = rdMolDescriptors.CalcMolFormula(mol)

    # --- Positive Mode Adducts ---
    pos_adducts = {}
    for adduct, mass_shift in ADDUCTS_POSITIVE.items():
        if "2+" in adduct:
            mz = round((exact_mass + 2 * mass_shift) / 2, 4)
        else:
            mz = round(exact_mass + mass_shift, 4)
        pos_adducts[adduct] = mz
    data["Positive Mode Adducts"] = pos_adducts

    # --- Negative Mode Adducts ---
    neg_adducts = {}
    for adduct, mass_shift in ADDUCTS_NEGATIVE.items():
        if "2-" in adduct:
            mz = round((exact_mass + 2 * mass_shift) / 2, 4)
        else:
            mz = round(exact_mass + mass_shift, 4)
        neg_adducts[adduct] = mz
    data["Negative Mode Adducts"] = neg_adducts

    # --- Predicted Fragment Ions (Neutral Losses from [M+H]+) ---
    precursor_mh = exact_mass + 1.007276
    fragments_pos = {}
    fragments_neg = {}

    # Check which neutral losses are relevant based on molecular structure
    relevant_losses = get_relevant_neutral_losses(mol)

    for loss_name, loss_mass in relevant_losses.items():
        frag_mz_pos = round(precursor_mh - loss_mass, 4)
        if frag_mz_pos > 50:  # Fragment must be reasonable
            fragments_pos[f"[M+H-{loss_name}]+"] = frag_mz_pos

    precursor_mh_neg = exact_mass - 1.007276
    for loss_name, loss_mass in relevant_losses.items():
        frag_mz_neg = round(precursor_mh_neg - loss_mass, 4)
        if frag_mz_neg > 50:
            fragments_neg[f"[M-H-{loss_name}]-"] = frag_mz_neg

    data["Predicted Fragments (Positive)"] = fragments_pos
    data["Predicted Fragments (Negative)"] = fragments_neg

    # --- Suggested MRM Transitions ---
    transitions = suggest_mrm_transitions(mol, data)
    data["Suggested MRM Transitions"] = transitions

    # --- Isotope Pattern Info ---
    data["Isotope Info"] = get_isotope_info(mol)

    return data


def get_relevant_neutral_losses(mol):
    """Determine which neutral losses are relevant based on molecular structure."""
    if mol is None:
        return {}

    relevant = {}
    smiles = Chem.MolToSmiles(mol)
    formula = rdMolDescriptors.CalcMolFormula(mol)

    # Always include H2O if oxygen present
    if "O" in formula:
        relevant["H2O"] = NEUTRAL_LOSSES["H2O"]

    # CO loss - if carbonyl groups present
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3]=O")):
        relevant["CO"] = NEUTRAL_LOSSES["CO"]

    # CO2 loss - if carboxyl groups
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=O)[OX2H1]")) or \
       mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=O)[OX2]")):
        relevant["CO2"] = NEUTRAL_LOSSES["CO2"]

    # NH3 loss - if amine groups
    if "N" in formula:
        relevant["NH3"] = NEUTRAL_LOSSES["NH3"]

    # HCN loss - if nitrogen in aromatic ring
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[nR]")):
        relevant["HCN"] = NEUTRAL_LOSSES["HCN"]

    # CH3OH loss - if methyl ester
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=O)[OX2][CH3]")):
        relevant["CH3OH"] = NEUTRAL_LOSSES["CH3OH"]

    # Formic acid loss
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=O)[OX2H1]")):
        relevant["HCOOH (Formic acid)"] = NEUTRAL_LOSSES["HCOOH (Formic acid)"]

    # Phosphoric acid loss
    if "P" in formula:
        relevant["H3PO4 (Phosphoric acid)"] = NEUTRAL_LOSSES["H3PO4 (Phosphoric acid)"]

    # SO3 loss - if sulfo groups
    if "S" in formula and mol.HasSubstructMatch(Chem.MolFromSmarts("[SX4](=O)(=O)")):
        relevant["SO3"] = NEUTRAL_LOSSES["SO3"]

    # Halogen losses
    if "F" in formula:
        relevant["HF"] = NEUTRAL_LOSSES["HF"]
    if "Cl" in formula:
        relevant["HCl"] = NEUTRAL_LOSSES["HCl"]
    if "Br" in formula:
        relevant["HBr"] = NEUTRAL_LOSSES["HBr"]

    # Glucose loss - for glycosides
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[OX2][CX4][CX4][CX4][CX4][CX4][OX2]")):
        relevant["C6H10O5 (Glucose)"] = NEUTRAL_LOSSES["C6H10O5 (Glucose)"]

    # If nothing specific found, add common losses
    if not relevant:
        relevant["H2O"] = NEUTRAL_LOSSES["H2O"]
        relevant["CO"] = NEUTRAL_LOSSES["CO"]

    return relevant


def suggest_mrm_transitions(mol, mrm_data):
    """Suggest top MRM transitions for method development."""
    transitions = []
    exact_mass = mrm_data["Monoisotopic Mass"]

    # Primary transition: [M+H]+ -> strongest fragment
    precursor_pos = round(exact_mass + 1.007276, 4)
    precursor_neg = round(exact_mass - 1.007276, 4)

    pos_fragments = mrm_data.get("Predicted Fragments (Positive)", {})
    neg_fragments = mrm_data.get("Predicted Fragments (Negative)", {})

    # Positive mode transitions
    for i, (frag_name, frag_mz) in enumerate(pos_fragments.items()):
        if i >= 3:  # Top 3 transitions
            break
        ce = estimate_collision_energy(precursor_pos, frag_mz, "positive")
        transitions.append({
            "Mode": "Positive",
            "Precursor (m/z)": precursor_pos,
            "Product (m/z)": frag_mz,
            "Transition": f"{precursor_pos} -> {frag_mz}",
            "Fragment": frag_name,
            "Est. CE (eV)": ce,
            "Use": "Quantifier" if i == 0 else "Qualifier"
        })

    # Negative mode transitions
    for i, (frag_name, frag_mz) in enumerate(neg_fragments.items()):
        if i >= 2:  # Top 2 transitions
            break
        ce = estimate_collision_energy(precursor_neg, frag_mz, "negative")
        transitions.append({
            "Mode": "Negative",
            "Precursor (m/z)": precursor_neg,
            "Product (m/z)": frag_mz,
            "Transition": f"{precursor_neg} -> {frag_mz}",
            "Fragment": frag_name,
            "Est. CE (eV)": ce,
            "Use": "Quantifier" if i == 0 else "Qualifier"
        })

    return transitions


def estimate_collision_energy(precursor_mz, product_mz, mode="positive"):
    """
    Estimate collision energy based on precursor mass.
    These are rough estimates - actual CE should be optimized experimentally.
    Based on general linear CE ramp: CE = slope * precursor_mz + intercept
    """
    if mode == "positive":
        # General positive mode estimation
        ce = round(0.03 * precursor_mz + 5, 0)
    else:
        # General negative mode estimation
        ce = round(0.025 * precursor_mz + 5, 0)

    return max(5, min(ce, 80))  # Clamp between 5 and 80 eV


def get_isotope_info(mol):
    """Get basic isotope pattern information."""
    formula = rdMolDescriptors.CalcMolFormula(mol)
    info = {}

    # Count elements for isotope considerations
    atom_counts = {}
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        atom_counts[symbol] = atom_counts.get(symbol, 0) + 1

    info["Contains Cl/Br"] = "Cl" in atom_counts or "Br" in atom_counts
    if "Cl" in atom_counts:
        info["Chlorine Atoms"] = atom_counts["Cl"]
        info["Isotope Pattern"] = "M, M+2 pattern (3:1 ratio per Cl)"
    elif "Br" in atom_counts:
        info["Bromine Atoms"] = atom_counts["Br"]
        info["Isotope Pattern"] = "M, M+2 pattern (1:1 ratio per Br)"
    else:
        info["Isotope Pattern"] = "Standard pattern (M, M+1, M+2)"

    if "S" in atom_counts:
        info["Sulfur Atoms"] = atom_counts["S"]
        info["Sulfur Note"] = "S-34 contributes ~4.4% to M+2 per S atom"

    return info
