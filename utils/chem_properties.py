# -*- coding: utf-8 -*-
"""
RDKit-based physicochemical property calculations.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen, rdMolDescriptors, QED
from rdkit.Chem import Draw, AllChem
from io import BytesIO
from PIL import Image
import base64


def get_mol_from_smiles(smiles: str):
    """Convert SMILES to RDKit Mol object."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        AllChem.Compute2DCoords(mol)
    return mol


def generate_structure_image(mol, size=(450, 350)):
    """Generate 2D structure image from RDKit Mol object."""
    if mol is None:
        return None
    img = Draw.MolToImage(mol, size=size)
    return img


def generate_structure_image_base64(mol, size=(450, 350)):
    """Generate 2D structure image as base64 string."""
    if mol is None:
        return None
    img = Draw.MolToImage(mol, size=size)
    buffered = BytesIO()
    img.save(buffered, format="PNG")
    return base64.b64encode(buffered.getvalue()).decode()


def calculate_physicochemical_properties(mol):
    """Calculate comprehensive physicochemical properties."""
    if mol is None:
        return None

    props = {}

    # --- Basic Molecular Properties ---
    props["Molecular Formula"] = rdMolDescriptors.CalcMolFormula(mol)
    props["Molecular Weight (g/mol)"] = round(Descriptors.MolWt(mol), 4)
    props["Exact Mass (Monoisotopic)"] = round(Descriptors.ExactMolWt(mol), 6)
    props["Heavy Atom Count"] = Lipinski.HeavyAtomCount(mol)
    props["Number of Atoms"] = mol.GetNumAtoms()
    props["Number of Bonds"] = mol.GetNumBonds()

    # --- Lipophilicity & Solubility ---
    props["ALogP (Wildman-Crippen)"] = round(Crippen.MolLogP(mol), 4)
    props["Molar Refractivity (MR)"] = round(Crippen.MolMR(mol), 4)

    # --- Polar Surface Area ---
    props["Topological PSA (A2)"] = round(Descriptors.TPSA(mol), 4)

    # --- Hydrogen Bond Donors & Acceptors ---
    props["H-Bond Donors (HBD)"] = Lipinski.NumHDonors(mol)
    props["H-Bond Acceptors (HBA)"] = Lipinski.NumHAcceptors(mol)

    # --- Rotatable Bonds ---
    props["Rotatable Bonds"] = Lipinski.NumRotatableBonds(mol)

    # --- Ring Information ---
    props["Number of Rings"] = Descriptors.RingCount(mol)
    props["Aromatic Rings"] = Descriptors.NumAromaticRings(mol)
    props["Aliphatic Rings"] = Descriptors.NumAliphaticRings(mol)
    props["Saturated Rings"] = Descriptors.NumSaturatedRings(mol)

    # --- Heteroatoms ---
    props["Number of Heteroatoms"] = Lipinski.NumHeteroatoms(mol)
    props["N+O Count"] = Lipinski.NOCount(mol)
    props["NH+OH Count"] = Lipinski.NHOHCount(mol)

    # --- Fraction & Complexity ---
    props["Fraction Csp3"] = round(Lipinski.FractionCSP3(mol), 4)
    props["Bertz CT (Complexity)"] = round(Descriptors.BertzCT(mol), 2)

    # --- QED (Drug-Likeness Score) ---
    try:
        props["QED (Drug-Likeness)"] = round(QED.qed(mol), 4)
    except Exception:
        props["QED (Drug-Likeness)"] = "N/A"

    # --- Formal Charge ---
    props["Formal Charge"] = Chem.GetFormalCharge(mol)

    # --- InChI & InChIKey ---
    try:
        from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey
        props["InChI"] = MolToInchi(mol)
        props["InChIKey"] = InchiToInchiKey(props["InChI"])
    except Exception:
        try:
            props["InChI"] = Chem.MolToInchi(mol)
            props["InChIKey"] = Chem.InchiToInchiKey(props["InChI"])
        except Exception:
            props["InChI"] = "N/A"
            props["InChIKey"] = "N/A"

    return props


def evaluate_drug_likeness(props):
    """Evaluate drug-likeness rules."""
    if props is None:
        return None

    rules = {}

    # Lipinski's Rule of Five
    mw = props.get("Molecular Weight (g/mol)", 0)
    logp = props.get("ALogP (Wildman-Crippen)", 0)
    hbd = props.get("H-Bond Donors (HBD)", 0)
    hba = props.get("H-Bond Acceptors (HBA)", 0)
    ro5_violations = sum([
        mw > 500,
        logp > 5,
        hbd > 5,
        hba > 10
    ])
    rules["Lipinski Rule of 5"] = {
        "MW <= 500": "PASS" if mw <= 500 else "FAIL",
        "LogP <= 5": "PASS" if logp <= 5 else "FAIL",
        "HBD <= 5": "PASS" if hbd <= 5 else "FAIL",
        "HBA <= 10": "PASS" if hba <= 10 else "FAIL",
        "Violations": ro5_violations,
        "Pass": ro5_violations <= 1
    }

    # Veber's Rules (Oral Bioavailability)
    psa = props.get("Topological PSA (A2)", 0)
    rtb = props.get("Rotatable Bonds", 0)
    rules["Veber Rules"] = {
        "PSA <= 140 A2": "PASS" if psa <= 140 else "FAIL",
        "Rotatable Bonds <= 10": "PASS" if rtb <= 10 else "FAIL",
        "Pass": psa <= 140 and rtb <= 10
    }

    # Ghose Filter
    ghose_pass = (160 <= mw <= 480) and (-0.4 <= logp <= 5.6)
    rules["Ghose Filter"] = {
        "160 <= MW <= 480": "PASS" if 160 <= mw <= 480 else "FAIL",
        "-0.4 <= LogP <= 5.6": "PASS" if -0.4 <= logp <= 5.6 else "FAIL",
        "Pass": ghose_pass
    }

    # Egan Filter (Oral Bioavailability)
    rules["Egan Filter"] = {
        "PSA <= 131.6 A2": "PASS" if psa <= 131.6 else "FAIL",
        "LogP <= 5.88": "PASS" if logp <= 5.88 else "FAIL",
        "Pass": psa <= 131.6 and logp <= 5.88
    }

    return rules


def get_functional_groups(mol):
    """Identify common functional groups present in the molecule."""
    if mol is None:
        return []

    groups = []
    smarts_patterns = {
        "Hydroxyl (-OH)": "[OX2H]",
        "Carboxyl (-COOH)": "[CX3](=O)[OX2H1]",
        "Amine (-NH2)": "[NX3;H2;!$(NC=O)]",
        "Secondary Amine (-NH-)": "[NX3;H1;!$(NC=O)]",
        "Tertiary Amine (-N<)": "[NX3;H0;!$(NC=O);!$(N=O)]",
        "Amide (-CONH-)": "[NX3][CX3](=[OX1])",
        "Ester (-COO-)": "[CX3](=O)[OX2H0]",
        "Ether (-O-)": "[OD2]([#6])[#6]",
        "Aldehyde (-CHO)": "[CX3H1](=O)[#6]",
        "Ketone (C=O)": "[#6][CX3](=O)[#6]",
        "Nitro (-NO2)": "[$([NX3](=O)=O),$([NX3+](=O)[O-])]",
        "Nitrile (-CN)": "[NX1]#[CX2]",
        "Thiol (-SH)": "[#16X2H]",
        "Sulfide (-S-)": "[#16X2H0]",
        "Sulfonamide (-SO2NH-)": "[SX4](=[OX1])(=[OX1])([NX3])",
        "Phosphate (-PO4)": "[PX4](=[OX1])([OX2])([OX2])",
        "Halogen (F)": "[F]",
        "Halogen (Cl)": "[Cl]",
        "Halogen (Br)": "[Br]",
        "Halogen (I)": "[I]",
        "Aromatic Ring": "c1ccccc1",
        "Phenol": "[OX2H][cX3]:[c]",
    }

    for name, smarts in smarts_patterns.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            matches = mol.GetSubstructMatches(pattern)
            groups.append({"Group": name, "Count": len(matches)})

    return groups
