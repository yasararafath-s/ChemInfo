"""
PubChem API integration for compound data retrieval.
"""

import re
import requests
import time
import json

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUBCHEM_VIEW = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view"

# Regex to detect CAS numbers (e.g. 50-78-2, 58-08-2, 103-90-2)
CAS_PATTERN = re.compile(r"^\d{2,7}-\d{2}-\d$")


def is_cas_number(text: str) -> bool:
    """Check if a string looks like a CAS Registry Number."""
    return bool(CAS_PATTERN.match(text.strip()))


def search_compound_by_name(name: str):
    """Search PubChem by compound name and return CID + basic info."""
    try:
        # Search for CID by name
        url = f"{PUBCHEM_BASE}/compound/name/{requests.utils.quote(name)}/JSON"
        resp = requests.get(url, timeout=15)

        if resp.status_code == 200:
            data = resp.json()
            compounds = data.get("PC_Compounds", [])
            if compounds:
                compound = compounds[0]
                cid = compound.get("id", {}).get("id", {}).get("cid")
                return extract_compound_data(compound, cid)
        elif resp.status_code == 404:
            return {"error": f"Compound '{name}' not found in PubChem."}

        return {"error": f"PubChem API error: {resp.status_code}"}

    except requests.exceptions.Timeout:
        return {"error": "PubChem request timed out. Please try again."}
    except Exception as e:
        return {"error": f"Error querying PubChem: {str(e)}"}


def get_compound_by_cid(cid: int):
    """Get compound data by PubChem CID."""
    try:
        url = f"{PUBCHEM_BASE}/compound/cid/{cid}/JSON"
        resp = requests.get(url, timeout=15)

        if resp.status_code == 200:
            data = resp.json()
            compounds = data.get("PC_Compounds", [])
            if compounds:
                return extract_compound_data(compounds[0], cid)

        return {"error": f"CID {cid} not found."}
    except Exception as e:
        return {"error": str(e)}


def extract_compound_data(compound, cid):
    """Extract useful data from PubChem compound JSON."""
    result = {"CID": cid}

    # Extract properties
    props = compound.get("props", [])
    for prop in props:
        urn = prop.get("urn", {})
        label = urn.get("label", "")
        name = urn.get("name", "")
        value = prop.get("value", {})

        if label == "IUPAC Name" and name == "Preferred":
            result["IUPAC Name"] = value.get("sval", "")
        elif label == "SMILES":
            # Accept any SMILES variant (Canonical, Isomeric, Absolute, Connectivity)
            sval = value.get("sval", "")
            if sval and "SMILES" not in result:
                result["SMILES"] = sval
        elif label == "InChI" and name == "Standard":
            result["InChI"] = value.get("sval", "")
        elif label == "InChIKey" and name == "Standard":
            result["InChIKey"] = value.get("sval", "")
        elif label == "Molecular Formula":
            result["Molecular Formula"] = value.get("sval", "")
        elif label == "Molecular Weight":
            result["PubChem MW"] = value.get("sval", "")
        elif label == "Weight" and name == "MonoIsotopic":
            result["PubChem Monoisotopic Mass"] = value.get("sval", "")

    return result


def get_compound_description(cid: int):
    """Get compound description/summary from PubChem."""
    try:
        url = f"{PUBCHEM_BASE}/compound/cid/{cid}/description/JSON"
        resp = requests.get(url, timeout=15)

        if resp.status_code == 200:
            data = resp.json()
            descriptions = data.get("InformationList", {}).get("Information", [])
            for desc in descriptions:
                description = desc.get("Description", "")
                if description and len(description) > 20:
                    return {
                        "description": description,
                        "source": desc.get("DescriptionSourceName", "PubChem")
                    }

        return {"description": "No description available.", "source": ""}
    except Exception:
        return {"description": "Could not fetch description.", "source": ""}


def get_compound_synonyms(cid: int, limit=15):
    """Get compound synonyms/alternative names."""
    try:
        url = f"{PUBCHEM_BASE}/compound/cid/{cid}/synonyms/JSON"
        resp = requests.get(url, timeout=15)

        if resp.status_code == 200:
            data = resp.json()
            info_list = data.get("InformationList", {}).get("Information", [])
            if info_list:
                synonyms = info_list[0].get("Synonym", [])
                return synonyms[:limit]

        return []
    except Exception:
        return []


def get_compound_safety(cid: int):
    """Get GHS safety information."""
    try:
        url = f"{PUBCHEM_VIEW}/data/compound/{cid}/JSON?heading=GHS+Classification"
        resp = requests.get(url, timeout=15)

        if resp.status_code == 200:
            data = resp.json()
            return {"available": True, "data": "GHS data available on PubChem"}
        return {"available": False}
    except Exception:
        return {"available": False}


def get_pubchem_links(cid: int):
    """Generate useful PubChem reference links."""
    if cid is None:
        return {}

    return {
        "PubChem Compound": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}",
        "PubChem 3D Viewer": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}#section=3D-Conformer",
        "PubChem Safety": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}#section=Safety-and-Hazards",
        "PubChem Spectra": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}#section=Spectra",
        "PubChem BioAssays": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}#section=Biological-Test-Results",
        "PubChem Literature": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}#section=Literature",
        "PubChem Patents": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}#section=Patents",
    }


def get_external_references(cid: int, inchikey: str = None, cas_number: str = None):
    """Generate external database reference links for a compound.

    Returns a dict of {database_name: url} for various chemistry databases.
    """
    refs = {}

    # PubChem
    if cid:
        refs["PubChem"] = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"

    # ChemSpider (search by InChIKey or name)
    if inchikey and inchikey != "N/A":
        refs["ChemSpider"] = f"https://www.chemspider.com/Search.aspx?q={inchikey}"

    # DrugBank (search)
    if cid:
        refs["DrugBank"] = f"https://go.drugbank.com/unearth/q?searcher=drugs&query={cid}"

    # CompTox EPA (search by CAS or InChIKey)
    if cas_number:
        refs["CompTox (EPA)"] = f"https://comptox.epa.gov/dashboard/chemical/details/{cas_number}"
    elif inchikey and inchikey != "N/A":
        refs["CompTox (EPA)"] = f"https://comptox.epa.gov/dashboard/chemical/details/{inchikey}"

    # ChEBI
    if inchikey and inchikey != "N/A":
        refs["ChEBI"] = f"https://www.ebi.ac.uk/chebi/advancedSearchFT.do?searchString={inchikey}"

    # Wikipedia
    if cas_number:
        refs["Wikipedia"] = f"https://en.wikipedia.org/wiki/Special:Search?search={cas_number}"

    # NIST Chemistry WebBook
    if cas_number:
        refs["NIST WebBook"] = f"https://webbook.nist.gov/cgi/cbook.cgi?ID={cas_number}&Units=SI"

    return refs


def get_similar_compounds(cid: int, threshold=90, max_results=5):
    """Find structurally similar compounds in PubChem."""
    try:
        url = (f"{PUBCHEM_BASE}/compound/fastsimilarity_2d/cid/{cid}/cids/JSON"
               f"?Threshold={threshold}&MaxRecords={max_results}")
        resp = requests.get(url, timeout=20)

        if resp.status_code == 200:
            data = resp.json()
            cids = data.get("IdentifierList", {}).get("CID", [])
            # Remove the query compound itself
            cids = [c for c in cids if c != cid][:max_results]
            return cids

        return []
    except Exception:
        return []


def search_compound_by_cas(cas_number: str):
    """Search PubChem by CAS Registry Number and return compound data."""
    try:
        # PubChem can resolve CAS numbers via the /name/ endpoint
        url = f"{PUBCHEM_BASE}/compound/name/{requests.utils.quote(cas_number.strip())}/JSON"
        resp = requests.get(url, timeout=15)

        if resp.status_code == 200:
            data = resp.json()
            compounds = data.get("PC_Compounds", [])
            if compounds:
                compound = compounds[0]
                cid = compound.get("id", {}).get("id", {}).get("cid")
                result = extract_compound_data(compound, cid)
                result["CAS Number"] = cas_number.strip()
                return result
        elif resp.status_code == 404:
            return {"error": f"CAS Number '{cas_number}' not found in PubChem."}

        return {"error": f"PubChem API error: {resp.status_code}"}

    except requests.exceptions.Timeout:
        return {"error": "PubChem request timed out. Please try again."}
    except Exception as e:
        return {"error": f"Error querying PubChem: {str(e)}"}


def get_cas_number(cid: int):
    """Retrieve CAS Registry Number for a compound from PubChem synonyms.

    CAS numbers are typically found in the synonyms list and follow the
    pattern: digits-digits-digit (e.g. 50-78-2 for aspirin).
    """
    try:
        synonyms = get_compound_synonyms(cid, limit=50)
        for syn in synonyms:
            if CAS_PATTERN.match(syn.strip()):
                return syn.strip()
        return None
    except Exception:
        return None
