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


# ============================================================
# COMPOUND CLASSIFICATION (Pharmacological & Chemical)
# ============================================================

def _walk_sections(section, target_heading):
    """Recursively walk PUG View section tree to find sections by heading."""
    results = []
    heading = section.get("TOCHeading", "")
    if heading.lower() == target_heading.lower():
        results.append(section)
    for child in section.get("Section", []):
        results.extend(_walk_sections(child, target_heading))
    return results


def _extract_string_values(section):
    """Extract all string values from a PUG View section's Information list."""
    values = []
    for info in section.get("Information", []):
        val = info.get("Value", {})
        for sv in val.get("StringWithMarkup", []):
            text = sv.get("String", "").strip()
            if text:
                values.append(text)
    return values


def _parse_fda_classes(info_list):
    """Parse FDA Pharmacological Classification entries by tag type.

    FDA entries have format: "Tag [CODE] - Value"
    Tags: [EPC] Established Pharmacologic Class, [MoA] Mechanism of Action,
          [PE] Physiologic Effect, [CS] Chemical Structure
    """
    epc = []
    moa = []
    pe = []
    cs = []
    seen = set()

    for info in info_list:
        if info.get("Name") != "Pharmacological Classes":
            continue
        for sv in info.get("Value", {}).get("StringWithMarkup", []):
            text = sv.get("String", "").strip()
            if not text or ";" in text:
                # Skip combined entries (they repeat individual ones)
                continue
            if text in seen:
                continue
            seen.add(text)

            if "[EPC]" in text:
                label = text.split(" - ", 1)[-1].strip() if " - " in text else text
                epc.append(label)
            elif "[MoA]" in text:
                label = text.split(" - ", 1)[-1].strip() if " - " in text else text
                moa.append(label)
            elif "[PE]" in text:
                label = text.split(" - ", 1)[-1].strip() if " - " in text else text
                pe.append(label)
            elif "[CS]" in text:
                label = text.split(" - ", 1)[-1].strip() if " - " in text else text
                cs.append(label)

    return epc, moa, pe, cs


def _parse_pharmacology_data(data):
    """Parse Pharmacology and Biochemistry heading data from PUG View."""
    result = {
        "pharmacological_group": "",
        "pharmacological_subgroups": [],
        "mechanism_of_action": "",
        "physiologic_effects": [],
        "atc_codes": [],
        "mesh_terms": [],
    }

    try:
        record = data.get("Record", {})
        sections = record.get("Section", [])

        for top_section in sections:
            # --- FDA Pharmacological Classification ---
            for s in _walk_sections(top_section, "FDA Pharmacological Classification"):
                info_list = s.get("Information", [])
                epc, moa, pe, cs = _parse_fda_classes(info_list)

                if epc:
                    result["pharmacological_group"] = epc[0]
                    if len(epc) > 1:
                        result["pharmacological_subgroups"].extend(epc[1:])
                if moa:
                    result["mechanism_of_action"] = moa[0]
                if pe:
                    result["physiologic_effects"] = pe
                if cs:
                    result["pharmacological_subgroups"].extend(cs)

            # --- MeSH Pharmacological Classification ---
            # Use the Info "Name" field (short label) not the description
            for s in _walk_sections(top_section, "MeSH Pharmacological Classification"):
                for info in s.get("Information", []):
                    name = info.get("Name", "").strip()
                    if name and name not in result["mesh_terms"]:
                        result["mesh_terms"].append(name)

            # If no FDA group, try MeSH terms as fallback
            if not result["pharmacological_group"] and result["mesh_terms"]:
                result["pharmacological_group"] = result["mesh_terms"][0]
                if len(result["mesh_terms"]) > 1:
                    result["pharmacological_subgroups"] = result["mesh_terms"][1:]

            # --- ATC Code ---
            atc_pattern = re.compile(r"^[A-Z]\d{2}[A-Z]{2}\d{2}$")
            for s in _walk_sections(top_section, "ATC Code"):
                for info in s.get("Information", []):
                    for sv in info.get("Value", {}).get("StringWithMarkup", []):
                        code = sv.get("String", "").strip()
                        if not code:
                            continue
                        # Keep actual ATC codes (e.g., N02BA01) and
                        # combined codes (e.g., "B01AC06; N02BA01")
                        # and descriptive codes (e.g., "N02BA01 - Acetylsalicylic acid")
                        if atc_pattern.match(code):
                            # Only add bare code if no descriptive version exists
                            pass  # Will be handled by descriptive version
                        elif " - " in code and atc_pattern.match(code.split(" - ")[0].strip()):
                            if code not in result["atc_codes"]:
                                result["atc_codes"].append(code)

            # --- Mechanism of Action (detailed text) ---
            if not result["mechanism_of_action"]:
                for s in _walk_sections(top_section, "Mechanism of Action"):
                    for info in s.get("Information", []):
                        for sv in info.get("Value", {}).get("StringWithMarkup", []):
                            text = sv.get("String", "").strip()
                            if text and len(text) > 20:
                                # Take first meaningful description
                                result["mechanism_of_action"] = text[:300]
                                break
                        if result["mechanism_of_action"]:
                            break

    except Exception:
        pass

    return result


def _parse_chemical_classes(data):
    """Parse Chemical Classification heading data from PUG View."""
    result = {
        "chemical_group": "",
        "chemical_subgroups": [],
    }

    try:
        record = data.get("Record", {})
        sections = record.get("Section", [])

        all_classes = []
        for top_section in sections:
            for s in _walk_sections(top_section, "Chemical Classes"):
                vals = _extract_string_values(s)
                all_classes.extend(vals)

            # Also check "Classification" heading
            if not all_classes:
                for s in _walk_sections(top_section, "Classification"):
                    vals = _extract_string_values(s)
                    all_classes.extend(vals)

            # Try ChEBI Ontology as fallback
            if not all_classes:
                for s in _walk_sections(top_section, "ChEBI Ontology"):
                    vals = _extract_string_values(s)
                    all_classes.extend(vals)

        if all_classes:
            result["chemical_group"] = all_classes[0]
            if len(all_classes) > 1:
                result["chemical_subgroups"] = all_classes[1:]

    except Exception:
        pass

    return result


def get_compound_classification(cid: int):
    """Get pharmacological and chemical classification for a compound.

    Makes two targeted PUG View API requests:
    1. Pharmacology and Biochemistry → FDA class, MeSH, ATC codes
    2. Chemical and Physical Properties → Chemical classes

    Returns a dict with classification data, or empty dict on failure.
    """
    classification = {
        "pharmacological_group": "",
        "pharmacological_subgroups": [],
        "chemical_group": "",
        "chemical_subgroups": [],
        "atc_codes": [],
        "mechanism_of_action": "",
        "physiologic_effects": [],
        "mesh_terms": [],
    }

    if not cid:
        return classification

    # Request 1: Pharmacology and Biochemistry
    try:
        url = f"{PUBCHEM_VIEW}/data/compound/{cid}/JSON?heading=Pharmacology+and+Biochemistry"
        resp = requests.get(url, timeout=10)
        if resp.status_code == 200:
            pharma_data = _parse_pharmacology_data(resp.json())
            classification.update(pharma_data)
    except Exception:
        pass

    # Request 2: Chemical and Physical Properties (for chemical classes)
    try:
        url = f"{PUBCHEM_VIEW}/data/compound/{cid}/JSON?heading=Chemical+and+Physical+Properties"
        resp = requests.get(url, timeout=10)
        if resp.status_code == 200:
            chem_data = _parse_chemical_classes(resp.json())
            classification.update(chem_data)
    except Exception:
        pass

    # If chemical classes still empty, try Names and Identifiers section
    if not classification["chemical_group"]:
        try:
            url = f"{PUBCHEM_VIEW}/data/compound/{cid}/JSON?heading=Names+and+Identifiers"
            resp = requests.get(url, timeout=10)
            if resp.status_code == 200:
                chem_data = _parse_chemical_classes(resp.json())
                if chem_data.get("chemical_group"):
                    classification["chemical_group"] = chem_data["chemical_group"]
                    classification["chemical_subgroups"] = chem_data.get("chemical_subgroups", [])
        except Exception:
            pass

    return classification
