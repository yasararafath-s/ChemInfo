# -*- coding: utf-8 -*-
"""
ChemInfo Dashboard - Comprehensive Chemical Compound Analyzer
=============================================================
"""

import streamlit as st
import pandas as pd
import time
import traceback
from io import BytesIO

try:
    # RDKit imports
    from rdkit import Chem
    from rdkit.Chem import Draw

    # Local utility modules
    from utils.chem_properties import (
        get_mol_from_smiles,
        generate_structure_image,
        calculate_physicochemical_properties,
        evaluate_drug_likeness,
        get_functional_groups,
    )
    from utils.mrm_calculator import calculate_mrm_data
    from utils.pubchem_api import (
        search_compound_by_name,
        search_compound_by_cas,
        get_compound_description,
        get_compound_synonyms,
        get_pubchem_links,
        get_external_references,
        get_similar_compounds,
        get_compound_by_cid,
        get_cas_number,
        is_cas_number,
    )
except Exception as e:
    st.error(f"Failed to import dependencies: {e}")
    st.code(traceback.format_exc())
    st.stop()


# ============================================================
# PAGE CONFIGURATION
# ============================================================
st.set_page_config(
    page_title="ChemInfo Dashboard",
    page_icon="\U0001f9ea",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ============================================================
# THEME TOGGLE
# ============================================================
if "theme" not in st.session_state:
    st.session_state.theme = "Dark"

# Theme selector in sidebar (placed before other sidebar content)
with st.sidebar:
    st.session_state.theme = st.selectbox(
        "Theme",
        ["Light", "Dark"],
        index=["Light", "Dark"].index(st.session_state.theme),
        key="theme_selector",
    )

is_dark = st.session_state.theme == "Dark"

# ============================================================
# CUSTOM CSS (adapts to selected theme)
# ============================================================
if is_dark:
    theme_css = """
    <style>
        /* ===== DARK THEME - Comprehensive Styling ===== */

        /* --- Base App --- */
        .stApp {
            background-color: #0d1117;
            color: #e6edf3;
        }

        /* --- Custom Headers --- */
        .main-header {
            font-size: 2.4rem;
            font-weight: 800;
            background: linear-gradient(135deg, #58a6ff 0%, #a371f7 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            text-align: center;
            margin-bottom: 0.3rem;
            letter-spacing: -0.5px;
        }
        .sub-header {
            font-size: 0.95rem;
            color: #8b949e;
            text-align: center;
            margin-bottom: 2rem;
            letter-spacing: 0.3px;
        }
        .section-header {
            font-size: 1.3rem;
            font-weight: 600;
            color: #79c0ff;
            border-bottom: 2px solid #1f6feb;
            padding-bottom: 6px;
            margin-top: 1.5rem;
        }

        /* --- Sidebar --- */
        section[data-testid="stSidebar"] {
            background-color: #161b22;
            border-right: 1px solid #21262d;
        }
        section[data-testid="stSidebar"] * { color: #e6edf3 !important; }
        section[data-testid="stSidebar"] .stSelectbox label,
        section[data-testid="stSidebar"] .stRadio label,
        section[data-testid="stSidebar"] .stTextInput label,
        section[data-testid="stSidebar"] h1, section[data-testid="stSidebar"] h2,
        section[data-testid="stSidebar"] h3, section[data-testid="stSidebar"] h4,
        section[data-testid="stSidebar"] p, section[data-testid="stSidebar"] li,
        section[data-testid="stSidebar"] span,
        section[data-testid="stSidebar"] .stMarkdown { color: #e6edf3 !important; }
        section[data-testid="stSidebar"] small,
        section[data-testid="stSidebar"] .stCaption,
        section[data-testid="stSidebar"] caption { color: #8b949e !important; }
        section[data-testid="stSidebar"] hr { border-color: #30363d !important; }

        /* --- Sidebar Inputs --- */
        section[data-testid="stSidebar"] .stTextInput input,
        section[data-testid="stSidebar"] .stSelectbox > div > div {
            background-color: #0d1117 !important;
            border-color: #30363d !important;
            color: #e6edf3 !important;
        }

        /* --- Main Content Inputs --- */
        .stTextInput input, .stTextArea textarea, .stNumberInput input {
            background-color: #161b22 !important;
            border: 1px solid #30363d !important;
            color: #e6edf3 !important;
            border-radius: 8px !important;
        }
        .stTextInput input:focus, .stTextArea textarea:focus {
            border-color: #1f6feb !important;
            box-shadow: 0 0 0 3px rgba(31, 111, 235, 0.2) !important;
        }
        .stSelectbox > div > div {
            background-color: #161b22 !important;
            border-color: #30363d !important;
            color: #e6edf3 !important;
        }

        /* --- Tabs --- */
        .stTabs [data-baseweb="tab-list"] {
            background-color: #161b22;
            border-radius: 10px;
            padding: 4px;
            gap: 4px;
            border: 1px solid #21262d;
        }
        .stTabs [data-baseweb="tab-list"] button {
            font-size: 1.0rem;
            font-weight: 500;
            color: #8b949e !important;
            background-color: transparent;
            border-radius: 8px;
            padding: 8px 16px;
            border: none !important;
        }
        .stTabs [data-baseweb="tab-list"] button[aria-selected="true"] {
            background-color: #1f6feb !important;
            color: #ffffff !important;
            font-weight: 600;
        }
        .stTabs [data-baseweb="tab-list"] button:hover {
            color: #e6edf3 !important;
            background-color: #21262d;
        }

        /* --- Metrics --- */
        div[data-testid="stMetric"] {
            background-color: #161b22;
            border: 1px solid #21262d;
            border-radius: 10px;
            padding: 12px 16px;
        }
        div[data-testid="stMetric"] label {
            color: #8b949e !important;
            font-size: 0.85rem !important;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }
        div[data-testid="stMetric"] div[data-testid="stMetricValue"] {
            color: #58a6ff !important;
            font-weight: 700 !important;
        }

        /* --- DataFrames / Tables --- */
        .stDataFrame, div[data-testid="stDataFrame"] {
            border-radius: 10px;
            overflow: hidden;
            border: 1px solid #21262d;
        }
        .stDataFrame [data-testid="glideDataEditor"],
        div[data-testid="stDataFrame"] > div {
            background-color: #0d1117 !important;
        }

        /* --- Expanders --- */
        div[data-testid="stExpander"] {
            background-color: #161b22;
            border: 1px solid #21262d !important;
            border-radius: 10px !important;
            margin-bottom: 8px;
        }
        div[data-testid="stExpander"] details summary {
            padding: 12px 16px;
        }
        div[data-testid="stExpander"] details summary p {
            font-size: 1.05rem;
            font-weight: 600;
            color: #e6edf3 !important;
        }
        div[data-testid="stExpander"] details > div {
            border-top: 1px solid #21262d;
        }

        /* --- Info / Success / Warning / Error Boxes --- */
        div[data-testid="stAlert"] {
            border-radius: 10px !important;
            border: none !important;
        }
        .stAlert > div[data-baseweb="notification"][kind="info"] {
            background-color: rgba(31, 111, 235, 0.12) !important;
            color: #79c0ff !important;
        }
        .stAlert > div[data-baseweb="notification"][kind="success"] {
            background-color: rgba(63, 185, 80, 0.12) !important;
            color: #56d364 !important;
        }
        .stAlert > div[data-baseweb="notification"][kind="warning"] {
            background-color: rgba(210, 153, 34, 0.12) !important;
            color: #e3b341 !important;
        }
        .stAlert > div[data-baseweb="notification"][kind="negative"] {
            background-color: rgba(248, 81, 73, 0.12) !important;
            color: #f85149 !important;
        }
        /* Streamlit info/success/warning/error direct styling */
        div[role="alert"] {
            border-radius: 10px !important;
        }

        /* --- Buttons --- */
        .stButton > button {
            background-color: #21262d !important;
            color: #e6edf3 !important;
            border: 1px solid #30363d !important;
            border-radius: 8px !important;
            font-weight: 500;
            transition: all 0.2s ease;
        }
        .stButton > button:hover {
            background-color: #30363d !important;
            border-color: #8b949e !important;
        }
        .stButton > button[kind="primary"],
        .stFormSubmitButton > button {
            background: linear-gradient(135deg, #1f6feb, #388bfd) !important;
            color: #ffffff !important;
            border: none !important;
            font-weight: 600;
        }
        .stFormSubmitButton > button:hover {
            background: linear-gradient(135deg, #388bfd, #58a6ff) !important;
        }

        /* --- Download Buttons --- */
        .stDownloadButton > button {
            background-color: #161b22 !important;
            color: #58a6ff !important;
            border: 1px solid #1f6feb !important;
            border-radius: 8px !important;
            font-weight: 500;
        }
        .stDownloadButton > button:hover {
            background-color: #1f6feb !important;
            color: #ffffff !important;
        }

        /* --- Code Blocks --- */
        .stCode, code, pre {
            background-color: #161b22 !important;
            border: 1px solid #21262d !important;
            border-radius: 8px !important;
            color: #e6edf3 !important;
        }

        /* --- Markdown Links --- */
        a { color: #58a6ff !important; }
        a:hover { color: #79c0ff !important; text-decoration: underline; }

        /* --- Dividers --- */
        hr { border-color: #21262d !important; }

        /* --- Progress Bar --- */
        .stProgress > div > div > div {
            background: linear-gradient(90deg, #1f6feb, #a371f7) !important;
            border-radius: 10px;
        }
        .stProgress > div > div {
            background-color: #21262d !important;
            border-radius: 10px;
        }

        /* --- File Uploader --- */
        div[data-testid="stFileUploader"] {
            background-color: #161b22;
            border: 2px dashed #30363d !important;
            border-radius: 10px;
            padding: 16px;
        }

        /* --- Markdown Tables (Welcome screen) --- */
        .stMarkdown table {
            border-collapse: collapse;
            width: 100%;
        }
        .stMarkdown th {
            background-color: #161b22 !important;
            color: #58a6ff !important;
            border: 1px solid #30363d !important;
            padding: 10px 14px !important;
            font-weight: 600;
        }
        .stMarkdown td {
            background-color: #0d1117 !important;
            color: #e6edf3 !important;
            border: 1px solid #21262d !important;
            padding: 8px 14px !important;
        }
        .stMarkdown tr:hover td {
            background-color: #161b22 !important;
        }

        /* --- Captions --- */
        .stCaption, caption, small {
            color: #8b949e !important;
        }

        /* --- Images (structure viewer) --- */
        .stImage {
            background-color: #ffffff;
            border-radius: 10px;
            padding: 8px;
            border: 1px solid #21262d;
        }

        /* --- JSON viewer --- */
        div[data-testid="stJson"] {
            background-color: #161b22 !important;
            border: 1px solid #21262d !important;
            border-radius: 10px !important;
        }

        /* --- Spinner --- */
        .stSpinner > div { color: #58a6ff !important; }

        /* --- Badges --- */
        .pass-badge {
            background: linear-gradient(135deg, #238636, #2ea043);
            color: white;
            padding: 4px 12px;
            border-radius: 20px;
            font-weight: 600;
            font-size: 0.85rem;
        }
        .fail-badge {
            background: linear-gradient(135deg, #da3633, #f85149);
            color: white;
            padding: 4px 12px;
            border-radius: 20px;
            font-weight: 600;
            font-size: 0.85rem;
        }

        /* --- Scrollbar --- */
        ::-webkit-scrollbar { width: 8px; height: 8px; }
        ::-webkit-scrollbar-track { background: #0d1117; }
        ::-webkit-scrollbar-thumb { background: #30363d; border-radius: 4px; }
        ::-webkit-scrollbar-thumb:hover { background: #484f58; }
    </style>
    """
else:
    theme_css = """
    <style>
        /* ===== LIGHT THEME - Comprehensive Styling ===== */

        /* --- Base App --- */
        .stApp {
            background-color: #ffffff;
            color: #1f2937;
        }

        /* --- Custom Headers --- */
        .main-header {
            font-size: 2.4rem;
            font-weight: 800;
            background: linear-gradient(135deg, #1E3A5F 0%, #2563eb 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            text-align: center;
            margin-bottom: 0.3rem;
            letter-spacing: -0.5px;
        }
        .sub-header {
            font-size: 0.95rem;
            color: #6b7280;
            text-align: center;
            margin-bottom: 2rem;
            letter-spacing: 0.3px;
        }
        .section-header {
            font-size: 1.3rem;
            font-weight: 600;
            color: #1E3A5F;
            border-bottom: 2px solid #3b82f6;
            padding-bottom: 6px;
            margin-top: 1.5rem;
        }

        /* --- Sidebar --- */
        section[data-testid="stSidebar"] {
            background-color: #f8fafc;
            border-right: 1px solid #e5e7eb;
        }
        section[data-testid="stSidebar"] * { color: #1f2937 !important; }
        section[data-testid="stSidebar"] .stSelectbox label,
        section[data-testid="stSidebar"] .stRadio label,
        section[data-testid="stSidebar"] .stTextInput label,
        section[data-testid="stSidebar"] h1, section[data-testid="stSidebar"] h2,
        section[data-testid="stSidebar"] h3, section[data-testid="stSidebar"] h4,
        section[data-testid="stSidebar"] p, section[data-testid="stSidebar"] li,
        section[data-testid="stSidebar"] span,
        section[data-testid="stSidebar"] .stMarkdown { color: #1f2937 !important; }
        section[data-testid="stSidebar"] small,
        section[data-testid="stSidebar"] .stCaption,
        section[data-testid="stSidebar"] caption { color: #6b7280 !important; }
        section[data-testid="stSidebar"] hr { border-color: #e5e7eb !important; }

        /* --- Tabs --- */
        .stTabs [data-baseweb="tab-list"] {
            background-color: #f1f5f9;
            border-radius: 10px;
            padding: 4px;
            gap: 4px;
            border: 1px solid #e2e8f0;
        }
        .stTabs [data-baseweb="tab-list"] button {
            font-size: 1.0rem;
            font-weight: 500;
            color: #64748b !important;
            background-color: transparent;
            border-radius: 8px;
            padding: 8px 16px;
            border: none !important;
        }
        .stTabs [data-baseweb="tab-list"] button[aria-selected="true"] {
            background-color: #2563eb !important;
            color: #ffffff !important;
            font-weight: 600;
        }
        .stTabs [data-baseweb="tab-list"] button:hover {
            color: #1f2937 !important;
            background-color: #e2e8f0;
        }

        /* --- Metrics --- */
        div[data-testid="stMetric"] {
            background-color: #f8fafc;
            border: 1px solid #e2e8f0;
            border-radius: 10px;
            padding: 12px 16px;
        }
        div[data-testid="stMetric"] label {
            color: #6b7280 !important;
            font-size: 0.85rem !important;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }
        div[data-testid="stMetric"] div[data-testid="stMetricValue"] {
            color: #1E3A5F !important;
            font-weight: 700 !important;
        }

        /* --- DataFrames / Tables --- */
        .stDataFrame, div[data-testid="stDataFrame"] {
            border-radius: 10px;
            overflow: hidden;
            border: 1px solid #e2e8f0;
        }

        /* --- Expanders --- */
        div[data-testid="stExpander"] {
            background-color: #f8fafc;
            border: 1px solid #e2e8f0 !important;
            border-radius: 10px !important;
            margin-bottom: 8px;
        }
        div[data-testid="stExpander"] details summary p {
            font-size: 1.05rem;
            font-weight: 600;
            color: #1f2937 !important;
        }
        div[data-testid="stExpander"] details > div {
            border-top: 1px solid #e2e8f0;
        }

        /* --- Buttons --- */
        .stFormSubmitButton > button {
            background: linear-gradient(135deg, #2563eb, #3b82f6) !important;
            color: #ffffff !important;
            border: none !important;
            border-radius: 8px !important;
            font-weight: 600;
        }
        .stFormSubmitButton > button:hover {
            background: linear-gradient(135deg, #3b82f6, #60a5fa) !important;
        }

        /* --- Download Buttons --- */
        .stDownloadButton > button {
            background-color: #f8fafc !important;
            color: #2563eb !important;
            border: 1px solid #3b82f6 !important;
            border-radius: 8px !important;
            font-weight: 500;
        }
        .stDownloadButton > button:hover {
            background-color: #2563eb !important;
            color: #ffffff !important;
        }

        /* --- Progress Bar --- */
        .stProgress > div > div > div {
            background: linear-gradient(90deg, #2563eb, #7c3aed) !important;
            border-radius: 10px;
        }

        /* --- Images (structure viewer) --- */
        .stImage {
            border-radius: 10px;
            border: 1px solid #e2e8f0;
        }

        /* --- Markdown Tables --- */
        .stMarkdown th {
            background-color: #f1f5f9 !important;
            color: #1E3A5F !important;
            border: 1px solid #e2e8f0 !important;
            padding: 10px 14px !important;
            font-weight: 600;
        }
        .stMarkdown td {
            border: 1px solid #e2e8f0 !important;
            padding: 8px 14px !important;
        }
        .stMarkdown tr:hover td {
            background-color: #f8fafc !important;
        }

        /* --- Badges --- */
        .pass-badge {
            background: linear-gradient(135deg, #16a34a, #22c55e);
            color: white;
            padding: 4px 12px;
            border-radius: 20px;
            font-weight: 600;
            font-size: 0.85rem;
        }
        .fail-badge {
            background: linear-gradient(135deg, #dc2626, #ef4444);
            color: white;
            padding: 4px 12px;
            border-radius: 20px;
            font-weight: 600;
            font-size: 0.85rem;
        }
    </style>
    """

st.markdown(theme_css, unsafe_allow_html=True)


# ============================================================
# HEADER
# ============================================================
st.markdown('<p class="main-header">ChemInfo Dashboard</p>', unsafe_allow_html=True)
st.markdown('<p class="sub-header">Comprehensive Chemical Compound Analyzer | Structure | Properties | MRM | References</p>', unsafe_allow_html=True)


# ============================================================
# SIDEBAR - INPUT
# ============================================================
with st.sidebar:
    st.header("Compound Input")

    input_mode = st.radio(
        "Input Method:",
        ["Chemical Name", "CAS Number", "SMILES String", "Batch Upload (CSV/Excel)"],
        index=0,
    )

    compound_name = None
    smiles_input = None
    cas_input = None
    batch_df = None
    analyze_btn = False

    if input_mode == "Chemical Name":
        with st.form(key="name_form"):
            compound_name = st.text_input(
                "Enter compound name:",
                placeholder="e.g., Aspirin, Caffeine, Ibuprofen",
            )
            analyze_btn = st.form_submit_button(">> Analyze", type="primary", use_container_width=True)
        st.caption("Type a name and press Enter or click Analyze")

    elif input_mode == "CAS Number":
        with st.form(key="cas_form"):
            cas_input = st.text_input(
                "Enter CAS Number:",
                placeholder="e.g., 50-78-2, 58-08-2, 103-90-2",
            )
            analyze_btn = st.form_submit_button(">> Analyze", type="primary", use_container_width=True)
        st.caption("Enter a CAS Registry Number (e.g., 50-78-2 for Aspirin)")

    elif input_mode == "SMILES String":
        with st.form(key="smiles_form"):
            smiles_input = st.text_input(
                "Enter SMILES:",
                placeholder="e.g., CC(=O)OC1=CC=CC=C1C(=O)O",
            )
            compound_name = st.text_input(
                "Compound name (optional):",
                placeholder="e.g., Aspirin",
            )
            analyze_btn = st.form_submit_button(">> Analyze", type="primary", use_container_width=True)
        st.caption("Type SMILES and press Enter or click Analyze")

    elif input_mode == "Batch Upload (CSV/Excel)":
        st.markdown("**Upload a file with compound names or SMILES**")
        uploaded_file = st.file_uploader(
            "Upload CSV or Excel",
            type=["csv", "xlsx", "xls"],
        )
        if uploaded_file:
            try:
                if uploaded_file.name.endswith(".csv"):
                    batch_df = pd.read_csv(uploaded_file)
                else:
                    batch_df = pd.read_excel(uploaded_file)
                st.success(f"Loaded {len(batch_df)} rows")
                st.dataframe(batch_df.head(), use_container_width=True)

                # Select column
                col_options = batch_df.columns.tolist()
                selected_col = st.selectbox(
                    "Select the column with compound names/SMILES:",
                    col_options,
                )
                input_type_col = st.radio(
                    "Column contains:",
                    ["Compound Names", "SMILES"],
                )
            except Exception as e:
                st.error(f"Error reading file: {e}")
                batch_df = None

        analyze_btn = st.button(">> Analyze", type="primary", use_container_width=True)

    st.divider()
    st.markdown("### About")
    st.markdown("""
    **ChemInfo Dashboard** provides:
    - 2D Structure Visualization
    - Physicochemical Properties
    - Drug-Likeness Rules
    - MRM Mass Spec Data
    - PubChem References
    - Batch Processing
    """)


# ============================================================
# HELPER FUNCTIONS
# ============================================================
def analyze_single_compound(name=None, smiles=None, cas=None):
    """Analyze a single compound and return all data."""
    result = {}

    # Step 1: Resolve compound via PubChem
    pubchem_data = None

    # --- CAS Number lookup ---
    if cas:
        with st.spinner(f"Searching PubChem for CAS '{cas}'..."):
            pubchem_data = search_compound_by_cas(cas)

        if "error" in pubchem_data:
            st.error(pubchem_data["error"])
            return None

        smiles = pubchem_data.get("SMILES") or pubchem_data.get("Canonical SMILES") or pubchem_data.get("Isomeric SMILES")
        if not smiles:
            st.error(f"Could not find SMILES for CAS '{cas}' in PubChem.")
            return None

        # Use IUPAC name as the display name
        name = pubchem_data.get("IUPAC Name", f"CAS {cas}")

    # --- Name lookup ---
    elif name and not smiles:
        with st.spinner(f"Searching PubChem for '{name}'..."):
            pubchem_data = search_compound_by_name(name)

        if "error" in pubchem_data:
            st.error(pubchem_data["error"])
            return None

        smiles = pubchem_data.get("SMILES") or pubchem_data.get("Canonical SMILES") or pubchem_data.get("Isomeric SMILES")
        if not smiles:
            st.error(f"Could not find SMILES for '{name}' in PubChem.")
            return None

    elif name and smiles:
        # If both provided, still try PubChem for extra data
        pubchem_data = search_compound_by_name(name)
        if "error" in pubchem_data:
            pubchem_data = None

    if not smiles:
        st.error("No SMILES string available. Please provide a compound name, CAS number, or SMILES.")
        return None

    # Step 2: Create RDKit Mol object
    mol = get_mol_from_smiles(smiles)
    if mol is None:
        st.error(f"Invalid SMILES: {smiles}")
        return None

    result["smiles"] = smiles
    result["mol"] = mol
    result["pubchem"] = pubchem_data
    result["name"] = name or (pubchem_data.get("IUPAC Name", "Unknown") if pubchem_data else "Unknown")

    # Step 3: Calculate properties
    result["properties"] = calculate_physicochemical_properties(mol)
    result["drug_likeness"] = evaluate_drug_likeness(result["properties"])
    result["functional_groups"] = get_functional_groups(mol)
    result["mrm"] = calculate_mrm_data(mol)

    # Step 4: PubChem extras + CAS number + external references
    if pubchem_data and pubchem_data.get("CID"):
        cid = pubchem_data["CID"]
        result["description"] = get_compound_description(cid)
        result["synonyms"] = get_compound_synonyms(cid)
        result["links"] = get_pubchem_links(cid)
        # Get CAS number if not already available
        if pubchem_data.get("CAS Number"):
            result["cas_number"] = pubchem_data["CAS Number"]
        else:
            result["cas_number"] = get_cas_number(cid)
        # External database references
        inchikey = result["properties"].get("InChIKey")
        result["external_refs"] = get_external_references(
            cid=cid,
            inchikey=inchikey,
            cas_number=result["cas_number"],
        )
    else:
        result["description"] = None
        result["synonyms"] = []
        result["links"] = {}
        result["cas_number"] = cas if cas else None
        result["external_refs"] = {}

    return result


def display_compound_report(result):
    """Display the full compound analysis report."""
    if result is None:
        return

    name = result.get("name", "Unknown")
    smiles = result.get("smiles", "")
    mol = result.get("mol")
    props = result.get("properties", {})
    drug_rules = result.get("drug_likeness", {})
    func_groups = result.get("functional_groups", [])
    mrm = result.get("mrm", {})
    pubchem = result.get("pubchem", {})
    description = result.get("description")
    synonyms = result.get("synonyms", [])
    links = result.get("links", {})
    cas_number = result.get("cas_number")

    # ---- Title & Structure ----
    st.markdown(f'<p class="section-header">{name}</p>', unsafe_allow_html=True)

    col_struct, col_info = st.columns([1, 1.5])

    with col_struct:
        img = generate_structure_image(mol, size=(500, 400))
        if img:
            st.image(img, caption=f"2D Structure of {name}", use_container_width=True)

        st.code(smiles, language=None)

        # Show CAS Number and PubChem CID
        id_parts = []
        if cas_number:
            id_parts.append(f"**CAS:** {cas_number}")
        if pubchem and pubchem.get("CID"):
            id_parts.append(f"**PubChem CID:** [{pubchem['CID']}](https://pubchem.ncbi.nlm.nih.gov/compound/{pubchem['CID']})")
        if id_parts:
            st.markdown(" &nbsp;|&nbsp; ".join(id_parts))

    with col_info:
        if description and description.get("description"):
            st.markdown("**Description:**")
            st.info(description["description"][:500])

        # Key metrics
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("MW", f"{props.get('Molecular Weight (g/mol)', 'N/A')}")
        m2.metric("LogP", f"{props.get('ALogP (Wildman-Crippen)', 'N/A')}")
        m3.metric("PSA", f"{props.get('Topological PSA (A2)', 'N/A')}")
        m4.metric("QED", f"{props.get('QED (Drug-Likeness)', 'N/A')}")

        m5, m6, m7, m8 = st.columns(4)
        m5.metric("HBD", f"{props.get('H-Bond Donors (HBD)', 'N/A')}")
        m6.metric("HBA", f"{props.get('H-Bond Acceptors (HBA)', 'N/A')}")
        m7.metric("RotBonds", f"{props.get('Rotatable Bonds', 'N/A')}")
        m8.metric("Rings", f"{props.get('Number of Rings', 'N/A')}")

        if synonyms:
            with st.expander(f"Synonyms ({len(synonyms)} names)"):
                st.write(", ".join(synonyms))

    # ---- TABS ----
    tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
        "Properties",
        "Drug-Likeness",
        "MRM / Mass Spec",
        "Functional Groups",
        "References",
        "Export",
    ])

    # ---- TAB 1: Properties ----
    with tab1:
        st.markdown('<p class="section-header">Physicochemical Properties</p>', unsafe_allow_html=True)

        # Build property list with CAS at the top
        prop_rows = []
        if cas_number:
            prop_rows.append({"Property": "CAS Registry Number", "Value": cas_number})
        prop_rows.extend(
            [{"Property": k, "Value": v} for k, v in props.items()
             if k not in ["InChI", "InChIKey"]]
        )
        prop_df = pd.DataFrame(prop_rows)
        st.dataframe(prop_df, use_container_width=True, hide_index=True, height=500)

        # InChI in expander (they're long)
        with st.expander("InChI & InChIKey"):
            st.text(f"InChI:    {props.get('InChI', 'N/A')}")
            st.text(f"InChIKey: {props.get('InChIKey', 'N/A')}")

    # ---- TAB 2: Drug-Likeness ----
    with tab2:
        st.markdown('<p class="section-header">Drug-Likeness Evaluation</p>', unsafe_allow_html=True)

        for rule_name, rule_data in drug_rules.items():
            is_pass = rule_data.get("Pass", False)

            with st.expander(f"{rule_name}  {'[PASS]' if is_pass else '[FAIL]'}", expanded=True):
                criteria = {k: v for k, v in rule_data.items() if k not in ["Pass", "Violations"]}
                criteria_df = pd.DataFrame(
                    [{"Criterion": k, "Result": v} for k, v in criteria.items()]
                )
                st.dataframe(criteria_df, use_container_width=True, hide_index=True)

                if "Violations" in rule_data:
                    st.markdown(f"**Violations:** {rule_data['Violations']}")

    # ---- TAB 3: MRM / Mass Spec ----
    with tab3:
        st.markdown('<p class="section-header">MRM Mass Spectrometry Data</p>', unsafe_allow_html=True)

        if mrm:
            # Basic mass info
            mc1, mc2, mc3 = st.columns(3)
            mc1.metric("Monoisotopic Mass", f"{mrm.get('Monoisotopic Mass', 'N/A')}")
            mc2.metric("Average Mass", f"{mrm.get('Average Mass', 'N/A')}")
            mc3.metric("Formula", f"{mrm.get('Molecular Formula', 'N/A')}")

            # Adducts
            st.markdown("#### Positive Mode Adducts")
            pos_adducts = mrm.get("Positive Mode Adducts", {})
            if pos_adducts:
                pos_df = pd.DataFrame(
                    [{"Adduct": k, "m/z": v} for k, v in pos_adducts.items()]
                )
                st.dataframe(pos_df, use_container_width=True, hide_index=True)

            st.markdown("#### Negative Mode Adducts")
            neg_adducts = mrm.get("Negative Mode Adducts", {})
            if neg_adducts:
                neg_df = pd.DataFrame(
                    [{"Adduct": k, "m/z": v} for k, v in neg_adducts.items()]
                )
                st.dataframe(neg_df, use_container_width=True, hide_index=True)

            # Predicted fragments
            st.markdown("#### Predicted Fragment Ions (Neutral Losses)")
            col_fp, col_fn = st.columns(2)

            with col_fp:
                st.markdown("**Positive Mode**")
                pos_frags = mrm.get("Predicted Fragments (Positive)", {})
                if pos_frags:
                    pf_df = pd.DataFrame(
                        [{"Fragment": k, "m/z": v} for k, v in pos_frags.items()]
                    )
                    st.dataframe(pf_df, use_container_width=True, hide_index=True)
                else:
                    st.info("No predicted fragments in positive mode.")

            with col_fn:
                st.markdown("**Negative Mode**")
                neg_frags = mrm.get("Predicted Fragments (Negative)", {})
                if neg_frags:
                    nf_df = pd.DataFrame(
                        [{"Fragment": k, "m/z": v} for k, v in neg_frags.items()]
                    )
                    st.dataframe(nf_df, use_container_width=True, hide_index=True)
                else:
                    st.info("No predicted fragments in negative mode.")

            # Suggested MRM transitions
            st.markdown("#### Suggested MRM Transitions")
            transitions = mrm.get("Suggested MRM Transitions", [])
            if transitions:
                trans_df = pd.DataFrame(transitions)
                st.dataframe(trans_df, use_container_width=True, hide_index=True)
                st.caption("Note: Collision energies are estimated. Optimize experimentally for your instrument.")
            else:
                st.info("No MRM transitions could be predicted.")

            # Isotope info
            st.markdown("#### Isotope Pattern Information")
            iso_info = mrm.get("Isotope Info", {})
            if iso_info:
                for k, v in iso_info.items():
                    st.markdown(f"- **{k}:** {v}")

    # ---- TAB 4: Functional Groups ----
    with tab4:
        st.markdown('<p class="section-header">Functional Groups Detected</p>', unsafe_allow_html=True)

        if func_groups:
            fg_df = pd.DataFrame(func_groups)
            st.dataframe(fg_df, use_container_width=True, hide_index=True)
        else:
            st.info("No common functional groups detected.")

    # ---- TAB 5: References ----
    with tab5:
        st.markdown('<p class="section-header">References & External Links</p>', unsafe_allow_html=True)

        external_refs = result.get("external_refs", {})

        # External Database Links (like ChemDash)
        if external_refs:
            col_ref1, col_ref2 = st.columns(2)
            with col_ref1:
                st.markdown("#### Database References")
                if cas_number:
                    st.metric("CAS Number", cas_number)
                if pubchem and pubchem.get("CID"):
                    st.metric("PubChem CID", pubchem["CID"])
                st.markdown("#### External Links")
                link_icons = {
                    "PubChem": "\U0001f310",
                    "ChemSpider": "\U0001f52c",
                    "DrugBank": "\U0001f48a",
                    "CompTox (EPA)": "\U0001f3ed",
                    "ChEBI": "\U0001f9ea",
                    "Wikipedia": "\U0001f4d6",
                    "NIST WebBook": "\U0001f4ca",
                }
                for ref_name, ref_url in external_refs.items():
                    icon = link_icons.get(ref_name, "\U0001f517")
                    st.markdown(f"{icon} [{ref_name}]({ref_url})")

            with col_ref2:
                st.markdown("#### Compound Identifiers")
                if pubchem and pubchem.get("IUPAC Name"):
                    st.text_area("IUPAC Name", pubchem["IUPAC Name"], height=80)
                if props.get("InChIKey") and props["InChIKey"] != "N/A":
                    st.code(f"InChIKey: {props['InChIKey']}", language=None)
                if synonyms:
                    with st.expander(f"Synonyms ({len(synonyms)} names)"):
                        st.write(", ".join(synonyms[:15]))
        else:
            st.info("No external references available (compound may not be in PubChem).")

        st.markdown("---")
        st.markdown("#### PubChem Section Links")
        if links:
            for link_name, url in links.items():
                st.markdown(f"- [{link_name}]({url})")
        else:
            st.info("No PubChem links available.")

        # PubChem data
        if pubchem and not pubchem.get("error"):
            with st.expander("PubChem Raw Data"):
                st.json(pubchem)

    # ---- TAB 6: Export ----
    with tab6:
        st.markdown('<p class="section-header">Export Data</p>', unsafe_allow_html=True)

        # Export all properties as CSV
        export_data = {}
        if cas_number:
            export_data["CAS Number"] = cas_number
        export_data.update(props)
        if mrm:
            export_data["Monoisotopic Mass"] = mrm.get("Monoisotopic Mass")
            for adduct, mz in mrm.get("Positive Mode Adducts", {}).items():
                export_data[f"MRM {adduct}"] = mz
            for adduct, mz in mrm.get("Negative Mode Adducts", {}).items():
                export_data[f"MRM {adduct}"] = mz

        export_df = pd.DataFrame([export_data])
        csv_data = export_df.to_csv(index=False)

        st.download_button(
            label="[Download] Properties as CSV",
            data=csv_data,
            file_name=f"{name.replace(' ', '_')}_properties.csv",
            mime="text/csv",
            use_container_width=True,
        )

        # Export MRM transitions
        transitions = mrm.get("Suggested MRM Transitions", []) if mrm else []
        if transitions:
            trans_df = pd.DataFrame(transitions)
            mrm_csv = trans_df.to_csv(index=False)
            st.download_button(
                label="[Download] MRM Transitions as CSV",
                data=mrm_csv,
                file_name=f"{name.replace(' ', '_')}_MRM_transitions.csv",
                mime="text/csv",
                use_container_width=True,
            )

        # Export structure image
        img = generate_structure_image(mol, size=(800, 600))
        if img:
            buf = BytesIO()
            img.save(buf, format="PNG")
            st.download_button(
                label="[Download] Structure Image (PNG)",
                data=buf.getvalue(),
                file_name=f"{name.replace(' ', '_')}_structure.png",
                mime="image/png",
                use_container_width=True,
            )


# ============================================================
# MAIN ANALYSIS LOGIC
# ============================================================
if analyze_btn:

    # ---- SINGLE COMPOUND (Name, CAS, or SMILES) ----
    if input_mode in ["Chemical Name", "SMILES String"]:
        if not compound_name and not smiles_input:
            st.warning("Please enter a compound name or SMILES string.")
        else:
            result = analyze_single_compound(
                name=compound_name if compound_name else None,
                smiles=smiles_input if smiles_input else None,
            )
            if result:
                display_compound_report(result)

    elif input_mode == "CAS Number":
        if not cas_input:
            st.warning("Please enter a CAS Registry Number.")
        else:
            result = analyze_single_compound(cas=cas_input.strip())
            if result:
                display_compound_report(result)

    # ---- BATCH PROCESSING ----
    elif input_mode == "Batch Upload (CSV/Excel)" and batch_df is not None:
        compounds = batch_df[selected_col].dropna().unique().tolist()

        if not compounds:
            st.warning("No compounds found in the selected column.")
        else:
            st.markdown(f'<p class="section-header">Batch Analysis: {len(compounds)} compounds</p>',
                        unsafe_allow_html=True)

            progress_bar = st.progress(0)
            status_text = st.empty()

            all_results = []
            all_mrm = []
            failed = []

            for i, compound in enumerate(compounds):
                status_text.text(f"Analyzing {i+1}/{len(compounds)}: {compound}")
                progress_bar.progress((i + 1) / len(compounds))

                try:
                    if input_type_col == "Compound Names":
                        result = analyze_single_compound(name=compound)
                    else:
                        result = analyze_single_compound(smiles=compound)

                    if result and result.get("properties"):
                        row = {"Compound": compound, "SMILES": result["smiles"]}
                        row.update(result["properties"])

                        # Add MRM data
                        mrm = result.get("mrm", {})
                        if mrm:
                            row["Monoisotopic Mass"] = mrm.get("Monoisotopic Mass")
                            pos_adducts = mrm.get("Positive Mode Adducts", {})
                            if "[M+H]+" in pos_adducts:
                                row["[M+H]+ (m/z)"] = pos_adducts["[M+H]+"]
                            if "[M+Na]+" in pos_adducts:
                                row["[M+Na]+ (m/z)"] = pos_adducts["[M+Na]+"]
                            neg_adducts = mrm.get("Negative Mode Adducts", {})
                            if "[M-H]-" in neg_adducts:
                                row["[M-H]- (m/z)"] = neg_adducts["[M-H]-"]

                            # Add top transitions
                            transitions = mrm.get("Suggested MRM Transitions", [])
                            for j, t in enumerate(transitions[:3]):
                                row[f"MRM Transition {j+1}"] = t.get("Transition", "")
                                row[f"CE {j+1} (eV)"] = t.get("Est. CE (eV)", "")

                        all_results.append(row)

                        # Collect MRM transitions separately
                        for t in mrm.get("Suggested MRM Transitions", []):
                            t_row = {"Compound": compound, **t}
                            all_mrm.append(t_row)

                    else:
                        failed.append(compound)

                except Exception as e:
                    failed.append(compound)

                # Rate limit PubChem requests
                time.sleep(0.3)

            progress_bar.empty()
            status_text.empty()

            # Display results
            if all_results:
                st.success(f"Successfully analyzed {len(all_results)} / {len(compounds)} compounds")

                results_df = pd.DataFrame(all_results)
                st.dataframe(results_df, use_container_width=True, hide_index=True)

                # Download buttons
                col_dl1, col_dl2 = st.columns(2)
                with col_dl1:
                    csv = results_df.to_csv(index=False)
                    st.download_button(
                        "[Download] All Properties (CSV)",
                        data=csv,
                        file_name="batch_compound_properties.csv",
                        mime="text/csv",
                        use_container_width=True,
                    )

                with col_dl2:
                    if all_mrm:
                        mrm_df = pd.DataFrame(all_mrm)
                        mrm_csv = mrm_df.to_csv(index=False)
                        st.download_button(
                            "[Download] All MRM Transitions (CSV)",
                            data=mrm_csv,
                            file_name="batch_MRM_transitions.csv",
                            mime="text/csv",
                            use_container_width=True,
                        )

                # Show individual reports (without outer expander to avoid nesting)
                st.markdown("---")
                st.markdown("### Detailed Reports")
                selected_compound = st.selectbox(
                    "Select compound for detailed view:",
                    compounds[:20],
                    index=0,
                    key="batch_detail_select",
                )
                if selected_compound:
                    if input_type_col == "Compound Names":
                        r = analyze_single_compound(name=selected_compound)
                    else:
                        r = analyze_single_compound(smiles=selected_compound)
                    if r:
                        display_compound_report(r)

            if failed:
                st.warning(f"Failed to analyze: {', '.join(failed)}")

else:
    # Welcome screen
    st.markdown("---")
    st.markdown("""
    ### Welcome! Get started:
    1. **Enter a compound name** (e.g., `Aspirin`, `Caffeine`, `Metformin`) in the sidebar
    2. Or **enter a SMILES string** (e.g., `CC(=O)OC1=CC=CC=C1C(=O)O`)
    3. Or **upload a CSV/Excel** file with multiple compounds
    4. Click **>> Analyze**

    ### Example compounds to try:
    | Name | Use Case |
    |------|----------|
    | Aspirin | Classic NSAID |
    | Caffeine | Stimulant |
    | Metformin | Diabetes drug |
    | Ibuprofen | Anti-inflammatory |
    | Paracetamol | Analgesic |
    | Quercetin | Natural flavonoid |
    | Curcumin | Turmeric compound |
    | Resveratrol | Antioxidant |
    """)
