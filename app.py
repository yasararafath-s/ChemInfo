# -*- coding: utf-8 -*-
"""
ChemInfo Dashboard - Comprehensive Chemical Compound Analyzer
=============================================================
"""

import streamlit as st
import pandas as pd
import time
from io import BytesIO

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
    get_compound_description,
    get_compound_synonyms,
    get_pubchem_links,
    get_similar_compounds,
    get_compound_by_cid,
)


# ============================================================
# PAGE CONFIGURATION
# ============================================================
st.set_page_config(
    page_title="ChemInfo Dashboard",
    page_icon="*",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ============================================================
# CUSTOM CSS
# ============================================================
st.markdown("""
<style>
    .main-header {
        font-size: 2.2rem;
        font-weight: 700;
        color: #3498db;
        text-align: center;
        margin-bottom: 0.5rem;
    }
    .sub-header {
        font-size: 1.0rem;
        opacity: 0.7;
        text-align: center;
        margin-bottom: 2rem;
    }
    .section-header {
        font-size: 1.3rem;
        font-weight: 600;
        color: #3498db;
        border-bottom: 2px solid #3498db;
        padding-bottom: 5px;
        margin-top: 1.5rem;
    }
    .pass-badge {
        background-color: #27ae60;
        color: white;
        padding: 3px 10px;
        border-radius: 12px;
        font-weight: 600;
    }
    .fail-badge {
        background-color: #e74c3c;
        color: white;
        padding: 3px 10px;
        border-radius: 12px;
        font-weight: 600;
    }
    .stTabs [data-baseweb="tab-list"] button {
        font-size: 1.05rem;
    }
    div[data-testid="stExpander"] details summary p {
        font-size: 1.1rem;
        font-weight: 600;
    }
</style>
""", unsafe_allow_html=True)


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
        ["Chemical Name", "SMILES String", "Batch Upload (CSV/Excel)"],
        index=0,
    )

    compound_name = None
    smiles_input = None
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
def analyze_single_compound(name=None, smiles=None):
    """Analyze a single compound and return all data."""
    result = {}

    # Step 1: Resolve compound via PubChem if name provided
    pubchem_data = None
    if name and not smiles:
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
        st.error("No SMILES string available. Please provide a compound name or SMILES.")
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

    # Step 4: PubChem extras
    if pubchem_data and pubchem_data.get("CID"):
        cid = pubchem_data["CID"]
        result["description"] = get_compound_description(cid)
        result["synonyms"] = get_compound_synonyms(cid)
        result["links"] = get_pubchem_links(cid)
    else:
        result["description"] = None
        result["synonyms"] = []
        result["links"] = {}

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

    # ---- Title & Structure ----
    st.markdown(f'<p class="section-header">{name}</p>', unsafe_allow_html=True)

    col_struct, col_info = st.columns([1, 1.5])

    with col_struct:
        img = generate_structure_image(mol, size=(500, 400))
        if img:
            st.image(img, caption=f"2D Structure of {name}", use_container_width=True)

        st.code(smiles, language=None)

        if pubchem and pubchem.get("CID"):
            st.markdown(f"**PubChem CID:** [{pubchem['CID']}](https://pubchem.ncbi.nlm.nih.gov/compound/{pubchem['CID']})")

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

        prop_df = pd.DataFrame(
            [{"Property": k, "Value": v} for k, v in props.items()
             if k not in ["InChI", "InChIKey"]]
        )
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

        if links:
            for link_name, url in links.items():
                st.markdown(f"- [{link_name}]({url})")
        else:
            st.info("No PubChem links available (compound may not be in PubChem).")

        # PubChem data
        if pubchem and not pubchem.get("error"):
            with st.expander("PubChem Raw Data"):
                st.json(pubchem)

    # ---- TAB 6: Export ----
    with tab6:
        st.markdown('<p class="section-header">Export Data</p>', unsafe_allow_html=True)

        # Export all properties as CSV
        export_data = {**props}
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

    # ---- SINGLE COMPOUND (Name or SMILES) ----
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

                # Show individual reports in expanders
                st.markdown("---")
                st.markdown("### Detailed Reports")
                for compound in compounds[:20]:  # Limit detailed view
                    if input_type_col == "Compound Names":
                        r = analyze_single_compound(name=compound)
                    else:
                        r = analyze_single_compound(smiles=compound)
                    if r:
                        with st.expander(f"[Detail] {compound}"):
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
