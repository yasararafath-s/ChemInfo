# ChemInfo Dashboard

**Live App:** [https://yasararafath-s-cheminfo-app-tpzusl.streamlit.app](https://yasararafath-s-cheminfo-app-tpzusl.streamlit.app/)

A comprehensive chemical compound analyzer built with **Streamlit** and **RDKit**. Enter a compound name, SMILES string, or upload a batch file to instantly generate physicochemical properties, 2D structures, MRM mass spectrometry data, drug-likeness evaluations, and PubChem references.

## Features

- **Compound Lookup** — Search by chemical name (via PubChem API) or SMILES string
- **2D Structure Visualization** — RDKit-generated molecular structure images
- **Physicochemical Properties** — 30+ descriptors including MW, LogP, PSA, HBD, HBA, QED, and more
- **Drug-Likeness Evaluation** — Lipinski Rule of 5, Veber Rules, Ghose Filter, Egan Filter
- **MRM Mass Spectrometry Data** — Precursor ions, adducts (positive/negative mode), predicted fragment ions, suggested MRM transitions with estimated collision energies
- **Functional Group Detection** — Identifies 22+ functional groups via SMARTS patterns
- **PubChem Integration** — Descriptions, synonyms, and direct links to PubChem resources
- **Batch Processing** — Upload CSV/Excel files to analyze multiple compounds at once
- **Data Export** — Download properties, MRM transitions, and structure images

## Screenshots

After launching, the dashboard provides a sidebar for input and a main area with tabbed results:

- **Properties Tab** — Full table of physicochemical descriptors
- **Drug-Likeness Tab** — Pass/fail evaluation for multiple drug-likeness rules
- **MRM / Mass Spec Tab** — Adducts, fragments, and suggested MRM transitions
- **Functional Groups Tab** — Detected chemical functional groups
- **References Tab** — PubChem links and raw data
- **Export Tab** — Download buttons for CSV and PNG files

## Try It Online (No Installation)

### Streamlit Community Cloud (Recommended)

You can deploy and use this app directly in the browser — no installation required:

1. Go to [share.streamlit.io](https://share.streamlit.io)
2. Sign in with your **GitHub** account
3. Click **New app** and select:
   - Repository: `ChemInfo`
   - Branch: `master`
   - Main file path: `app.py`
4. Click **Deploy** — your app will be live in a few minutes

### GitHub Codespaces

Run the full development environment in your browser:

1. Go to the [ChemInfo repository](https://github.com/yasararafath-s/ChemInfo)
2. Click the green **Code** button → **Codespaces** → **Create codespace on master**
3. Once the environment loads, run in the terminal:
   ```bash
   pip install -r requirements.txt
   streamlit run app.py
   ```
4. Codespaces will show a popup to open the app in the browser

---

## Local Installation

### Prerequisites

- Python 3.10 or higher
- pip (Python package manager)

### Quick Setup (Windows)

1. Clone the repository:
   ```bash
   git clone https://github.com/yasararafath-s/ChemInfo.git
   cd ChemInfo
   ```

2. Run the installation script:
   ```bash
   INSTALL.bat
   ```

3. Launch the dashboard:
   ```bash
   START.bat
   ```

4. Open your browser at: **http://localhost:8501**

### Manual Setup (All Platforms)

1. Clone the repository:
   ```bash
   git clone https://github.com/yasararafath-s/ChemInfo.git
   cd ChemInfo
   ```

2. Create and activate a virtual environment:
   ```bash
   python -m venv venv

   # Windows
   venv\Scripts\activate

   # macOS/Linux
   source venv/bin/activate
   ```

3. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

4. Run the application:
   ```bash
   streamlit run app.py --server.port 8501
   ```

## Usage

### Single Compound Analysis

1. Select **Chemical Name** or **SMILES String** from the sidebar
2. Enter the compound name (e.g., `Aspirin`, `Caffeine`, `Metformin`) or SMILES
3. Press **Enter** or click **>> Analyze**

### Batch Analysis

1. Select **Batch Upload (CSV/Excel)** from the sidebar
2. Upload a CSV or Excel file with a column of compound names or SMILES
3. Select the appropriate column and input type
4. Click **>> Analyze**

A sample file (`sample_compounds.csv`) is included with 15 example compounds.

## Project Structure

```
ChemInfo/
├── app.py                  # Main Streamlit application
├── requirements.txt        # Python dependencies
├── packages.txt            # System dependencies (Streamlit Cloud)
├── sample_compounds.csv    # Example compounds for batch testing
├── INSTALL.bat             # Windows installation script
├── START.bat               # Windows launch script
├── .gitignore              # Git ignore rules
├── LICENSE                 # MIT License
├── README.md               # This file
├── .streamlit/
│   └── config.toml         # Streamlit theme & server config
└── utils/
    ├── __init__.py
    ├── chem_properties.py  # RDKit property calculations & drug-likeness
    ├── mrm_calculator.py   # MRM mass spectrometry calculations
    └── pubchem_api.py      # PubChem API integration
```

## Dependencies

| Package | Purpose |
|---------|---------|
| streamlit | Web dashboard framework |
| rdkit | Cheminformatics toolkit |
| pandas | Data manipulation |
| requests | HTTP requests (PubChem API) |
| Pillow | Image processing |
| openpyxl | Excel file support |

## API Usage

This application uses the **PubChem REST API** (PUG REST) for compound data retrieval. No API key is required. Please respect PubChem's usage policies and avoid excessive automated requests.

## Caution

> **This project was created using AI-assisted vibe coding.** While the code has been tested and verified to work, please use it with caution. Always validate results against trusted reference sources before using them for research, publication, or clinical purposes. The authors are not responsible for any inaccuracies in the calculated data.

## Notes

- MRM collision energies are **estimated** using linear ramp models. Always optimize experimentally for your specific instrument.
- Drug-likeness rules are guidelines, not absolute predictors of drug behavior.
- The application does not store any user data or compound information.

## License

MIT License
