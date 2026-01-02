<div align="center">
  <table border="0" cellpadding="0" cellspacing="0" style="border: none;">
    <tr>
      <td style="border: none; vertical-align: middle;">
        <img src="assets/RinPy_128x128.png" alt="RinPy logo">
      </td>
      <td style="border: none; vertical-align: middle; padding-left: 25px;">
        <h2 style="margin-top: -12px; margin-left: -40px">
          RinPy – Residue Interaction Network for Protein Structures
        </h2>
      </td>
    </tr>
  </table>
</div>

📖 Description
------------
The **RinPy** is designed for constructing, visualizing, and analyzing Residue Interaction Networks (RINs). RIN describes a protein as a network of nodes interconnected by weighted edges. In this network, each node represents a residue, nucleotide or a ligand at the average coordinate of its atomic coordinates. The edge weight between two nodes is defined by the local interaction strength between the two residues. The average coordinates of the residues and/or nucleotides are placed at the Cα atom or P atom, respectively, for the protein-RNA/DNA complexes. Each node is annotated with attributes such as Chain ID, Residue Number, Insertion Code, and its Cartesian coordinates.

✨ Features
-------------

- Converts protein complexes to residue interaction networks (RINs).
- Supports weighted edges based on the local interaction strength or affinity.
- Annotates nodes with Chain ID, Residue Number, Insertion Code, and Cartesian coordinates.
- Can be integrated with molecular visualization system such as PyMOL.

🖥️ RinPy GUI
--------------
You can download the standalone graphical user interface (GUI) version of **RinPy** from [here](https://drive.google.com/drive/folders/1GlLva31y7Ebpmpd8Dk6uQmGHCem2vWfO?usp=drive_link).

⚙️ Installation
-----------------

### 📌 Prerequisites (Important)

- **Python ≥ 3.10** is required.
- RinPy depends on **NetworkX ≥ 3.4**, which requires Python 3.10 or newer.

To ensure reliable management of the Python environment and scientific dependencies, we strongly recommend using *
*Miniconda** or **Anaconda**.

- **Miniconda** (lightweight, recommended):  
  https://docs.conda.io/en/latest/miniconda.html

- **Anaconda** (full distribution):  
  https://www.anaconda.com/products/distribution

---

### 🚀 Installation via PyPI (Recommended)

RinPy is available on **PyPI** ([RinPy](https://pypi.org/project/rinpy/)) and can be installed directly using **pip**.

The following steps demonstrate how to create and activate a conda virtual environment, install RinPy, verify the installation, and run the program from the command line:

```bash
# Create a conda virtual environment
conda create -n rinpy_env python=3.10 -y

# Activate the environment
conda activate rinpy_env or source activate rinpy_env

# Install RinPy
pip install rinpy

# Check installation
rinpy --help

# Run RinPy
rinpy rinpy --input_path INPUT_PATH --output_path OUTPUT_PATH --calculation_option_file path/to/calculation_options.json
```

#### Parameters from the terminal

- `--input_path`: Input directory including PDB files
- `--output_path`: Output directory
- `--calculation_option_file`: JSON file containing parameters. To download, click [calculation_options.json](https://github.com/zehrasarica/rinpy/tree/main/src/rinpy/calculation_options.json).

### 🔧 Installation from Source (Alternative)

If you prefer to install RinPy from source, follow the steps below:

1. Clone the repository:
   git clone https://github.com/zehrasarica/rinpy.git

2. Navigate to the project folder:
   cd rinpy

3. Create a Python virtual environment (optional but recommended):
   ```bash
   conda create -n rinpy_env python=3.10
   source activate rinpy_env or conda activate rinpy_env

4. Install dependencies:
   pip install -r requirements.txt

🚀 Usage
---------
RinPy can be used programmatically via the `RINProcess` API within your Python scripts. Create a Python file named `main.py`, insert the content given below, and execute the script via the terminal or an equivalent environment.

### Basic Example

```python
from rinpy import RINProcess

# Define calculation options as a JSON-like dictionary
calculation_options = {
    'remove_hydrogen': {
        'is_checked': True,
        'value': 0
    },
    'betweenness': {
        'is_checked': True,
        'value': 5
    },
    'closeness': {
        'is_checked': True,
        'value': 100
    },
    'degree': {
        'is_checked': True,
        'value': 100
    },
    'cluster_number': {
        'is_checked': True,
        'value': 3
    },
    'cutoff': {
        'is_checked': True,
        'value': 4.5
    }
}

# Initialize RINProcess, 
#
# output_path is mandatory.
# For input, provide ONLY ONE of the following:
#   - input_path
#   - pdb_ids
#   - trajectory_file
# The system checks inputs in the following order:
# input_path → pdb_ids → trajectory_file
#
# If using pdb_ids, set input_path and trajectory_file to None.

rin = RINProcess(
    input_path="path/to/input/files",
    output_path="path/to/output",
    pdb_ids=None,  # list of PDB IDs to process if download from protein data bank such ["4OBE", "4DSN"].
    ligand_dict=None,  # optional ligand information
    calculation_options=calculation_options,
    trajectory_file="path/to/input/files",
    stride=1  # The default is 1. This parameter is used in conjunction with trajectory_file.
)

# Start the process
rin.start_process()
```

---------------------------------------------------------

📝 Notes
----------
**RinPy** requires output_path and only one of the following: input_path, pdb_ids, or trajectory_file. Providing multiple
inputs is not allowed; if more than one is given, input_path will take precedence. Input processing order: input_path →
pdb_ids → trajectory_file.

- **output_path**: Folder where processed RINs and results will be saved
- **input_path**: Folder containing your input PDB files
- **pdb_ids**: List of specific PDB IDs to process
- **ligand_dict**: Optional dictionary with ligand information
- **calculation_options**: JSON-like dictionary specifying which calculations to run and their parameters
- **trajectory_file**: The MD trajectory file (pdb format) which contains multiple snaphots from MD.
- **stride**: Default is **1**. This parameter is used in conjunction with trajectory_file parameter.


📄 License
------------
MIT License. See LICENSE file for details.

📘How to Cite
---------------

If you use this repository, please cite this study as follows:

```bibtex
@article{rinpy,
  title={RinPy, a Python Package for Residue Interaction Network Model to Analyze Protein Structures and Predict Ligand Binding Sites},
  author={Sarica, Z.; Sungur, F. A.; Kurkcuoglu, O.},
  journal={},
  volume={},
  year={},
  publisher={}
}
```

📬 Contact
------------
For questions, contact: zehraacar559@gmail.com, sarica16@itu.edu.tr
