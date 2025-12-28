RINPY – Residue Interaction Network for 3D Protein Structures
======================================

📖 Description
------------
The **RinPy** is designed to predict allosteric sites in 3D protein structures
by converting them into the network. In this network, each node represents a
residue, where a residue composed of multiple atoms is simplified to a single
alpha carbon (CA) atom positioned at the average coordinates of all atoms in that residue.
Each node is annotated with attributes such as Chain ID, Residue Number, Insertion Code, and its
3D coordinates (X, Y, Z). The edges between nodes are weighted based on the local
interaction strength or affinity between residue atoms.

✨ Features
-------------

- Converts 3D protein structures to residue interaction networks (RINs)
- Supports edge weighting based on the affinity
- Annotates nodes with chain, residue number, insertion code, and coordinates
- Can be integrated with visualization tools such as PyMOL

🖥️ RINPY GUI
--------------
You can download the standalone graphical user interface (GUI) version of **RinPy** from [here](YOUR_GOOGLE_DRIVE_LINK).

⚙️ Installation
-----------------
### 📌 Prerequisites (Important)

- **Python ≥ 3.10** is required.
- RinPy depends on **NetworkX ≥ 3.4**, which requires Python 3.10 or newer.

To ensure reliable management of the Python environment and scientific dependencies, we strongly recommend using **Miniconda** or **Anaconda**.

- **Miniconda** (lightweight, recommended):  
  https://docs.conda.io/en/latest/miniconda.html

- **Anaconda** (full distribution):  
  https://www.anaconda.com/products/distribution

---

### 🚀 Installation via PyPI (Recommended)

RinPy is available on **PyPI** and can be installed directly using `pip`:

```bash
pip install rinpy
```

### 🔧 Installation from Source (Alternative)
If you prefer to install RinPy from source, follow the steps below:
1. Clone the repository:
   git clone https://github.com/zehrasarica/rinpy.git

2. Navigate to the project folder:
   cd rinpy

3. Create a Python virtual environment (optional but recommended):
   ```bash
   conda create -n rinpy python=3.10
   source activate rinpy or conda activate rinpy

4. Install dependencies:
   pip install -r requirements.txt

🚀 Usage
---------
RinPy can be used programmatically via the `RINProcess` API within your Python scripts.

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
        'value': 10
    },
    'closeness': {
        'is_checked': False,
        'value': 100
    },
    'degree': {
        'is_checked': False,
        'value': 100
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
    pdb_ids=None,  # list of PDB IDs to process if download from protein data bank such ["1ABC", "2DEF"].
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
RinPy requires output_path and only one of the following: input_path, pdb_ids, or trajectory_file. Providing multiple
inputs is not allowed; if more than one is given, input_path will take precedence. Input processing order: input_path →
pdb_ids → trajectory_file.

- **output_path**: Folder where processed RINs and results will be saved
- **input_path**: Folder containing your input PDB files
- **pdb_ids**: List of specific PDB IDs to process
- **ligand_dict**: Optional dictionary with ligand information
- **calculation_options**: JSON-like dictionary specifying which calculations to run and their parameters
- **trajectory_file**: The MD trajectory file (pdb format) which contains multiple snaphots from MD.
- **stride**: Default is **1**. This parameter is used in conjunction with trajectory_file parameter.

🤝 Contributing
-----------------
Contributions, bug reports, and feature requests are welcome! Please follow
standard GitHub workflow: fork → branch → pull request.

🪪 License
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
