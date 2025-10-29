# -*- coding: utf-8 -*-
"""
    __author__ = 'Zehra Sarica'
    __credits__ = ''
    __email__ = ['sarica16@itu.edu.tr','zehraacar559@gmail.com']
"""

RESIDUE_EDGES_FILE: str = "residue_edges.csv"
RESIDUE_AVERAGE_COORDINATES_FILE: str = "residue_average_coordinates.csv"
RESIDUE_AVERAGE_PDB_FILE: str = "residue_average.pdb"
CENTRALITY_CSV_TEMPLATE: str = "centrality_{type}.csv"
CENTRALITY_PDB_TEMPLATE: str = "centrality_{type}.pdb"

HIGH_PERCENTAGE_TEMPLATE: str = "{type}_high_percentage_residues.csv"
FREQUENCY_HIGH_PERCENTAGE_TEMPLATE: str = "frequency_{type}_high_percentage"

TOP_PERCENTAGE_TILE: str = "top_percentage_quantile"
SOURCE_RESIDUE: str = 'source_residue'
TARGET_RESIDUE: str = 'target_residue'
GRAPH_NAME: str = "rinpy"

VALUE: str = "value"

X: str = "x"
Y: str = "y"
Z: str = "z"

REMOVE_HYDROGEN: str = "remove_hydrogen"
IS_CHECKED: str = "is_checked"
CUTOFF: str = "cutoff"

ELEMENT_SYMBOL: str = 'element_symbol'
HYDROGEN_SYMBOL: str = 'H'
PREPROCESS_SAVE_EXT: str = "_pre.pdb"
WITH_HET_ATOM_SAVE_EXT: str = "_with_het_atom.pdb"

ATOM_NUMBER: str = "atom_number"
ATOM_NAME: str = "atom_name"
CENTRALITY_SCORE: str = "centrality_score"
B_FACTOR: str = "b_factor"
OCCUPANCY: str = 'occupancy'
RESIDUE_NUMBER: str = "residue_number"
MODE: str = "mode"
RESIDUE_NAME: str = "residue_name"
CHAIN_ID: str = "chain_id"
INSERTION: str = "insertion"
PDB_ID: str = "pdb_id"
TWO_DECIMAL_FLOAT_FORMATTER: str = '{:.2f}'
X_COORD: str = "x_coord"
Y_COORD: str = "y_coord"
Z_COORD: str = "z_coord"
RESIDUE_INDEX: str = "residue_index"
BLANK_4: str = "blank_4"

CA_ATOM_NAME: str = "CA"
P_ATOM_NAME: str = "P"
PDB_EXT: str = ".pdb"
STAR_PRINT_COUNT: int = 30

CALCULATION_OPTION_DEFAULT_JSON_PATH: str = 'calculation_options.json'

PDB_COLUMNS = ['record_name', 'atom_number', 'blank_1', 'atom_name', 'alt_loc',
               'residue_name', 'blank_2', 'chain_id', 'residue_number', 'insertion',
               'blank_3', 'x_coord', 'y_coord', 'z_coord', 'occupancy', 'b_factor',
               'blank_4', 'segment_id', 'element_symbol', 'charge', 'line_idx']
