# -*- coding: utf-8 -*-
"""
    __author__ = 'Zehra Sarica'
    __credits__ = ''
    __email__ = ['sarica16@itu.edu.tr','zehraacar559@gmail.com']
"""

import logging
import os
import re
import shutil
from enum import Enum

import pandas as pd
from biopandas.pdb import PandasPdb

from rinpy.constants import *


class PdbRecord(Enum):
    ATOM = 'ATOM',
    HETATM = 'HETATM',
    ANISOU = 'ANISOU',
    OTHERS = 'OTHERS'


class CentralityType(Enum):
    BET = 'BET',
    CLOS = 'CLOS',
    DEG = 'DEG'

    @property
    def display_name(self):
        names = {
            CentralityType.BET: "betweenness",
            CentralityType.CLOS: "closeness",
            CentralityType.DEG: "degree"
        }
        return names[self]

    def display_name_capitalized(self):
        return self.display_name.capitalize()


def create_empty_ppdb_atom_df(atom_df):
    ppdb = PandasPdb()
    ppdb.df[PdbRecord.ATOM.name] = atom_df
    return ppdb


def create_ppdb_atom_df(atom_df):
    df = pd.DataFrame(columns=PDB_COLUMNS)
    for idx, row in atom_df.iterrows():
        if 'record_name' in atom_df.columns:
            df.at[idx, 'record_name'] = row['record_name']
        else:
            df.at[idx, 'record_name'] = ''
        if 'atom_number' in atom_df.columns:
            df.at[idx, 'atom_number'] = row['atom_number']
        else:
            df.at[idx, 'atom_number'] = 0
        if 'blank_1' in atom_df.columns:
            df.at[idx, 'blank_1'] = row['blank_1']
        else:
            df.at[idx, 'blank_1'] = ''
        if 'atom_name' in atom_df.columns:
            df.at[idx, 'atom_name'] = row['atom_name']
        else:
            df.at[idx, 'atom_name'] = ''
        if 'alt_loc' in atom_df.columns:
            df.at[idx, 'alt_loc'] = row['alt_loc']
        else:
            df.at[idx, 'alt_loc'] = ''
        if 'residue_name' in atom_df.columns:
            df.at[idx, 'residue_name'] = row['residue_name']
        else:
            df.at[idx, 'residue_name'] = ''
        if 'blank_2' in atom_df.columns:
            df.at[idx, 'blank_2'] = row['blank_2']
        else:
            df.at[idx, 'blank_2'] = ''
        if 'chain_id' in atom_df.columns:
            df.at[idx, 'chain_id'] = row['chain_id']
        else:
            df.at[idx, 'chain_id'] = ''
        if 'residue_number' in atom_df.columns:
            df.at[idx, 'residue_number'] = row['residue_number']
        else:
            df.at[idx, 'residue_number'] = 0
        if 'insertion' in atom_df.columns:
            df.at[idx, 'insertion'] = row['insertion']
        else:
            df.at[idx, 'insertion'] = ''
        if 'blank_3' in atom_df.columns:
            df.at[idx, 'blank_3'] = row['blank_3']
        else:
            df.at[idx, 'blank_3'] = ''
        if 'x_coord' in atom_df.columns:
            df.at[idx, 'x_coord'] = row['x_coord']
        if 'y_coord' in atom_df.columns:
            df.at[idx, 'y_coord'] = row['y_coord']
        if 'z_coord' in atom_df.columns:
            df.at[idx, 'z_coord'] = row['z_coord']
        if 'occupancy' in atom_df.columns:
            df.at[idx, 'occupancy'] = row['occupancy']
        else:
            df.at[idx, 'occupancy'] = 0
        if 'b_factor' in atom_df.columns:
            df.at[idx, 'b_factor'] = row['b_factor']
        if 'blank_4' in atom_df.columns:
            df.at[idx, 'blank_4'] = row['blank_4']
        else:
            df.at[idx, 'blank_4'] = ''
        if 'segment_id' in atom_df.columns:
            df.at[idx, 'segment_id'] = row['segment_id']
        else:
            df.at[idx, 'segment_id'] = ''
        if 'element_symbol' in atom_df.columns:
            df.at[idx, 'element_symbol'] = row['element_symbol']
        else:
            df.at[idx, 'element_symbol'] = ''
        if 'charge' in atom_df.columns:
            df.at[idx, 'charge'] = row['charge']
    return df


def create_folder_not_exists(file_path=None):
    if not os.path.exists(file_path):
        os.makedirs(file_path)
        logging.info("The '{}' directory is created!".format(file_path))


def convert_pdb_to_pandas_pdb(pdb_path=None):
    ppdb = PandasPdb()
    ppdb.read_pdb(pdb_path)
    return ppdb


def clear_existing_file(file_path=None):
    if os.path.exists(file_path):
        with open(file_path, "r+") as file_object:
            file_object.truncate()
        file_object.close()


def read_text_file_to_list(input_file):
    df_list = []
    with open(input_file) as f:
        for line in f:
            line = line.strip()
            columns = re.split('\s+', line, maxsplit=4)
            df_list.append(columns)
    return df_list


def write_ppdb_to_pdb_file(ppdb, out_filename):
    ppdb.to_pdb(path=out_filename,
                records=[PdbRecord.ATOM.name],
                gz=False,
                append_newline=True)


def write_ppdb_to_pdb_file_atom_and_hetatom(ppdb, out_filename):
    ppdb.to_pdb(path=out_filename,
                records=[PdbRecord.ATOM.name, PdbRecord.HETATM.name],
                gz=False,
                append_newline=True)


def write_ppdb_to_pdb_file_all(ppdb, out_filename):
    ppdb.to_pdb(path=out_filename,
                records=None,
                gz=False,
                append_newline=True)


def convert_pdb_df_to_atom_ppdb(pdb_df):
    ppdb = PandasPdb()
    ppdb.df[PdbRecord.ATOM.name] = pdb_df
    return ppdb


def write_centrality_scores_to_file(filename_path=None, score_dict=None, formatter='%d %9.7f\n'):
    if score_dict is None:
        score_dict = dict()
    score_dict = dict(sorted(score_dict.items()))
    with open(filename_path, 'w') as f:
        for key, value in score_dict.items():
            f.write(formatter % (key, value))


def delete_folders_in_path(path):
    """ Deletes only folders and their contents for a given path.
    If the path is a file, no deletion occurs.

    Parameters:
        path (str): The path to check for folders to delete.
    """
    if not os.path.exists(path):
        logging.error(f"The given path does not exist: {path}")
        return

    if os.path.isdir(path):
        for item in os.listdir(path):
            item_path = os.path.join(path, item)
            if os.path.isdir(item_path):
                shutil.rmtree(item_path)


def get_df(output_file_path=None, columns: list = None, dtypes: dict = None,
           sort_keys: list = None) -> pd.DataFrame:
    """ Convert the given text file into a pandas dataframe and return it.

    Parameters:
    -----------
    output_file_path : str, optional
        The file path of the output file to be converted into dataframe.

    columns : list of str, optional
        A list specifying the column names to be used in the resulting DataFrame.

    dtypes : dict, optional
        A dictionary specifying data types for the columns.
    sort_keys : list of str, optional
        A list of column names to sort the DataFrame by.
    Returns:
    --------
    pd.DataFrame
        The resulting DataFrame loaded from the given file.
    """
    if len(columns) != len(dtypes):
        raise ValueError("The number of columns must match the number of dtype specifications.")

    missing_columns = [col for col in columns if col not in dtypes]

    if missing_columns:
        raise ValueError(
            f"The following columns are not present in the dtype dictionary: {', '.join(missing_columns)}")

    df = pd.DataFrame(read_text_file_to_list(output_file_path), columns=columns)

    for col, dtype in dtypes.items():
        if col in df.columns:
            df[col] = df[col].astype(dtype)

    if sort_keys is not None:
        df = df.sort_values(by=sort_keys)

    return df


def get_atom_pdb_df(ppdb=None):
    """ return the dataframe for ATOM entries in the PPDB file

    Parameters:
    -----------
    ppdb : Biopandas file
        if None, raise an error
    """
    if ppdb is None:
        raise ValueError("ppdb must not be None. Please provide a valid biopandas file.")

    df = ppdb.df[PdbRecord.ATOM.name]
    df[ATOM_NUMBER] = df[ATOM_NUMBER].astype(int)
    df = df.sort_values(by=[ATOM_NUMBER], ascending=True)
    return df


def get_ppdb(pdb_file_path):
    """Converts the given PDB file to the biopandas PandasPdb object"""
    return convert_pdb_to_pandas_pdb(pdb_file_path)


def get_residue_id(node):
    res = str(node[RESIDUE_NUMBER])
    ins = node.get(INSERTION, '')
    return f"{res}{ins}" if ins not in [None, '', "''"] else res


def get_residue_id_by_tuple(residue_tuple: tuple):
    """ residue_tuple: tuple like ('A', 270, '', 'A') or ('A', 270, 'A', 'G')
    returns: str like '270', '270A', '270B' etc.
    """
    res_num = residue_tuple[1]
    ins_code = residue_tuple[2] if len(residue_tuple) > 2 else ''
    if ins_code and ins_code not in ["", "''", None]:
        return f"{res_num}{ins_code}"
    else:
        return str(res_num)


def get_base_pdb_df(residue_to, destination_output_path, pdb_name):
    base_pdb_df = convert_pdb_to_pandas_pdb(os.path.join(destination_output_path, pdb_name, f'{pdb_name}.pdb'))
    atom_df = base_pdb_df.df[PdbRecord.ATOM.name]
    hetatm_df = base_pdb_df.df[PdbRecord.HETATM.name]

    atom_keys = list(zip(atom_df[RESIDUE_NUMBER], atom_df[CHAIN_ID], atom_df[INSERTION]))
    hetatm_keys = list(zip(hetatm_df[RESIDUE_NUMBER], hetatm_df[CHAIN_ID], hetatm_df[INSERTION]))

    base_pdb_df.df[PdbRecord.ATOM.name][B_FACTOR] = [residue_to.get(key, 0.0) for key in atom_keys]
    base_pdb_df.df[PdbRecord.HETATM.name][B_FACTOR] = [residue_to.get(key, 0.0) for key in hetatm_keys]

    return base_pdb_df
