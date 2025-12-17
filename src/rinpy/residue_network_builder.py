# -*- coding: utf-8 -*-
"""
    __author__ = 'Zehra Sarica'
    __credits__ = ''
    __email__ = ['sarica16@itu.edu.tr','zehraacar559@gmail.com']
"""

import logging
import multiprocessing as mp
import os
from math import trunc
from typing import cast, Tuple

import numpy as np
import pandas as pd

from rinpy import utils
from rinpy.constants import *
from rinpy.log_messages import LOG_MESSAGE
from rinpy.log_util import log_elapsed_time
from rinpy.utils import create_folder_not_exists, PdbRecord, clear_existing_file, write_ppdb_to_pdb_file, \
    create_empty_ppdb_atom_df, create_ppdb_atom_df

SOURCE = 'source'
TARGET = 'target'
AFFINITY = 'affinity'


class ResidueGraphBuilder:
    """ The Residue Graph Builder module generates output data for constructing residue interaction networks based on
        affinity calculations between α-carbon atoms. It computes atom-to-atom distances between source and target
        residues to identify potential interactions.

        Attributes
        ----------
        pdb_name : str, default: None
            PDB ID which must be same with file name, such as 3t0t (3t0t.pdb).

        pdb_path : str, default: None
            Location of the PDB file that will be used for contact calculation (path/to/3t0t.pdb).

        use_preprocess : bool, default: False
            Whether to preprocess the PDB before analysis. The given PDB file contains additional information
            such as {'ATOM', 'HETATM', 'ANISOU', 'OTHERS'}.If this is True, then only ATOM records from
            the PDB file will be extracted as preprocess and will be saved as 3t0t_pre.pdb in the given output of 3t0t.

        cutoff : float, default: 4.5 Armstrong (Å)
            The distance between atoms will be considered if it is equal or less than the given cutoff.

        destination_output_path : str, default: output, this path must be provided.
            Location of the output results in which all files for each PDB file process are stored in this path.

        het_atom_list: list, default: None
            List of hetero atoms to be included as a node in the generated graph.

        remove_hydrogen: bool, default: False
            Remove the hydrogen atoms from the given PDB file.

        num_workers: int, default: None The number of CPU cores will be used for the heavy atom-atom calculation.
        If None, available CPU cores will be used.
    """

    def __init__(self, pdb_name: str = None,
                 pdb_path: str = None,
                 use_preprocess: bool = False,
                 cutoff: float = 4.5,
                 destination_output_path: str = None,
                 het_atom_list: list = None,
                 remove_hydrogen: bool = True,
                 num_workers=None):
        if destination_output_path is None:
            raise ValueError('You must provide an output path to proceed.')
        self.pdb_name = pdb_name
        self.pdb_path = pdb_path
        self.het_atom_list = het_atom_list
        self.remove_hydrogen = remove_hydrogen
        self.ppdb = utils.get_ppdb(pdb_file_path=self.pdb_path)
        self.cutoff = cutoff
        logging.info(f"Cutoff: {self.cutoff}")
        self.num_workers = num_workers
        workers_msg = (str(self.num_workers) if self.num_workers is not None else "Detected CPU cores will be used.")
        logging.info(f"CPU workers requested: {workers_msg}")
        self.destination_output_path = destination_output_path
        self.actual_residue_number_map = dict()

        if use_preprocess:
            self._preprocess()

        if het_atom_list is not None:
            self._append_het_atoms()

        self.output_path = os.path.join(self.destination_output_path, self.pdb_name)

        create_folder_not_exists(self.output_path)

        self.residue_edges_file_path = os.path.join(self.destination_output_path, self.pdb_name,
                                                    f'{self.pdb_name}_{RESIDUE_EDGES_FILE}')

        self.residue_average_coordinates_file_path = os.path.join(self.destination_output_path, self.pdb_name,
                                                                  f'{self.pdb_name}_{RESIDUE_AVERAGE_COORDINATES_FILE}')

        self.residue_average_pdb_file_path = os.path.join(self.destination_output_path, self.pdb_name,
                                                          f'{self.pdb_name}_{RESIDUE_AVERAGE_PDB_FILE}')

    def get_number_of_amino_acid(self):
        """ Returns the total number of amino acid for the given PDB file """
        return self._get_pdb_df()[PdbRecord.ATOM.name].groupby([CHAIN_ID, RESIDUE_NUMBER, INSERTION]).ngroups

    def get_actual_residue_number_map(self):
        """ Returns the actual residue map for the given PDB file. This file contains Chain ID and residue number
        corresponding residue index.
        """
        return self.actual_residue_number_map

    def _get_pdb_df(self):
        """ Returns dictionary storing pandas DataFrames for PDB record section.
        The dictionary keys are {'ATOM', 'HETATM', 'ANISOU', 'OTHERS'}
        """
        return self.ppdb.df

    def _preprocess(self):
        """ Only gets the ATOM record section from the dictionary that stores {'ATOM', 'HETATM', 'ANISOU', 'OTHERS'}"""
        if self.remove_hydrogen:
            logging.info(f'Removing hydrogen...')
            self.ppdb.df[PdbRecord.ATOM.name] = self.ppdb.df[PdbRecord.ATOM.name][
                self.ppdb.df[PdbRecord.ATOM.name][ELEMENT_SYMBOL] != HYDROGEN_SYMBOL]

            self.ppdb.to_pdb(
                path=os.path.join(os.path.dirname(self.pdb_path), f'{self.pdb_name}{PREPROCESS_SAVE_EXT}'),
                records=[PdbRecord.ATOM.name],
                gz=False,
                append_newline=True)

            logging.info(f"{self.pdb_name}{PREPROCESS_SAVE_EXT} {LOG_MESSAGE['PREPROCESS_SAVED']}")

    def _append_het_atoms(self):
        """ Append HETERO ATOMS to the ATOM section of the pdb file which will be included in the graph generated."""
        if self.het_atom_list is not None:
            logging.info(f'Appending HETERO ATOMS ...')
            self.ppdb.df[PdbRecord.ATOM.name] = pd.concat(
                [self.ppdb.df[PdbRecord.ATOM.name], self._get_selected_het_atoms()], ignore_index=True)

            self.ppdb.to_pdb(
                path=os.path.join(os.path.dirname(self.pdb_path), f'{self.pdb_name}{WITH_HET_ATOM_SAVE_EXT}'),
                records=[PdbRecord.ATOM.name],
                gz=False,
                append_newline=True)

            logging.info(f"{self.pdb_name}{WITH_HET_ATOM_SAVE_EXT} {LOG_MESSAGE['PREPROCESS_SAVED']}")

    def clear_files(self):
        """ Clears created and stored files for the given PDB file."""
        clear_existing_file(file_path=self.residue_average_pdb_file_path)
        clear_existing_file(file_path=self.residue_average_coordinates_file_path)
        clear_existing_file(file_path=self.residue_edges_file_path)

    def _get_selected_het_atoms(self):
        """ het_atom is tuple which contains (chain_id, residue_number, residue_name, insertion) """
        matching_rows_df_list = []
        if self.het_atom_list is not None:
            het_atom_df = self._get_pdb_df()[PdbRecord.HETATM.name]
            for het_atom in self.het_atom_list:
                matched_rows = het_atom_df[(het_atom_df[CHAIN_ID] == het_atom[0]) &
                                           (het_atom_df[RESIDUE_NUMBER] == het_atom[1]) &
                                           (het_atom_df[RESIDUE_NAME] == het_atom[2]) &
                                           (het_atom_df[INSERTION] == het_atom[3])]
                if not matched_rows.empty:
                    matching_rows_df_list.append(matched_rows)

        matching_rows_df = pd.concat(matching_rows_df_list,
                                     ignore_index=True) if matching_rows_df_list else pd.DataFrame(columns=PDB_COLUMNS)
        return matching_rows_df

    def _calculate_affinity(self, index_pair, residue_number_list, residue_coord_dict):
        """ Calculates the local interaction strength between two residues based on the spatial distance
        between their atom coordinates.

        Parameters
        ----------
        index_pair : tuple[int, int]
            A tuple containing the indices (i, j) of the source and target residues within the residue_number_list.

        residue_number_list : list[tuple[str, int, str]]
            A list of residue identifiers where each item is a tuple of (chain ID, residue number, insertion code).

        residue_coord_dict : dict[tuple[str, int, str], list[dict[str, float]]]
            A dictionary mapping residue identifiers to their 3D coordinates.
            Each value is a list containing one dictionary with keys: 'x_coord', 'y_coord', and 'z_coord'.

        Returns
        -------
        dict or None
            A dictionary containing the 1-based indices of the source and target residues
            and the calculated affinity value, or None if the affinity could not be computed.

            Example:
            {
                'source': 1,
                'target': 2,
                'affinity': 0.12345678
            }
        """
        i, j = index_pair
        source_list = residue_coord_dict[residue_number_list[i]]
        target_list = residue_coord_dict[residue_number_list[j]]
        if not source_list or not target_list:
            return None
        source_coords = np.array([(source[X_COORD], source[Y_COORD], source[Z_COORD]) for source in source_list])
        target_coords = np.array([(target[X_COORD], target[Y_COORD], target[Z_COORD]) for target in target_list])
        distances = np.linalg.norm(source_coords[:, np.newaxis] - target_coords, axis=2)
        affinity_count = np.sum(distances <= self.cutoff)
        if affinity_count != 0:
            affinity = affinity_count / (np.sqrt(len(source_list)) * np.sqrt(len(target_list)))
            return {
                SOURCE: i + 1,
                TARGET: j + 1,
                AFFINITY: trunc(affinity * np.power(10, 8)) / np.power(10, 8)
            }
        else:
            return None

    @staticmethod
    def _get_default_atom_name(atom_df):
        if (atom_df[ATOM_NAME] == P_ATOM_NAME).any():
            return P_ATOM_NAME
        return CA_ATOM_NAME

    def calculate_contact(self):
        """Calculate the atom-atom distance between residues."""
        residue_average_pdb_list = []
        residue_average_coords_list = []
        atoms_df = self._get_pdb_df()[PdbRecord.ATOM.name]

        residue_chain_dict = {(str(chain_id), int(residue_number), str(insertion)): rows for
                              (chain_id, residue_number, insertion), rows in
                              atoms_df.groupby([CHAIN_ID, RESIDUE_NUMBER, INSERTION])}
        residue_chain_dict = {key: residue_chain_dict[key] for key in
                              sorted(residue_chain_dict.keys(), key=lambda x: (x[0], x[1], x[2]))}
        residue_coord_dict = {}
        residue_number_list = []

        atom_name = self._get_default_atom_name(next(iter(residue_chain_dict.values())))
        logging.info(f'Default atom name for {self.pdb_name} is {atom_name}')

        residue_index = 1
        for res_num, res_atom_df in residue_chain_dict.items():
            residue_number_list.append(res_num)
            average_values = res_atom_df[[X_COORD, Y_COORD, Z_COORD, B_FACTOR]].mean()
            x_coord_average = average_values[X_COORD]
            y_coord_average = average_values[Y_COORD]
            z_coord_average = average_values[Z_COORD]
            b_factor_average = average_values[B_FACTOR]
            temp_coord_list = []

            for _, row in res_atom_df.iterrows():
                temp_coord_list.append({X_COORD: row[X_COORD],
                                        Y_COORD: row[Y_COORD],
                                        Z_COORD: row[Z_COORD]
                                        })

            residue_coord_dict[res_num] = temp_coord_list

            chain_id, residue_number, insertion = cast(Tuple[str, int, str], res_num)

            self.actual_residue_number_map[residue_index] = tuple(
                list(res_num) + [res_atom_df[RESIDUE_NAME].iloc[-1]])
            residue_average_pdb_list.append({RECORD_NAME: PdbRecord.ATOM.name,
                                             ATOM_NUMBER: residue_index,
                                             ATOM_NAME: atom_name,
                                             RESIDUE_NAME: res_atom_df[RESIDUE_NAME].iloc[-1],
                                             CHAIN_ID: chain_id,
                                             RESIDUE_NUMBER: residue_number,
                                             INSERTION: insertion,
                                             X_COORD: x_coord_average,
                                             Y_COORD: y_coord_average,
                                             Z_COORD: z_coord_average,
                                             B_FACTOR: b_factor_average
                                             })
            residue_average_coords_list.append({RESIDUE_NUMBER: residue_index,
                                                X: THREE_DECIMAL_FLOAT_FORMATTER.format(x_coord_average),
                                                Y: THREE_DECIMAL_FLOAT_FORMATTER.format(y_coord_average),
                                                Z: THREE_DECIMAL_FLOAT_FORMATTER.format(z_coord_average)
                                                })

            residue_index += 1

        residue_average_coords_df = pd.DataFrame(residue_average_coords_list)

        residue_average_coords_df.to_csv(
            path_or_buf=str(self.residue_average_coordinates_file_path),
            sep='\t',
            index=False,
            header=False)

        write_ppdb_to_pdb_file(create_empty_ppdb_atom_df(create_ppdb_atom_df(pd.DataFrame(residue_average_pdb_list))),
                               self.residue_average_pdb_file_path)

        residue_edges_list = self._process_affinity_calculation(residue_number_list, residue_coord_dict)

        residue_edges_df = pd.DataFrame(residue_edges_list, columns=[SOURCE, TARGET, AFFINITY])

        residue_edges_df.to_csv(str(self.residue_edges_file_path),
                                sep='\t',
                                index=False,
                                header=False)

    @log_elapsed_time(message="CALCULATE AFFINITY")
    def _process_affinity_calculation(self, residue_number_list, residue_coord_dict):
        """ Calculates all affinity calculation for residues.

        Parameters
        ----------
        residue_number_list : list[tuple[str, int, str]]
            A list of residue identifiers where each item is a tuple of (chain ID, residue number, insertion code).

        residue_coord_dict : dict[tuple[str, int, str], list[dict[str, float]]]
            A dictionary mapping residue identifiers to their 3D coordinates.
            Each value is a list containing one dictionary with keys: 'x_coord', 'y_coord', and 'z_coord'.

        Returns
        -------
        list of dict
            List of a dictionary containing the 1-based indices of the source and target residues
            and the calculated affinity value.

            Example:
            [{
                'source': 1,
                'target': 2,
                'affinity': 0.12345678
            },
            {
                'source': 1,
                'target': 3,
                'affinity': 2.1558478
            },
            ]
        """

        total_residue_number = len(residue_number_list)
        index_pairs = [(i, j) for i in range(total_residue_number - 1) for j in range(i + 1, total_residue_number)]
        cpu_count = mp.cpu_count()
        if self.num_workers is not None:
            cpu_count = self.num_workers
        logging.info(f"CPU cores detected: {mp.cpu_count()} | CPU workers requested: {cpu_count}")
        with mp.Pool(cpu_count) as pool:
            results = pool.starmap(self._calculate_affinity,
                                   [(index_pair, residue_number_list, residue_coord_dict) for index_pair in
                                    index_pairs])

        residue_edges_list = [result for result in results if result is not None]
        residue_edges_list = sorted(residue_edges_list, key=lambda x: (x[SOURCE], x[TARGET]))
        return residue_edges_list


def main():
    pass


if __name__ == '__main__':
    main()
