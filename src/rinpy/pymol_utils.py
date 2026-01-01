# -*- coding: utf-8 -*-
"""
    __author__ = 'Zehra Sarica'
    __email__ = ['sarica16@itu.edu.tr','zehraacar559@gmail.com']
"""

import logging
import os

from rinpy.constants import CHAIN_ID, RESIDUE_NUMBER, INSERTION
from rinpy.utils import CentralityType


class PymolUtils:
    def __init__(self, pdb_name, destination_output_path, actual_residue_number_map):
        self.pdb_name = pdb_name
        self.destination_output_path = destination_output_path
        self.actual_residue_number_map = actual_residue_number_map

    def get_full_save_path(self, filename: str, extension: str = "png"):
        return os.path.join(self.destination_output_path, self.pdb_name, f"{self.pdb_name}_{filename}.{extension}")

    def generate_pymol_query_string(self, out_df, centrality_type: CentralityType):

        full_path = self.get_full_save_path(filename=f'{centrality_type.display_name}_pymol_queries', extension="txt")

        query_str_df = out_df.groupby(CHAIN_ID).apply(
            lambda df: list(zip(df[RESIDUE_NUMBER], df[INSERTION]))).reset_index(name="residues")

        query_strings = []

        for _, row in query_str_df.iterrows():
            residues = row["residues"]
            formatted_residues = '+'.join(
                f"{res_num}{ins}" if ins not in [None, '', "''"] else str(res_num)
                for res_num, ins in residues
            )
            query = f"select {self.pdb_name} and chain {row[CHAIN_ID]} and resi {formatted_residues}"
            query_strings.append(query)

        with open(full_path, "w") as file:
            for query in query_strings:
                file.write(query + "\n")

        logging.info(f'PyMOL query string has been saved to {full_path}.')

    def export_pymol_script(self, full_path_to_pdb, residues, filename=None, is_actual_residue=False, full_path=None):
        """ Generate a PyMOL script to visualize residues.

        Parameters:
        - residues: List of residue numbers (int)
        """
        if full_path is None:
            full_path = self.get_full_save_path(filename=filename, extension="pml")

        with open(full_path, 'w') as file:
            file.write(f"load ./{os.path.basename(full_path_to_pdb)}\n")
            file.write("hide everything\n")
            file.write("show cartoon\n")
            file.write("color gray80\n\n")
            if is_actual_residue:
                selected_residues = "+".join(str(residue_index) for residue_index in residues)
            else:
                selected_residues = "+".join(
                    f"{res_info[1]}{res_info[2]}" if res_info[2] not in [None, '', "''"] else str(res_info[1])
                    for residue_index in residues
                    if (res_info := self.actual_residue_number_map[residue_index])
                )

            file.write(f"select residues, resi {selected_residues}\n")
            file.write("show cartoon, residues\n")
            file.write("color warmpink, residues\n")
            file.write("set stick_radius, 0.25, residues\n\n")

            file.write("set label_size, 16\n")
            file.write("set label_font_id, 7\n")
            file.write("label residues and name CA, resn + resi\n")
            file.write("set label_position, [0.0, -0.3, -0.5], residues\n")
            file.write("set label_color, white\n\n")

            file.write("\nzoom residues, 10\n")
            file.write("deselect\n")

    def export_pymol_script_hinge(self, full_path_to_pdb, residues, filename=None, is_actual_residue=False,
                                  full_path=None, high_percentage_residues=None):
        """ Generate a PyMOL script to visualize residues.

        Parameters:
        - residues: List of residue numbers (int)
        """
        if full_path is None:
            full_path = self.get_full_save_path(filename=filename, extension="pml")

        with open(full_path, 'w') as file:
            file.write(f"load ./{os.path.basename(full_path_to_pdb)}\n")
            file.write("hide everything\n")
            file.write("show cartoon\n")
            file.write("color gray80\n\n")
            if is_actual_residue:
                selected_residues = "+".join(str(residue_index) for residue_index in residues)
            else:
                selected_residues = "+".join(
                    f"{res_info[1]}{res_info[2]}" if res_info[2] not in [None, '', "''"] else str(res_info[1])
                    for residue_index in residues
                    if (res_info := self.actual_residue_number_map[residue_index])
                )

            file.write(f"select hinge_residues, resi {selected_residues}\n")
            file.write("show cartoon, hinge_residues\n")
            file.write("color yellow, hinge_residues\n")
            file.write("set stick_radius, 0.25, hinge_residues\n\n")

            if high_percentage_residues is not None:
                selected_high_percentage_residues = "+".join(
                    f"{res_num}{ins}" if ins not in [None, '', "''"] else f"{res_num}"
                    for res_num, ins in high_percentage_residues
                )
                file.write(f"select hub_residues, resi {selected_high_percentage_residues}\n")
                file.write("show cartoon, hub_residues\n")
                file.write("color warmpink, hub_residues\n")
                file.write("set stick_radius, 0.25, hub_residues\n\n")

            file.write("set label_size, 16\n")
            file.write("set label_font_id, 7\n")

            file.write("label hinge_residues and name CA, resn + resi\n")
            file.write("set label_position, [0.0, -0.3, -0.5], hinge_residues\n")

            if high_percentage_residues is not None:
                file.write("label hub_residues and name CA, resn + resi\n")
                file.write("set label_position, [0.0, -0.3, -0.5], hub_residues\n")

            file.write("set label_color, white\n\n")

            file.write("\nzoom hinge_residues or hub_residues, 10\n")
            file.write("deselect\n")


def main():
    pass


if __name__ == '__main__':
    main()
