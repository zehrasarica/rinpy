# -*- coding: utf-8 -*-
"""
    __author__ = 'Zehra Sarica'
    __credits__ = ''
    __email__ = ['sarica16@itu.edu.tr','zehraacar559@gmail.com']
"""

import logging
import multiprocessing as mp
import os

import networkx as nx
import numpy as np

from rinpy import utils
from rinpy.constants import CENTRALITY_PDB_TEMPLATE, HIGH_PERCENTAGE_TEMPLATE, B_FACTOR, \
    RESIDUE_NUMBER, VALUE, CENTRALITY_CSV_TEMPLATE, RESIDUE_EDGES_FILE, RESIDUE_AVERAGE_COORDINATES_FILE, IS_CHECKED, \
    RESIDUE_NAME, CHAIN_ID, INSERTION, X_COORD, Z_COORD, Y_COORD, RESIDUE_INDEX, TOP_PERCENTAGE_TILE, SOURCE_RESIDUE, \
    TARGET_RESIDUE, X, Y, Z, GRAPH_NAME
from rinpy.graph_plotter import GraphPlotter
from rinpy.log_util import log_details
from rinpy.pymol_utils import PymolUtils
from rinpy.utils import CentralityType

CENTRALITY_TYPE_KEY = 'centrality_type'
QUANTILE_PERCENTAGE = "quantile_percentage"
OUTFILE = "outfile"
PDB_FILE = "pdb_file"


class CentralityAnalyzer:
    """ The Centrality Analyzer module generates output data for constructing residue interaction networks based on
        affinity calculations between α-carbon atoms. It computes atom-to-atom distances between source and target
        residues to identify potential contacts.

        Attributes
        ----------
        pdb_name : str, default: None
            PDB ID which must be same with file name, such as 3t0t (3t0t.pdb).

        number_of_amino_acid : int, default: None
            The number of residue in the given PDB file.

        destination_output_path : str, default: output, this path must be provided.
            Location of the output results in which all files for each PDB file process are stored in this path.

        calculation_options: dict, default: None
            The calculation options as the input for the whole process.
            Example:
            {
                'remove_hydrogen': {
                    'is_checked': True,
                    'value': 0
                },
                'betweenness': {
                    'is_checked': True,
                    'value': 5
                },
                'closeness': {
                    'is_checked': False,
                    'value': 10
                },
                'degree': {
                    'is_checked': False,
                    'value': 10
                },
                'cluster_number': {
                    'is_checked': None,
                    'value': 4
                },
                'cutoff': {
                    'is_checked': True,
                    'value': 4.5
                },
            }
        actual_residue_number_map: dict, default: None
        """

    def __init__(self, pdb_name, number_of_amino_acid, destination_output_path=None,
                 calculation_options=None, actual_residue_number_map=None):
        if destination_output_path is None:
            raise ValueError('You must provide an output path to proceed.')
        self.pdb_name = pdb_name
        self.calculation_options = calculation_options
        self.number_of_amino_acid = number_of_amino_acid
        logging.info(f"Number of amino acid: {self.number_of_amino_acid}")
        self.destination_output_path = destination_output_path
        utils.create_folder_not_exists(os.path.join(self.destination_output_path))
        self.actual_residue_number_map = actual_residue_number_map
        self.edges = None
        self.weights = None
        self.number_of_nodes = 0
        self.graph = None
        self.graph_plotter = GraphPlotter(pdb_name=self.pdb_name,
                                          destination_output_path=self.destination_output_path,
                                          actual_residue_number_map=self.actual_residue_number_map)

        self.pymol_utils = PymolUtils(pdb_name=pdb_name,
                                      destination_output_path=destination_output_path,
                                      actual_residue_number_map=actual_residue_number_map)

    def _compute_top_percentage_quantile(self, centrality_type: CentralityType = CentralityType.DEG,
                                         centrality_pdb_file: str = CENTRALITY_PDB_TEMPLATE.format(
                                             type=CentralityType.BET.display_name),
                                         high_percentage_file: str = HIGH_PERCENTAGE_TEMPLATE.format(
                                             type=CentralityType.BET.display_name),
                                         quantile_percentage=5):
        """ Creates quantile files ({pdb_id}_betweenness_high_percentage_residues.csv,
        {pdb_id}_closeness_high_percentage_residues.csv, and {pdb_id}_degree_high_percentage_residues.csv)
        for the given input file ({pdb_id}_centrality_betweenness.csv or {pdb_id}_centrality_closeness.csv or
        {pdb_id}_centrality_degree.csv) and its corresponding hub file with the given quantile percentage.
        Also, it generates PyMOL query strings to visualizes hub residues on the protein 3D structure.

        Parameters:
        -----------
        centrality_type : CentralityType, default CentralityType.DEG
            The centrality type for centrality analysis.
            strength.
        centrality_pdb_file : str, default centrality_betweennes.pdb
            The PDB file contains the calculated centrality score (betweenness, closeness, and degree) in
            the b-factor column of the PDB.
        high_percentage_file : str, default betweennes_high_percentage_residues.csv
            The file will store the top quantile scores of the given input.
        quantile_percentage : int, default 5
            The threshold value to obtain the top quantile of the given percentage.
        """

        log_detail_name = f'{centrality_type.display_name} - {self.pdb_name}'
        logging.info(f'Given quantile percentage for {log_detail_name} : {quantile_percentage}')
        centrality_pdb_file_path = os.path.join(os.path.join(str(self.destination_output_path), self.pdb_name),
                                                centrality_pdb_file)
        hub_df = utils.get_atom_pdb_df(ppdb=utils.convert_pdb_to_pandas_pdb(centrality_pdb_file_path))
        b_factors = hub_df[B_FACTOR].to_numpy(dtype=np.float64)
        threshold = np.percentile(b_factors, 100 - quantile_percentage)
        logging.info(f'Calculated threshold value for {log_detail_name} : {threshold}')
        out_df = hub_df[hub_df[B_FACTOR] >= threshold]
        out_path = os.path.join(self.destination_output_path, self.pdb_name, high_percentage_file)
        out_df.loc[:, INSERTION] = out_df[INSERTION].apply(lambda x: "''" if x == "" else x)

        out_df.to_csv(path_or_buf=out_path,
                      sep=';',
                      index=False,
                      header=False,
                      float_format="%.5f",
                      columns=[RESIDUE_NAME, CHAIN_ID, RESIDUE_NUMBER, INSERTION, B_FACTOR])

        self._create_pymol_data(out_df=out_df, centrality_type=centrality_type)

    def _create_pymol_data(self, out_df, centrality_type: CentralityType):
        """ Creates pymol query strings to be saved in a txt file for later use and also create pymol command to visualize
        residues with top quantile percentage centrality scores
        Parameters:
        -----------
        out_df: pd.DataFrame
            Data frame containing the data to be saved
        centrality_type: CentralityType
            can be the following: CentralityType.BET, CentralityType.CLOS, and CentralityType.DEG
        """
        self.pymol_utils.generate_pymol_query_string(out_df=out_df, centrality_type=centrality_type)

        original_pdb = f'{self.pdb_name}_original_{CENTRALITY_PDB_TEMPLATE.format(type=centrality_type.display_name)}'

        residue_list = [f"{res_num}{ins}" if ins not in [None, '', "''"] else str(res_num)
                        for res_num, ins in zip(out_df[RESIDUE_NUMBER], out_df[INSERTION])
                        ]

        self.pymol_utils.export_pymol_script(
            full_path_to_pdb=os.path.join(self.destination_output_path, self.pdb_name, original_pdb),
            residues=residue_list,
            filename=f'original_{centrality_type.display_name}_{TOP_PERCENTAGE_TILE}',
            is_actual_residue=True)

        self.pymol_utils.export_pymol_script(
            full_path_to_pdb=os.path.join(self.destination_output_path, self.pdb_name, f'{self.pdb_name}.pdb'),
            residues=residue_list,
            filename=f'{centrality_type.display_name}_{TOP_PERCENTAGE_TILE}',
            is_actual_residue=True)

        self.pymol_utils.export_pymol_script(
            full_path_to_pdb=os.path.join(self.destination_output_path, self.pdb_name, f'{self.pdb_name}_pre.pdb'),
            residues=residue_list,
            filename=f'{centrality_type.display_name}_{TOP_PERCENTAGE_TILE}_pre',
            is_actual_residue=True)

    def _generate_quantile_setting_dict(self):
        """ Generate a quantile setting dictionary for quantile calculation.

        Returns
        -------
        Settings for each centrality type
        """
        return {
            CentralityType.BET: {
                PDB_FILE: self._get_pdb_file(centrality_type=CentralityType.BET),
                OUTFILE: self._get_outfile(centrality_type=CentralityType.BET),
                QUANTILE_PERCENTAGE: self.calculation_options[CentralityType.BET.display_name][VALUE]
            },
            CentralityType.CLOS: {
                PDB_FILE: self._get_pdb_file(centrality_type=CentralityType.CLOS),
                OUTFILE: self._get_outfile(centrality_type=CentralityType.CLOS),
                QUANTILE_PERCENTAGE: self.calculation_options[CentralityType.CLOS.display_name][VALUE]
            },
            CentralityType.DEG: {
                PDB_FILE: self._get_pdb_file(centrality_type=CentralityType.DEG),
                OUTFILE: self._get_outfile(centrality_type=CentralityType.DEG),
                QUANTILE_PERCENTAGE: self.calculation_options[CentralityType.DEG.display_name][VALUE]
            }
        }

    def _get_pdb_file(self, centrality_type: CentralityType) -> str:
        """returns formatted pdb file name based on centrality type"""
        return f'{self.pdb_name}_{CENTRALITY_PDB_TEMPLATE.format(type=centrality_type.display_name)}'

    def _get_outfile(self, centrality_type: CentralityType) -> str:
        """returns formatted outfile name for high percentage based on centrality type"""
        return f'{self.pdb_name}_{HIGH_PERCENTAGE_TEMPLATE.format(type=centrality_type.display_name)}'

    def _get_cvs_filename(self, centrality_type: CentralityType) -> str:
        """returns formatted CSV filename based on centrality type"""
        return f'{self.pdb_name}_{CENTRALITY_CSV_TEMPLATE.format(type=centrality_type.display_name)}'

    def _calculate_quantile(self, centrality_type: CentralityType):
        """ Calculates the quantile value for the given centrality_type

        Parameters:
        -----------
        centrality_type : CentralityType
            centrality_type can be CentralityType.CLOS or CentralityType.DEG or CentralityType.BET
        """
        quantile_setting_dict = self._generate_quantile_setting_dict()
        self._compute_top_percentage_quantile(
            centrality_type=centrality_type,
            centrality_pdb_file=quantile_setting_dict[centrality_type][PDB_FILE],
            high_percentage_file=quantile_setting_dict[centrality_type][OUTFILE],
            quantile_percentage=quantile_setting_dict[centrality_type][QUANTILE_PERCENTAGE])

    def _get_input_files(self) -> list:
        input_files = []
        if self._is_bet_checked():
            input_files.append(CentralityType.BET)

        if self._is_clos_checked():
            input_files.append(CentralityType.CLOS)

        if self._is_deg_checked():
            input_files.append(CentralityType.DEG)

        return input_files

    def calculate_quantiles(self, use_parallel: bool = False):
        """ Calculates the quantile values for the three centrality types

        Parameters:
        -----------
        use_parallel : bool, default False
            if it is true, it uses multiprocessing to calculate the quantiles for each centrality type.
            Otherwise, it calculates the quantiles for each centrality sequentially.
            The quantile results for each centrality type are saved in files.
        """
        if use_parallel:
            input_files = self._get_input_files()
            with mp.Pool() as pool:
                pool.map(self._calculate_quantile, input_files)
        else:
            if self._is_bet_checked():
                self._compute_top_percentage_quantile(
                    centrality_type=CentralityType.BET,
                    centrality_pdb_file=self._get_pdb_file(centrality_type=CentralityType.BET),
                    high_percentage_file=self._get_outfile(centrality_type=CentralityType.BET),
                    quantile_percentage=self.calculation_options[CentralityType.BET.display_name][VALUE])

            if self._is_clos_checked():
                self._compute_top_percentage_quantile(
                    centrality_type=CentralityType.CLOS,
                    centrality_pdb_file=self._get_pdb_file(centrality_type=CentralityType.CLOS),
                    high_percentage_file=self._get_outfile(centrality_type=CentralityType.CLOS),
                    quantile_percentage=self.calculation_options[CentralityType.CLOS.display_name][VALUE])

            if self._is_deg_checked():
                self._compute_top_percentage_quantile(
                    centrality_type=CentralityType.DEG,
                    centrality_pdb_file=self._get_pdb_file(centrality_type=CentralityType.DEG),
                    high_percentage_file=self._get_outfile(centrality_type=CentralityType.DEG),
                    quantile_percentage=self.calculation_options[CentralityType.DEG.display_name][VALUE])

    def _is_deg_checked(self) -> bool:
        """ returns true if the centrality type degree is checked"""
        return self.calculation_options[CentralityType.DEG.display_name][IS_CHECKED]

    def _is_clos_checked(self) -> bool:
        """ returns true if the centrality type closeness is checked"""
        return self.calculation_options[CentralityType.CLOS.display_name][IS_CHECKED]

    def _is_bet_checked(self) -> bool:
        """ returns true if the centrality type betweenness is checked"""
        return self.calculation_options[CentralityType.BET.display_name][IS_CHECKED]

    def calculate_all_scores(self, use_parallel=True):
        """ Calculates the centrality scores for betweenness, closeness, and degree of the generated graph based on
        the local interaction strength.

        Parameters:
        -----------
        use_parallel : bool, default False
            if it is true, it uses multiprocessing to calculate the centrality scores
            for betweenness, closeness, and degree. Otherwise, it calculates the centrality score for each sequentially.
        """
        if use_parallel:
            centrality_types = self._get_input_files()
            with mp.Pool() as pool:
                pool.map(self._calculate_centrality_score, centrality_types)
        else:
            if self._is_bet_checked():
                self._calculate_centrality_score(centrality_type=CentralityType.BET)
            if self._is_clos_checked():
                self._calculate_centrality_score(centrality_type=CentralityType.CLOS)
            if self._is_deg_checked():
                self._calculate_centrality_score(centrality_type=CentralityType.DEG)

    @staticmethod
    def _centrality_message(*args, **kwargs):
        """ Generate the log message for _calculate_centrality_score method of the CentralityCalculation class based on
        given args and kwargs.
        Parameters:
        -----------
        *args : tuple
        **kwargs : dict

        Return:
        -----------
        The generated log message based on the centrality type.
        """
        centrality_type = kwargs.get(CENTRALITY_TYPE_KEY, None)
        if centrality_type is None and len(args) > 1:
            centrality_type = args[1]
        if centrality_type is None:
            centrality_type = CentralityType.BET
        if CentralityType.BET == centrality_type:
            centrality_type_name = CentralityType.BET.display_name_capitalized()
        elif CentralityType.CLOS == centrality_type:
            centrality_type_name = CentralityType.CLOS.display_name_capitalized()
        else:
            centrality_type_name = CentralityType.DEG.display_name_capitalized()

        return f"Centrality score calculation for {centrality_type_name}"

    @log_details(_centrality_message)
    def _calculate_centrality_score(self, centrality_type=CentralityType.BET):
        r""" Calculate the centrality score for the given centrality type.

        Parameters:
        -----------
        centrality_type : CentralityType, default CentralityType.BET
        The centrality type to calculate the score from the generated graph.
        Also, necessary figures are generated from the centrality score.
        These figures are:
            2D plot of centrality score on the graph
            2D plot of centrality score for each residue index
            Histogram of centrality score for each centrality type
        In addition, the degree scores are plotted on the additional calculation like laplacian.
        Check out the details of plot_graph_2d_degree and plot_graph_2d_degree_with_coords in GraphPlotter class.
        """
        centrality_root_path = os.path.join(str(self.destination_output_path), str(self.pdb_name))
        utils.create_folder_not_exists(centrality_root_path)
        if CentralityType.BET == centrality_type:
            filename = self._get_cvs_filename(centrality_type=centrality_type)
            filename_path = os.path.join(centrality_root_path, filename)
            scores_dict = nx.betweenness_centrality(G=self.graph, weight="weight")
            utils.write_centrality_scores_to_file(
                filename_path=filename_path,
                score_dict=scores_dict)

            logging.info(f"{filename_path} has been saved to the disk!")

            self.graph_plotter.plot_graph_2d(graph=self.graph,
                                             filename="betweenness_graph",
                                             centrality_scores=dict(sorted(scores_dict.items())),
                                             title="Betweenness Centrality Scores",
                                             cmap='jet',
                                             node_size=60,
                                             alpha=1,
                                             font_size=6,
                                             font_color='black')

            self.graph_plotter.plot_figures(scores_dict=scores_dict,
                                            centrality_type_name=CentralityType.BET.display_name)

            self.graph_plotter.plot_histogram(scores_dict=scores_dict,
                                              centrality_type_name=CentralityType.BET.display_name)

        elif CentralityType.CLOS == centrality_type:
            filename = self._get_cvs_filename(centrality_type=centrality_type)
            filename_path = os.path.join(centrality_root_path, filename)
            scores_dict = nx.closeness_centrality(G=self.graph, distance='weight')
            utils.write_centrality_scores_to_file(
                filename_path=filename_path,
                score_dict=scores_dict)

            logging.info(f"{filename_path} has been saved to the disk!")

            self.graph_plotter.plot_graph_2d(graph=self.graph,
                                             filename="closeness_graph",
                                             centrality_scores=dict(sorted(scores_dict.items())),
                                             title="Closeness Centrality Scores",
                                             cmap='jet',
                                             node_size=60,
                                             alpha=1,
                                             font_size=6,
                                             font_color='black')

            self.graph_plotter.plot_figures(scores_dict=scores_dict,
                                            centrality_type_name=CentralityType.CLOS.display_name)

            self.graph_plotter.plot_histogram(scores_dict=scores_dict,
                                              centrality_type_name=CentralityType.CLOS.display_name)
        else:
            filename = self._get_cvs_filename(centrality_type=CentralityType.DEG)
            filename_path = os.path.join(centrality_root_path, filename)
            scores_dict = dict(self.graph.degree())
            utils.write_centrality_scores_to_file(
                filename_path=filename_path,
                score_dict=scores_dict,
                formatter='%d %8.5f\n')

            logging.info(f"{filename_path} has been saved to the disk!")

            self.graph_plotter.plot_graph_2d_degree(graph=self.graph, title="Degree")

            self.graph_plotter.plot_graph_2d_degree_with_coords(graph=self.graph,
                                                                title="Degree Centrality Scores")

            self.graph_plotter.plot_figures(scores_dict=dict(self.graph.degree()),
                                            centrality_type_name=CentralityType.DEG.display_name,
                                            use_decimal=False)

            self.graph_plotter.plot_histogram(scores_dict=dict(self.graph.degree()),
                                              centrality_type_name=CentralityType.DEG.display_name)

    def build_graph(self):
        """ Builds a graph using networkx from contact and coordinate data.

        This method constructs a graph where each node represents a residue in a protein structure.
        Nodes are added with various attributes such as coordinates (x, y, z), chain id, residue number,
        and insertion information. Edges are added based on the provided contact data, with weights representing
        the interaction strength between residues.

        Returns
        -------
        None
            This method modifies the `self.graph` object in place, adding nodes and edges.

        Process:
        --------
        1. Loads contact data (source, target residues, and interaction value).
        2. Loads coordinates for each residue (x, y, z) and additional information such as chain ID, residue number,
           and insertion code.
        3. Adds nodes with their respective attributes to the graph.
        4. Adds weighted edges based on the contact data.

        Notes
        -----
        - If any residue is missing in `actual_residue_number_map`, it is skipped with a warning.
        - The contact data edges are weighted using the inverse of the interaction value from the contact dataframe.
        """
        node_attrs = [X, Y, Z, CHAIN_ID, RESIDUE_NUMBER, INSERTION, RESIDUE_NAME]
        residue_edges_file_path = os.path.join(self.destination_output_path, self.pdb_name,
                                               f'{self.pdb_name}_{RESIDUE_EDGES_FILE}')
        residue_average_coords_file_path = os.path.join(self.destination_output_path, self.pdb_name,
                                                        f'{self.pdb_name}_{RESIDUE_AVERAGE_COORDINATES_FILE}')

        residue_edges_df = utils.get_df(output_file_path=residue_edges_file_path,
                                        columns=[SOURCE_RESIDUE, TARGET_RESIDUE, VALUE],
                                        dtypes={SOURCE_RESIDUE: int, TARGET_RESIDUE: int, VALUE: float},
                                        sort_keys=[SOURCE_RESIDUE, TARGET_RESIDUE])

        residue_average_coords_df = utils.get_df(residue_average_coords_file_path,
                                                 columns=[RESIDUE_INDEX, X_COORD, Y_COORD, Z_COORD],
                                                 dtypes={RESIDUE_INDEX: int,
                                                         X_COORD: float,
                                                         Y_COORD: float,
                                                         Z_COORD: float},
                                                 sort_keys=[RESIDUE_INDEX])

        sources = residue_edges_df[SOURCE_RESIDUE].to_numpy(dtype=np.int32)
        targets = residue_edges_df[TARGET_RESIDUE].to_numpy(dtype=np.int32)
        weights = 1. / residue_edges_df[VALUE].to_numpy(dtype=np.float64)

        self.graph = nx.Graph(name=GRAPH_NAME)

        missing_residues = 0

        for residue_index, x_coord, y_coord, z_coord in zip(residue_average_coords_df[RESIDUE_INDEX],
                                                            residue_average_coords_df[X_COORD],
                                                            residue_average_coords_df[Y_COORD],
                                                            residue_average_coords_df[Z_COORD]):

            if residue_index not in self.actual_residue_number_map:
                logging.warning(f"Warning: residue {residue_index} not found in actual_residue_number_map.")
                missing_residues += 1
                continue

            chain_id, residue_number, insertion, residue_name = self.actual_residue_number_map[residue_index]
            attrs = {key: value for key, value in
                     zip(node_attrs, [x_coord, y_coord, z_coord, chain_id, residue_number, insertion, residue_name])}

            self.graph.add_node(residue_index, **attrs)

        if missing_residues:
            logging.warning(f"Skipped {missing_residues} residues due to missing metadata.")

        self.graph.add_weighted_edges_from(zip(sources, targets, weights))

        self.number_of_nodes = self.graph.number_of_nodes()

        self.graph_plotter.plot_graph_2d(graph=self.graph,
                                         filename="graph",
                                         title="Residue Interaction Network",
                                         equal_axis=True)

        self.graph_plotter.plot_graph_3d(graph=self.graph,
                                         filename="graph_3D",
                                         title="Residue Interaction Network 3D")
        return self.graph


def main():
    pass


if __name__ == '__main__':
    main()
