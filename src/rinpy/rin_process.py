# -*- coding: utf-8 -*-
"""
    __author__ = 'Zehra Sarica'
    __email__ = ['sarica16@itu.edu.tr','zehraacar559@gmail.com']
"""

import argparse
import glob
import json
from datetime import datetime
from pathlib import Path

from rinpy import log_util, logging_config
from rinpy.centrality_analyzer import CentralityAnalyzer
from rinpy.centrality_pdb_mapper import CentralityPdbMapper
from rinpy.font_utils import load_fonts_once
from rinpy.hinge_analyzer import HingeAnalyzer
from rinpy.log_util import *
from rinpy.quantile_analyzer import QuantileAnalyzer
from rinpy.residue_network_builder import ResidueNetworkBuilder
from rinpy.trajectory import Trajectory
from rinpy.utils import *

PDB_ID_LENGTH = 4
CLUSTER_NUMBER = 'cluster_number'


class RINProcess:
    _fonts_initialized = False

    def __init__(self, output_path=None, input_path=None, pdb_ids=None, trajectory_file=None, stride=1,
                 calculation_options=None, ligand_dict=None, num_workers=None):

        if not RINProcess._fonts_initialized:
            load_fonts_once()
            RINProcess._fonts_initialized = True

        self._validate_input_options(output_path=output_path,
                                     input_path=input_path,
                                     pdb_ids=pdb_ids,
                                     trajectory_file=trajectory_file)

        self.output_path = output_path
        self.input_path = input_path
        self.pdb_ids = pdb_ids
        self.trajectory_file = trajectory_file
        self.stride = stride
        self.ligand_dict = ligand_dict
        self.num_workers = num_workers

        self.calculation_options = calculation_options

        if self.calculation_options is None or not isinstance(calculation_options, dict):
            self.calculation_options = self._get_default_calculation_options()

    @staticmethod
    def _validate_input_options(output_path, input_path, pdb_ids, trajectory_file):
        """Validates the input options provided. The input can be one of the following:
            input_path, pdb_ids, trajectory_file
        """
        if output_path is None:
            raise ValueError('You must provide an output path to proceed.')

        inputs_provided = sum(x is not None for x in [input_path, pdb_ids, trajectory_file])
        if inputs_provided == 0:
            raise ValueError("You must provide exactly one input: 'input_path', 'pdb_ids', or 'trajectory_file'.")
        elif inputs_provided > 1:
            raise ValueError(
                "Only one input can be provided. Please provide only one of 'input_path', 'pdb_ids', or 'trajectory_file'.")

    @staticmethod
    def _get_default_calculation_options():
        """Returns the default calculation options from the calculation_options.json file
        in the current working directory."""
        json_file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), CALCULATION_OPTION_DEFAULT_JSON_PATH)
        try:
            with open(json_file_path, 'r') as json_file:
                data = json.load(json_file)
            if data is None:
                raise ValueError("Data in the JSON file is None. Please provide valid data.")
            return data
        except FileNotFoundError:
            print(f"File not found: {json_file_path}")
        except json.JSONDecodeError:
            print(f"Error decoding JSON in {json_file_path}")
        except ValueError as e:
            print(e)
        except Exception as e:
            print(f"An unexpected error occurred: {e}")

    @staticmethod
    @log_time("Fetching PDB IDs from Protein Data Bank (RCSB)")
    def _fetch_pdb_ids(pdb_ids, output_path=None):
        """Fetches PDB IDs from the Protein Data Bank (RCSB)."""
        total_start_time = datetime.strptime(datetime.now().strftime("%H:%M:%S"), "%H:%M:%S")
        fetched_pdb_files = []
        output_path_pdb_ids = []
        for pdb_id in pdb_ids:
            if len(pdb_id) == PDB_ID_LENGTH:
                destination_path = os.path.join(output_path, pdb_id)
                create_folder_not_exists(destination_path)
                pdb_id_file = os.path.join(str(destination_path), pdb_id + PDB_EXT)
                if not os.path.exists(pdb_id_file):
                    ppdb = PandasPdb().fetch_pdb(pdb_id)
                    write_ppdb_to_pdb_file_all(ppdb, pdb_id_file)
                    fetched_pdb_files.append(pdb_id)
                    output_path_pdb_ids.append(pdb_id_file)
                else:
                    logging.warning(f"{pdb_id} has already fetched to {pdb_id_file}.")
                    output_path_pdb_ids.append(pdb_id_file)
            else:
                logging.error(f"{pdb_id} is NOT VALID PDB ID. It must be 4 character like 3t05!")
        if fetched_pdb_files:
            logging.info(
                f"PDB IDS: {fetched_pdb_files} have been fetched from the protein data bank and saved to the path: {output_path}")
            total_end_time = datetime.strptime(datetime.now().strftime("%H:%M:%S"), "%H:%M:%S")
            delta = total_end_time - total_start_time
            logging.info(f"Total fetch duration is {delta.total_seconds()} seconds.")

        return output_path_pdb_ids

    @log_with_stars("Residue Network Builder")
    @log_time("Residue Network Builder")
    def _process_residue_network_builder(self, pdb_name=None, pdb_path=None, use_preprocess=True,
                                         destination_output_path=None,
                                         calculation_options: dict = None) -> ResidueNetworkBuilder:
        """Process Residue Network Builder to prepare necessary input files for Residue Interaction Network."""
        het_atom_list = None

        if self.ligand_dict is not None and pdb_name is not None and pdb_name in self.ligand_dict:
            het_atom_list = self.ligand_dict[pdb_name]

        cutoff_value = DEFAULT_CUTOFF

        if CUTOFF in calculation_options:
            option = calculation_options[CUTOFF]
            if option.get(IS_CHECKED):
                cutoff_value = option[VALUE]

        rgb = ResidueNetworkBuilder(pdb_name=pdb_name,
                                    pdb_path=pdb_path,
                                    use_preprocess=use_preprocess,
                                    cutoff=cutoff_value,
                                    destination_output_path=destination_output_path,
                                    het_atom_list=het_atom_list,
                                    remove_hydrogen=calculation_options[REMOVE_HYDROGEN][IS_CHECKED],
                                    num_workers=self.num_workers)

        rgb.calculate_contact()

        return rgb

    @log_with_stars("Centrality Analyzer")
    @log_time("Centrality Analyzer")
    def _process_centrality_analyzer(self, pdb_name: str = None, number_of_amino_acid=None,
                                     destination_output_path=None,
                                     calculation_options=None,
                                     actual_residue_number_map=None) -> CentralityAnalyzer:
        """Process centrality analysis from the generated network based on atom-atom distance of a given PDB."""

        ca = CentralityAnalyzer(pdb_name=pdb_name,
                                number_of_amino_acid=number_of_amino_acid,
                                destination_output_path=destination_output_path,
                                calculation_options=calculation_options,
                                actual_residue_number_map=actual_residue_number_map)
        ca.build_network()
        ca.calculate_all_scores(use_parallel=True)
        return ca

    @log_with_stars("Hinge Residues Analyzer")
    @log_time("Hinge Residues Analyzer")
    def _process_hinge_analyzer(self, graph, pdb_name, destination_output_path,
                                actual_residue_number_map, calculation_options) -> HingeAnalyzer:
        """Process hinges from the generate network."""
        hc = HingeAnalyzer(graph=graph,
                           pdb_name=pdb_name,
                           destination_output_path=destination_output_path,
                           actual_residue_number_map=actual_residue_number_map)
        nums_cluster = None
        if CLUSTER_NUMBER in calculation_options and calculation_options[CLUSTER_NUMBER][IS_CHECKED]:
            nums_cluster = calculation_options[CLUSTER_NUMBER][VALUE]
        hc.compute_hinge_residues_with_sign(num_modes=nums_cluster)
        return hc

    @log_with_stars("Centrality PDB Mapper")
    @log_time("Centrality PDB Mapper")
    def _process_centrality_pdb_mapper(self, graph=None, pdb_name: str = None, destination_output_path: str = None,
                                       calculation_options: dict = None) -> CentralityPdbMapper:
        """Process creating PDBs to visualize by mapping centrality scores into the given PDB structure."""
        cpm = CentralityPdbMapper(graph=graph, pdb_name=pdb_name, destination_output_path=destination_output_path,
                                  calculation_options=calculation_options)
        cpm.process()
        return cpm

    @log_with_stars("Quantile Analyzer")
    @log_time("Quantile Analyzer")
    def _process_quantiles(self, ca: CentralityAnalyzer = None, use_parallel: bool = False):
        """Process quantiles to calculate the quantile values for the three centrality types."""
        ca.calculate_quantiles(use_parallel=use_parallel)

    @log_time("Quantile Analyzer")
    def _process_quantile_centrality(self, high_percentage_dict=None, centrality_type=CentralityType.BET,
                                     destination_output_path=None) -> QuantileAnalyzer:
        """Process quantile calculations to find common residues with high centrality scores among the given PDBs"""
        qa = QuantileAnalyzer(high_percentage_dict=high_percentage_dict, centrality_type=centrality_type,
                              destination_output_path=destination_output_path)
        qa.run_analysis()

        return qa

    @log_with_stars("Processing Quantile Outputs")
    @log_time("Processing Quantile Outputs")
    def _process_quantile_output(self, bet_high_percentage_dict=None, clos_high_percentage_dict=None,
                                 deg_high_percentage_dict=None, destination_output_path=None, calculation_options=None):
        if calculation_options[CentralityType.BET.display_name][IS_CHECKED] and bet_high_percentage_dict:
            logging.info(f"bet_high_percentage_dict: {bet_high_percentage_dict}")
            self._process_quantile_centrality(high_percentage_dict=bet_high_percentage_dict,
                                              centrality_type=CentralityType.BET,
                                              destination_output_path=destination_output_path)

        if calculation_options[CentralityType.CLOS.display_name][IS_CHECKED] and clos_high_percentage_dict:
            logging.info(f"clos_high_percentage_dict: {clos_high_percentage_dict}")
            self._process_quantile_centrality(high_percentage_dict=clos_high_percentage_dict,
                                              centrality_type=CentralityType.CLOS,
                                              destination_output_path=destination_output_path)
        if calculation_options[CentralityType.DEG.display_name][IS_CHECKED] and deg_high_percentage_dict:
            logging.info(f"deg_high_percentage_dict: {deg_high_percentage_dict}")
            self._process_quantile_centrality(high_percentage_dict=deg_high_percentage_dict,
                                              centrality_type=CentralityType.DEG,
                                              destination_output_path=destination_output_path)

    @log_time("Extracting PDB IDs to be processed.")
    def _extract_output_path_pdb_ids(self, output_path, input_path, pdb_ids, trajectory_file) -> list:
        """Extracts PDB IDs from the given input option."""
        output_path_pdb_ids = []
        if input_path is not None and input_path != "":
            logging.info(f"input_path: {input_path}")
            input_pdb_files = [pdb for pdb in glob.glob(os.path.join(input_path, "*.pdb"))]
            for input_pdb in input_pdb_files:
                input_pdb_path = Path(input_pdb)
                pdb_code = input_pdb_path.stem
                destination_path = os.path.join(output_path, pdb_code)
                create_folder_not_exists(destination_path)
                shutil.copy(input_pdb, destination_path)
                output_path_pdb_ids.append(os.path.join(destination_path, input_pdb_path.name))
        elif pdb_ids is not None and pdb_ids:
            logging.info(f"fetch_pdb_ids: {self.pdb_ids}")
            fetch_pdb_ids = [x.strip() for x in self.pdb_ids.split(',')]
            output_path_pdb_ids = self._fetch_pdb_ids(fetch_pdb_ids, output_path)
        elif trajectory_file is not None and trajectory_file != "":
            logging.info(f"trajectory_file: {trajectory_file}")
            trajectory = Trajectory(input_trajectory_file=self.trajectory_file, output_path=output_path,
                                    stride=self.stride)
            output_path_pdb_ids = trajectory.parse_to_pdb()
        else:
            logging.warning(f"There is no proper input selection to be processed! Please check your input!")
            raise Exception("Please Select Only One Input Option!!!")

        return output_path_pdb_ids

    def start_process(self) -> bool:
        # clear the debug.log file before starting the processes.
        logging_config.clear_logs()
        output_path = self.output_path
        logging.info(f"output_path: {output_path}")

        output_path_pdb_ids = self._extract_output_path_pdb_ids(output_path=output_path,
                                                                input_path=self.input_path,
                                                                pdb_ids=self.pdb_ids,
                                                                trajectory_file=self.trajectory_file)
        bet_high_percentage_dict = {}
        clos_high_percentage_dict = {}
        deg_high_percentage_dict = {}
        total_amino_acids_map = {}
        pdb_ids = []

        if output_path_pdb_ids:
            logging.info(f"Processing PDB files: {output_path_pdb_ids}")
            total_start_time = datetime.strptime(datetime.now().strftime("%H:%M:%S"), "%H:%M:%S")
            for pdb_file_path in output_path_pdb_ids:
                pdb_name = Path(str(pdb_file_path)).stem
                pdb_ids.append(pdb_name)
                start_time = datetime.strptime(datetime.now().strftime("%H:%M:%S"), "%H:%M:%S")

                log_util.log_star_message(pdb_name.upper())
                logging.info(f'Start time: {start_time.time()}')

                residue_network_builder = self._process_residue_network_builder(pdb_name=pdb_name,
                                                                                pdb_path=pdb_file_path,
                                                                                use_preprocess=True,
                                                                                destination_output_path=output_path,
                                                                                calculation_options=self.calculation_options)

                number_of_amino_acid = residue_network_builder.get_number_of_amino_acid()
                total_amino_acids_map[pdb_name] = number_of_amino_acid

                ca = self._process_centrality_analyzer(pdb_name=pdb_name,
                                                       number_of_amino_acid=number_of_amino_acid,
                                                       destination_output_path=output_path,
                                                       calculation_options=self.calculation_options,
                                                       actual_residue_number_map=residue_network_builder.get_actual_residue_number_map())

                bet_high_percentage_filename = f'{pdb_name}_{HIGH_PERCENTAGE_TEMPLATE.format(type=CentralityType.BET.display_name)}'
                clos_high_percentage_filename = f'{pdb_name}_{HIGH_PERCENTAGE_TEMPLATE.format(type=CentralityType.CLOS.display_name)}'
                deg_high_percentage_filename = f'{pdb_name}_{HIGH_PERCENTAGE_TEMPLATE.format(type=CentralityType.DEG.display_name)}'

                bet_high_percentage_dict[pdb_name] = os.path.join(output_path, pdb_name, bet_high_percentage_filename)
                clos_high_percentage_dict[pdb_name] = os.path.join(output_path, pdb_name, clos_high_percentage_filename)
                deg_high_percentage_dict[pdb_name] = os.path.join(output_path, pdb_name, deg_high_percentage_filename)

                cpm = self._process_centrality_pdb_mapper(graph=ca.graph,
                                                          pdb_name=pdb_name,
                                                          destination_output_path=output_path,
                                                          calculation_options=self.calculation_options)

                self._process_quantiles(ca=ca, use_parallel=True)

                hc = self._process_hinge_analyzer(graph=ca.graph,
                                                  pdb_name=pdb_name,
                                                  destination_output_path=output_path,
                                                  actual_residue_number_map=residue_network_builder.get_actual_residue_number_map(),
                                                  calculation_options=self.calculation_options)

                end_time = datetime.strptime(datetime.now().strftime("%H:%M:%S"), "%H:%M:%S")
                logging.info(f'End time: {end_time.time()}')
                delta = end_time - start_time
                logging.info(f"Elapsed duration for {pdb_name.upper()} is {delta.total_seconds()} seconds")
                log_util.log_star_message(pdb_name.upper())

            if len(output_path_pdb_ids) > 1:
                self._process_quantile_output(bet_high_percentage_dict=bet_high_percentage_dict,
                                              clos_high_percentage_dict=clos_high_percentage_dict,
                                              deg_high_percentage_dict=deg_high_percentage_dict,
                                              destination_output_path=output_path,
                                              calculation_options=self.calculation_options)

            logging.info(
                f"Total PDBs: {len(total_amino_acids_map)}, Total amino acids: {sum(total_amino_acids_map.values())}")

            total_end_time = datetime.strptime(datetime.now().strftime("%H:%M:%S"), "%H:%M:%S")
            delta = total_end_time - total_start_time
            total_seconds = int(delta.total_seconds())
            hours, remainder = divmod(total_seconds, 3600)
            minutes, seconds = divmod(remainder, 60)

            logging.info(f"Total elapsed duration is {delta.total_seconds()} seconds ({hours}h {minutes}m {seconds}s).")
            logging_config.save_logs(output_dir=output_path)

            return True
        else:
            logging.warning(f"There is no PDB file to be processed! Please check your input!")
            raise Exception("NO PDB FILES TO BE PROCESSED!")


def _load_dict_file(file_path: Path):
    if not file_path.exists():
        print(f"Warning: ligand_dict file not found at '{file_path}'. Skipping.")
        return None
    try:
        with file_path.open('r') as f:
            return json.load(f)
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON: {e}")
        return None


def build_parser():
    parser = argparse.ArgumentParser(
        description="RinPy: Residue Interaction Network construction and analysis from protein structures and trajectories."
    )
    parser.add_argument("--output_path",
                        type=str,
                        help="Specify output path which will store the results.",
                        required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--input_path",
        type=str,
        help="Specify the path of PDB files which will be analyzed."
    )
    group.add_argument(
        "--pdb_ids",
        type=str,
        help="Load a PDB file from the Protein Data Bank via its unique 4-letter ID, e.g., '6epl, 6epm, 6epn'"
    )
    group.add_argument(
        "--trajectory_file",
        type=str,
        help="Specify the trajectory file to be analyzed."
    )
    parser.add_argument("--stride",
                        type=int,
                        required=False,
                        help="Parsing the trajectory file every this stride number.",
                        default=1)
    parser.add_argument("--calculation_option_file",
                        type=Path,
                        required=True,
                        help="Required: Path to JSON file containing the calculation options in the dictionary.",
                        default=None)
    parser.add_argument("--ligand_dict_file",
                        type=Path,
                        required=False,
                        help="Optional: Path to JSON file containing the ligand dictionary.",
                        default=None)
    parser.add_argument("--num_workers",
                        type=int,
                        default=None,
                        help="Number of CPU worker processes to use (default: use all detected CPU cores.)")
    args = parser.parse_args()
    return args


def main():
    args = build_parser()

    ligand_dict = None

    if args and args.ligand_dict_file is not None:
        ligand_dict = _load_dict_file(args.ligand_dict_file)

    calculation_options = _load_dict_file(args.calculation_option_file)

    rin = RINProcess(input_path=args.input_path,
                     output_path=args.output_path,
                     pdb_ids=args.pdb_ids,
                     ligand_dict=ligand_dict,
                     calculation_options=calculation_options,
                     trajectory_file=args.trajectory_file,
                     stride=args.stride,
                     num_workers=args.num_workers)

    rin.start_process()


if __name__ == '__main__':
    main()
