# -*- coding: utf-8 -*-
"""
    __author__ = 'Zehra Sarica'
    __email__ = ['sarica16@itu.edu.tr','zehraacar559@gmail.com']
"""

from rinpy.utils import *

_END = "END"
_ATOM = "ATOM"
_OCCUPANCY = 'occupancy'


class Trajectory:

    def __init__(self, input_trajectory_file, output_path, stride=10):
        if output_path is None:
            raise ValueError('You must provide an output path to proceed.')

        self.input_trajectory_file = input_trajectory_file
        self.output_path = output_path
        self.stride = stride

    def parse_to_pdb(self):
        """ All extracted PDBs absolute path list
        :return: return list: ['output/0001/0001.pdb', 'output/0002/0002.pdb', 'output/0003/0003.pdb']
        """
        delete_folders_in_path(self.output_path)

        snapshots = []
        current_snapshot = []

        lines = self.read_lines()

        for line in lines:
            if line.startswith(_END):
                snapshots.append(current_snapshot)
                current_snapshot = []
            elif line.startswith(_ATOM):
                current_snapshot.append(line)

        num_snapshots = len(snapshots)
        output_path_pdb_ids = []

        for i in range(0, num_snapshots, self.stride):
            snapshot = snapshots[i]
            filename = f"{i + 1:04d}"
            parent_dir = os.path.join(self.output_path, filename)
            create_folder_not_exists(parent_dir)
            output_file = os.path.join(parent_dir, f"{filename}.pdb")

            ppdb = PandasPdb()
            ppdb.read_pdb_from_list(snapshot)
            ppdb.df[_ATOM].loc[ppdb.df[_ATOM][
                _OCCUPANCY].isna(), _OCCUPANCY] = 0  # if occupancy is NAN value, 0 value is set as default.
            write_ppdb_to_pdb_file(ppdb, output_file)
            output_path_pdb_ids.append(output_file)
            logging.info(f"Snapshot {i + 1} saved as '{output_file}'")

        logging.info(f"{output_path_pdb_ids}")
        return output_path_pdb_ids

    def read_lines(self):
        with open(self.input_trajectory_file, 'r') as file:
            lines = file.readlines()
        return lines


def main():
    pass


if __name__ == '__main__':
    main()
