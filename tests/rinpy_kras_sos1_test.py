from rinpy import RINProcess


def start_rinpy_test():
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
        },
    }

    rinpy = RINProcess(input_path="kras_sos1_input",
                       output_path="output",
                       pdb_ids=None,
                       ligand_dict=None,
                       calculation_options=calculation_options,
                       trajectory_file=None,
                       stride=1)
    rinpy.start_process()


def main():
    start_rinpy_test()


if __name__ == '__main__':
    main()
