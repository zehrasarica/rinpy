# -*- coding: utf-8 -*-
"""
    __author__ = 'Zehra Sarica'
    __credits__ = ''
    __email__ = ['sarica16@itu.edu.tr','zehraacar559@gmail.com']
"""

import logging
import os
import time
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.offline as pyo

from rinpy import log_util
from rinpy import utils
from rinpy.constants import RESIDUE_NUMBER, MODE, CHAIN_ID, INSERTION, HIGH_PERCENTAGE_TEMPLATE, RESIDUE_NAME, \
    CENTRALITY_SCORE, TXT_EXT, PDB_EXT, PML_EXT, CSV_EXT, HTML_EXT, PNG_EXT, X, Y, Z
from rinpy.log_util import log_details
from rinpy.pymol_utils import PymolUtils
from rinpy.style_config import FONT_STYLES, FONT_FAMILY, EDGE_COLOR, COLOR_PALETTE
from rinpy.utils import CentralityType

HINGE_MODES_DIR = "hinge_modes"


class HingeAnalyzer:
    def __init__(self, graph, pdb_name, destination_output_path, actual_residue_number_map):
        self.graph = graph
        self.pdb_name = pdb_name
        self.destination_output_path = Path(destination_output_path)
        self.actual_residue_number_map = actual_residue_number_map
        mpl.rcParams['font.family'] = 'Times New Roman'
        mpl.rcParams['font.size'] = 12
        self.pymol_utils = PymolUtils(pdb_name=pdb_name,
                                      destination_output_path=destination_output_path,
                                      actual_residue_number_map=actual_residue_number_map)

        utils.create_folder_not_exists(self.destination_output_path / self.pdb_name / HINGE_MODES_DIR)

    def compute_laplacian_modes(self, num_modes):
        nodes = list(self.graph.nodes())
        nodes.sort()
        laplacian = nx.laplacian_matrix(self.graph, nodelist=nodes).astype(float)
        eigen_values, eigen_vectors = np.linalg.eigh(laplacian.toarray())
        fiedler_vector = eigen_vectors[:, 1].real
        filtered_data = eigen_vectors[:, 1:num_modes + 1]  # Fiedler +  3 more

        start = time.time()
        self.plot_fiedler(fiedler_vector=fiedler_vector)
        log_util.log_elapsed_time1("plot_fiedler: ", start, time.time())

        return eigen_values, eigen_vectors, fiedler_vector, filtered_data

    @staticmethod
    def find_hinge_residues_with_sign(mode_data):
        hinge_residues = set()
        for i in range(1, len(mode_data)):
            if np.sign(mode_data[i]) != np.sign(mode_data[i - 1]):
                hinge_residues.add(i + 1)
        return sorted(hinge_residues)

    def save_and_plot_mode_to_file(self, mode_data, mode_name, hinge_residues):
        mode_output_path = os.path.join(self.destination_output_path, self.pdb_name, HINGE_MODES_DIR)
        utils.create_folder_not_exists(mode_output_path)
        nodes = self.graph.nodes()
        actual_hinge_residues_tuple = self.get_actual_hinge_residues_tuple(nodes=nodes,
                                                                           hinge_residues=hinge_residues)

        full_path = os.path.join(str(mode_output_path), f"{mode_name}_hinge_residues{TXT_EXT}")
        self.write_to_file(actual_hinge_residues=actual_hinge_residues_tuple, full_path=full_path)

        residue_keys = [(nodes[n][RESIDUE_NUMBER], nodes[n][CHAIN_ID], nodes[n][INSERTION]) for n in nodes]

        mode_df = pd.DataFrame({
            RESIDUE_NUMBER: residue_keys,
            'mode': mode_data
        })

        residue_to_mode = mode_df.set_index(RESIDUE_NUMBER)[MODE].to_dict()

        base_pdb_df = utils.get_base_pdb_df(residue_to=residue_to_mode,
                                            destination_output_path=self.destination_output_path,
                                            pdb_name=self.pdb_name)

        out_filename = os.path.join(str(mode_output_path), f"{mode_name}{PDB_EXT}")
        utils.write_ppdb_to_pdb_file_atom_and_hetatom(ppdb=base_pdb_df, out_filename=out_filename)

        full_path = os.path.join(str(mode_output_path), f"{mode_name}{PML_EXT}")

        high_percentage_residues_df = pd.read_csv(os.path.join(str(self.destination_output_path),
                                                               self.pdb_name,
                                                               f'{self.pdb_name}_{HIGH_PERCENTAGE_TEMPLATE.format(type=CentralityType.BET.display_name)}'),
                                                  sep=";",
                                                  header=None,
                                                  names=[RESIDUE_NAME, CHAIN_ID, RESIDUE_NUMBER, INSERTION,
                                                         CENTRALITY_SCORE])
        self.pymol_utils.export_pymol_script_hinge(
            full_path_to_pdb=out_filename,
            residues=hinge_residues,
            full_path=full_path,
            high_percentage_residues=list(
                zip(high_percentage_residues_df[RESIDUE_NUMBER], high_percentage_residues_df[INSERTION])))

    def _save_selected_eigen_vectors(self, eigen_vectors):
        filtered_data_df = pd.DataFrame(eigen_vectors)
        num_columns = filtered_data_df.shape[1]
        headers = list(range(1, num_columns + 1))
        filtered_data_df.columns = headers

        residue_ids = [
            f"{self.actual_residue_number_map[i + 1][1]}"
            + (f"{self.actual_residue_number_map[i + 1][2]}" if len(self.actual_residue_number_map[i + 1]) > 2 and
                                                                self.actual_residue_number_map[i + 1][2] else "")
            if (i + 1) in self.actual_residue_number_map else str(i + 1)
            for i in range(filtered_data_df.shape[0])
        ]

        chain_ids = [
            f"{self.actual_residue_number_map[i + 1][0]}"
            if (i + 1) in self.actual_residue_number_map else str(i + 1)
            for i in range(filtered_data_df.shape[0])
        ]

        filtered_data_df.insert(0, "Residue Number", residue_ids)
        filtered_data_df.insert(1, "Chain Id", chain_ids)

        filtered_data_path = self.destination_output_path / self.pdb_name / HINGE_MODES_DIR / f'{self.pdb_name}_eigenvectors{CSV_EXT}'
        filtered_data_df.to_csv(filtered_data_path, index=False)

    def compute_hinge_residues_with_sign(self, num_modes=None):
        if num_modes is None:
            num_modes = 4

        eigen_values, eigen_vectors, fiedler_vector, filtered_data = self.compute_laplacian_modes(
            num_modes=num_modes)

        self._save_selected_eigen_vectors(eigen_vectors=filtered_data)

        for i in range(filtered_data.shape[1]):
            mode_data = filtered_data[:, i]
            hinge_residues = self.find_hinge_residues_with_sign(mode_data=mode_data)
            if hinge_residues and len(hinge_residues) > 0:
                mode_data = np.where(mode_data < 0, -1, 1)
                self.save_and_plot_mode_to_file(mode_data=mode_data,
                                                mode_name=f'{self.pdb_name}_laplacian_mode_{(i + 1)}',
                                                hinge_residues=hinge_residues)

                node_ids = sorted(list(self.graph.nodes()))
                cluster_labels = {node_ids[i]: int(mode_data[i]) for i in range(len(node_ids))}
                full_path = os.path.join(self.destination_output_path, self.pdb_name, HINGE_MODES_DIR,
                                         f'{self.pdb_name}_laplacian_mode_{(i + 1)}_hinge_interactive_clusters_3d{HTML_EXT}')
                self.plot_graph_interactive_clusters_3d(self.graph, cluster_labels, hinge_residues, full_path)
                start = time.time()

                p = os.path.join(self.destination_output_path, self.pdb_name, HINGE_MODES_DIR,
                                 f"{self.pdb_name}_laplacian_mode_{(i + 1)}_graph_clusters_2d{PNG_EXT}")
                self.plot_graph_clusters_2d(self.graph, cluster_labels, hinge_residues, full_path=p)
                log_util.log_elapsed_time1("plot_graph_clusters_2d: ", start, time.time())

    @staticmethod
    def get_actual_hinge_residues_tuple(nodes, hinge_residues):
        return [(
            f"{nodes[hr][RESIDUE_NAME]};{nodes[hr][CHAIN_ID]};{nodes[hr][RESIDUE_NUMBER]};'{nodes[hr][INSERTION]}'")
            for hr in hinge_residues
        ]

    def get_full_save_path(self, filename: str, extension: str = "png"):
        return os.path.join(self.destination_output_path, self.pdb_name, f"{self.pdb_name}_{filename}.{extension}")

    def write_to_file(self, actual_hinge_residues, full_path=None):
        if full_path is None:
            full_path = self.get_full_save_path(filename="hinge_residues", extension="txt")

        with open(full_path, 'w') as f:
            for i, hinge_residue in enumerate(actual_hinge_residues):
                if i < len(actual_hinge_residues) - 1:
                    f.write(f"{hinge_residue}\n")
                else:
                    f.write(f"{hinge_residue}")

    @log_details("Plotting Fiedler vector")
    def plot_fiedler(self, fiedler_vector, sort=False):
        """ Plots and optionally saves the Fiedler vector.

        Parameters:
        - fiedler_vector (np.ndarray): The Fiedler vector to plot.
        - save_path (str, optional): Path to save the figure (without extension). If None, just shows the plot.
        - sort (bool): Whether to plot the sorted Fiedler vector (for smoother visualization).
        - file_format (str): Image format for saving (e.g., 'png', 'pdf').
        - dpi (int): Resolution of the saved image.
        """
        full_path = self.get_full_save_path(filename="fiedler_vector")
        vector = fiedler_vector.copy()
        if sort:
            vector = np.sort(vector)

        fig, ax = plt.subplots(figsize=(10, 6))

        next_round = int(np.ceil(len(vector) / 100.0)) * 100
        ax.set_xlim(0, next_round)

        residue_labels = [
            f"{i + 1} ({self.actual_residue_number_map[i + 1][0]}, {utils.get_residue_id_by_tuple(self.actual_residue_number_map[i + 1])})"
            for i in range(len(vector))
            if (i + 1) in self.actual_residue_number_map
        ]
        x_ticks = ax.get_xticks().astype(int)
        x_ticks = x_ticks[x_ticks < len(vector)]

        if (len(vector) - 1) not in x_ticks:
            x_ticks = np.append(x_ticks, len(vector) - 1)

        ax.plot(vector, marker='o', linestyle='-', color='blue', label='Fiedler Vector')
        ax.axhline(y=0, color='red', linestyle='--', linewidth=1, label='Zero Line')
        ax.set_title(f"Fiedler Vector - {self.pdb_name}", fontdict=FONT_STYLES["title"])
        ax.set_xticks(x_ticks)
        ax.set_xticklabels([f"{residue_labels[i]}" for i in x_ticks], rotation=45, **{
            "fontsize": FONT_STYLES["xtick"]["labelsize"],
            "fontfamily": FONT_STYLES["xtick"]["fontfamily"]
        })
        ax.set_xlabel("Residue Index (Chain Id, Residue Number)", fontdict=FONT_STYLES["xlabel"])

        for label in ax.get_yticklabels():
            label.set_fontsize(FONT_STYLES["ytick"]["labelsize"])
            label.set_fontfamily(FONT_STYLES["ytick"]["fontfamily"])

        ax.set_ylabel("Value", fontdict=FONT_STYLES["ylabel"])
        ax.grid(True, linestyle=':', alpha=0.7)
        ax.legend(prop={"family": FONT_FAMILY, "size": FONT_STYLES["legend"]["fontsize"]})

        next_round = int(np.ceil(len(vector) / 100.0)) * 100
        ax.set_xlim(0, next_round)

        fig.tight_layout()

        fig.savefig(full_path, dpi=300)
        logging.info(f"Saving Fiedler vector to: {os.path.abspath(full_path)}")
        plt.close(fig)

    def plot_graph_clusters_2d(self, graph, cluster_labels, hinge_residues, full_path=None):
        if full_path is None:
            full_path = self.get_full_save_path(filename="graph_clusters_2d")
        unique_clusters = sorted(set(cluster_labels.values()))
        color_map = {cluster_id: COLOR_PALETTE[i % len(COLOR_PALETTE)] for i, cluster_id in
                     enumerate(unique_clusters)}

        fig, ax = plt.subplots()

        for cluster_id in unique_clusters:
            cluster_nodes = [n for n in graph.nodes() if cluster_labels[n] == cluster_id]
            xs = [graph.nodes[n]['x'] for n in cluster_nodes]
            ys = [graph.nodes[n]['y'] for n in cluster_nodes]
            ax.scatter(xs, ys, s=100, label=f"Cluster ({cluster_id})", color=color_map[cluster_id],
                       edgecolor=color_map[cluster_id], linewidth=1, alpha=0.6)

            for n in cluster_nodes:
                x, y = graph.nodes[n]['x'], graph.nodes[n]['y']
                ax.text(x, y, str(n), fontsize=6, ha='center', va='center', color='black', alpha=0.6)

        for u, v in graph.edges():
            x = [graph.nodes[u]['x'], graph.nodes[v]['x']]
            y = [graph.nodes[u]['y'], graph.nodes[v]['y']]
            ax.plot(x, y, color=EDGE_COLOR, alpha=0.5, linewidth=0.5)

        hinge_x = [graph.nodes[n]['x'] for n in hinge_residues]
        hinge_y = [graph.nodes[n]['y'] for n in hinge_residues]
        ax.scatter(hinge_x, hinge_y, s=120, c='yellow', edgecolors='black', marker='*',
                   label='Hinge Residues', zorder=10)

        for n in hinge_residues:
            x, y = graph.nodes[n]['x'], graph.nodes[n]['y']
            ax.text(x, y, str(n), fontsize=7, ha='center', va='center', color='black', fontweight='bold')

        ax.set_xlabel("X Coordinate", fontdict=FONT_STYLES["xlabel"])
        ax.set_ylabel("Y Coordinate", fontdict=FONT_STYLES["ylabel"])
        ax.set_title("2D Spectral Clustering of Residue Interaction Network", fontdict=FONT_STYLES["title"])
        ax.legend(prop={"family": FONT_FAMILY, "size": FONT_STYLES["legend"]["fontsize"]})
        ax.axis('equal')
        ax.grid(True, which='both', axis='both', color='gray', linestyle='-', linewidth=0.5, alpha=0.3)

        self.update_tick_label_style(ax.get_xticklabels())
        self.update_tick_label_style(ax.get_yticklabels())

        fig.tight_layout()
        fig.savefig(full_path, dpi=300)
        plt.close(fig)

    @staticmethod
    def update_tick_label_style(labels):
        for label in labels:
            label.set_fontsize(FONT_STYLES["ytick"]["labelsize"])
            label.set_fontfamily(FONT_STYLES["ytick"]["fontfamily"])

    def plot_graph_clusters_3d(self, graph, cluster_labels, hinge_residues):
        """ 3D plot (X-Y-Z) with manually assigned high-contrast colors."""
        full_path = self.get_full_save_path(filename="graph_clusters_3d")

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        unique_clusters = sorted(set(cluster_labels.values()))
        color_map = {cluster_id: COLOR_PALETTE[i % len(COLOR_PALETTE)] for i, cluster_id in
                     enumerate(unique_clusters)}

        legend_handles = []

        for cluster_id in unique_clusters:
            cluster_nodes = [n for n in graph.nodes() if cluster_labels[n] == cluster_id]
            xs = [graph.nodes[n]['x'] for n in cluster_nodes]
            ys = [graph.nodes[n]['y'] for n in cluster_nodes]
            zs = [graph.nodes[n]['z'] for n in cluster_nodes]
            scatter = ax.scatter(xs, ys, zs, s=100, label=f"Cluster ({cluster_id})", color=color_map[cluster_id])
            legend_handles.append(scatter)

        for u, v in graph.edges():
            x = [graph.nodes[u][X], graph.nodes[v][X]]
            y = [graph.nodes[u][Y], graph.nodes[v][Y]]
            z = [graph.nodes[u][Z], graph.nodes[v][Z]]
            ax.plot(x, y, z, color=EDGE_COLOR, alpha=0.5, linewidth=0.5)

        hinge_scatter = []
        for n in hinge_residues:
            if n in graph.nodes():
                x, y, z = graph.nodes[n][X], graph.nodes[n][Y], graph.nodes[n][Z]
                ax.scatter(x, y, z, s=150, c='yellow', edgecolors='black', marker='*')
                ax.text(x, y, z, str(n), fontsize=7, ha='center', va='center', color='black', fontweight='bold')
                hinge_scatter.append(ax.scatter([], [], s=150, c='yellow', edgecolors='black', marker='*'))

        legend_handles.append(hinge_scatter[0])
        ax.legend(handles=legend_handles,
                  labels=[f"Cluster ({i})" for i in unique_clusters] + ["Hinge Residues"],
                  prop={"family": FONT_FAMILY, "size": FONT_STYLES["legend"]["fontsize"]}, loc='upper right',
                  bbox_to_anchor=(1.2, 1.01))

        ax.set_xlabel("X Coordinate", fontdict=FONT_STYLES["xlabel"])
        ax.set_ylabel("Y Coordinate", fontdict=FONT_STYLES["ylabel"])
        ax.set_zlabel("Z Coordinate", fontdict=FONT_STYLES["zlabel"])
        plt.title("3D Spectral Clustering of Residue Interaction Network", fontdict=FONT_STYLES["title"])

        self.update_tick_label_style(ax.get_xticklabels())
        self.update_tick_label_style(ax.get_yticklabels())
        self.update_tick_label_style(ax.get_zticklabels())

        plt.tight_layout()
        fig.savefig(full_path, dpi=300)
        plt.close(fig)

    # @staticmethod
    @staticmethod
    def plot_graph_interactive_clusters_3d(graph, cluster_labels, hinge_residues, full_path):
        unique_clusters = sorted(set(cluster_labels.values()))
        color_map = {
            cluster_id: COLOR_PALETTE[i % len(COLOR_PALETTE)]
            for i, cluster_id in enumerate(unique_clusters)
        }
        node_traces = []
        for cluster_id in unique_clusters:
            cluster_nodes = [n for n in graph.nodes if cluster_labels[n] == cluster_id]
            x = [graph.nodes[n][X] for n in cluster_nodes]
            y = [graph.nodes[n][Y] for n in cluster_nodes]
            z = [graph.nodes[n][Z] for n in cluster_nodes]

            labels = [
                (
                    f"{graph.nodes[n][RESIDUE_NAME]} "
                    f"(Residue Number: {utils.get_residue_id(graph.nodes[n])}, Chain: {graph.nodes[n][CHAIN_ID]}), "
                    f"Cluster ({cluster_id})"
                )
                for n in cluster_nodes
            ]

            trace = go.Scatter3d(
                x=x,
                y=y,
                z=z,
                mode='markers',
                marker=dict(size=5, color=color_map[cluster_id]),
                text=labels,
                hoverinfo='text',
                name=f'Cluster {cluster_id}'
            )
            node_traces.append(trace)

        # Create edge trace
        edge_x, edge_y, edge_z = [], [], []
        for u, v in graph.edges():
            edge_x += [graph.nodes[u][X], graph.nodes[v][X], None]
            edge_y += [graph.nodes[u][Y], graph.nodes[v][Y], None]
            edge_z += [graph.nodes[u][Z], graph.nodes[v][Z], None]

        edge_trace = go.Scatter3d(
            x=edge_x,
            y=edge_y,
            z=edge_z,
            mode='lines',
            line=dict(color='#1f4e79', width=1),
            hoverinfo='none',
            name='Edges'
        )

        hinge_trace = go.Scatter3d(
            x=[graph.nodes[n][X] for n in hinge_residues],
            y=[graph.nodes[n][Y] for n in hinge_residues],
            z=[graph.nodes[n][Z] for n in hinge_residues],
            mode='markers+text',
            marker=dict(size=7, color='yellow', symbol='diamond', line=dict(color='black', width=1)),
            # text=[f"Hinge {n}" for n in hinge_residues],
            # hoverinfo='text',
            name='Hinge Residues',
            hoverinfo='text',  # hides x, y, z and only shows hovertext
            hovertext=[
                f"<b>Residue Name: {graph.nodes[n][RESIDUE_NAME]}</b><br>"
                f"<b>Chain: {graph.nodes[n][CHAIN_ID]}</b><br>"
                f"<b>Hinge Residue Number: {utils.get_residue_id(graph.nodes[n])}</b>"
                f"<br>────────────────────<br>"
                f"x={graph.nodes[n]['x']:.3f}<br>"
                f"y={graph.nodes[n]['y']:.3f}<br>"
                f"z={graph.nodes[n]['z']:.3f}"
                for n in hinge_residues],
            hoverlabel=dict(
                bgcolor='yellow',
                font=dict(color='black')
            )
        )

        fig = go.Figure(data=[edge_trace, hinge_trace] + node_traces)
        fig.update_layout(
            scene=dict(
                xaxis_title='X Coordinate',
                yaxis_title='Y Coordinate',
                zaxis_title='Z Coordinate'
            ),
            margin=dict(l=0, r=0, b=0, t=40),
            title='Interactive 3D Spectral Clustering of Residue Interaction Network'
        )

        pyo.plot(fig, filename=full_path, auto_open=False)
        logging.info(f"Interactive plot saved to: {full_path}")

    def export_pymol_script(self, hinge_residues):
        """ Generate a PyMOL script to visualize hinge residues.
        Parameters:
            - hinge_residues: List of hinge residue numbers (int)
        """

        full_path = self.get_full_save_path(filename="highlight_hinge_residues", extension="pml")
        with open(full_path, 'w') as file:
            file.write(f"load {self.pdb_name}.pdb\n")
            file.write("hide everything\n")
            file.write("show cartoon\n")
            file.write("color gray80\n\n")

            hinge_sel = "+".join(
                str(self.actual_residue_number_map[residue_index][1]) for residue_index in hinge_residues)

            file.write(f"select hinge_residues, resi {hinge_sel}\n")
            file.write("show sticks, hinge_residues\n")
            file.write("color yellow, hinge_residues\n")
            file.write("set stick_radius, 0.25, hinge_residues\n\n")

            file.write("set label_size, 18\n")
            file.write("set label_font_id, 7\n")
            file.write("label hinge_residues and name CA, resn + resi\n")
            file.write("set label_position, [0.0, -0.3, -0.5], hinge_residues\n")
            file.write("set label_color, white\n\n")

            file.write("\nzoom hinge_residues, 10\n")


def main():
    pass


if __name__ == '__main__':
    main()
