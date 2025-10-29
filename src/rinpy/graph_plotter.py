# -*- coding: utf-8 -*-
"""
    __author__ = 'Zehra Sarica'
    __credits__ = ''
    __email__ = ['sarica16@itu.edu.tr','zehraacar559@gmail.com']
"""

import logging
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import plotly.graph_objects as go
import seaborn as sns
from matplotlib import ticker
from scipy.sparse.linalg import eigs

from rinpy import utils
from rinpy.constants import CHAIN_ID, RESIDUE_NAME
from rinpy.style_config import FONT_STYLES, FONT_FAMILY, EDGE_COLOR, COLOR_PALETTE

X_COORD_LABEL = "X Coordinate"
Y_COORD_LABEL = "Y Coordinate"
Z_COORD_LABEL = "Z Coordinate"


class GraphPlotter:
    def __init__(self, pdb_name, destination_output_path, actual_residue_number_map):
        self.pdb_name = pdb_name
        self.destination_output_path = destination_output_path
        self.actual_residue_number_map = actual_residue_number_map
        mpl.rcParams['font.family'] = FONT_FAMILY

    def get_full_save_path(self, filename: str, extension: str = "png"):
        return os.path.join(self.destination_output_path, self.pdb_name, f"{self.pdb_name}_{filename}.{extension}")

    @staticmethod
    def _apply_colorbar_style(cbar, cbar_label=None):
        """Apply consistent style to colorbar."""
        cbar.ax.tick_params(labelsize=FONT_STYLES["colorbar"]["tick_params"]["labelsize"])
        if cbar_label is not None:
            cbar.set_label(
                "Colorbar Label",
                fontsize=FONT_STYLES["colorbar"]["label"]["fontsize"],
                fontfamily=FONT_STYLES["colorbar"]["label"]["fontfamily"]
            )

    def plot_graph_2d(self, graph, filename="graph", centrality_scores: dict = None, title=None, cmap='jet',
                      node_size=40, alpha=0.8, font_size=6, font_color='black', equal_axis=False):
        """ Creates a 2D plot of a graph, optionally with centrality-based coloring of nodes.

           Parameters:
           - graph: The networkx graph object to be plotted.
           - filename: The name of the file to save the plot.
           - centrality_scores: Dictionary with node IDs as keys and centrality scores as values.
           - title: The title to be displayed on the plot.
           - cmap: Colormap for node color when centrality_scores is provided.
           - node_size: Size of nodes in the plot.
           - alpha: Transparency level for nodes and edges.
           - font_size: Font size for node labels.
           - font_color: Font color for node labels.
           - equal_axis: If True, sets the axes to be equal in scale.
       """

        full_path = self.get_full_save_path(filename=filename)

        plt.figure()

        pos = {n: (graph.nodes[n]['x'], graph.nodes[n]['y']) for n in graph.nodes}

        node_colors = [centrality_scores.get(n, 0) for n in graph.nodes] if centrality_scores else None

        nodes = nx.draw_networkx_nodes(graph,
                                       pos,
                                       node_color=node_colors,
                                       cmap=plt.get_cmap(cmap) if centrality_scores is not None else None,
                                       node_size=node_size,
                                       alpha=alpha)
        nx.draw_networkx_labels(graph,
                                pos,
                                labels={n: str(n) for n in graph.nodes()},
                                font_color=font_color,
                                font_size=font_size)

        nx.draw_networkx_edges(graph,
                               pos,
                               edge_color='gray',
                               alpha=alpha,
                               width=0.5)

        if node_colors:
            cbar = plt.colorbar(nodes)
            self._apply_colorbar_style(cbar)

        ax = plt.gca()
        ax.set_xlabel(X_COORD_LABEL, fontdict=FONT_STYLES["xlabel"])
        ax.set_ylabel(Y_COORD_LABEL, fontdict=FONT_STYLES["ylabel"])
        ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True,
                       labelsize=FONT_STYLES["tick"]["labelsize"])
        ax.grid(True, alpha=0.3)

        if title is not None:
            plt.title(f"{title} - {self.pdb_name}", fontdict=FONT_STYLES["title"])

        if equal_axis:
            ax.axis('equal')

        plt.tight_layout()
        plt.savefig(full_path, dpi=300)
        plt.close()

    def plot_graph_3d(self, graph, filename="graph_3D", title="3D Graph Visualization"):
        """ Creates a 3D plot of a graph using Plotly and saves it as an interactive HTML file. """

        full_path = self.get_full_save_path(filename=filename, extension="html")

        pos = {n: (graph.nodes[n]['x'], graph.nodes[n]['y'], graph.nodes[n]['z']) for n in graph.nodes}
        edge_x, edge_y, edge_z = [], [], []
        for u, v in graph.edges():
            x0, y0, z0 = pos[u]
            x1, y1, z1 = pos[v]
            edge_x += [x0, x1, None]
            edge_y += [y0, y1, None]
            edge_z += [z0, z1, None]

        edge_trace = go.Scatter3d(
            x=edge_x,
            y=edge_y,
            z=edge_z,
            mode='lines',
            line=dict(width=1, color='gray'),
            hoverinfo='none',
            showlegend=True,
            name="Edge"
        )

        chain_to_nodes = {}
        for node_id in graph.nodes():
            current_node = graph.nodes[node_id]
            chain_id = current_node['chain_id']
            chain_to_nodes.setdefault(chain_id, []).append(node_id)

        node_traces = []
        for i, (chain_id, node_ids) in enumerate(chain_to_nodes.items()):
            color = COLOR_PALETTE[i] if i < len(COLOR_PALETTE) else None
            node_x = [graph.nodes[node_id]['x'] for node_id in node_ids]
            node_y = [graph.nodes[node_id]['y'] for node_id in node_ids]
            node_z = [graph.nodes[node_id]['z'] for node_id in node_ids]
            hovertext = [
                (f"Chain: {graph.nodes[node_id][CHAIN_ID]}<br>"
                 f"Residue Number: {utils.get_residue_id(graph.nodes[node_id])}<br>"
                 f"Residue Name: {graph.nodes[node_id][RESIDUE_NAME]}<br>"
                 f"X: {graph.nodes[node_id]['x']:.3f}<br>"
                 f"Y: {graph.nodes[node_id]['y']:.3f}<br>Z: {graph.nodes[node_id]['z']:.3f}")
                for node_id in node_ids
            ]

            labels = [utils.get_residue_id(graph.nodes[n]) for n in node_ids]

            trace = go.Scatter3d(
                x=node_x, y=node_y, z=node_z,
                mode='markers+text',
                marker=dict(size=6, color=color),
                text=labels,
                textposition='top center',
                hoverinfo='text',
                hovertext=hovertext,
                name=f"Chain {chain_id}",
                legendgroup=f"Chain {chain_id}",
                showlegend=True
            )
            node_traces.append(trace)

        fig = go.Figure(data=[edge_trace] + node_traces)
        fig.update_layout(
            title=f"{title} - {self.pdb_name}",
            showlegend=True,
            margin=dict(l=0, r=0, b=0, t=40),
            scene=dict(
                xaxis=dict(title=X_COORD_LABEL),
                yaxis=dict(title=Y_COORD_LABEL),
                zaxis=dict(title=Z_COORD_LABEL),
                aspectmode='data'
            )
        )

        fig.write_html(full_path)

    def plot_graph_2d_degree(self, graph, filename="degree_graph", title=None, node_circle_size=80, equal_axis=False):

        full_path = self.get_full_save_path(filename=filename, extension="png")
        laplacian = nx.laplacian_matrix(graph).todense()
        eigenvalues, eigenvectors = eigs(laplacian, k=4, which='SM')
        w = eigenvectors[:, 0].real
        pos = {residue_index: (graph.nodes[residue_index]['x'], graph.nodes[residue_index]['y']) for residue_index in
               graph.nodes}
        plt.figure()
        node_color = ['r' if w[i] >= 0 else 'k' for i in range(len(w))]

        nx.draw_networkx_nodes(graph,
                               pos,
                               node_color=node_color,
                               node_size=node_circle_size)
        nx.draw_networkx_edges(graph,
                               pos,
                               edge_color=EDGE_COLOR,
                               width=0.5)

        ax = plt.gca()
        ax.set_xlabel(X_COORD_LABEL, fontdict=FONT_STYLES["xlabel"])
        ax.set_ylabel(Y_COORD_LABEL, fontdict=FONT_STYLES["ylabel"])
        ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True,
                       labelsize=FONT_STYLES["tick"]["labelsize"])
        ax.grid(True, alpha=0.3)

        if title is not None:
            plt.title(f"{title} - {self.pdb_name}", fontdict=FONT_STYLES["title"])

        if equal_axis:
            ax.axis('equal')

        plt.savefig(full_path, dpi=300)
        plt.close()

    def plot_graph_2d_degree_with_coords(self, graph, filename="degree_coords_graph", cmap='jet', title=None,
                                         equal_axis=False):
        full_path = self.get_full_save_path(filename=filename, extension="png")
        pos = {n: (graph.nodes[n]['x'], graph.nodes[n]['y']) for n in graph.nodes}

        plt.figure()

        deg = np.array([degree for node, degree in graph.degree()])

        node_sizes = 20 * np.sqrt(deg - np.min(deg) + 0.2)

        node_colors = deg

        nodes = nx.draw_networkx_nodes(graph,
                                       pos,
                                       cmap=cmap,
                                       node_color=node_colors,
                                       node_size=node_sizes)

        nx.draw_networkx_edges(graph,
                               pos,
                               edge_color=EDGE_COLOR,
                               width=0.5,
                               alpha=0.5)
        cbar = plt.colorbar(nodes)
        self._apply_colorbar_style(cbar)

        ax = plt.gca()
        ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True,
                       labelsize=FONT_STYLES["tick"]["labelsize"])
        ax.axis('on')

        ax.set_xlabel(X_COORD_LABEL, fontdict=FONT_STYLES["xlabel"])
        ax.set_ylabel(Y_COORD_LABEL, fontdict=FONT_STYLES["ylabel"])

        if title is not None:
            plt.title(f"{title} - {self.pdb_name}", fontdict=FONT_STYLES["title"])

        if equal_axis:
            ax.axis('equal')

        plt.savefig(full_path, dpi=300)
        plt.close()

    def plot_figures(self, scores_dict: dict, centrality_type_name: str, use_decimal: bool = True) -> None:
        """ Plots centrality scores for residues using actual residue numbers on the x-axis.

        Parameters:
        - scores_dict (dict): Mapping of residue index to centrality score.
        - centrality_type_name (str): Type of centrality (e.g., 'betweenness', 'closeness').
        - use_decimal (bool): Whether to format y-axis with decimal values.
        """
        full_path = self.get_full_save_path(filename=f"{centrality_type_name}_scores", extension="png")
        sorted_data = dict(sorted(scores_dict.items()))
        keys = list(sorted_data.keys())
        values = list(sorted_data.values())

        plt.figure(figsize=(10, 6))
        plt.plot(keys, values, linestyle='-')

        plt.ylim([0, max(values) * 1.05])

        if use_decimal:
            plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.2f}'))
        else:
            plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x)}'))

        margin = 0.05 * (max(keys) - min(keys))
        plt.xlim([min(keys) - margin, max(keys) + margin])

        ax = plt.gca()
        ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
        auto_x_ticks = ax.get_xticks().astype(int)
        visible_x_ticks = [tick for tick in auto_x_ticks if min(keys) <= tick <= max(keys)]
        if min(keys) not in visible_x_ticks:
            visible_x_ticks.insert(0, min(keys))
        if max(keys) not in visible_x_ticks:
            visible_x_ticks.append(max(keys))

        xtick_labels = [
            f"{i} ({self.actual_residue_number_map[i][0]}, {utils.get_residue_id_by_tuple(self.actual_residue_number_map[i])})"
            if i in self.actual_residue_number_map else i
            for i in visible_x_ticks
        ]
        plt.xticks(ticks=visible_x_ticks, labels=xtick_labels, rotation=45, fontsize=FONT_STYLES['xtick']['labelsize'],
                   fontname=FONT_FAMILY)
        plt.yticks(fontsize=FONT_STYLES['ytick']['labelsize'], fontname=FONT_FAMILY)

        plt.xlabel('Residue Index (Chain Id, Residue Number)', fontdict=FONT_STYLES["xlabel"])
        plt.ylabel(f'{centrality_type_name.capitalize()} Centrality Score', fontdict=FONT_STYLES["ylabel"])

        plt.title(f"PDB ID: {self.pdb_name}", fontdict=FONT_STYLES["title"])

        plt.grid(alpha=0.5)
        plt.tight_layout()

        plt.savefig(full_path, dpi=300)
        logging.info(f"{centrality_type_name.capitalize()} scores have been saved to {full_path}.")

    def plot_histogram(self, scores_dict: dict, centrality_type_name: str) -> None:
        """ Plots and saves a histogram of centrality scores.

            Parameters:
            -----------
            scores_dict : dict
                A dictionary where keys are residue indices and values are centrality scores.

            centrality_type_name : str
                The name of the centrality type (e.g., 'closeness', 'betweenness') to be used in axis labels and file naming.

            Description:
            ------------
            This method generates a histogram showing the frequency distribution of centrality scores
            provided in `scores_dict`. The histogram is saved as a PNG file with a filename based on
            `centrality_type_name`. Axis labels, tick formatting, and layout are adjusted for clarity.
        """

        full_path = self.get_full_save_path(filename=f"{centrality_type_name}_scores_histogram", extension="png")

        sorted_data = dict(sorted(scores_dict.items()))

        values = list(sorted_data.values())

        plt.figure(figsize=(10, 6))

        hist_plot = sns.histplot(values, kde=False, bins=10)

        plt.xlabel(f'{centrality_type_name.capitalize()} Centrality', fontdict=FONT_STYLES["xlabel"])
        plt.ylabel('Frequency', fontdict=FONT_STYLES["xlabel"])

        plt.xlim(left=0, right=np.ceil(max(values) * 10) / 10)
        max_y_value = hist_plot.get_ylim()[1]
        plt.ylim(bottom=0, top=np.ceil(max_y_value))
        plt.title(f"PDB ID: {self.pdb_name}", fontdict=FONT_STYLES["title"])
        plt.tight_layout()
        plt.savefig(full_path, dpi=300)
        plt.close()
        logging.info(f"{centrality_type_name.capitalize()} scores has been saved to {full_path}.")


def main():
    pass


if __name__ == '__main__':
    main()
