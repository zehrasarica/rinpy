# -*- coding: utf-8 -*-
"""
    __author__ = 'Zehra Sarica'
    __credits__ = ''
    __email__ = ['sarica16@itu.edu.tr','zehraacar559@gmail.com']
"""

import os.path
import re

import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import seaborn as sns

from rinpy.constants import FREQUENCY_HIGH_PERCENTAGE_TEMPLATE, CHAIN_ID, RESIDUE_NUMBER, RESIDUE_NAME, PDB_ID, \
    CENTRALITY_SCORE, INSERTION, EXCEL_EXT, CSV_EXT, HTML_EXT
from rinpy.style_config import FONT_STYLES, FONT_FAMILY
from rinpy.utils import CentralityType

RESIDUE_ID = "residue_id"
OCCURRENCES = "occurrences"
PDBS = "pdbs"
FREQUENCY = "frequency"
PDB_IDS = 'pdb_ids'
COUNT = 'count'

_COLONS = [PDB_ID, RESIDUE_NAME, CHAIN_ID, RESIDUE_NUMBER, INSERTION, CENTRALITY_SCORE]


class QuantileAnalyzer:
    def __init__(self, high_percentage_dict=None, centrality_type=CentralityType.BET, destination_output_path=None):
        self.high_percentage_dict = high_percentage_dict
        self.centrality_type = centrality_type
        self.destination_output_path = destination_output_path
        self.outfile_df = pd.DataFrame(columns=_COLONS)
        self.merge_files_and_convert_to_df()

    def merge_files_and_convert_to_df(self):
        file_row_list = []
        for pdb_id, out_file_path in self.high_percentage_dict.items():
            file_out_df = pd.read_csv(out_file_path,
                                      sep=";",
                                      names=[RESIDUE_NAME, CHAIN_ID, RESIDUE_NUMBER, INSERTION, CENTRALITY_SCORE])
            for _, row in file_out_df.iterrows():
                file_row_list.append({
                    PDB_ID: pdb_id,
                    RESIDUE_NAME: row[RESIDUE_NAME],
                    CHAIN_ID: row[CHAIN_ID],
                    RESIDUE_NUMBER: row[RESIDUE_NUMBER],
                    INSERTION: str(row[INSERTION]).strip() if pd.notna(row[INSERTION]) else "",
                    CENTRALITY_SCORE: row[CENTRALITY_SCORE]
                })

        new_data_df = pd.DataFrame(file_row_list)
        self.outfile_df = pd.concat([self.outfile_df, new_data_df], ignore_index=True)

    def _find_frequency(self):
        """Finds the frequency of each residue along with its name, number, chain id, and insertion code"""
        return self.outfile_df.groupby([RESIDUE_NAME, RESIDUE_NUMBER, CHAIN_ID, INSERTION]).size()

    def run_analysis(self):
        """Runs the analysis to find frequencies and save the results to a csv file and plot bar plot and heatmaps"""
        fre_df = self._find_frequency()
        sorted_out_file_list = self._write_to_file(fre_df)
        self._save_bar_plot(sorted_out_file_list)
        self._save_heatmap(sorted_out_file_list)
        self._save_heatmap_interactive(sorted_out_file_list)

    def _write_to_file(self, fre_df):
        fre_df = fre_df.sort_values(ascending=False)
        out_file_list = []
        for index, value in fre_df.items():
            row_df = self.outfile_df[
                (self.outfile_df[RESIDUE_NAME] == index[0]) & (self.outfile_df[RESIDUE_NUMBER] == index[1]) & (
                        self.outfile_df[CHAIN_ID] == index[2]) & (
                        self.outfile_df[INSERTION] == index[3])]
            pdb_ids = []
            for r in row_df.iterrows():
                pdb_ids.append(r[1][PDB_ID])

            out_file_list.append({
                RESIDUE_NAME: index[0],
                RESIDUE_NUMBER: index[1],
                CHAIN_ID: index[2],
                INSERTION: index[3],
                COUNT: value,
                PDB_IDS: pdb_ids
            })

        sorted_out_file_list = sorted(out_file_list, key=lambda x: (x[CHAIN_ID], x[RESIDUE_NAME], x[COUNT]))

        out_file_df = pd.DataFrame(sorted_out_file_list)
        out_file_df[PDB_IDS] = out_file_df[PDB_IDS].astype(str).str.replace("'", '')
        output_filename = os.path.join(self.destination_output_path, FREQUENCY_HIGH_PERCENTAGE_TEMPLATE.format(
            type=self.centrality_type.display_name))

        self._export_to_frequency_with_pdb_ids_to_csv_file(out_file_df=out_file_df, output_filename=output_filename)

        out_file_df[PDB_IDS] = out_file_df[PDB_IDS].astype(str).str.replace('[\]\[]', '')

        self._write_to_excel(out_file_df=out_file_df, output_filename=output_filename)

        self._export_to_frequency_without_pdb_ids_to_csv_file(fre_df=fre_df, output_filename=output_filename)

        return sorted_out_file_list

    @staticmethod
    def _write_to_excel(out_file_df, output_filename):
        """Writes dataframe to the Excel file."""
        df_copy = out_file_df.copy()
        headers = [RESIDUE_NAME, RESIDUE_NUMBER, CHAIN_ID, INSERTION, FREQUENCY, PDBS]
        df_copy.columns = headers
        out_file_df.to_excel(f"{output_filename}{EXCEL_EXT}", index=False)

    @staticmethod
    def _export_to_frequency_with_pdb_ids_to_csv_file(out_file_df, output_filename):
        max_residue_name_len = out_file_df[RESIDUE_NAME].astype(str).map(len).max()
        max_residue_num_len = out_file_df[RESIDUE_NUMBER].astype(str).map(len).max()
        max_chain_id_len = out_file_df[CHAIN_ID].astype(str).map(len).max()
        max_insertion_len = out_file_df[INSERTION].astype(str).map(len).max()
        max_count_len = out_file_df[COUNT].astype(str).map(len).max()

        out_file_df[RESIDUE_NAME] = out_file_df[RESIDUE_NAME].astype(str).map(
            lambda x: x.ljust(max_residue_name_len))
        out_file_df[RESIDUE_NUMBER] = out_file_df[RESIDUE_NUMBER].astype(str).map(
            lambda x: x.rjust(max_residue_num_len))
        out_file_df[CHAIN_ID] = out_file_df[CHAIN_ID].astype(str).map(lambda x: x.ljust(max_chain_id_len))
        out_file_df[INSERTION] = out_file_df[INSERTION].astype(str).map(lambda x: x.ljust(max_insertion_len))
        out_file_df[COUNT] = out_file_df[COUNT].astype(str).map(lambda x: x.rjust(max_count_len))

        out_file_df.to_csv(
            path_or_buf=f"{output_filename}{CSV_EXT}",
            sep='\t',
            header=False,
            index=False)

    @staticmethod
    def _export_to_frequency_without_pdb_ids_to_csv_file(fre_df, output_filename):
        index_df = fre_df.index.to_frame(index=False)

        max_residue_name_len = index_df.iloc[:, 0].astype(str).map(len).max()
        max_residue_num_len = index_df.iloc[:, 1].astype(str).map(len).max()
        max_chain_id_len = index_df.iloc[:, 2].astype(str).map(len).max()
        max_insertion_len = index_df.iloc[:, 3].astype(str).map(len).max()
        max_count_len = fre_df.astype(str).map(len).max()

        index_df.iloc[:, 0] = index_df.iloc[:, 0].astype(str).str.ljust(max_residue_name_len)
        index_df.iloc[:, 1] = index_df.iloc[:, 1].astype(str).str.rjust(max_residue_num_len)
        index_df.iloc[:, 2] = index_df.iloc[:, 2].astype(str).str.ljust(max_chain_id_len)
        index_df.iloc[:, 3] = index_df.iloc[:, 3].astype(str).str.ljust(max_insertion_len)
        counts_padded = fre_df.astype(str).str.rjust(max_count_len)

        final_df = pd.concat([index_df, counts_padded.reset_index(drop=True)], axis=1)

        final_df.to_csv(
            path_or_buf=f"{output_filename}_sorted_without_pdb_ids.csv",
            sep='\t',
            header=False,
            index=False)

    def _save_bar_plot(self, sorted_out_file_list):

        unique_chain_ids = sorted(set(entry[CHAIN_ID] for entry in sorted_out_file_list))
        colors = plt.get_cmap('jet')(np.linspace(0, 1, len(unique_chain_ids)))
        residue_counts = {chain_id: {} for chain_id in unique_chain_ids}
        for chain_id in unique_chain_ids:
            chain_entries = [entry for entry in sorted_out_file_list if entry[CHAIN_ID] == chain_id]
            chain_entries = sorted(chain_entries, key=lambda x: (x[RESIDUE_NUMBER], x[INSERTION]))

            for entry in chain_entries:
                residue_number = entry[RESIDUE_NUMBER]
                insertion = entry[INSERTION]
                count = entry[COUNT]

                key = (residue_number, insertion)

                if key not in residue_counts[chain_id]:
                    residue_counts[chain_id][key] = 0

                residue_counts[chain_id][key] += count

        unique_residue_numbers = set()
        for counts in residue_counts.values():
            unique_residue_numbers.update(counts.keys())

        unique_residue_numbers = sorted(unique_residue_numbers)

        x_positions = np.arange(len(unique_residue_numbers))

        bar_width = 0.2

        fig, ax = plt.subplots(figsize=(12, 6))

        for i, (chain_id, color) in enumerate(zip(unique_chain_ids, colors)):
            counts = [residue_counts[chain_id].get(r, 0) for r in unique_residue_numbers]
            x_positions_adjusted = x_positions + i * bar_width
            ax.bar(x_positions_adjusted, counts, width=bar_width, label=f"Chain {chain_id}", color=color)

        x_labels = [
            f"{num}{ins}" if ins != "''" else str(num)
            for num, ins in unique_residue_numbers
        ]

        ax.set_xticks(x_positions)
        ax.set_xticklabels(x_labels, rotation=90, **{
            "fontsize": FONT_STYLES["xtick_small"]["labelsize"],
            "fontfamily": FONT_STYLES["xtick"]["fontfamily"]
        })

        ax.set_xlabel("Residue Number", fontdict=FONT_STYLES["xlabel"])
        ax.set_ylabel("Frequency of High Percentage Residues", fontdict=FONT_STYLES["xlabel"])
        ax.legend(loc="upper right")
        ax.legend(prop={"family": FONT_FAMILY, "size": FONT_STYLES["legend"]["fontsize"]})

        plt.tight_layout()

        plot_output_filename = os.path.join(self.destination_output_path,
                                            f"frequency_high_percentage_{self.centrality_type.display_name}_bar_plot")
        plt.savefig(plot_output_filename, dpi=300)

    def _save_heatmap(self, sorted_out_file_list):

        total_pdb_ids = len(set(pdb_id for entry in sorted_out_file_list for pdb_id in entry[PDB_IDS]))

        data = [[entry[RESIDUE_NAME], entry[RESIDUE_NUMBER], entry[CHAIN_ID], entry[INSERTION], entry[COUNT],
                 entry[PDB_IDS]]
                for entry in sorted_out_file_list]

        df = pd.DataFrame(data, columns=[RESIDUE_NAME, RESIDUE_NUMBER, CHAIN_ID, INSERTION, FREQUENCY, OCCURRENCES])
        df[RESIDUE_NUMBER] = df[RESIDUE_NUMBER].astype(int)
        df[FREQUENCY] = df[FREQUENCY].astype(int)
        df[RESIDUE_ID] = df.apply(
            lambda row: f"{row[RESIDUE_NUMBER]}{row[INSERTION]}" if row[INSERTION] else str(row[RESIDUE_NUMBER]),
            axis=1)

        heatmap_data = df.pivot_table(index=CHAIN_ID, columns=RESIDUE_ID, values=FREQUENCY, aggfunc='sum', fill_value=0)

        heatmap_data = heatmap_data.sort_index(ascending=False)
        heatmap_data = heatmap_data[sorted(heatmap_data.columns, key=lambda x: (
            int(''.join(filter(str.isdigit, x))), ''.join(filter(str.isalpha, x))))]

        heatmap_data = heatmap_data.apply(pd.to_numeric, errors='coerce')

        max_freq_data = int(np.nanmax(heatmap_data.values))

        max_freq = max(max_freq_data, total_pdb_ids)

        if max_freq <= 10:
            tick_vals = np.arange(0, max_freq + 1)
        else:
            step = max(1, round(max_freq / 10))
            tick_vals = list(range(0, max_freq + 1, step))
            if max_freq not in tick_vals:
                tick_vals.append(max_freq)

        custom_cmap = sns.diverging_palette(240, 10, as_cmap=True)

        plt.figure(figsize=(14, 8))
        sns_heatmap = sns.heatmap(
            heatmap_data,
            cmap=custom_cmap,
            annot=True,
            fmt="d",
            annot_kws={"size": 8},
            linewidths=0.5,
            linecolor="white",
            cbar_kws={
                "ticks": tick_vals,
                "label": "Frequency"
            },
            vmin=0,
            vmax=max_freq
        )
        for text in sns_heatmap.texts:
            text.set_rotation(90)

        sns_heatmap.set_xticks(np.arange(len(heatmap_data.columns)) + 0.5)

        formatted_columns = [
            f"{res_num}{ins}" if ins else str(res_num)
            for col in heatmap_data.columns
            for res_num, ins in [self._parse_residue_string(col)]
        ]
        sns_heatmap.set_xticklabels(
            formatted_columns,
            rotation=90,
            fontsize=FONT_STYLES['xtick_small_small']['labelsize'],
            fontname=FONT_FAMILY
        )

        cbar = sns_heatmap.collections[0].colorbar
        self.apply_colorbar_style(cbar, cbar_label="Frequency")

        plt.title(f"Residue Frequency Heatmap ({self.centrality_type.display_name_capitalized()})",
                  fontdict=FONT_STYLES["title"])
        plt.xlabel("Residue Number", fontdict=FONT_STYLES["xlabel"])
        plt.ylabel("Chain ID", fontdict=FONT_STYLES["ylabel"])
        plt.xticks(rotation=90, fontsize=FONT_STYLES['xtick_medium']['labelsize'], fontname=FONT_FAMILY)
        plt.yticks(fontsize=FONT_STYLES['xtick']['labelsize'], fontname=FONT_FAMILY)
        plt.tight_layout()
        plot_output_filename = os.path.join(self.destination_output_path,
                                            f"frequency_high_percentage_{self.centrality_type.display_name}_heatmap")
        plt.savefig(plot_output_filename, dpi=300)

    @staticmethod
    def _parse_residue_string(s: str):
        s = s.strip()
        if s.endswith("''"):
            residue_number = int(s[:-2])
            insertion = ''
        else:
            match = re.fullmatch(r"(\d+)([A-Za-z]?)", s)
            if not match:
                raise ValueError(f"Invalid residue string format: {s}")
            residue_number = int(match.group(1))
            insertion = match.group(2)
        return residue_number, insertion

    def _save_heatmap_interactive(self, sorted_out_file_list):

        total_pdb_ids = len(set(pdb_id for entry in sorted_out_file_list for pdb_id in entry[PDB_IDS]))

        data = [[entry[RESIDUE_NAME], entry[RESIDUE_NUMBER], entry[CHAIN_ID], entry[INSERTION], entry[COUNT],
                 entry[PDB_IDS]]
                for entry in sorted_out_file_list]

        df = pd.DataFrame(data, columns=[RESIDUE_NAME, RESIDUE_NUMBER, CHAIN_ID, INSERTION, FREQUENCY, OCCURRENCES])
        df[RESIDUE_NUMBER] = df[RESIDUE_NUMBER].astype(str)
        df[FREQUENCY] = df[FREQUENCY].astype(int)
        df[RESIDUE_ID] = df.apply(
            lambda row: f"{row[RESIDUE_NUMBER]}{row[INSERTION]}" if row[INSERTION] else str(row[RESIDUE_NUMBER]),
            axis=1)

        heatmap_data = df.pivot_table(index=CHAIN_ID, columns=RESIDUE_ID, values=FREQUENCY, aggfunc='sum', fill_value=0)

        heatmap_data = heatmap_data.sort_index(ascending=False)
        heatmap_data = heatmap_data[sorted(heatmap_data.columns, key=lambda x: (
            int(''.join(filter(str.isdigit, x))), ''.join(filter(str.isalpha, x))))]

        agg_df = df.groupby([CHAIN_ID, RESIDUE_NUMBER, INSERTION]).agg({
            RESIDUE_NAME: lambda x: ', '.join(sorted(set(x))),
            OCCURRENCES: lambda x: ', '.join(sorted(set(
                pdb.strip()
                for sublist in x
                for pdb in str(sublist).strip("[]").replace('"', '').split(',')
            )))
        }).reset_index()

        agg_df[RESIDUE_ID] = agg_df.apply(
            lambda row: f"{row[RESIDUE_NUMBER]}{row[INSERTION]}" if row[INSERTION] else str(row[RESIDUE_NUMBER]),
            axis=1)

        residue_matrix = agg_df.pivot(index=CHAIN_ID, columns=RESIDUE_ID, values=RESIDUE_NAME).reindex(
            index=heatmap_data.index, columns=heatmap_data.columns, fill_value="")

        pdb_matrix = agg_df.pivot(index=CHAIN_ID, columns=RESIDUE_ID, values=OCCURRENCES).reindex(
            index=heatmap_data.index, columns=heatmap_data.columns, fill_value="")

        def wrap_pdb_ids(pdb_str, wrap_after=8):
            if not isinstance(pdb_str, str):
                return "None"
            ids = [p.strip() for p in pdb_str.strip("[]").replace('"', '').split(",") if p.strip()]
            if not ids:
                return "None"
            chunks = [', '.join(ids[i:i + wrap_after]) for i in range(0, len(ids), wrap_after)]
            return "<br>".join(chunks)

        custom_hovertext = heatmap_data.copy().astype(str)
        for row in heatmap_data.index:
            for col in heatmap_data.columns:
                resname = residue_matrix.loc[row, col]
                pdbs_raw = pdb_matrix.loc[row, col]
                freq = heatmap_data.loc[row, col]

                resname = resname if pd.notna(resname) else "None"
                pdbs_wrapped = wrap_pdb_ids(pdbs_raw) if pd.notna(pdbs_raw) else "None"
                residue_number, insertion = self._parse_residue_string(col)
                custom_hovertext.loc[row, col] = (
                    f"<b>Residue Number</b>: {residue_number}{insertion}<br>"
                    f"<b>Chain ID</b>: {row}<br>"
                    f"<b>Residue Name(s)</b>: {resname}<br>"
                    f"<b>Frequency</b>: {freq}<br>"
                    f"<b>PDB IDs</b>: {pdbs_wrapped}"
                )

        custom_cmap = sns.diverging_palette(240, 10, as_cmap=True)
        custom_plotly_colorscale = self.matplotlib_to_plotly(custom_cmap)
        formatted_columns = [
            f"{res_num}{ins}" if ins else str(res_num)
            for col in heatmap_data.columns
            for res_num, ins in [self._parse_residue_string(col)]
        ]
        fig = px.imshow(
            heatmap_data,
            x=formatted_columns,
            color_continuous_scale=custom_plotly_colorscale,
            aspect='auto',
            labels=dict(color='Frequency'),
            title=f"Residue Frequency Heatmap ({self.centrality_type.display_name_capitalized()})",
            text_auto=True,
        )

        fig.update_layout(
            xaxis_title="Residue Number",
            yaxis_title="Chain ID",
            xaxis=dict(
                tickangle=45,
                tickmode='linear',
                dtick=1,
                tickfont=dict(size=10),
            ),
            font=dict(family=FONT_FAMILY, size=14),
        )

        max_freq_data = int(df["frequency"].max())

        max_freq = max(max_freq_data, total_pdb_ids)

        if max_freq <= 10:
            tick_vals = list(range(0, max_freq + 1))
        else:
            step = max(1, round(max_freq / 10))
            tick_vals = list(range(0, max_freq + 1, step))
            if max_freq not in tick_vals:
                tick_vals.append(max_freq)

        fig.update_coloraxes(
            colorbar=dict(
                tickmode='array',
                tickvals=tick_vals,
                ticktext=[str(i) for i in tick_vals],
                title=dict(
                    text='Frequency',
                    side='right',
                    font=dict(size=FONT_STYLES['legend']['title_fontsize'],
                              family=FONT_STYLES['legend']['fontfamily']),
                ),
            ),
            cmax=max_freq,
            cmin=0
        )

        fig.update_traces(hovertemplate="%{customdata}", customdata=custom_hovertext.values)

        plot_output_filename = os.path.join(
            self.destination_output_path,
            f"frequency_high_percentage_{self.centrality_type.display_name}_heatmap{HTML_EXT}"
        )
        config = {
            "scrollZoom": True,
            "doubleClick": "reset",
            "displayModeBar": True,
            "modeBarButtonsToAdd": ["pan2d", "zoomIn2d", "zoomOut2d"]
        }
        fig.update_layout(dragmode='pan')

        fig.write_html(plot_output_filename, config=config)

    @staticmethod
    def matplotlib_to_plotly(cmap, n_colors=256):
        norm = mcolors.Normalize(vmin=0, vmax=1)
        scalar_map = cm.ScalarMappable(norm=norm, cmap=cmap)
        return [
            [i / (n_colors - 1), mcolors.to_hex(scalar_map.to_rgba(i / (n_colors - 1)))]
            for i in range(n_colors)
        ]

    @staticmethod
    def apply_colorbar_style(cbar, cbar_label=None):
        """Apply consistent style to colorbar ticks and label."""
        styles = FONT_STYLES["colorbar"]
        cbar.ax.tick_params(labelsize=styles["tick_params"]["labelsize"])
        for tick in cbar.ax.get_yticklabels():
            tick.set_fontfamily(styles["tick_params"].get("fontfamily", "sans-serif"))
        if cbar_label is not None:
            cbar.set_label(
                cbar_label,
                fontsize=styles["label"]["fontsize"],
                fontfamily=styles["label"]["fontfamily"])


def main():
    pass


if __name__ == '__main__':
    main()
