# Admissible pre-release markers:
#   X.YaN   # Alpha release
#   X.YbN   # Beta release
#   X.YrcN  # Release Candidate
#   X.Y     # Final release
#
# Dev branch marker is: 'X.Y.dev' or 'X.Y.devN' where N is an integer.
# 'X.Y.dev0' is the canonical version of 'X.Y.dev'
#
import logging
import os
from importlib import resources

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt

__version__ = "0.0.1dev"
__author__ = "Zehra Sarica <zehraacar559@gmail.com, sarica16@itu.edu.tr>"

import rinpy.logging_config

from rinpy.rin_process import RINProcess

__all__ = ["RINProcess"]


def _load_fonts():
    """Load packaged Times New Roman fonts dynamically."""
    try:
        with resources.path("rinpy.fonts", "times.ttf") as f:
            font_dir = os.path.dirname(f)
            for font_file in os.listdir(font_dir):
                if font_file.endswith(".ttf"):
                    fm.fontManager.addfont(os.path.join(font_dir, font_file))

        plt.rcParams["font.family"] = "Times New Roman"
        logging.info(f"Active font: {plt.rcParams['font.family'][0]}")
    except Exception as e:
        print(f"Warning: Could not load Times New Roman fonts ({e}). Using default font.")
        plt.rcParams["font.family"] = "DejaVu Serif"
        logging.warning(f"Could not load Times font ({e}). Using default {plt.rcParams['font.family'][0]}.")


_load_fonts()
