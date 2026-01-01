# -*- coding: utf-8 -*-
"""
    __author__ = 'Zehra Sarica'
    __email__ = ['sarica16@itu.edu.tr','zehraacar559@gmail.com']
"""

import logging
import os
from importlib import resources

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt

_fonts_loaded = False


def load_fonts_once():
    """ Load packaged Times New Roman fonts exactly once per interpreter. """
    global _fonts_loaded
    if _fonts_loaded:
        return

    try:
        with resources.as_file(
                resources.files("rinpy") / "fonts"
        ) as font_dir:
            for font_file in os.listdir(font_dir):
                if font_file.endswith(".ttf"):
                    fm.fontManager.addfont(
                        os.path.join(font_dir, font_file)
                    )

        plt.rcParams["font.family"] = "Times New Roman"
        logging.info(f"Active font: {plt.rcParams['font.family'][0]}")

    except Exception as e:
        plt.rcParams["font.family"] = "DejaVu Serif"
        logging.warning(f"Could not load Times font ({e}). Using default {plt.rcParams['font.family'][0]}.")

    _fonts_loaded = True
