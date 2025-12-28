# Admissible pre-release markers:
#   X.YaN   # Alpha release
#   X.YbN   # Beta release
#   X.YrcN  # Release Candidate
#   X.Y     # Final release
#
# Dev branch marker is: 'X.Y.dev' or 'X.Y.devN' where N is an integer.
# 'X.Y.dev0' is the canonical version of 'X.Y.dev'
#

__version__ = "0.0.1dev"
__author__ = "Zehra Sarica <zehraacar559@gmail.com, sarica16@itu.edu.tr>"

import rinpy.logging_config

from rinpy.rin_process import RINProcess

__all__ = ["RINProcess"]
