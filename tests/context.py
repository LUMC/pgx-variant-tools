import os
import sys
sys.path.insert(1, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from variant_tools import alignment, variant_calling, mutalyzer, utils, vcf_tools, str_search
