import pkg_resources
import os

# Constants for intactness
# In format start, end, maximum_allowed_deletion
# HXB2 coordinates, must be converted for other subtypes

DEFAULT_FORWARD_ORFs = [(790, 2292, 30), (2085, 5096, 30), (6225, 8795, 100)]
DEFAULT_REVERSE_ORFS = []
DEFAULT_ORF_LENGTH = 1000
DEFAULT_ERROR_BAR = 100

# Packaging Signal Location -- HXB2 coordinates only
DEFAULT_PSI_LOCUS = (681, 810)
PSI_ERROR_TOLERANCE = 8 #8 Nucleotide deletion tolerance

# Rev Response Element Location -- HXB2 coordinates only
DEFAULT_RRE_LOCUS = (7756, 8021)
RRE_ERROR_TOLERANCE = 20 

