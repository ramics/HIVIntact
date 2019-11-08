import pkg_resources
import os

# Constants for intactness
# In format start, end, maximum_allowed_deletion
# HXB2 coordinates, must be converted for other subtypes

# All coordinates zero-numbered


DEFAULT_FORWARD_ORFs = [(790, 2291, 30), (2085, 5095, 30), (6225, 8793, 100)]
DEFAULT_REVERSE_ORFS = []
DEFAULT_ORF_LENGTH = 1000
DEFAULT_ERROR_BAR = 0

# Major Splice Donor Site Location -- HXB2 coordinates only
DEFAULT_MSD_SITE_LOCUS = 743
DEFAULT_MSD_SEQUENCE = "GT"

# Packaging Signal Location -- HXB2 coordinates only
DEFAULT_PSI_LOCUS = (680, 809)
PSI_ERROR_TOLERANCE = 8 #8 Nucleotide deletion tolerance

# Rev Response Element Location -- HXB2 coordinates only
DEFAULT_RRE_LOCUS = (7755, 8020)
RRE_ERROR_TOLERANCE = 20 


