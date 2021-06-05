import pkg_resources
import os

# Constants for intactness
# In format start, end, maximum_allowed_deletion
# HXB2 coordinates, must be converted for other subtypes

# All coordinates zero-numbered


# (start, end, error_bar, breaks_intactness)
DEFAULT_FORWARD_ORFs = [("gag", 790, 2291, 30), 
                        ("pol", 2085, 5095, 30), 
                        ("env", 6225, 8793, 100)]
DEFAULT_REVERSE_ORFS = []
DEFAULT_SMALL_FORWARD_ORFS = [("vif", 5041, 5619, 30),   
                      ("vpr", 5559, 5850, 30), 
                      ("tat_exon1", 5831, 6045, 30),
                      ("rev_exon1", 5970, 6045, 30),
                      ("vpu", 6062, 6310, 30),
                      ("tat_exon2", 8379, 8469, 30),
                      ("rev_exon2", 8379, 8653, 30),
                      ("nef", 8797, 9417, 30)]
DEFAULT_SMALL_REVERSE_ORFS = []
                       
DEFAULT_ORF_LENGTH = 1000
DEFAULT_SMALL_ORF_LENGTH = 100
DEFAULT_ERROR_BAR = 1

# Major Splice Donor Site Location -- HXB2 coordinates only
DEFAULT_MSD_SITE_LOCUS = 743
DEFAULT_MSD_SEQUENCE = "GT"

# Packaging Signal Location -- HXB2 coordinates only
DEFAULT_PSI_LOCUS = (680, 809)
PSI_ERROR_TOLERANCE = 10 #10 Nucleotide deletion tolerance

# Rev Response Element Location -- HXB2 coordinates only
DEFAULT_RRE_LOCUS = (7755, 8020)
RRE_ERROR_TOLERANCE = 20 


