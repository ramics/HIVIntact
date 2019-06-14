import os
import pkg_resources
from Bio import Seq, SeqIO, SeqRecord

import util.wrappers as wrappers


REFERENCE_DIR = pkg_resources.resource_filename('util', 'subtype_alignments/')

def HXB2():
    """
    Return the sequence of HXB2, the standard HIV reference sequence
    """
    return SeqIO.read(os.path.join(REFERENCE_DIR, "HXB2.fasta"), "fasta")

def subtypes():
	"""
    List all currently available HIV subtypes
    """
	return [f.replace(".fasta", "") for f in os.listdir(REFERENCE_DIR)]

def alignment_file(subtype):
	"""
    Return an alignment file associated with an HIV subtype.

    Args:
        subtype: folder in which to put temporary files.
    """
	return os.path.join(REFERENCE_DIR, subtype + ".fasta")

def subtype_sequence(subtype):
	"""
    Return an example sequence associated with an HIV subtype.

    Args:
        subtype: folder in which to put temporary files.
    """
	alignment = list(SeqIO.parse(os.path.join(REFERENCE_DIR, subtype + ".fasta"), "fasta"))
	return SeqRecord.SeqRecord(
        Seq.Seq(str(alignment[0].seq).replace("-","").replace("\n", "")),
        id = alignment[0].id,
        name = alignment[0].name
        )

def convert_from_hxb2_to_subtype(working_dir, position, subtype):
    """
    Convert a position number in HXB2 to the equivalent in another subtype.

    Args:
        working_dir: working folder in which to place temporary files
        position: hxb2 coordinate position to convert
        subtype: subtype position to convert to
    """

    sequences = [HXB2(), subtype_sequence(subtype)]

    alignment = wrappers.mafft(working_dir, sequences)

    hxb2_pos = 0
    subtype_pos = 0
    for i in range(len(alignment[0])):
        if hxb2_pos == position:
            return subtype_pos
        if alignment[0][i] != "-":
            hxb2_pos += 1
        if alignment[1][i] != "-":
            subtype_pos += 1


