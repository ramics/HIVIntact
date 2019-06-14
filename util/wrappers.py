import os
import subprocess
import sys
import uuid

from Bio import SeqIO, AlignIO

def mafft(working_dir, sequences):
    '''
    Call mafft on a set of sequences and return the resulting alignment.

    Args:
        sequences: Sequences to be aligned

    Returns:
        An AlignIO object
    '''

    alignment_input = os.path.join(working_dir, str(uuid.uuid4()) + ".fasta")
    alignment_output = os.path.join(working_dir, str(uuid.uuid4()) + ".fasta")

    SeqIO.write(sequences, alignment_input, "fasta")

    command = "mafft " + alignment_input + " > " + alignment_output

    # TODO: figure out how to make this a logging call - not currently
    # behaving, hence the pipe above.
    devnull = open(os.devnull, 'w')
    subprocess.call(command, shell=True, stderr=devnull)
    devnull.close()

    alignment = AlignIO.read(alignment_output, "fasta")

    os.unlink(alignment_input)
    os.unlink(alignment_output)

    return alignment