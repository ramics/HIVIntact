# BioPython and pysam, while excellent for a lot of things, aren't great
# when you know your data is well-formed and want to do something very simple
# with it: they're very clunky.  Hence, much simpler iterators

try:
    from itertools import zip_longest as zip_longest
except:
    from itertools import izip_longest as zip_longest

def parse_fastq(handle):
    """
    Convert a file handle to a simple 4-lines-at-a-time fastq iterator.

    Args:
        handle: The file handle.

    Returns:
        The iterator.
    """
    return zip_longest(*[handle]*4)

# order of elements in the list from the sam specification
QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, \
RNEXT, PNEXT, TLEN, SEQ, QUAL = range(11)

def parse_sam(handle):
    """
    Convert a file handle to a list-like sam iterator.

    Args:
        handle: The file handle.

    Returns:
        The iterator.
    """
    return (h.split('\t') for h in handle if not h.startswith('@'))