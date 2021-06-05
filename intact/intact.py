import json
import os
import pkg_resources
import re
import subprocess
import sys
import uuid
from Bio import AlignIO, pairwise2, Seq, SeqIO, SeqRecord

import util.constants as const
import util.subtypes as st
import util.wrappers as wrappers


START_CODON = 'ATG'

WRONGORFNUMBER_ERROR = "WrongORFNumber"
MISPLACEDORF_ERROR   = "MisplacedORF"
LONGDELETION_ERROR   = "LongDeletion"
DELETIONINORF_ERROR  = "DeletionInOrf"
FRAMESHIFTINORF_ERROR  = "FrameshiftInOrf"
MSDMUTATED_ERROR = "MajorSpliceDonorSiteMutated"
PSIDELETION_ERROR    = "PackagingSignalDeletion"
PSINOTFOUND_ERROR    = "PackagingSignalNotComplete"
RREDELETION_ERROR    = "RevResponseElementDeletion"

class IntactnessError:
    def __init__(self, sequence_name, error, message):
        self.sequence_name = sequence_name
        self.error = error
        self.message = message

class ORF:
    def __init__(self, orientation, start, end):
        self.orientation = orientation
        self.start = start
        self.end = end

def has_mutated_major_splice_donor_site(alignment, 
                                        splice_donor_start_pos, 
                                        splice_donor_end_pos,
                                        splice_donor_sequence):
    """
    Determines whether the major splice donor site is mutated.
    Keyword Args:
        
        alignment -- multiple sequence alignment object containing the 
                     reference and query sequence.
        splice_donor_start_pos -- first position of splice donor site
        splice_donor_end_pos -- last position of splice donor site
        splice_donor_sequence - sequence of splice donor site
    """

    sd_begin = [m.start() for m in re.finditer(r"[^-]",
                str(alignment[0].seq))][splice_donor_start_pos]
    
    sd_end = [m.start() for m in re.finditer(r"[^-]",
              str(alignment[0].seq))][splice_donor_end_pos]
    
    sd = alignment[1].seq[sd_begin:(sd_end + 1)]

    # splice donor site is missing from sequence
    if all([x == "-" for x in sd]):
        return IntactnessError(
                alignment[1].id, MSDMUTATED_ERROR,
                "Query sequence has a missing splice donor site, " 
                + "".join(sd.upper()) + "."
                )

    if sd.upper() != splice_donor_sequence.upper():

        return IntactnessError(
                alignment[1].id, MSDMUTATED_ERROR,
                "Query sequence has a mutated splice donor site, " 
                + "".join(sd.upper()) + "."
                )

    return None

    

    
    

def has_packaging_signal(alignment, psi_locus, psi_tolerance):
    """
    Determines presence and possible intactness of HIV 
    Packaging Signal Region.
    
    
    Keyword Args:
        
        alignment -- multiple sequence alignment object containing the 
                     reference and query sequence.
        psi_locus -- tuple containing start and end coordinates of PSI wrt
                     the reference being used.
        psi_tolerance -- number of deletions in query tolerated to be intact
        
        
    Internal Args:
        
        packaging_begin -- Aligned PSI start position.
        packaging_end -- Aligned PSI end position.
        query_start -- beginning position of query sequence.
        query_psi -- extracted PSI region from query
        query_psi_deletions -- number of deletions in query PSI region
        
        
    Return:
        
        PSINOTFOUND_ERROR -- IntactnessError denoting query does not encompass
                             the complete PSI region.
        PSIDELETION_ERROR -- IntactnessError denoting query likely contains
                             defective PSI.
        None -- Denotes intact PSI.    
    """
    packaging_begin = [m.start() for m in re.finditer(r"[^-]",
                       str(alignment[0].seq))][psi_locus[0]]
    query_start = [m.start() for m in re.finditer(r"[^-]",
                   str(alignment[1].seq))][0]
    packaging_end = [m.start() for m in re.finditer(r"[^-]",
                     str(alignment[0].seq))][psi_locus[1]]
    # if query_start > packaging_begin:
    #     return IntactnessError(
    #             alignment[1].id, PSINOTFOUND_ERROR,
    #             "Query Start at reference position " + str(query_start)
    #             + ". Does not encompass PSI at positions "
    #             + str(packaging_begin) + " to " + str(packaging_end) + "."
    #             )
    # #/end if
    query_psi = str(alignment[1].seq[packaging_begin:packaging_end])
    query_psi_deletions = len(re.findall(r"-", query_psi))
    if query_psi_deletions > psi_tolerance:
        return IntactnessError(
                alignment[1].id, PSIDELETION_ERROR,
                "Query Sequence exceeds maximum deletion tolerance in PSI. " + 
                "Contains " + str(query_psi_deletions) + " deletions with max "
                + "tolerance of " + str(psi_tolerance) + " deletions."
                )
    #/end if
    return None
#/end def has_packaging_signal

def has_rev_response_element(alignment, rre_locus, rre_tolerance):
    """
    Determines presence and possible intactness of HIV 
    Packaging Signal Region.
    
    
    Keyword Args:
        
        alignment -- multiple sequence alignment object containing the 
                     reference and query sequence.
        rre_locus -- tuple containing start and end coordinates of RRE wrt
                     the reference being used.
        RRE_tolerance -- number of deletions in query tolerated to be intact
        
        
    Internal Args:
        
        rre_begin -- Aligned RRE start position.
        rre_end -- Aligned RRE end position.
        query_rre -- extracted RRE region from query
        query_rre_deletions -- number of deletions in query RRE region
        
        
    Return:
        
        RREDELETION_ERROR -- IntactnessError denoting query likely contains
                             defective RRE.
        None -- Denotes intact RRE.    
    """
    rre_begin = [m.start() for m in re.finditer(r"[^-]",
                       str(alignment[0].seq))][rre_locus[0]]
    rre_end = [m.start() for m in re.finditer(r"[^-]",
                     str(alignment[0].seq))][rre_locus[1]]
    query_rre = str(alignment[1].seq[rre_begin:rre_end])
    query_rre_deletions = len(re.findall(r"-", query_rre))
    if query_rre_deletions > rre_tolerance:
        return IntactnessError(
                alignment[1].id, RREDELETION_ERROR,
                "Query Sequence exceeds maximum deletion tolerance in RRE. " + 
                "Contains " + str(query_rre_deletions) + " deletions with max "
                + "tolerance of " + str(rre_tolerance) + " deletions."
                )
    #/end if
    return None
#/end def has_rev_response_element

def reading_frames_single_stranded(alignment, sequence, length):
    """
    Find all reading frames longer than length in the forward strand
    of the given sequence.

    Args:
        sequence: the sequence to check.
        length: the minimum nucleotide length of a reading frame

    Returns:
        A list of tuples of (frame_start, frame_end)
    """

    # figure out where the query starts w.r.t HXB2 in case full
    # genome consensus is not being used
    # offset = re.search(r'[^-]', str(alignment[1].seq)).start()

    # for each position in the query, figure out how many inserts and
    # deletes we've seen with respect to the reference
    delete_offset = []
    insert_offset = []
    delete_count = 0
    insert_count = 0

    for i in range(len(alignment[0])):
        if alignment[1][i] == "-":
            delete_count += 1
            continue
        if alignment[0][i] == "-":
            insert_count += 1
        delete_offset.append(delete_count)
        insert_offset.append(insert_count)

    long_frames = []

    for frame in range(0, 3):

        for_translation = sequence.seq[frame:]
        for _ in range(3 - len(sequence.seq[frame:]) % 3):
            for_translation += 'N'

        protein = Seq.translate(for_translation)

        current_start = 0
        current_len = 0
        

        for i, elem in enumerate(protein):
            if elem == "*":
                fs = current_start * 3 + 1 + frame
                frame_start = fs + delete_offset[fs] - insert_offset[fs]
                fe = i * 3 + 1 + frame
                frame_end = fe + delete_offset[fe] - insert_offset[fe]

                if current_len * 3 >= length:
                    long_frames.append(
                        (frame_start, frame_end + 1, delete_offset[fe] - delete_offset[fs], insert_offset[fe] - insert_offset[fs])
                    )
                current_len = 0
                continue
            elif current_len == 0:
                current_start = i

            current_len += 1

    long_frames.sort(key=lambda x: x[0])            

    return long_frames

def alignment_score(alignment):
    """
    Simple score for an alignment (just to find out if it's forward or
    reverse - absolutely not a true genetic distance)
    """

    return sum([a==b for a, b in zip(alignment[0].seq, alignment[1].seq)])

def small_frames(
    alignment, sequence, length, 
    expected, error_bar, reverse = False
):
    """
    Check for presence of small reading frames
    """
    frames = reading_frames_single_stranded(
                           alignment,
                           sequence, length)
    f_type = "forward"
    if reverse:
        tmp_reference = SeqRecord.SeqRecord(Seq.reverse_complement(alignment[0].seq),
                                        id = alignment[0].id,
                                        name = alignment[0].name
                                        )
        tmp_subtype = SeqRecord.SeqRecord(Seq.reverse_complement(alignment[1].seq),
                                        id = alignment[1].id,
                                        name = alignment[1].name
                                        )
        tmp_sequence = SeqRecord.SeqRecord(Seq.reverse_complement(sequence.seq),
                                        id = sequence.id,
                                        name = sequence.name
                                        )

        reverse_alignment = [tmp_reference, tmp_subtype]
        frames = reading_frames_single_stranded(
                            reverse_alignment,
                            tmp_sequence,
                            length)
        f_type = "reverse"

    if len(frames) == 0:
        return [IntactnessError(
            sequence.id, WRONGORFNUMBER_ERROR,
            "No ORFs >" + str(length) + " bases found.")]

    errors = []
    for e in expected:
        best_match = (0, 0)
        best_match_delta = 10000000
        found_inside_match = False
        for f in frames:
            inside = f[0] < e[1] and f[1] > e[2]
            delta = abs(e[1] - f[0]) + abs(e[2] - f[1])
            # only compare inside matches to the best inside match
            if inside == found_inside_match and delta < best_match_delta:
                best_match = f
                best_match_delta = delta
            # if we find an inside match, erase all non-inside matches
            if inside and not found_inside_match:
                found_inside_match = True
                best_match = f
                best_match_delta = delta

        # ORF lengths and locations are incorrect
        if (best_match[0] - e[1] > error_bar \
        or e[2] - best_match[1] > error_bar): 
            errors.append(IntactnessError(
                sequence.id, MISPLACEDORF_ERROR,
                "Expected a smaller ORF, " + str(e[0]) + ", at " + str(e[1]) 
                + "-" + str(e[2]) 
                + " in the " + f_type + " strand, got closest match " 
                + str(best_match[0]) 
                + "-" + str(best_match[1])
            )) 
        else:
            insertions = len(re.findall(r"-", str(alignment[0].seq[e[1]:e[2]])))
            deletions = len(re.findall(r"-", str(alignment[1].seq[e[1]:e[2]])))

            # Max deletion allowed in ORF exceeded
            if deletions > e[3]:

                errors.append(IntactnessError(
                    sequence.id, DELETIONINORF_ERROR,
                    "Smaller ORF " + str(e[0]) + " at " + str(e[1]) 
                    + "-" + str(e[2]) 
                    + " can have maximum deletions "
                    + str(e[3]) + ", got " 
                    + str(deletions)
                ))

            # Check for frameshift in ORF
            if (deletions - insertions) % 3 != 0:

                errors.append(IntactnessError(
                    sequence.id, FRAMESHIFTINORF_ERROR,
                    "Smaller ORF " + str(e[0]) + " at " + str(e[1]) 
                    + "-" + str(e[2]) 
                    + " contains an out of frame indel: insertions " + str(insertions)
                    + " deletions " + str(deletions) + "."
                ))

        
    return errors
       

def has_reading_frames(
    alignment, reference,
    sequence, length, 
    forward_expected, reverse_expected, error_bar):
    """
    Check that a sequences has the appropriate number of reading frames
    longer than a certain length in the appropriate strands and
    positions.

    Args:
        sequence: the sequence to check.
        length: the minimum nucleotide length of a reading frame.

    Returns:
        A list of tuples of (frame_start, frame_end)
    """
    
    forward_frames = reading_frames_single_stranded(
                           alignment,
                           sequence, length)

    tmp_reference = SeqRecord.SeqRecord(Seq.reverse_complement(alignment[0].seq),
                                       id = alignment[0].id,
                                       name = alignment[0].name
                                       )
    tmp_subtype = SeqRecord.SeqRecord(Seq.reverse_complement(alignment[1].seq),
                                       id = alignment[1].id,
                                       name = alignment[1].name
                                       )
    tmp_sequence = SeqRecord.SeqRecord(Seq.reverse_complement(sequence.seq),
                                       id = sequence.id,
                                       name = sequence.name
                                       )

    reverse_alignment = [tmp_reference, tmp_subtype]
    reverse_frames = reading_frames_single_stranded(
                           reverse_alignment,
                           tmp_sequence,
                           length)


    orfs = []
    for f_type, got, expected in [
                                ("forward", forward_frames, forward_expected),
                                ("reverse", reverse_frames, reverse_expected)
                                 ]:
        for got_elem in got:
            orfs.append(ORF(f_type, got_elem[0], got_elem[1]))


    for f_type, got, expected in [
                                ("forward", forward_frames, forward_expected),
                                ("reverse", reverse_frames, reverse_expected)
                                 ]:


        if len(got) != len(expected) and len(expected) > 0:
            return orfs, [IntactnessError(
                sequence.id, WRONGORFNUMBER_ERROR,
                "Expected " + str(len(expected)) 
                + " " + f_type 
                + " ORFs, got " + str(len(got))
            )]
    
    errors = []
    for f_type, got, expected in [
                                ("forward", forward_frames, forward_expected),
                                ("reverse", reverse_frames, reverse_expected)
                                 ]:
        for got_elem, expected_elem in zip(got, expected):

            # ORF lengths and locations are incorrect
            if got_elem[0] - expected_elem[1] > error_bar \
            or expected_elem[2] - got_elem[1] > error_bar: 

                errors.append(IntactnessError(
                    sequence.id, MISPLACEDORF_ERROR,
                    "Expected an ORF, " + str(expected_elem[0]) + ", at " + str(expected_elem[1]) 
                    + "-" + str(expected_elem[2]) 
                    + " in the " + f_type + " strand, got " 
                    + str(got_elem[0]) 
                    + "-" + str(got_elem[1])
                ))

            # Max deletion allowed in ORF exceeded
            if got_elem[2] > expected_elem[3]:

                errors.append(IntactnessError(
                    sequence.id, DELETIONINORF_ERROR,
                    "ORF " + str(expected_elem[0]) + " at " + str(got_elem[0]) 
                    + "-" + str(got_elem[1]) 
                    + " can have maximum deletions "
                    + str(expected_elem[3]) + ", got " 
                    + str(got_elem[2])
                ))

            # Check for frameshift deletion in ORF
            if (got_elem[2] - got_elem[3]) % 3 != 0:

                errors.append(IntactnessError(
                    sequence.id, FRAMESHIFTINORF_ERROR,
                    "ORF " + str(expected_elem[0]) + " at " + str(got_elem[0]) 
                    + "-" + str(got_elem[1]) 
                    + " contains an out of frame indel, deletions " 
                    + str(got_elem[2]) + " insertions " + str(got_elem[3]) + "."
                ))

            


    return orfs, errors

def intact( working_dir,
            input_file,
            subtype,
            include_packaging_signal,
            include_rre,
            check_major_splice_donor_site,
            include_small_orfs,
            hxb2_forward_orfs = const.DEFAULT_FORWARD_ORFs,
            hxb2_reverse_orfs = const.DEFAULT_REVERSE_ORFS,
            hxb2_small_orfs = const.DEFAULT_SMALL_FORWARD_ORFS,
            hxb2_psi_locus = const.DEFAULT_PSI_LOCUS,
            hxb2_rre_locus = const.DEFAULT_RRE_LOCUS,
            hxb2_msd_site_locus = const.DEFAULT_MSD_SITE_LOCUS,
            min_orf_length = const.DEFAULT_ORF_LENGTH,
            error_bar = const.DEFAULT_ERROR_BAR):
    """
    Check if a set of consensus sequences in a FASTA file is intact.

    Args:
        input_folder: folder of files from NGS machine.

    Returns:
        Name of a file containing all consensus sequences.
    """

    intact_sequences = []
    non_intact_sequences = []
    orfs = {}
    errors = {}

    # convert ORF positions to appropriate subtype
    forward_orfs, reverse_orfs, small_orfs = [
    [                       
        (
            n,
            st.convert_from_hxb2_to_subtype(working_dir, s, subtype), 
            st.convert_from_hxb2_to_subtype(working_dir, e, subtype), 
            delta
        ) \
        for (n, s, e, delta) in orfs
    ] \
    for orfs in [hxb2_forward_orfs, hxb2_reverse_orfs, hxb2_small_orfs]
    ]

    # convert PSI locus and RRE locus to appropriate subtype
    psi_locus = [st.convert_from_hxb2_to_subtype(working_dir, x, subtype) for x in hxb2_psi_locus]
    rre_locus = [st.convert_from_hxb2_to_subtype(working_dir, x, subtype) for x in hxb2_rre_locus]

    reference = st.subtype_sequence(subtype)

    with open(input_file, 'r') as in_handle:
        for sequence in SeqIO.parse(in_handle, "fasta"):
            
            sequence_errors = []

            reverse_sequence = SeqRecord.SeqRecord(Seq.reverse_complement(sequence.seq),
                                       id = sequence.id + " [REVERSED]",
                                       name = sequence.name
                                       )

            
            alignment = wrappers.mafft(working_dir, [reference, sequence])  
            reverse_alignment = wrappers.mafft(working_dir, [reference, reverse_sequence])

            forward_score = alignment_score(alignment)
            reverse_score = alignment_score(reverse_alignment)
            if alignment_score(reverse_alignment) > alignment_score(alignment):
                print("Reversing sequence " + sequence.id + "; forward score " 
                        + str(forward_score) + "; reverse score " + str(reverse_score))
                alignment = reverse_alignment
                sequence = reverse_sequence

            sequence_orfs, orf_errors = has_reading_frames(
                alignment,
                reference, sequence, min_orf_length,
                forward_orfs, reverse_orfs, error_bar)
            sequence_errors.extend(orf_errors)

            small_orf_errors = small_frames(
                alignment, sequence, 100, 
                small_orfs, error_bar, reverse = False)
            if include_small_orfs:
                sequence_errors.extend(small_orf_errors)

            hxb2_found_orfs = [ORF(
                                    o.orientation,
                                    st.convert_from_subtype_to_hxb2(working_dir, o.start, o.orientation, subtype),
                                    st.convert_from_subtype_to_hxb2(working_dir, o.end, o.orientation, subtype)
                              ) for o in sequence_orfs]
            
            if include_packaging_signal:
                missing_psi_locus = has_packaging_signal(alignment,
                                                     psi_locus,
                                                     const.PSI_ERROR_TOLERANCE)
                if missing_psi_locus is not None:
                    sequence_errors.append(missing_psi_locus)

            if include_rre:
                missing_rre_locus = has_rev_response_element(alignment,
                                                         rre_locus,
                                                         const.RRE_ERROR_TOLERANCE
                                                         )
                if missing_rre_locus is not None:
                    sequence_errors.append(missing_rre_locus)
            
            if check_major_splice_donor_site:
                mutated_splice_donor_site = has_mutated_major_splice_donor_site(alignment,
                                                st.convert_from_hxb2_to_subtype(working_dir, hxb2_msd_site_locus, subtype),
                                                st.convert_from_hxb2_to_subtype(working_dir, hxb2_msd_site_locus + 1, subtype),
                                                const.DEFAULT_MSD_SEQUENCE)
                if mutated_splice_donor_site is not None:
                    sequence_errors.append(mutated_splice_donor_site)

            orfs[sequence.id] = hxb2_found_orfs
            if len(sequence_errors) == 0:
                intact_sequences.append(sequence)
            else:
                non_intact_sequences.append(sequence)
                errors[sequence.id] = sequence_errors
            
            # add the small orf errors after the intactness check if not included
            if not include_small_orfs:
                if sequence.id in errors:
                    errors[sequence.id].extend(small_orf_errors)
                else:
                    errors[sequence.id] = small_orf_errors

    intact_file = os.path.join(os.getcwd(), "intact.fasta")
    with open(intact_file, 'w') as f:
       SeqIO.write(intact_sequences, f, "fasta")

    non_intact_file = os.path.join(os.getcwd(), "nonintact.fasta")
    with open(non_intact_file, 'w') as f:
        SeqIO.write(non_intact_sequences, f, "fasta")
    
    orf_file = os.path.join(os.getcwd(), "orfs.json")
    with open(orf_file, 'w') as f:
        f.write(json.dumps({seq: [x.__dict__ for x in sorfs] \
                            for seq, sorfs in orfs.items()},
                            indent=4))

    error_file = os.path.join(os.getcwd(), "errors.json")
    with open(error_file, 'w') as f:
        f.write(json.dumps({seq: [x.__dict__ for x in serrors] \
                            for seq,serrors in errors.items()}, 
                            indent=4))

    return intact_file, non_intact_file, orf_file, error_file
#/end def intact
#/end intact.py