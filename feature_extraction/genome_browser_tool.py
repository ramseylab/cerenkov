# This file is copyright 2002 Jim Kent, but license is hereby
# granted for all use - public, private or commercial.

# Bin indexing system used in UCSC Genome Browser
#   See http://genomewiki.ucsc.edu/index.php/Bin_indexing_system

# Note that `bin` is NOT a index column. Its ability to accelerate queries is limited.

binOffsets = [512+64+8+1, 64+8+1, 8+1, 1, 0]
binOffsetsExtended = [4096+512+64+8+1, 512+64+8+1, 64+8+1, 8+1, 1, 0]

_binFirstShift = 17  # How much to shift to get to finest bin.
_binNextShift = 3  # How much to shift to get to next larger bin.
_binOffsetOldToExtended = 4681  # From binRange.h


def __bin_from_range_standard(start, end):
    """
    Given start,end in chromosome coordinates, assign it a bin.
    There's a bin for each 128k segment, for each 1M segment, for each 8M segment, for each 64M segment,
    and for each chromosome (which is assumed to be less than 512M.)
    A range goes into the smallest bin it will fit in./
    """

    start_bin = start
    end_bin = end-1
    start_bin >>= _binFirstShift
    end_bin >>= _binFirstShift

    for i in range(0, len(binOffsets)):
        if start_bin == end_bin:
            return binOffsets[i] + start_bin
        start_bin >>= _binNextShift
        end_bin >>= _binNextShift

    raise ValueError("start {}, end {} out of range in findBin (max is 512M)".format(start, end))


# Add one new level to get coverage past chrom sizes of 512 Mb.
# Effective limit is now the size of an integer since chrom start and
#     end coordinates are always being used in int's == 2Gb-1
def __bin_from_range_extended(start, end):
    """
    Given start,end in chromosome coordinates, assign it a bin.
    There's a bin for each 128k segment, for each 1M segment, for each 8M segment, for each 64M segment,
    for each 512M segment, and one top level bin for 4Gb.

    Note, since start and end are int's, the practical limit is up to 2Gb-1, and thus,
    only four result bins on the second level.

    A range goes into the smallest bin it will fit in.
    """

    start_bin = start
    end_bin = end-1
    start_bin >>= _binFirstShift
    end_bin >>= _binFirstShift
    for i in range(0, len(binOffsetsExtended)):
        if start_bin == end_bin:
            return _binOffsetOldToExtended + binOffsetsExtended[i] + start_bin
        start_bin >>= _binNextShift
        end_bin >>= _binNextShift

    raise ValueError("start {}, end {} out of range in findBin (max is 2Gb)".format(start, end))


def bin_from_range(start, end):
    # Initial implementation is used when `chromEnd` is less than or equal to 536,870,912 = 2^29
    # Extended implementation is used when `chromEnd` is greater than 536,870,912 = 2^29 and
    #   less than 2,147,483,647 = 2^31 - 1
    if end <= 2**29:
        return __bin_from_range_standard(start, end)
    else:
        return __bin_from_range_extended(start, end)
