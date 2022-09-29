# standard library
import sys
import logging
# outside dependencies
import pysam

logging.basicConfig(
    level=logging.CRITICAL,
    format='[%(asctime)s] {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s',
    datefmt='%H:%M:%S',
    stream=sys.stdout
)

def add_read_group_and_tags(bampath_in, bampath_out, genome_path,
                   genome_index_path, id_length, insertion_length):
    """
    :param bampath_in: path to the sorted, indexed bam file with UMI IDs
                       in the ID entry
    :param bampath_out: the same bam, but with the UMI ID added as a RG: tag
    :param genome_path:
    :param genome_index_path:
    :param id_length: length of the barcode which is appended, after an _,
                      to the read id

    Look at each read, extract the ID, add it as a tag (tag prefix RG) and
    write out to another bam. Output is the same bam, but with RG:<index>
    added to each line
    """

    # TODO: Right now, this loops over the file twice to create the header,
    #       and then to add the tags to each read. Need to figure out how to
    #       add the header without re-writing the whole bam so that only loops
    #       one time.

    # open files ---------------------------------------------------------------
    # including the index allows random access on disc
    genome = pysam.FastaFile(genome_path, genome_index_path)
    # open the input bam
    input_bamfile = pysam.AlignmentFile(bampath_in, "rb")

    # Get the set of unique barcodes in the input_bamfile ----------------------
    # extract current header
    new_header = input_bamfile.header.copy().to_dict()
    # instantiate a set object
    barcode_set = set()

    # loop through the reads in the input bam, extract unique barcodes, add to
    # barcode set
    for read in input_bamfile.fetch():
        # extract the ID from the ID line
        barcode = read.query_name[-id_length:]
        barcode_set.add(barcode)
    # Create new read group header. Note: this is used below in the tagged_bam
    new_header['RG'] = [{'ID': barcode} for barcode in barcode_set]

    # Add RG, XI and XZ tags to each line in the bam file ----------------------
    input_bamfile = pysam.AlignmentFile(bampath_in, "rb")
    # open a file to which to write. Note: the new header
    tagged_bam = pysam.AlignmentFile(bampath_out, "wb", header = new_header)
    for read in input_bamfile.fetch():
        # this is here to serve as a breakpoint in testing
        # if(read.query_name == "NB501801:571:HNM3KAFX3:4:21409:11423:19883_GAATCAATTCACTACGTCAACATTTTTCTATCGA"):
        #     print("here!")

        # Extract RG, XI and XZ tags -------------------------------------------
        tag_dict = dict()

        # set Read Group
        tag_dict['RG'] = read.query_name[-id_length:]

        # (using the bitwise operator) check if the read is unmapped,
        # if so, set the region_dict start and end to *, indicating that there is
        # no alignment, and so there is no start and end region for the alignment
        if read.flag & 0x4:
            tag_dict['XI'] = "*"
            tag_dict['XZ'] = "*"
        # if the bit flag 0x10 is set, the read reverse strand. Handle accordingly
        elif read.flag & 0x10:

            # A cigartuple looks like [(0,4), (2,2), (1,6),..,(4,68)] if read
            # is reverse complement. If it is forward, it would have the (4,68),
            # in this case, in the first position.
            # The first entry in the tuple is the cigar operation and the
            # second is the length. Note that pysam does order the tuples in the
            # reverse order from the sam cigar specs, so cigar 30M would be
            # (0,30). 4 is cigar S or BAM_CSOFT_CLIP. The list operation below
            # extracts the length of cigar operation 4 and returns a integer.
            # if 4 DNE, then soft_clip_length is 0.
            try:
                soft_clip_length = read.cigartuples[-1][1] \
                    if read.cigartuples[-1][0] == 4 \
                    else 0
            except TypeError:
                sys.exit(f"In bamfile {bampath_in}, for read {read.query_name}, "
                f"cigar string {read.cigartuples} is not parse-able")

            # The insertion point is at the end of the alignment
            # note that this is -1 because per the docs
            # reference_end points to one past the last aligned residue.
            read_3_prime = (read.reference_end-1)+soft_clip_length

            # this is the `insert_length` number bases which precede the
            # read (after adjusting for soft clipping)
            try:
                # if the soft-clip adjustment put the 3 prime end beyond the
                # end of the chrom, set XI to *
                if(read_3_prime > genome.get_reference_length(read.reference_name)):
                    tag_dict['XI'] = "*"
                    tag_dict['XZ'] = "*"
                # if the endpoint of the insertion sequence is off the end of
                # the chrom, set XZ to *
                elif(read_3_prime+1+insertion_length >=
                genome.get_reference_length(read.reference_name)):
                    tag_dict['XI'] = read_3_prime
                    tag_dict['XZ'] = "*"
                else:
                    # This is the first base -- adjusted for soft clipping -- in the
                    # read which cover the genome
                    tag_dict['XI'] = read_3_prime
                    tag_dict['XZ'] = genome.fetch(read.reference_name,
                                            read_3_prime+1,
                                            read_3_prime+1+insertion_length).upper()
            except ValueError:
                sys.exit(f"In bamfile {bampath_in}, for read {read.query_name}, "
                f"insert region {read.reference_name}:{read_3_prime+1}-"\
                    f"{read_3_prime+1+insertion_length} is out of bounds")

        # else, Read is in the forward orientation. Note that a single end
        # forward strand read with no other flags will have flag 0
        else:

            # see if clause for lengthy explanation. This examines the first
            # operation in the cigar string. If it is a soft clip (code 4),
            # the length of the soft clipping is stored. Else there is 0 soft
            # clipping
            try:
                soft_clip_length = read.cigartuples[0][1] \
                    if read.cigartuples[0][0] == 4 \
                    else 0
            except TypeError as exc:
                raise TypeError(f"In bamfile {bampath_in}, for read {read.query_name}, "\
                    f"cigar string {read.cigartuples} is not parse-able") from exc
            # extract insert position
            read_3_prime = read.reference_start - soft_clip_length

            # this is the `insert_length` number bases which precede the
            # read (after adjusting for soft clipping)
            try:
                # if the 3 prime end, after soft clipping, is less than 0, set
                # XI to *
                if(read_3_prime < 0):
                    tag_dict['XI'] = "*"
                    tag_dict['XZ'] = "*"
                # if the insertion sequence extends beyond the beginning of the
                # chrom, set to *
                elif(read_3_prime-insertion_length < 0):
                    tag_dict['XI'] = read_3_prime
                    tag_dict['XZ'] = "*"
                else:
                    # This is the first base -- adjusted for soft clipping -- in the
                    # read which cover the genome
                    tag_dict['XI'] = read_3_prime
                    tag_dict['XZ'] = genome.fetch(read.reference_name,
                                            read_3_prime-insertion_length,
                                            read_3_prime).upper()
            except ValueError as exc:
                raise ValueError(f"In bamfile {bampath_in}, for read {read.query_name}, "
                f"insert region {read.reference_name}:{read_3_prime-insertion_length}-"\
                    f"{read_3_prime} is out of bounds") from exc

        # Set tags -------------------------------------------------------------
        # the list comprehension outputs [None, None, None, ...]. Out catches
        # this so it isn't printed to std out
        # how else to do this?
        out = [read.set_tag(tag,tag_str) for tag, tag_str in tag_dict.items()]
        # Write to file --------------------------------------------------------
        tagged_bam.write(read)

    # Close files --------------------------------------------------------------
    genome.close()
    tagged_bam.close()
    input_bamfile.close()

    # Re-index the output ------------------------------------------------------
    # This is only here only to prevent the warning message:
    # bam timestamp and read timestamp are different
    # that you sometimes get when the bam is modified after the index is created
    pysam.index(bampath_out)