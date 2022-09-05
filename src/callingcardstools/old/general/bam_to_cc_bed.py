# standard library
import sys
from copy import deepcopy
import logging
# outside dependencies
import pysam
import pandas as pd

logging.basicConfig(
    level=logging.CRITICAL,
    format='[%(asctime)s] {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s',
    datefmt='%H:%M:%S',
    stream=sys.stdout
)


def bam_to_cc_bed(bampath, require_exact_map, require_exact_length,
               mapq_filter, bed_6_3_colnames):
    """Iterate through the bam file and extract to bed custom format reads which
    meet quality/alignment thresholds

    :param bampath: path to the sorted, indexed bam file with RG, XZ and XI tags
    :param mapq_filter: exclude reads less than or equal to this value
    :param require_exact_map: True to accept and count only those reads which
      map to the reference exactly -- no clipping, no indels
    :param require_exact_length: require that the full length of the read map.
      Note that this automatically sets require_exact_map to True
    :param bed_6_3_colnames: names of the 6 + 3 bed columns (see official bed
    format)
    :param dry_run: default False. Set to True for testing. Purpose is to be
        able to set a breakpoint before the write statement. Setting to True
        prevents writing to file.

    :return: none. Write the custom bed file to file as tsv. output name is
             the input name with .bam replaced by .bed

    """

    # create dictionaries to hold the paritions of read1s
    passing_reads = {k:[] for k in bed_6_3_colnames}
    failing_reads = deepcopy(passing_reads)

    # open the bame file and begin looping over reads
    input_bamfile = pysam.AlignmentFile(bampath, "rb")
    for read in input_bamfile.fetch():
        # this is here to serve as a breakpoint in the tests
        # if(read.query_name == 'MN00200:583:000H3LYJT:1:22107:3155:20247_TAACGTCAATTTTACGCAGACTATCTTTCTAGGGTTAA'):
        #     print("here!")

        # only consider read 1
        if read.is_read1 or not read.is_paired:

            count_read = False

            # Below are further criteria to determine if read should be counted.
            # In either case -- paired or not -- throw out the read if the
            # mapping quality is too low

            if not read.is_paired:
                # do not count read if unmapped, secondary, failing qc
                # or is supplementary alignment
                count_read = True if \
                    sum([int(read.flag) & x for x in [0x4,0x100,0x200,0x800]]) == 0 \
                    and int(read.mapping_quality) >= mapq_filter \
                    else False
            else:
                # 83 and 99 are the flags for correctly paired R1 on reverse and
                # forward strand respectively
                count_read = True if int(read.flag) in [83,99] and \
                    int(read.mapping_quality) >= mapq_filter \
                    else False

            # the NM tag stores the edit distance. Edit distance 0 signifies
            # an exact match between the query and ref seqs
            # require that read. Note that this does not test for any
            # clipping, only that the aligned portion of the read aligns
            # exactly
            # note that this would eliminate reads which align fully, but with
            # SNPs
            if require_exact_map:
                count_read = False if read.get_tag('NM') != 0 else count_read

            # require the alignment length to be exactly the reference length.
            # note that this does allow indels. set require_exact_map to
            # require exact match along the entire query
            if require_exact_length:
                count_read = False if read.query_length != \
                    read.query_alignment_length else count_read

            # if the 3 prime end of the read extends beyond the end of the
            # chromosome, then tag XI is set to *. Exclude these
            if read.get_tag("XI") == "*":
                count_read = False

            if read.is_unmapped:
                strand = "*"
                start = "*"
                end = "*"
            else:
                strand = '-' if read.is_reverse else "+"
                try:
                    start = read.get_tag("XI")+1 if strand == "-" \
                        else read.get_tag("XI") - len(read.get_tag("XZ"))
                except ValueError:
                    start = "*"
                except TypeError:
                    start = "*"
                try:
                    end = start + len(read.get_tag("XZ")) if strand == '-' \
                        else read.get_tag('XI')
                except ValueError:
                    end = "*"
                except TypeError:
                    end = "*"

            chrom = "*" if read.is_unmapped \
                else input_bamfile.get_reference_name(read.reference_id)

            # add read to correct dict -- this partitions read1s
            if count_read:
                try:
                    passing_reads['chrom'].append(chrom)
                    passing_reads['chromStart'].append(start)
                    passing_reads['chromEnd'].append(end)
                    passing_reads['barcode'].append(read.get_tag("RG"))
                    passing_reads['aln_mapq'].append(read.mapping_quality)
                    passing_reads['strand'].append(strand)
                    passing_reads['reads'].append(1)
                    passing_reads['insert_seq'].append(read.get_tag("XZ"))
                    passing_reads['aln_flag'].append(read.flag)
                except KeyError as err:
                    msg = "A expected key in the passing_reads dict does not match " +\
                        "the bed fields. If you see this error, you should open an " +\
                            "issue report on github. Please post this error:\n%s " %err
                    logging.critical(msg, exc_info=(sys.exc_info()))
                    raise
            else:
                try:
                    failing_reads['chrom'].append(chrom)
                    failing_reads['chromStart'].append(start)
                    failing_reads['chromEnd'].append(end)
                    failing_reads['barcode'].append(read.get_tag("RG"))
                    failing_reads['aln_mapq'].append(read.mapping_quality)
                    failing_reads['strand'].append(strand)
                    failing_reads['reads'].append(1)
                    failing_reads['insert_seq'].append(read.get_tag("XZ"))
                    failing_reads['aln_flag'].append(read.flag)
                except KeyError as err:
                    msg = "A expected key in the failing_reads dict does not match " +\
                        "the bed fields. If you see this error, you should open an " +\
                            "issue report on github. Please post this error:\n%s " %err
                    logging.critical(msg, exc_info=(sys.exc_info()))
                    raise

    # close the bamfile
    input_bamfile.close()

    # convert read dicts to dataframes
    passing_reads_df = pd.DataFrame(data=passing_reads)
    failing_reads_df = pd.DataFrame(data=failing_reads)
   
    return (passing_reads_df,failing_reads_df)