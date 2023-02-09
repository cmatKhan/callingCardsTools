# pylint:disable=W1203
# standard library
import os
import tempfile
import logging
# outside dependencies
import pysam
import pandas as pd
# from memory_profiler import profile
# local dependencies
from callingcardstools.Alignment.AlignmentTagger import AlignmentTagger
from callingcardstools.QcStatusCoding.create_status_coder import create_status_coder # noqa

__all__ = ['tag_bam']

logging.getLogger(__name__).addHandler(logging.NullHandler())


def tag_bam(bampath, fasta_path, barcode_details_json,
            mapq_threshold=0, out_suffix="_tagged.bam", nthreads=5):
    """Iterate over a bam file, set tags and output updated bam with
     read groups added to the header, tags added to the reads.
     Also output a summary of the reads

    Args:
        fasta_path (str): Path to a fasta file
        bampath (str): path to the alignment file (bam)
        insertion_length (int): Expected length of the insertion sequence
        barcode_details_json (str): Path to the barcode details json file
        mapq_threshold (int, optional): mapq threshold below which to label a
         read as failing. Defaults to None.
        out_suffix (str, optional): suffix to append to the augmented bam file
         output. Defaults to "_tagged.bam".
        nthreads (int): Number of threads which pysam.AlignmentFile may use to
         decompress lines

    Returns:
        int: 0 if successful
    """
    print("tagging reads...")
    # temp_dir is automatically cleaned when context ends
    with tempfile.TemporaryDirectory() as temp_dir:
        # create temp file in temp dir -- note that temp dir is destroyed
        # when this context ends
        bampath_tmp = os.path.join(temp_dir, "tmp_tagged.bam")
        # create the path to store the (permanent) output bam
        bampath_out = os.path.splitext(os.path.basename(bampath))[0] + out_suffix  # noqa

        # open files
        # open the input bam
        input_bamfile = pysam.AlignmentFile(  # pylint:disable=E1101
            bampath, "rb",
            require_index=True,
            threads=nthreads)

        tmp_tagged_bam = pysam.AlignmentFile(  # pylint:disable=E1101
            bampath_tmp,
            "wb",
            header=input_bamfile.header)

        at = AlignmentTagger(barcode_details_json, fasta_path)

        status_coder = create_status_coder(
            mapq_threshold=mapq_threshold, 
            check_5_prime_clip=True)

        read_group_set = set()
        read_summary = []
        # until_eof will include unmapped reads, also
        for read in input_bamfile.fetch(until_eof=True):
            tagged_read = at.tag_read(read)

            status_code = status_coder(tagged_read)

            summary_record = {"id": tagged_read.get('read').query_name,
                              "status": status_code,
                              "mapq": tagged_read.get('read').mapping_quality,
                              "flag": tagged_read.get('read').flag,
                              "chr": tagged_read.get('read').reference_name,
                              "strand": "*" if tagged_read.get('read').is_unmapped  # noqa
                              else "-" if tagged_read.get('read').is_reverse else '+',  # noqa
                              "five_prime": tagged_read.get('read').get_tag("XS"),  # noqa
                              "insert_start": tagged_read.get('read').get_tag("XI"),  # noqa
                              "insert_stop": tagged_read.get('read').get_tag("XE"),  # noqa
                              "insert_seq": tagged_read.get('read').get_tag("XZ")}  # noqa

            # add the additional tagged elements, defined in
            # the barcode_details json
            for k, v in at.tagged_components.items():
                summary_record[k] = tagged_read.get('read').get_tag(v)

            read_summary.append(summary_record)

            tmp_tagged_bam.write(tagged_read.get('read'))
            # read_obj_list.append(tagged_read)

        # close the write handle so we can create a read handle
        tmp_tagged_bam.close()
        pysam.index(bampath_tmp)  # pylint:disable=E1101
        # copy alignments from the tmp file to the actual output so that
        # we can include the RG headers. It is frustrating that this
        # seems like the only way to do this in pysam.
        # TODO find a way to just add the header rather than having to
        # iterate over
        # the reads
        new_header = input_bamfile.header.to_dict()
        # Create new read group header. Note: this is used below in
        # the tagged_bam
        new_header['RG'] = [{'ID': rg} for rg in read_group_set]
        # open the tmp_tagged_bam for reading
        tmp_tagged_bam = pysam.AlignmentFile(bampath_tmp, "rb")   # pylint:disable=E1101 # noqa
        # open the final bam output path and add the updated header
        tagged_bam_output = pysam.AlignmentFile(  # pylint:disable=E1101
            bampath_out, 'wb', header=new_header)
        # iterate over the reads to re-write
        print("re-writing bam with updated header...")
        count = 0
        for read in tmp_tagged_bam.fetch():
            # for read in read_obj_list:
            tagged_bam_output.write(read)
            count += 1
        print(f"finished writing {count} lines to bam...")
        # close the temp bampath. Note that the whole temp directory will be
        # deleted when we leave the with TempDirectory as ... clause
        tmp_tagged_bam.close()

    # Close files
    tagged_bam_output.close()
    input_bamfile.close()

    # Re-index the output
    # This is only here only to prevent the warning message:
    # bam timestamp and read timestamp are different
    # that you sometimes get when the bam is modified
    # after the index is created
    print("indexing updated bam...")
    pysam.index(bampath_out)  # pylint:disable=E1101
    logging.info(f'{bampath_out} complete')
    return pd.DataFrame(read_summary)
