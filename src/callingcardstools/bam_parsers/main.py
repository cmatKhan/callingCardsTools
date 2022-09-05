# standard library
import sys
import os
import tempfile
import logging
# outside dependencies
import pysam
# local dependencies
from .ReadTagger import ReadTagger
from .BarcodeParser import BarcodeParser
from .StatusFlags import StatusFlags

logging.basicConfig(
    level=logging.CRITICAL,
    format='[%(asctime)s] {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s',
    datefmt='%H:%M:%S',
    stream=sys.stdout
)

def main(genome_path,genome_index_path,bampath,
         barcode_length,insertion_length,barcode_details_json,mapq_threshold = 0, 
         out_suffix = "_tagged.bam"):
    """Iterate over a bam file, set tags and output updated bam with read groups 
    added to the header, tags added to the reads. Also output a summary of the 
    reads

    Args:
        genome_path (str): Path to a fasta file
        genome_index_path (str): Path to the index file for the fasta (fai)
        bampath (str): path to the alignment file (bam)
        barcode_length (int): Expected length of the barcode
        insertion_length (int): Expected length of the insertion sequence
        barcode_details_json (str): Path to the barcode details json file
        mapq_threshold (int, optional): mapq threshold below which to label a 
        read as failing. Defaults to 0.
        out_suffix (str, optional): suffix to append to the augmented bam file 
        output. Defaults to "_tagged.bam".

    Returns:
        int: 0 if successful
    """
    # temp_dir is automatically cleaned when context ends
    with tempfile.TemporaryDirectory() as temp_dir:
        # create temp file in temp dir -- note that temp dir is destroyed 
        # when this context ends
        bampath_tmp = os.path.join(temp_dir, "tmp_tagged.bam")
        # create the path to store the (permenant) output bam
        bampath_out = os.path.join(os.path.splitext(bampath)[0], out_suffix)

        # open files
        # including the index allows random access on disc
        genome = pysam.FastaFile(genome_path, genome_index_path)
        # open the input bam
        input_bamfile = pysam.AlignmentFile(bampath, "rb")

        tmp_tagged_bam = pysam.AlignmentFile(bampath_tmp, "wb")

        rt = ReadTagger(genome_path, genome_index_path, 
                        barcode_length, insertion_length)

        bp = BarcodeParser(barcode_details_json)

        read_group_set = set()
        read_summary = []
        for read in input_bamfile.fetch():
            tagged_read = rt.tag_read(read)
            barcode_details = bp.barcode_check(tagged_read.get_tag("RG"))
            tagged_read.set_tag("XF", barcode_details['tf'])
            read_group_set.add(tagged_read.get_tag("RG"))

            status_code = 0
            if not barcode_details['pass']:
                status_code += StatusFlags['BARCODE']
            if tagged_read.mapping_quality < mapq_threshold:
                status_code += StatusFlags['MAPQ']
            try:
                if tagged_read.get_tag("XZ") not in bp.get_insert_seqs():
                    status_code += StatusFlags['INSERT_SEQ']
            except AttributeError as exc:
                logging.debug(f"insert sequence not found in Barcode Parser. {exc}")

            read_summary += {"id": tagged_read.query_name(),
                            "bc": tagged_read.get_tag("RG"),
                            "status": status_code,
                            "mapq": tagged_read.mapping_quality,
                            "three_prime": tagged_read.get("XT"),
                            "flag": tagged_read.flag,
                            "insert_start": tagged_read.get_tag("XI"),
                            "insert_stop": tagged_read.get_tag("XE"),
                            "insert_seq": tagged_read.get_tag("XZ")}

            tmp_tagged_bam.write(tagged_read)

        # close the write handle so we can create a read handle
        tmp_tagged_bam.close()

        # copy alignments from the tmp file to the actual output so that we can 
        # include the RG headers. It is frustrating that this seems like the 
        # only way to do this in pysam.
        # TODO find a way to just add the header rather than having to iterate over 
        # the reads
        new_header = input_bamfile.header.copy().to_dict()
        # Create new read group header. Note: this is used below in the tagged_bam
        new_header['RG'] = [{'ID': rg} for rg in read_group_set]
        # open the tmp_tagged_bam for reading
        tmp_tagged_bam = pysam.AlignmentFile(bampath_tmp, "rb")
        # open the final bam output path and add the updated header
        tagged_bam_output = \
            pysam.AlignmentFile(bampath_out, 'wb', header=new_header)
        # iterate over the reads to re-write
        for read in tmp_tagged_bam.fetch():
            tagged_bam_output.write(read)

        # close the temp bampath. Note that the whole temp directory will be 
        # deleted when we leave the with TempDirectory as ... clause
        tmp_tagged_bam.close()

    # Close files
    genome.close()
    tagged_bam_output.close()
    input_bamfile.close()

    # Re-index the output
    # This is only here only to prevent the warning message:
    # bam timestamp and read timestamp are different
    # that you sometimes get when the bam is modified 
    # after the index is created
    pysam.index(bampath_out)

    return 0
