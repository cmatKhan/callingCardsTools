# standard library
#import cProfile
import sys
import os
import argparse
import tempfile
import logging
# outside dependencies
import pysam
import pandas as pd
# from memory_profiler import profile
# local dependencies
from callingcardstools.AlignmentTagger import AlignmentTagger
from callingcardstools.StatusFlags import StatusFlags

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s',
    datefmt='%H:%M:%S',
    stream=sys.stdout
)

def run(bampath, fasta_path, barcode_details_json,
        mapq_threshold = 0, out_suffix = "_tagged.bam", nthreads=5):
    """Iterate over a bam file, set tags and output updated bam with read groups 
    added to the header, tags added to the reads. Also output a summary of the 
    reads

    Args:
        fasta_path (str): Path to a fasta file
        bampath (str): path to the alignment file (bam)
        barcode_length (int): Expected length of the barcode
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
        # create the path to store the (permenant) output bam
        bampath_out = os.path.splitext(os.path.basename(bampath))[0] + out_suffix
        summary_out = os.path.splitext(bampath_out)[0] + "_summary.csv"

        # open files
        # open the input bam
        input_bamfile = pysam.AlignmentFile(
            bampath, "rb", 
            require_index=True,
            threads = nthreads)

        tmp_tagged_bam = pysam.AlignmentFile(
            bampath_tmp, 
            "wb", 
            header = input_bamfile.header)

        at = AlignmentTagger(barcode_details_json,fasta_path)
        
        read_group_set = set()
        read_summary = []
        read_obj_list = []
        # until_eof will include unmapped reads, also
        for read in input_bamfile.fetch(until_eof=True):
            tagged_read = at.tag_read(read)

            status_code = 0
            # if the read is unmapped, add the flag, but don't check 
            # other alignment metrics
            if tagged_read.is_unmapped:
                status_code += StatusFlags.UNMAPPED.flag()
            else:
                if tagged_read.is_qcfail:
                    status_code += StatusFlags.ALIGNER_QC_FAIL.flag()
                if tagged_read.is_secondary or tagged_read.is_supplementary:
                    status_code += StatusFlags.NOT_PRIMARY.flag()
                if tagged_read.mapping_quality < mapq_threshold:
                    status_code +=  StatusFlags.MAPQ.flag()
                # if the read is clipped on the 5' end, flag
                if (read.is_forward and tagged_read.query_alignment_start != 0) \
                    or (read.is_reverse and tagged_read.query_alignment_end != read.infer_query_length()):
                    status_code += StatusFlags.FIVE_PRIME_CLIP.flag()
            # check the insert sequence
            try:
                if not at.insert_seqs == ["*"]:
                    if tagged_read.get_tag("XZ") not in at.insert_seqs:
                            status_code += StatusFlags.INSERT_SEQ.flag()
            except AttributeError as exc:
                logging.debug(f"insert sequence not found in Barcode Parser. {exc}")
            
            summary_record = {"id": tagged_read.query_name,
                            "status": status_code,
                            "mapq": tagged_read.mapping_quality,
                            "flag": tagged_read.flag,
                            "chr": tagged_read.reference_name,
                            "strand": "*" if tagged_read.is_unmapped \
                                else "-" if tagged_read.is_reverse else '+',
                            "five_prime": tagged_read.get_tag("XS"),
                            "insert_start": tagged_read.get_tag("XI"),
                            "insert_stop": tagged_read.get_tag("XE"),
                            "insert_seq": tagged_read.get_tag("XZ")}
            
            # add the additional tagged elements, defined in the barcode_details json
            for k,v in at.tagged_components.items():
                summary_record[k] = tagged_read.get_tag(v)

            read_summary.append(summary_record)

            tmp_tagged_bam.write(tagged_read)
            #read_obj_list.append(tagged_read)

        # close the write handle so we can create a read handle
        tmp_tagged_bam.close()
        pysam.index(bampath_tmp)
        # copy alignments from the tmp file to the actual output so that we can 
        # include the RG headers. It is frustrating that this seems like the 
        # only way to do this in pysam.
        # TODO find a way to just add the header rather than having to iterate over 
        # the reads
        new_header = input_bamfile.header.to_dict()
        # Create new read group header. Note: this is used below in the tagged_bam
        new_header['RG'] = [{'ID': rg} for rg in read_group_set]
        # open the tmp_tagged_bam for reading
        tmp_tagged_bam = pysam.AlignmentFile(bampath_tmp, "rb")
        # open the final bam output path and add the updated header
        tagged_bam_output = \
            pysam.AlignmentFile(bampath_out, 'wb', header=new_header)
        # iterate over the reads to re-write
        print("re-writing bam with updated header...")
        count = 0
        for read in tmp_tagged_bam.fetch():
        #for read in read_obj_list:
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
    pysam.index(bampath_out)

    # write out summary
    print("writing summary...")
    pd.DataFrame(read_summary).to_csv(summary_out, index = False)
    
    print("done!")

def parse_args(args=None):
    """parse command line arguments for tagBam

    Args:
        args (_type_, optional): _description_. Defaults to None.

    Returns:
        _type_: _description_
    """
    Description = ""
    Epilog = ""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("--bam",
                         help="path to the input bam file")
    parser.add_argument("--fasta",
                        help="Note that an index .fai file must exist in the same path",
                        default='False')
    parser.add_argument("--barcode_details",
                         help = "")
    parser.add_argument("--mapq_threshold",
                         help = "", required=False)
    # parser.add_argument("--out_suffix",
    #                      help = "")

    return parser.parse_args(args)

def main(args=None):
    """_summary_

    Args:
        args (_type_, optional): _description_. Defaults to None.

    Raises:
        FileNotFoundError: _description_
    """
   args = parse_args(args)

    # Check inputs
    input_path_list = [args.bam,
                       args.fasta,
                       args.fasta + '.fai',
                       args.barcode_details]
    for input_path in input_path_list:
        if not os.path.exists(input_path):
            raise FileNotFoundError("Input file DNE: %s" %input_path)

    # loop over the reads in the bam file and add the read group (header and tag)
    # and the XI and XZ tags
    # with cProfile.Profile() as pr:
    run(args.bam,
        args.fasta,
        args.barcode_details,
        int(args.mapq_threshold))
    # pr.print_stats()

if __name__ == "__main__":
    sys.exit(main())
