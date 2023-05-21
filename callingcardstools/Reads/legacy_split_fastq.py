"""
cc_filter_reads.py
written 9/1/16 by RDM

Modified 4/3/16 to remove /n line from split files.

Modified 12/17/18 to update to include restriction enzyme QC
Modified 12/21/18 to add quality control path information and raw folder

Read #1
The first 5 base pairs will be the primer barcode (index 0 to 4), followed by
17bp of transposon sequence (index 5 to 21):  AATTCACTACGTCAACA, after which
is genomic sequence

Read #2
The first 8bp will be the transposon barcode, followed by the restriction 
enzyme bank and then genomic DNA.  So it will look like
NNNNNNNN TCGA GCGC CCGG where the genomic sequence could start at any 
of the three restriction sites.

The barcode file should be of the following format: 
#expt name \t primer barcode \t transposon barcode 

usage 
python cc_filter_reads.py -r1 <read1 file> -r2 <read2 file> 
 -b<barcode file> -o <output path> -t <path to temp folder>
-rp <path to raw folder> 
--hammp <hamming distance for primer barcode>
--hammt <hamming distance for transposon barcode>


    Required
    -r1 {read 1 filename (full path)}
    -r2 {read 2 filename (full path)}

    Not Required
    -b {barcode file = ../scripts/barcodes.txt}
    -t {path to temp folder = ../temp/}
    -rp {path to raw folder = ../raw/} 
    -hp {hamming distance for primer bc =0}
    -tp {hamming distance for transposon bc = 0}
    -o {output path} -t {path to temp folder} -rp {path to raw folder}  
    --hammp {hamming distance for primer barcode}
    --hammt {hamming distance for transposon barcode}
    
    Description:

    1. Reads barcodes and corresponding experiments into a dictionary.
    2. Opens the read 1 and checks for transposon sequence.
    3. If the tranposon sequence is present, it checks to see if the primer 
    barcode matches the transposon barcode.
    4. If both filters are passed, it prints the reads to a file of the format: 
    exptname_primerbc_transposonbc_R1.fasta (or R2 or I2) in the raw directory.  
    Undetermined reads are outputted to the temp directory.  
    5. It also prints a master file of all of the reads file in the temp
    directory.
    6. The program then outputs a brief QC that lists the fraction of reads 
    that have a transposon match, the fraction of reads that have matching
    primer and transposon barcodes and the number of reads for each experiment 
    and the total number of reads analyzed.  It also computes, for each 
    experiment,the number of reads generated for each restriction enzyme.  
    This qc file is outputted to the temp directory, and is used by 
    perform_QC.py.
"""
import argparse
import logging
import csv
from Bio import SeqIO,Seq


__all__ = ['parse_args', 'legacy_split_fastq']

logger = logging.getLogger(__name__)


def read_barcode_file(barcode_filename):
    # This function reads in the experiment name, the primer barcode, and the transposon barcodes
    # from a file which is in the following format:
    # expt name \t primer barcode \t transposon barcode 1,transposon barcode 2, transposon barcode 3 etc.
    # It then returns a dictionary with key = a tuple (primer barcode, transposon barcode), value = expt_name
    # The last line of the file should not have a return character
    reader = csv.reader(open(barcode_filename, 'r'), delimiter='\t')
    d = {}
    for row in reader:
        exname, b1, b2 = row
        for transposon_barcode in b2.split(","):
            d[(b1, transposon_barcode)] = exname
    return d


def hamming_distance(s1, s2):
    # Return the Hamming distance between equal-length sequences
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def filter_reads(read1fn, read2fn, barcodefn, output, hammp, hammt):
    # This function does the following:
    # 1.  Reads barcodes and corresponding experiments into a dictionary
    # 2.  Opens the read 1 and checks for transposon sequence
    # 3.  If the tranposon sequence is present, it checks to see if the primer
    # barcode matches the transposon barcode
    # 4.  If both filters are passed, it prints the reads to a file of the format:
    # exptname_primerbc_transposonbc_R1.fasta (and/or R2)
    # 5.  It also prints a master file of all of the reads file
    # 6.  The program then outputs a brief QC that lists the fraction of reads that have a
    # transposon match, the fraction of reads that have matching primer and transposon barcodes
    # and the number of reads for each experiment and the total number of reads analyzed.

    # Define some infrequently changed variables

    FILTER_SEQUENCE = "AATTCACTACGTCAACA"
    PRIMER_BARCODE_START = 0
    PRIMER_BARCODE_END = 4
    TRANSPOSON_BARCODE_START = 0
    TRANSPOSON_BARCODE_END = 7
    RESTRICTION_BANK_SEQ = "TCGAGCGCCCGG"
    # MATCH BARCODES TO EXPERIMENT
    barcode_dict = read_barcode_file(barcodefn)
    print("I have read in the experiment barcodes.")

    # Put all filenames in this file for later
    filelist_filehandle = open(output+"cc_filelist.txt", 'w')

    # Make  a dictionary of filehandles for each barcode pair and undetermined barcodes
    # dictionary that contains the handles for each barcode pair
    r1_bcp_filehandle_dict = {}
    r2_bcp_filehandle_dict = {}  # same as above, but for read 2

    for key in barcode_dict.keys():
        r1_filename = output+barcode_dict[key] + \
            "_"+key[0]+"_"+key[1]+"_R1.fastq"
        # print the filename minus the _R1.fasta suffix
        # This makes the next step easier
        print(output+barcode_dict[key]+"_"+key[0] +
              "_"+key[1], file=filelist_filehandle)
        r1_bcp_filehandle_dict[key] = open(r1_filename, 'w')
        if read2fn:
            r2_filename = output + \
                barcode_dict[key]+"_"+key[0]+"_"+key[1]+"_R2.fastq"
            r2_bcp_filehandle_dict[key] = open(r2_filename, 'w')
    filelist_filehandle.close()

    # Make a filehandle to dump undetermined reads
    r1Undet_filehandle = open(output+"undetermined_R1.fastq", 'w')
    r2Undet_filehandle = open(output+"undetermined_R2.fastq", 'w')

    # Make a filehandle for the QC information
    qc_filehandle = open(output+"qc_filter_reads.txt", 'w')

    # get handles for read files
    r1Handle = open(read1fn, "rU")  # open read1 file
    if read2fn:
        r2Handle = open(read2fn, "rU")  # open read2 file

    # make iterators for index read and r2
    if read2fn:
        read2_record_iter = SeqIO.parse(r2Handle, "fastq")

    # initialize QC counters
    total_reads = 0
    reads_with_transposon_seq = 0
    matched_reads = 0
    expt_dict = {}
    rest_enzyme_dict = {}
    for expt in barcode_dict.values():
        expt_dict[expt] = 0
        rest_enzyme_dict[expt] = {"TaqAI": 0,
                                  "HinP1I": 0, "HpaII": 0, "Undet": 0}
    # intialize the experiment dictionary and restriction enzyme dictionaries

    # LOOP THROUGH READ 1 (and READ 2)
    for read1_record in SeqIO.parse(r1Handle, "fastq"):
        if read2fn:
            read2_record = next(read2_record_iter)
        total_reads += 1  # advance reads counter
        if FILTER_SEQUENCE in read1_record.seq[0:7+len(FILTER_SEQUENCE)]:
            reads_with_transposon_seq = reads_with_transposon_seq + 1
            primerbc = str(
                read1_record.seq[PRIMER_BARCODE_START:PRIMER_BARCODE_END+1])
            transbc = str(
                read2_record.seq[TRANSPOSON_BARCODE_START:TRANSPOSON_BARCODE_END+1])
            # try to correct barcodes if the hamming distance cutoff is g.t. zero
            if ((hammp > 0) or (hammt > 0)):
                for key in barcode_dict.keys():
                    if hamming_distance(primerbc, key[0]) <= hammp:
                        if hamming_distance(transbc, key[1]) <= hammt:
                            primerbc = key[0]
                            transbc = key[1]
            # print primerbc,transbc
            # is primer_bc transposon barcode pair in dictionary?
            if (primerbc, transbc) in barcode_dict:
                # if so, increment matched reads
                matched_reads += 1
                # update reads for the experiment
                expt_name = barcode_dict[(primerbc, transbc)]
                expt_dict[expt_name] += 1
                # check to see which restriction enzyme was cut in read 2 of this
                # read

                if RESTRICTION_BANK_SEQ in read2_record.seq[8:20]:
                    rest_enzyme_dict[expt_name]["HpaII"] += 1
                elif RESTRICTION_BANK_SEQ[0:8] in read2_record.seq[8:20]:
                    rest_enzyme_dict[expt_name]["HinP1I"] += 1
                elif RESTRICTION_BANK_SEQ[0:4] in read2_record.seq[8:20]:
                    rest_enzyme_dict[expt_name]["TaqAI"] += 1
                else:
                    rest_enzyme_dict[expt_name]["Undet"] += 1

                # output reads to the correct place
                print(read1_record.format("fastq")[0:len(read1_record.format(
                    "fastq"))-1], file=r1_bcp_filehandle_dict[(primerbc, transbc)])
                if read2fn:
                    print(read2_record.format("fastq")[0:len(read2_record.format(
                        "fastq"))-1], file=r2_bcp_filehandle_dict[(primerbc, transbc)])
            else:  # if there is no match, print reads to undetermined file
                print(read1_record.format("fastq")[
                      0:len(read1_record.format("fastq"))-1], file=r1Undet_filehandle)
                print(read2_record.format("fastq")[
                      0:len(read2_record.format("fastq"))-1], file=r2Undet_filehandle)
        else:  # if there is no match, print reads to undetermined file
            print(read1_record.format("fastq")[
                  0:len(read1_record.format("fastq"))-1], file=r1Undet_filehandle)
            print(read2_record.format("fastq")[
                  0:len(read2_record.format("fastq"))-1], file=r2Undet_filehandle)
    # print QC values to file
    print("There were "+str(total_reads)+" total reads", file=qc_filehandle)
    print(str(reads_with_transposon_seq/float(total_reads)) +
          " of the reads had a transposon sequence.", file=qc_filehandle)
    print(str(matched_reads/float(total_reads)) +
          " of the total reads also had matched barcodes.", file=qc_filehandle)
    print("\n", file=qc_filehandle)
    print("Experiment\tMatched Reads\tHpaII\tHinP1I\tTaqAI\tUndet", file=qc_filehandle)
    for key in expt_dict.keys():
        print(str(key)+"\t"+str(expt_dict[key]) +
              "\t"+str(rest_enzyme_dict[key]["HpaII"])+"\t" +
              str(rest_enzyme_dict[key]["HinP1I"])+"\t" +
              str(rest_enzyme_dict[key]["TaqAI"])+"\t" +
              str(rest_enzyme_dict[key]["Undet"]), file=qc_filehandle)
    qc_filehandle.close()

    # Close all filehandles
    r1Handle.close()
    r1Undet_filehandle.close()
    for key in r1_bcp_filehandle_dict.keys():
        r1_bcp_filehandle_dict[key].close()
    if read2fn:
        r2Handle.close()
        r2Undet_filehandle.close()
        for key in r2_bcp_filehandle_dict.keys():
            r2_bcp_filehandle_dict[key].close()


def parse_args(
        subparser: argparse.ArgumentParser,
        script_desc: str,
        common_args: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """This is intended to be used as a subparser for a parent parser passed 
    from __main__.py. This is the 'legacy' split fastq method from the 
    original Mitra lab pipeline. This method is not recommended for use.

    Args:
        subparser (argparse.ArgumentParser): See __main__.py -- this is the 
        subparser for the parent parser in __main__.py
        script_desc (str): Description of this script, which is set in 
        __main__.py. The description is set in __main__.py so that all of 
        the script descriptions are together in one spot and it is easier to 
        write a unified cmd line interface
        common_args (argparse.ArgumentParser): These are the common arguments 
        for all scripts in callingCardsTools, for instance logging level

    Returns:
        argparse.ArgumentParser: The subparser with the this additional 
        cmd line tool added to it -- intended to be gathered in __main__.py 
        to create a unified cmd line interface for the package
    """

    parser = subparser.add_parser(
        'legacy_split_fastq',
        help=script_desc,
        prog='legacy_split_fastq',
        parents=[common_args]
    )

    parser.set_defaults(func=legacy_split_fastq)

    parser.add_argument(
        '-r1',
        '--read1',
        help='Read 1 filename (full path)',
        required=True)

    parser.add_argument(
        '-r2',
        '--read2',
        help='Read2 filename (full path)',
        required=True)

    parser.add_argument(
        '-b',
        '--barcodefile',
        help='barcode filename (full path)',
        required=False,
        default='../scripts/barcodes.txt')

    parser.add_argument(
        '-o',
        '--output',
        help='path to output directory',
        required=False,
        default='../raw/')

    parser.add_argument(
        '--hammp',
        help='Primer barcode hamming distance',
        required=False,
        default=0)

    parser.add_argument(
        '--hammt',
        help='Transposon barcode hamming distance',
        required=False,
        default=0)

    return subparser


def legacy_split_fastq(args: argparse.Namespace) -> None:
    if not args.output[-1] == "/":
        args.output = args.output+"/"
    if not args.output[-1] == "/":
        args.output = args.output+"/"
    filter_reads(args.read1, args.read2, args.barcodefile,
                 args.output, int(args.hammp), int(args.hammt))
