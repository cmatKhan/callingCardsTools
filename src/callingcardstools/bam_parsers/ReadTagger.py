"""_summary_

Raises:
    FileNotFoundError: _description_
    ValueError: _description_
    FileNotFoundError: _description_
    TypeError: _description_
    ValueError: _description_

Returns:
    _type_: _description_
"""
# standard library
import os
import sys
import logging
from xml.dom.minidom import Attr
# outside dependencies
import pysam

logging.basicConfig(
    level=logging.CRITICAL,
    format='[%(asctime)s] {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s',
    datefmt='%H:%M:%S',
    stream=sys.stdout
)

class ReadTagger:
    """Given an indexed fasta file (genome), id length and insertion length, 
    this object can returned a read tagged with the RG, XT and XZ tags"""
    def __init__(self, fasta_path, fasta_index_path, barcode_length, insertion_length):
        # check genome and index paths
        for input_path in [fasta_path,fasta_index_path]:
            if not os.path.exists(input_path):
                raise FileNotFoundError(f"Input file DNE: {input_path}")
        # check types
        for input_val in [barcode_length, insertion_length]:
            if not isinstance(input_val, int):
                raise ValueError(f"Input value {input_val} must be an integer" \
                    f"currently a {type(input_val)}")
        # set attributes
        self.fasta = fasta_path
        self.fasta_fai = fasta_index_path
        self.barcode_length = barcode_length
        self.insertion_length = insertion_length
        # open genome
        self.genome = pysam.FastaFile(self.fasta, self.fasta_fai)
    
    def __del__(self):
        """ensure that the genome file is closed when deleted"""
        self.close()

    def close(self):
        """close the genome file"""
        self.genome.close()
    
    def open(self):
        """open the genome file and set the self.genome attribute"""
        self.genome =  pysam.FastaFile(self.fasta, self.fasta_fa)
    
    def set_fasta(self, fasta_path, fasta_index_path):
        """set path to new fasta file. Also closes old genome obj and opens 
        a new one with the updated fasta

        Args:
            fasta_path (str): filepath to fasta (genome) file
            fasta_index_path (str): path to fai index 

        Raises:
            FileNotFoundError: Raised if path to either fasta or fai DNE
        """
        # check genome and index paths
        for input_path in [fasta_path,fasta_index_path]:
            if not os.path.exists(input_path):
                raise FileNotFoundError(f"File DNE: {input_path}")
        # update paths
        self.fasta = fasta_path
        self.fasta_fai = fasta_index_path
        # close old genome, open new one
        self.close()
        self.open()
        
    
    def is_open(self):
        """check if genome file is open"""
        return self.genome.is_open()
    
    def tag_read(self, read):
        """given a AlignedSegment object, add RG, XT and XZ tags

        Args:
            read (AlignedSegment): An aligned segment object -- eg returned 
              in a for loop by interating over bam.fetch() object from pysam

        Raises:
            TypeError: Raised with the cigarstring is not parse-able in a given read
            ValueError: Raised when the insertion sequence indicies are out of bounds

        Returns:
            AlignedSegment: The same read, but with RG, XT and XZ tags added
        """
        # Extract RG, XT and XZ tags -------------------------------------------
        tag_dict = dict()

        # set Read Group
        tag_dict['RG'] = read.query_name[-self.barcode_length:]

        # (using the bitwise operator) check if the read is unmapped,
        # if so, set the region_dict start and end to *, indicating that there is
        # no alignment, and so there is no start and end region for the alignment
        if read.flag & 0x4:
            tag_dict['XT'] = "*"
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
                sys.exit(f"Read {read.query_name}, "
                f"cigar string {read.cigartuples} is not parse-able")

            # The insertion point is at the end of the alignment
            # note that this is -1 because per the docs
            # reference_end points to one past the last aligned residue.
            read_3_prime = (read.reference_end-1)+soft_clip_length

            # this is the `insert_length` number bases which precede the
            # read (after adjusting for soft clipping)
            try:
                # if the soft-clip adjustment put the 3 prime end beyond the
                # end of the chrom, set XT to *
                if(read_3_prime >
                   self.genome.get_reference_length(read.reference_name)):
                    tag_dict['XT'] = "*"
                    tag_dict['XI'] = "*"
                    tag_dict['XE'] = "*"
                    tag_dict['XZ'] = "*"
                # if the endpoint of the insertion sequence is off the end of
                # the chrom, set XZ to *
                elif(read_3_prime+1+self.insertion_length >=
                self.genome.get_reference_length(read.reference_name)):
                    tag_dict['XT'] = read_3_prime
                    tag_dict['XI'] = "*"
                    tag_dict['XE'] = "*"
                    tag_dict['XZ'] = "*"
                else:
                    # This is the first base -- adjusted for soft clipping -- 
                    # in the read which cover the genome
                    tag_dict['XT'] = read_3_prime
                    tag_dict['XI'] = read_3_prime - self.insertion_length
                    tag_dict['XE'] = read_3_prime
                    tag_dict['XZ'] = self.genome.fetch(read.reference_name,
                                            read_3_prime+1,
                                            read_3_prime+1 +
                                            self.insertion_length).upper()
            except ValueError:
                sys.exit(f"Read {read.query_name}, "
                f"insert region {read.reference_name}:{read_3_prime+1}-"\
                    f"{read_3_prime+1+self.insertion_length} is out of bounds")

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
                raise TypeError(f"Read {read.query_name}, "\
                    f"cigar string {read.cigartuples} is not parse-able") \
                        from exc
            # extract insert position
            read_3_prime = read.reference_start - soft_clip_length

            # this is the `insert_length` number bases which precede the
            # read (after adjusting for soft clipping)
            try:
                # if the 3 prime end, after soft clipping, is less than 0, set
                # XT to *
                if(read_3_prime < 0):
                    tag_dict['XT'] = "*"
                    tag_dict['XI'] = "*"
                    tag_dict['XE'] = "*"
                    tag_dict['XZ'] = "*"
                # if the insertion sequence extends beyond the beginning of the
                # chrom, set to *
                elif(read_3_prime-self.insertion_length < 0):
                    tag_dict['XT'] = read_3_prime
                    tag_dict['XI'] = "*"
                    tag_dict['XE'] = "*"
                    tag_dict['XZ'] = "*"
                else:
                    # This is the first base -- adjusted for soft clipping -- 
                    # in the read which cover the genome
                    tag_dict['XT'] = read_3_prime
                    tag_dict['XI'] = read_3_prime + 1
                    tag_dict['XE'] = read_3_prime + 1 + self.insertion_length
                    tag_dict['XZ'] = self.genome.fetch(read.reference_name,
                                            read_3_prime-self.insertion_length,
                                            read_3_prime).upper()
            except ValueError as exc:
                raise ValueError(f"Read {read.query_name}, "
                f"insert region "\
                    f"{read.reference_name}:{read_3_prime-self.insertion_length}-"\
                    f"{read_3_prime} is out of bounds") from exc

        # Set tags -------------------------------------------------------------
        for tag, tag_str in tag_dict.items():
            read.set_tag(tag, tag_str)
        
        return read
    