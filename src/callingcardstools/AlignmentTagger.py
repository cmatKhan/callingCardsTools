"""An object to facilitate adding tags to alignments in a bam file"""
# standard library
import os
import sys
import logging 
# local dependendecies
from .BarcodeParser import BarcodeParser
# outside dependencies
import pysam
from pysam import FastaFile

logging.getLogger(__name__).addHandler(logging.NullHandler())

class AlignmentTagger(BarcodeParser):
    """Given an indexed fasta file (genome), id length and insertion length, 
    this object can returned a read tagged with the RG, XS and XZ tags"""

    _fasta = ""
    _genome = ""
    _insertion_length = 1

    def __init__(self, fasta_path:str, barcode_details_json:str, insertion_length:int = 1) -> None:
        """_summary_

        Args:
            fasta_path (str): Path to a genome fasta file. Note that a fai index file created by samtools faidx must 
            exist at the same location.
            barcode_details_json (str): Path to the barcode_details json file
            insertion_length (int): The length of the sequence which is expected to fall directly before the 5' end of 
            the read. Defaults to 1.
        """
        super().__init__(barcode_details_json)
        self.fasta(fasta_path)
        self.insertion_length(insertion_length)
        self.fasta(fasta_path)
        # open the fasta file as a pysam.FastaFile obj
        self.open()
    
    def __del__(self):
        """ensure that the genome file is closed when deleted"""
        del self.genome
    
    @property
    def fasta(self) -> str:
        """path to the fasta file. The index .fai file MUST BE in the same directory"""
        return self._fasta
    @fasta.setter
    def fasta(self, new_fasta: str) -> None:
        if not os.path.exists(new_fasta+'.fai'):
            raise FileNotFoundError(f"Genome index not found for {new_fasta}. "\
                f"The index .fai file must exist in same path. "\
                    f"Use samtools faidx to create an index if one DNE")
        self._fasta = new_fasta
    
    @property
    def genome(self):
        """pysam FastaFile object"""
        return self._genome
    @genome.setter
    def genome(self, new_genome_obj:pysam.FastaFile):
        self._genome = new_genome_obj
    @genome.deleter
    def genome(self):
        try:
            self._genome.close()
        except AttributeError:
            pass

    @property
    def insertion_length(self):
        """length of the sequence expected to immediately preceed the 5' end of a read"""
        return self._insertion_length
    @insertion_length.setter
    def insertion_length(self,new_insertion_length):
        if not isinstance(new_insertion_length, int):
            raise ValueError(f"Input value {new_insertion_length} must be an integer" \
                f"currently a {type(new_insertion_length)}")
        self._insertion_length = new_insertion_length

    def open(self):
        """open the genome file and set the self.genome attribute"""
        self.genome(pysam.FastaFile(self.fasta, self.fasta+'.fai'))
    
    def is_open(self):
        """check if genome file is open"""
        return self.genome.is_open()
    
    def close(self):
        """close the genome file"""
        del self.genome
        
    def tag_read(self, read):
        """given a AlignedSegment object, add RG, XS and XZ tags

        Args:
            read (AlignedSegment): An aligned segment object -- eg returned 
              in a for loop by interating over bam.fetch() object from pysam

        Raises:
            TypeError: Raised with the cigarstring is not parse-able in a given read
            ValueError: Raised when the insertion sequence indicies are out of bounds

        Returns:
            AlignedSegment: The same read, but with RG, XS and XZ tags added
        """
        # Extract RG, XS and XZ tags -------------------------------------------
        tag_dict = dict()

        # set Read Group
        tag_dict['RG'] = read.query_name[-self.barcode_length:]

        # (using the bitwise operator) check if the read is unmapped,
        # if so, set the region_dict start and end to *, indicating that there is
        # no alignment, and so there is no start and end region for the alignment
        if read.flag & 0x4:
            tag_dict['XS'] = "*"
            tag_dict['XI'] = "*"
            tag_dict['XE'] = "*"
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
            read_5_prime = (read.reference_end-1)+soft_clip_length

            # this is the `insert_length` number bases which precede the
            # read (after adjusting for soft clipping)
            try:
                # if the soft-clip adjustment put the 3 prime end beyond the
                # end of the chrom, set XS to *
                if(read_5_prime >
                   self.genome.get_reference_length(read.reference_name)):
                    tag_dict['XS'] = "*"
                    tag_dict['XI'] = "*"
                    tag_dict['XE'] = "*"
                    tag_dict['XZ'] = "*"
                # if the endpoint of the insertion sequence is off the end of
                # the chrom, set XZ to *
                elif(read_5_prime+1+self.insertion_length >=
                self.genome.get_reference_length(read.reference_name)):
                    tag_dict['XS'] = read_5_prime
                    tag_dict['XI'] = "*"
                    tag_dict['XE'] = "*"
                    tag_dict['XZ'] = "*"
                else:
                    # This is the first base -- adjusted for soft clipping -- 
                    # in the read which cover the genome
                    tag_dict['XS'] = read_5_prime
                    tag_dict['XI'] = read_5_prime + 1
                    tag_dict['XE'] = read_5_prime + 1 + self.insertion_length
                    tag_dict['XZ'] = self.genome.fetch(read.reference_name,
                                            read_5_prime+1,
                                            read_5_prime+1 +
                                            self.insertion_length).upper()
            except ValueError:
                sys.exit(f"Read {read.query_name}, "
                f"insert region {read.reference_name}:{read_5_prime+1}-"\
                    f"{read_5_prime+1+self.insertion_length} is out of bounds")

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
            read_5_prime = read.reference_start - soft_clip_length

            # this is the `insert_length` number bases which precede the
            # read (after adjusting for soft clipping)
            try:
                # if the 5 prime end, after soft clipping, is less than 0, set
                # XS to *
                if(read_5_prime < 0):
                    tag_dict['XS'] = "*"
                    tag_dict['XI'] = "*"
                    tag_dict['XE'] = "*"
                    tag_dict['XZ'] = "*"
                # if the insertion sequence extends beyond the beginning of the
                # chrom, set to *
                elif(read_5_prime-self.insertion_length < 0):
                    tag_dict['XS'] = read_5_prime
                    tag_dict['XI'] = "*"
                    tag_dict['XE'] = "*"
                    tag_dict['XZ'] = "*"
                else:
                    # This is the first base -- adjusted for soft clipping -- 
                    # in the read which cover the genome
                    tag_dict['XS'] = read_5_prime
                    tag_dict['XI'] = read_5_prime - self.insertion_length
                    tag_dict['XE'] = read_5_prime
                    tag_dict['XZ'] = self.genome.fetch(read.reference_name,
                                            read_5_prime-self.insertion_length,
                                            read_5_prime).upper()
            except ValueError as exc:
                raise ValueError(f"Read {read.query_name}, "
                f"insert region "\
                    f"{read.reference_name}:{read_5_prime-self.insertion_length}-"\
                    f"{read_5_prime} is out of bounds") from exc

        # Set tags -------------------------------------------------------------
        for tag, tag_str in tag_dict.items():
            read.set_tag(tag, tag_str)
        
        return read
    