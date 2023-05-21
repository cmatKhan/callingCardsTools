import os
import gzip
import logging

from Bio import SeqIO

from callingcardstools.BarcodeParser.BarcodeParser import BarcodeParser

logger = logging.getLogger(__name__)

__all__ = ['ReadParser']


class ReadParser(BarcodeParser):
    """Given either single or paired end reads, use the provided barcode
     details json to examine expected read components

    Depending on the entries in the barcode details json, this class will parse
     the read(s), return the assembled components (this could include both what
     could be construed as a barcode as well as any other components), and the
     reads trimmed for the components labelled for trimming.
    
    Attributes:
        r1_path (str): File path to Read 1.
        r2_path (str): File path to Read 2.
        r1_handle (_io.TextIOWrapper): File handle for Read 1.
        r2_handle (_io.TextIOWrapper): File handle for Read 2.
        cur_r1 (SeqRecord): Current read record for Read 1.
        cur_r2 (SeqRecord): Current read record for Read 2.

    Example:
    ```
    >>> rb = ReadParser('/path/to/barcode_details.json')
    >>> r1 = next(SeqIO.parse('/path/to/r1.fq', format="fastq"))
    >>> r2 = next(SeqIO.parse('/path/to/r2.fq', format="fastq"))
    >>> read_dict = rb.process_read(r1,r2)
    ```
    """
    _r1_path = ""
    _r2_path = ""
    _r1_handle = ""
    _r2_handle = ""
    _cur_r1 = ""
    _cur_r2 = ""

    def __init__(self, barcode_details_json: str = "", r1: str = "", r2: str = "") -> None:  # noqa
        """Initializes the ReadParser with barcode details and optional 
        read files.

        Args:
            barcode_details_json (str): The JSON file path containing 
                barcode details.
            r1 (str): The path to the Read 1 file.
            r2 (str): The path to the Read 2 file.
        """
        if barcode_details_json:
            super().__init__(barcode_details_json)
        if r1:
            self.r1_path = r1
        if r2:
            self.r2_path = r2

    def __del__(self):
        self.close()

    @property
    def r1_path(self):
        """filepath to read 1"""
        return self._r1_path

    @r1_path.setter
    def r1_path(self, r1_path):
        try:
            self.fastq_path_parser(r1_path)
        except FileNotFoundError:
            raise
        except IOError:  # pylint:disable=W0706
            raise
        self._r1_path = r1_path

    @property
    def r2_path(self):
        """filepath to read 2"""
        return self._r2_path

    @r2_path.setter
    def r2_path(self, r2_path):
        try:
            self.fastq_path_parser(r2_path)
        except FileNotFoundError:
            raise
        except IOError:  # pylint:disable=W0706
            raise
        self._r2_path = r2_path

    @property
    def r1_handle(self):
        """open SeqIO file handle to read 1"""
        return self._r1_handle

    @r1_handle.setter
    def r1_handle(self, r1_handle):
        self._r1_handle = r1_handle

    @property
    def r2_handle(self):
        """open SeqIO file handle to read 2"""
        return self._r2_handle

    @r2_handle.setter
    def r2_handle(self, r2_handle):
        self._r2_handle = r2_handle

    @property
    def cur_r1(self):
        """Current read SeqRecord for read 1"""
        return self._cur_r1

    @cur_r1.setter
    def cur_r1(self, r1_seqrecord):
        self._cur_r1 = r1_seqrecord

    @property
    def cur_r2(self):
        """Current read SeqRecord for read 2"""
        return self._cur_r2

    @cur_r2.setter
    def cur_r2(self, r2_seqrecord):
        self._cur_r2 = r2_seqrecord

    def fastq_path_parser(self, fq_path: str) -> bool:
        """Checks if the FastQ file path is valid.

        Args:
            fq_path (str): The path to the FastQ file.

        Raises:
            FileNotFoundError: If the FastQ file does not exist.
            IOError: If the FastQ file extension is not .fastq, .fq, or .gz.

        Returns:
            bool: True if the file path is valid, otherwise False.
        """
        error_msg = 'fastq extension must be either .fastq or .fq. ' +\
            'it may be gzipped, eg .fq.gz'
        fq_extensions = {'.fastq', '.fq'}

        if not os.path.exists(fq_path):
            raise FileNotFoundError(f'{fq_path} Does Not Exist!')
        elif os.path.splitext(fq_path)[1] not in fq_extensions:
            if os.path.splitext(fq_path)[1] == ".gz":
                if not os.path.splitext(os.path.splitext(fq_path)[0])[1] in fq_extensions:  # noqa
                    raise IOError(error_msg)
            else:
                raise IOError(error_msg)

    def open(self) -> None:
        """Opens the read file(s).

        If the file is gzipped, it is opened with gzip, otherwise it's opened normally.
        In case of paired-end reads, both files are opened.

        Raises:
            AttributeError: If the Read 1 file path is not set.
        """
        if not self.r1_path:
            raise AttributeError('R1 Not Set')

        elif os.path.splitext(self.r1_path)[1] == ".gz":
            self.r1_handle = SeqIO.parse(
                gzip.open(self.r1_path, "rt"),
                format='fastq')
        else:
            self.r1_handle = SeqIO.parse(
                self.r1_path,
                format='fastq')

        if not self.r2_path:
            print('Only R1 set -- single end mode')
        else:
            print('Opening R1 and R2 in paired end mode')
            if os.path.splitext(self.r2_path)[1] == ".gz":
                self.r2_handle = SeqIO.parse(
                    gzip.open(self.r2_path, "rt"),
                    format='fastq')
            else:
                self.r2_handle = SeqIO.parse(
                    self.r2_path,
                    format='fastq')

    def close(self) -> None:
        """close file objects, if they are set"""
        if self.r1_handle:
            self.r1_handle = ""
        if self.r2_handle:
            self.r2_handle = ""

    def next(self) -> None:
        """Iterate the reads and set cur_r1 and cur_r2 (if paired end) to the
         next read in hte file

        Raises:
            AttributeError: If R1 is not open
            StopIteration: If the end of file is reached in either
                R1 or R2 (if it is set)
            ValueError: If after advancing the read ids for
                R1 and r2 (if paired end) do not match
        """
        r1_stop = False
        r2_stop = False
        if not self.r1_handle:
            raise AttributeError('R1 is not open -- cannot advance')
        elif not self.r2_handle:
            self.cur_r1 = next(self.r1_handle)
        else:
            try:
                self.cur_r1 = next(self.r1_handle)
            except StopIteration:
                r1_stop = True
            try:
                self.cur_r2 = next(self.r2_handle)
            except StopIteration:
                r2_stop = True
            if r1_stop != r2_stop:
                error_msg = 'The length of R1 and R2 is not the same'
                logger.warning(error_msg)  # pylint:disable=E1102
                raise IOError(error_msg)
            elif r1_stop or r2_stop:
                raise StopIteration
            if self.cur_r1.id != self.cur_r2.id:
                error_msg = f"r1 ID: {self.cur_r1.id} does not match r2 "\
                    f"fID: {self.cur_r2.id}"
                logger.critical(error_msg)  # pylint:disable=E1102
                raise ValueError(error_msg)

    def parse(self, set_name_to_empty: bool = True) -> dict:  # noqa
        """Using the barcode details, process a given read

        Args:
            set_name_to_empty (bool, optional): Set the name attribute to 
                empty. SeqIO sets to ID if a name DNE. Defaults to True.

        Returns:
            dict: A dictionary with the r1 and r2 SeqRecord objects, the 
                barcode and the unadulterated read ID
        """
        # extract this in case it is augmented with the barcode later on
        read_id = self.cur_r1.id
        # create the beginnigns out of the output dictionary
        read_dict = {
            'r1': self.cur_r1,
            'r2': self.cur_r2
        }
        # instantiate a dict to hold the sequences corresponding to
        # barcode components
        components_dict = {}
        # a dictionary to assist with trimming the reads
        offset_dict = {
            'r1': 0,
            'r2': 0
        }
        for end, read in read_dict.items():
            if read:
                for k, v in self.barcode_dict[end].items():
                    components_dict["_".join([end, k])] = str(
                        read.seq[v['index'][0]:v['index'][1]])
                    # adjust offset for trimming
                    if v['trim']:
                        left_index = v['index'][0] - (offset_dict[end])
                        right_index = v['index'][1] - (offset_dict[end])
                        read_dict[end] = \
                            read_dict[end][:left_index] + \
                            read_dict[end][right_index:]
                        # adjust offset
                        offset_dict[end] += v['index'][1]
        # set name to empty string if flag set
        if set_name_to_empty:
            read_dict['r1'].name = ""
            try:
                read_dict['r2'].name = ""
            except NameError:
                pass
        # add barcode to read_dict
        read_dict['components'] = components_dict
        read_dict['status'] = self.component_check(components_dict)
        # add parse-able bam_tag line to the read id
        append_to_read_id = []
        for k, v in read_dict['status']['details'].items():
            if v.get('bam_tag', None):
                tag_value = v['name'] if v['dist'] == 0 \
                    else "/".join([v['name'], str(v['dist'])])
                append_to_read_id.append("-".join([v['bam_tag'], tag_value]))
        append_to_read_id = ";".join(append_to_read_id)
        if append_to_read_id:
            augmented_read_id = "_".join([read_id, append_to_read_id])
            read_dict['r1'].id = augmented_read_id
            read_dict['r2'].id = augmented_read_id
        read_dict['id'] = read_id

        return read_dict
