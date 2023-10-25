# pylint:disable=W0622,C0103
# standard library
import os
import logging
import pickle
from collections import defaultdict
from typing import Iterable, DefaultDict
from functools import partial
import csv
import re
# outside dependencies
import pandas as pd

from callingcardstools.QC.StatusFlags import StatusFlags

__all__ = ['Qbed']

logger = logging.getLogger(__name__)


class InnerDefaultDict(defaultdict):
    """A nested defaultdict class.

    :param defaultdict: a nested defaultdict class
    :type defaultdict: defaultdict
    """

    def __init__(self, data_type=int):
        super().__init__(data_type)


class MiddleDefaultDict1(defaultdict):
    """A nested defaultdict class.

    :param defaultdict: a nested defaultdict class
    :type defaultdict: defaultdict
    """

    def __init__(self, data_type=int):
        super().__init__(partial(InnerDefaultDict, data_type))


class MiddleDefaultDict2(defaultdict):
    """A nested defaultdict class.

    :param defaultdict: a nested defaultdict class
    :type defaultdict: defaultdict
    """

    def __init__(self, data_type=int):
        super().__init__(partial(MiddleDefaultDict1, data_type))


class MiddleDefaultDict3(defaultdict):
    """A nested defaultdict class.

    :param defaultdict: a nested defaultdict class
    :type defaultdict: defaultdict
    """

    def __init__(self, data_type=int):
        super().__init__(partial(MiddleDefaultDict2, data_type))


class OuterDefaultDict(defaultdict):
    """A nested defaultdict class.

    :param defaultdict: a nested defaultdict class
    :type defaultdict: defaultdict
    """

    def __init__(self, data_type=int):
        super().__init__(partial(MiddleDefaultDict3, data_type))


class Qbed():
    """
    An object to write records from a tagged_read_dict to qbed file and
        qc files. See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8150125/
        for more details.

    Attributes:
        qbed_fields (list): List of strings. Values in list are column names
            for qbed file.
        qbed (OuterDefaultDict): A nested defaultdict object. The keys are the
            qbed_fields. The values are the counts of each record.
        status_dict (DefaultDict): A defaultdict object. The keys are the
            status flags. The values are the counts of each status flag.    
    """
    _qbed_fields: list
    _qbed: DefaultDict
    _status_dict: DefaultDict

    def __init__(self, pickle_path: str = None) -> None:
        """
        Create a ReadRecords object. This object will write records to
            a qbed file and a qc file.

        Args:
            pickle_path: Path to a pickle file to load. If None, then
                initialize a new ReadRecords object.
        """
        if pickle_path:
            if not os.path.exists(pickle_path):
                msg = f"Path to pickle file {pickle_path} does not exist"
                raise FileNotFoundError(msg)
            self.load(pickle_path)
        else:
            # set qbed fields
            self.qbed_fields = ['chr', 'start', 'end', 'depth',
                                'strand', 'annotation']
            self.qbed = OuterDefaultDict(int)
            self.status_dict = DefaultDict(int)

    @property
    def qbed_fields(self):
        """Get the qbed fields"""
        return self._qbed_fields

    @qbed_fields.setter
    def qbed_fields(self, value: list):
        """Set the qbed fields"""
        if not len(value) == 6:
            raise ValueError('qbed_fields must have 6 values')
        self._qbed_fields = value

    @property
    def qbed(self):
        """Get the qbed"""
        return self._qbed

    @qbed.setter
    def qbed(self, value: DefaultDict):
        """Set the qbed"""
        self._qbed = value

    @property
    def status_dict(self):
        """Get the status_dict"""
        return self._status_dict

    @status_dict.setter
    def status_dict(self, value: DefaultDict):
        """Set the status_dict"""
        self._status_dict = value

    # private methods ---------------------------------------------------------
    def _combine(self, other):
        """Combines the qbed property of two Qbed objects.

        Args:
            other (Qbed): Another Qbed object.

        Raises:
            ValueError: If the other object is not a Qbed object.
        """
        if not isinstance(other, Qbed):
            raise ValueError('The other object must be a Qbed object.')

        # Combine qbed property
        for chr, value1 in other.qbed.items():
            for start, value2 in value1.items():
                for end, value3 in value2.items():
                    for strand, value4 in value3.items():
                        for srt_seq, count in value4.items():
                            (self.qbed[chr]
                             [start]
                             [end]
                             [strand]
                             [srt_seq]) += count

        # Combine status_dict property
        for status, count in other.status_dict.items():
            self.status_dict[status] += count

    def _srt_writer(self, output_path: str,
                    single_srt_count: int,
                    multi_srt_count: int) -> None:
        # Open a TSV file for writing.
        with open(output_path, "w", newline="", encoding='utf-8') as tsvfile:
            fieldnames = ['srt_type', 'count']

            # Create a DictWriter instance with a tab delimiter.
            writer = csv.DictWriter(tsvfile,
                                    fieldnames=fieldnames,
                                    delimiter='\t')

            # Write
            writer.writeheader()
            writer.writerow({'srt_type': 'single_srt',
                            'count': single_srt_count})
            writer.writerow({'srt_type': 'multi_srt',
                            'count': multi_srt_count})

    # public methods ----------------------------------------------------------

    def load(self, file_path: str) -> None:
        """Load a BarcodeQcCounter object from a file using Pickle.

        Args:
            file_path (str): The file path where the object is stored.
        """
        logger.info("loading Qbed object from %s", file_path)
        with open(file_path, "rb") as file:
            file_data = pickle.load(file)
            if not isinstance(file_data, Qbed):
                raise TypeError(
                    f"{file_path} is not a Qbed object")
            # copy the data from the loaded object to the current instance
            self.qbed_fields = file_data.qbed_fields
            self.qbed = file_data.qbed
            self.status_dict = file_data.status_dict

    @classmethod
    def combine(
            cls, counters: Iterable["Qbed"]) -> "Qbed":
        """Combine multiple Qbed objects into a single object.

        Args:
            counters (Iterable[Qbed]): An iterable of
                Qbed objects.

        Returns:
            Qbed: A new Qbed object formed from the list of input Qbeds.
        """
        result = Qbed()

        for counter in counters:
            result._combine(counter)

        return result

    def __add__(self, other: "Qbed") -> "Qbed":
        """Add two Qbed objects together with the + operator."""
        if not isinstance(other, Qbed):
            raise TypeError("Both objects must be of type 'Qbed'")

        result = Qbed()
        return result.combine([self, other])

    def update(self,
               tagged_read: dict,
               status: int,
               insert_offset: int = 1,
               srt_tag: str = 'ST') -> None:
        """write records to both the raw qbed tmpfile and raw qc tmpfile.
         Note that these tempfiles will be destroyed when the object is
         destroyed.

        Args:
            tagged_read (dict): A pysam.AlignedSegment object which has been
                tagged with the appropriate calling cards tags based on the
                BarcodeParser object used to create the object.
            status (int): A value which reflects how the read performs
                based on pre-defined quality metrics. A status of 0 is
                considered a pass. A status of greater than 0 is a read which
                fails at least 1 quality metric
            insert_offset (int): number to add to tag XI value to calculate
                the end coordinate. For instance, if the start coord is the
                first T in TTAA, then the offset would be 4.
            srt_tag (str): The tag which corresponds to the SRT sequence
                of a given read. This will be included in the annotation
                column of the mammals qbed file.
        """
        if len({'read', 'barcode_details'}-tagged_read.keys()) > 0:
            raise KeyError('tagged_read must have keys '
                           '{"reads","barcode_details"}')

        if status == 0:
            # for mammals, the SRT tag is expected. This will raise a KeyError
            # if the SRT tag is not present
            try:
                srt_with_edit_dist = tagged_read['read'].get_tag(srt_tag)
                srt = re.sub(r'\/\d+', '', srt_with_edit_dist)
            except KeyError as exc:
                raise f"tagged_read must have SRT key {srt_tag}" from exc
            chr = tagged_read['read'].reference_name
            start = tagged_read['read'].get_tag('XI')
            end = tagged_read['read'].get_tag('XI') + insert_offset
            strand = '+' if tagged_read['read'].is_forward else '-'

            self.qbed[chr][start][end][strand][srt] += 1

        self.status_dict[status] += 1

    def write(self,
              filename: str,
              suffix: str = "",
              raw: bool = False) -> pd.DataFrame:
        """Translate the qbed object and status_dict to DataFrames and
        write to either a pickle or tsv file.

        Args:
            filename (str): The name of the file to write to.
            raw (bool): If True, write to a raw file. Otherwise,
                write to a tsv file.

        Returns:
            [pd.DataFrame, pd.DataFrame]: The qbed and status DataFrames
        """
        # create qbed DataFrame
        qbed_df = pd.DataFrame(columns=self.qbed_fields)
        concat_list = []  # List to hold the Series before concatenating
        single_srt_counter = 0
        multi_srt_counter = 0

        for chr, value1 in self.qbed.items():
            for start, value2 in value1.items():
                for end, value3 in value2.items():
                    for strand, value4 in value3.items():
                        locus_srt_set = set()
                        for srt_seq, count in value4.items():
                            locus_srt_set.add(srt_seq)
                            # Prepare a new series to add to qbed DataFrame
                            new_row = pd.Series([chr, start, end, count,
                                                strand, srt_seq],
                                                index=self.qbed_fields)
                            concat_list.append(new_row)
                            # add a hop record to the qbed DataFrame

                        # count single/multi srt as appropriate
                        if len(locus_srt_set) > 1:
                            multi_srt_counter += 1
                        else:
                            single_srt_counter += 1

        qbed_df = pd.concat([qbed_df,
                             pd.DataFrame(concat_list,
                                          columns=self.qbed_fields)],
                            ignore_index=True)

        # create status DataFrame
        status_df = pd.DataFrame(columns=['status', 'count'])
        concat_list = []  # List to hold the Series before concatenating

        for status, count in self.status_dict.items():
            status_str = ",".join(StatusFlags.decompose(status))
            # Prepare a new series to add to status DataFrame
            new_row = pd.Series([status_str, count], index=['status', 'count'])
            concat_list.append(new_row)

        status_df = pd.concat([status_df,
                               pd.DataFrame(concat_list,
                                            columns=['status', 'count'])],
                              ignore_index=True)

        # write to file
        if raw:
            pickle_output_file = filename + '_' + suffix + '_qbed.pkl' \
                if suffix else filename + '_qbed.pkl'
            logger.info("writing Qbed object to %s", pickle_output_file)
            with open(pickle_output_file, 'wb') as file:
                pickle.dump(self, file)
        else:
            qbed_output_file = filename + '_' + suffix + '.qbed' \
                if suffix else filename + '.qbed'
            logger.info("writing qbed to %s", qbed_output_file)
            qbed_df.to_csv(qbed_output_file,
                           sep='\t',
                           index=False,
                           header=False)

            qc_output_file = filename + '_' + suffix + '_aln_summary.tsv' \
                if suffix else filename + '_aln_summary.tsv'
            logger.info("writing qc summary to %s", qc_output_file)
            status_df.to_csv(qc_output_file,
                             sep='\t',
                             index=False,
                             header=False)

            srt_output_file = filename + '_' + suffix + '_srt_count.tsv' \
                if suffix else filename + '_srt_count.tsv'
            logger.info("writing srt summary to %s", srt_output_file)
            self._srt_writer(srt_output_file,
                             single_srt_counter,
                             multi_srt_counter)
