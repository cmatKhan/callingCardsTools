"""
This module contains the BarcodeQcCounter class, which is used for analyzing
barcode quality control data. The class can process and update barcode metrics,
combine multiple objects, and write the results to output files.

Author: Chase Mateusiak
Date: 2023-05-16
"""
from collections import defaultdict
from typing import DefaultDict, Iterable
import logging
import pickle
import os
import csv
from functools import partial

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


class OuterDefaultDict(defaultdict):
    """A nested defaultdict class.

    :param defaultdict: a nested defaultdict class
    :type defaultdict: defaultdict
    """

    def __init__(self, data_type=int):
        super().__init__(partial(MiddleDefaultDict2, data_type))


class BarcodeQcCounter:
    """A class for counting and processing barcode quality control data.

    Attributes:
        metrics (DefaultDict): A nested defaultdict containing the
            barcode metrics.
        ltr1_seq_dict (DefaultDict): A defaultdict storing the R1
            transposon sequences.

    """
    _metrics: DefaultDict
    _bc_status: DefaultDict

    def __init__(self, pickle_path: str = None) -> None:
        """Initializes a BarcodeQcCounter instance.

        Args:
            pickle_path (str, optional): Path to a pickled BarcodeQcCounter
                object. If provided, loads the object from the file.
                Defaults to None.

        Raises:
            FileNotFoundError: If the provided pickle path does not exist.
        """
        if pickle_path:
            if not os.path.exists(pickle_path):
                msg = f"Path to pickle file {pickle_path} does not exist"
                raise FileNotFoundError(msg)
            self.load(pickle_path)
        else:
            self.metrics = OuterDefaultDict(int)
            self.bc_status = OuterDefaultDict(bool)

    @property
    def metrics(self) -> defaultdict:
        """Returns the _metrics attribute.

        Returns:
            defaultdict: The _metrics attribute.
        """
        return self._metrics

    @metrics.setter
    def metrics(self, value: defaultdict) -> None:
        """Sets the _metrics attribute"""
        self._metrics = value

    @property
    def bc_status(self) -> defaultdict:
        """Returns the _bc_status attribute.

        Returns:
            defaultdict: The _bc_status attribute.
        """
        return self._bc_status

    @bc_status.setter
    def bc_status(self, value: defaultdict) -> None:
        """Sets the _bc_status attribute"""
        self._bc_status = value

    # private methods ---------------------------------------------------------
    def _combine(self, other: "BarcodeQcCounter") -> None:
        """Combine the metrics from another BarcodeQcCounter object.

        Args:
            other (BarcodeQcCounter): Another BarcodeQcCounter object
                whose metrics will be combined with this object.
        """
        def combine_dicts_additive(d1, d2):
            """Recursive function to combine two nested dictionaries."""
            for k, v in d2.items():
                if isinstance(v, dict):
                    d1[k] = combine_dicts_additive(d1.get(k, {}), v)
                else:
                    d1[k] = d1.get(k, 0) + v
            return d1

        def combine_dicts_bool(d1, d2):
            """Recursive function to combine two nested dictionaries."""
            for k, v in d2.items():
                if isinstance(v, dict):
                    d1[k] = combine_dicts_bool(d1.get(k, {}), v)
                else:
                    d1[k] = v
            return d1

        # combine the metrics dictionaries
        self._metrics = combine_dicts_additive(self._metrics, other.metrics)
        # combine the bc_status dictionaries
        self._bc_status = combine_dicts_bool(self._bc_status, other.bc_status)

    # public methods ----------------------------------------------------------

    def load(self, file_path: str) -> None:
        """Load a BarcodeQcCounter object from a file using Pickle.

        Args:
            file_path (str): The file path where the object is stored.
        """
        logger.info("loading BarcodeQcCounter object from %s", file_path)
        with open(file_path, "rb") as file:
            file_data = pickle.load(file)
            if not isinstance(file_data, BarcodeQcCounter):
                raise TypeError(
                    f"{file_path} is not a BarcodeQcCounter object")
            # copy the data from the loaded object to the current instance
            self.metrics = file_data.metrics
            self.bc_status = file_data.bc_status

    @classmethod
    def combine(
            cls, counters: Iterable["BarcodeQcCounter"]) -> "BarcodeQcCounter":
        """Combine multiple BarcodeQcCounter objects into a single object.

        Args:
            counters (Iterable[BarcodeQcCounter]): An iterable of
                BarcodeQcCounter objects.

        Returns:
            BarcodeQcCounter: A new BarcodeQcCounter object with the
                combined metrics.
        """
        result = BarcodeQcCounter()

        for counter in counters:
            result._combine(counter)

        return result

    def __add__(self, other: "BarcodeQcCounter") -> "BarcodeQcCounter":
        """Add two BarcodeQcCounter objects together with the + operator."""
        if not isinstance(other, BarcodeQcCounter):
            raise TypeError("Both objects must be of type 'BarcodeQcCounter'")

        result = BarcodeQcCounter()
        return result.combine([self, other])

    def update(self,
               pb_seq: str,
               ltr1_seq: str,
               ltr2_seq: str,
               srt_seq: str,
               bc_status: bool) -> None:
        """Updates the metrics with given component and deviation tuples.

        Args:
            pb_seq (str): The primer binding sequence.
            ltr1_seq (str): The left transposon sequence.
            ltr2_seq (str): The right transposon sequence.
            srt_seq (str): The sample barcode sequence.
            bc_status (bool): The barcode status.
        """

        (self._metrics
         [pb_seq]
         [ltr1_seq]
         [ltr2_seq]
         [srt_seq]) += 1

        (self._bc_status
         [pb_seq]
         [ltr1_seq]
         [ltr2_seq]
         [srt_seq]) = bc_status

    def write(self,
              filename: str,
              suffix: str = "",
              raw: bool = False) -> None:
        """Write a pickle and/or a comma-delimited file summarizing the
        barcode QC metrics.

        Args:
            filename (str, optional): The base filename for the output files.
                Defaults to "barcode_qc".
            suffix (str, optional): A suffix to be appended to the base
                filename. Defaults to an empty string.
            raw (bool, optional): If True, pickles the object.
                Defaults to False.
        """
        # if raw is true, then pickle the object
        if raw:
            pickle_path = filename + '_' + suffix + '_barcode_qc.pkl'\
                if suffix else filename + '_barcode_qc.pkl'
            logger.info("pickling barcode_qc object to %s{pick_path}")
            with open(pickle_path, "wb") as pickle_file:
                pickle.dump(self, pickle_file)

        else:
            # write the barcode qc metrics to a csv file
            tsv_path = filename + '_' + suffix + "_barcode_qc.tsv" \
                if suffix else filename + '_barcode_qc.tsv'
            logger.info("writing barcode qc metrics to %s", tsv_path)
            with open(tsv_path, "w", encoding='utf-8') as tsv_file:
                csv_writer = csv.writer(tsv_file, delimiter='\t')
                csv_writer.writerow([
                    "pb_seq",
                    "ltr1_seq",
                    "ltr2_seq",
                    "srt_seq",
                    "count",
                    "barcode_status"
                ])

                for pb_seq, ltr1_dict in self._metrics.items():
                    for ltr1_seq, ltr2_dict in ltr1_dict.items():
                        for ltr2_seq, srt_dict in ltr2_dict.items():
                            for srt_seq, count in srt_dict.items():
                                bc_status = ("pass" if
                                             (self._bc_status[pb_seq]
                                              [ltr1_seq]
                                              [ltr2_seq]
                                              [srt_seq])
                                             else "false")
                                csv_writer.writerow([
                                    pb_seq,
                                    ltr1_seq,
                                    ltr2_seq,
                                    srt_seq,
                                    count,
                                    bc_status
                                ])
