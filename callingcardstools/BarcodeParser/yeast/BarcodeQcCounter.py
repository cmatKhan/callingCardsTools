"""
This module contains the BarcodeQcCounter class, which is used for analyzing
barcode quality control data. The class can process and update barcode metrics,
combine multiple objects, and write the results to output files.

Author: Chase Mateusiak
Date: 2023-05-01
"""
import logging
import os
import pickle
from collections import defaultdict
from functools import partial
from typing import DefaultDict, Iterable

import pandas as pd  # pylint:disable=E0401
from edlib import align  # pylint:disable=E0401,E0611

logger = logging.getLogger(__name__)

# TODO: use this to re-write the barcodeparser object. reduce the barcode
# parse to common attr and methods for any organism, and then extend it
# for mammals/yeast. use this as the yeast basis. Do not align the barcode
# components more than once in the iteration of the fastq file.


class InnerDefaultDict(defaultdict):
    def __init__(self, data_type=int):
        super().__init__(data_type)


class MiddleDefaultDict1(defaultdict):
    def __init__(self, data_type=int):
        super().__init__(partial(InnerDefaultDict, data_type))


class MiddleDefaultDict2(defaultdict):
    def __init__(self, data_type=int):
        super().__init__(partial(MiddleDefaultDict1, data_type))


class OuterDefaultDict(defaultdict):
    def __init__(self, data_type=int):
        super().__init__(partial(MiddleDefaultDict2, data_type))


class BarcodeQcCounter:
    """A class for counting and processing barcode quality control data.

    Attributes:
        metrics (DefaultDict): A nested defaultdict containing the 
            barcode metrics.
        r1_transposon_seq_dict (DefaultDict): A defaultdict storing the R1 
            transposon sequences.

    """
    _metrics: DefaultDict
    _r1_transposon_seq_dict: DefaultDict

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
            self._metrics = OuterDefaultDict(int)
            self._r1_transposon_seq_dict = defaultdict(set)

    @property
    def metrics(self) -> defaultdict:
        """Returns the _metrics attribute.

        Returns:
            defaultdict: The _metrics attribute.
        """
        return self._metrics

    @property
    def r1_transposon_dict(self) -> defaultdict:
        """Returns the _r1_transposon_seq_dict attribute.

        Returns:
            defaultdict: the _r1_transposon_seq_dict attribute.
        """
        return self._r1_transposon_seq_dict

    # private methods ---------------------------------------------------------

    def _combine(self, other: "BarcodeQcCounter") -> None:
        """Combine the metrics from another BarcodeQcCounter object.

        Args:
            other (BarcodeQcCounter): Another BarcodeQcCounter object
                whose metrics will be combined with this object.
        """
        # Combine _metrics dictionaries
        for r1_transposon_edit_dist, r1_primer_dict in other.metrics.items():
            for r1_primer_seq, r2_transposon_dict in r1_primer_dict.items():
                for r2_transposon_seq, r2_restriction_enzyme_dict in \
                        r2_transposon_dict.items():
                    for r2_restriction_enzyme_name, count in \
                            r2_restriction_enzyme_dict.items():
                        (self._metrics
                         [r1_transposon_edit_dist]
                         [r1_primer_seq]
                         [r2_transposon_seq]
                         [r2_restriction_enzyme_name]) += count

        # Combine _r1_transposon_seq_dict dictionaries
        for r1_transposon_edit_dist, r1_transposon_seq_set in \
                other.r1_transposon_dict.items():
            self._r1_transposon_seq_dict[r1_transposon_edit_dist]\
                .update(r1_transposon_seq_set)

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
            self._metrics = file_data._metrics
            self._r1_transposon_seq_dict = file_data._r1_transposon_seq_dict

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
        result.combine(self)
        result.combine(other)
        return result

    def update(self,
               component_tuple: tuple,
               r1_transposon_edit_dist: int,
               r2_restriction_enzyme_name: str) -> None:
        """Updates the metrics with given component and deviation tuples.

        Args:
            component_tuple (tuple): A tuple containing R1 primer,
                R1 transposon, and R2 transposon sequences.
            r1_transposon_edit_dist (int): The edit distance between the
                R1 transposon sequence and the expected R1 transposon
            r2_restriction_enzyme_name (str): The R2 restriction enzyme name.
        """
        (r1_primer_seq,
         r1_transposon_seq,
         r2_transposon_seq) = component_tuple

        (self._metrics
         [r1_transposon_edit_dist]
         [r1_primer_seq]
         [r2_transposon_seq]
         [r2_restriction_enzyme_name]) += 1

        self._r1_transposon_seq_dict[r1_transposon_edit_dist]\
            .add(r1_transposon_seq)

    def _summarize_by_tf(self, component_dict: dict) -> None:
        """Summarizes the metrics by transcription factor (TF).

        Args:
            component_dict (dict): A dictionary containing keys for 
                'tf', 'r1_primers', and 'r2_transposons', and their 
                respective lists of values.

        Returns:
            tuple: A tuple containing R1 primer summary and R2 
                transposon summary.
        """
        #
        r1_primer_summary = []
        r2_transposon_summary = []
        # only iterate over those reads which had an r1 transposon seq
        # edit distance of n or less

        # r1_for_given_r2_dict = defaultdict(lambda: defaultdict(set))
        r1_for_given_r2_dict = MiddleDefaultDict1(set)
        for i, r1_transposon_dict in self._metrics.items():
            # first level of iteration is over the r1 primer keys.
            # The dictionary is a nested dictionary with the keys being
            # r2_transposon sequences and values another dicitonary with
            # the restriciton enzyme and count
            for r1_primer_seq, r1_primer_dict in r1_transposon_dict.items():
                # if the r1 primer sequence is the expected sequence
                # for a given tf, then iterate over the r2 transposon
                # entries and record the results
                for r2_transposon_seq, r2_transposon_seq_dict in \
                        r1_primer_dict.items():
                    # if the r2_transposon_seq is recognized, then save the
                    # r1_primer_seq. structure of the dict is:
                    # {'valid_r2_trans_seq': set(r1_primer_seq1, ...)}
                    if r2_transposon_seq in \
                            component_dict['r2_transposon']:
                        (r1_for_given_r2_dict
                         [i]
                         [r2_transposon_seq]
                         .add(r1_primer_seq))
                    # if the r1_primer_seq is an expected sequence, then
                    # iterate over the r2_transposon_seq_dict and record the
                    # results
                    if r1_primer_seq in component_dict['r1_primer']:
                        r1_primer_index = \
                            component_dict['r1_primer'].index(r1_primer_seq)
                        r2_transposon_target_seq = \
                            component_dict['r2_transposon'][r1_primer_index]
                        edit_dist = \
                            align(
                                r2_transposon_seq,
                                r2_transposon_target_seq)
                        r1_primer_record = {
                            "tf": component_dict['tf'][r1_primer_index],
                            "r1_primer_seq":
                            component_dict['r1_primer'][r1_primer_index],
                            "r1_transposon_edit_dist": i,
                            "r2_transposon_seq": r2_transposon_seq,
                            "r2_transposon_edit_dist":
                            edit_dist.get("editDistance")}
                        for restriction_enzyme, count in \
                                r2_transposon_seq_dict.items():
                            record_copy = r1_primer_record.copy()
                            record_copy.update({
                                'restriction_enzyme':
                                restriction_enzyme,
                                'count': count})
                            r1_primer_summary.append(record_copy)

        # in the second iteration, iterate over only those r1_primer_seqs with
        # a valid r2_transposon_seq
        for r1_transposon_ed, r1_transposon_ed_dict in \
                r1_for_given_r2_dict.items():
            for r2_transposon_seq, r1_primer_seq_set in \
                    r1_transposon_ed_dict.items():
                # extract the TF and expected r1_primer sequence for this
                # r2_transposon_seq and TF
                index = component_dict['r2_transposon']\
                    .index(r2_transposon_seq)
                tf = component_dict['tf'][index]
                r1_primer_expected = component_dict["r1_primer"][index]
                # iterate over all of the `r1_primer_seq` for this
                # r2_transposon_seq
                for r1_primer_query in r1_primer_seq_set:
                    # align the r1_primer to the expected r1_primer for this
                    # r2_transposon_seq and TF
                    edit_dist = align(
                        r1_primer_query,
                        r1_primer_expected)
                    # create the base record
                    r2_transposon_record = {
                        "tf": tf,
                        "r2_transposon_seq": r2_transposon_seq,
                        "r1_transposon_edit_dist": r1_transposon_ed,
                        "r1_primer_seq": r1_primer_query,
                        "r1_primer_edit_dist":
                        edit_dist.get("editDistance")}
                    for restriction_enzyme, count in \
                        (self._metrics
                         [r1_transposon_ed]
                         [r1_primer_query]
                         [r2_transposon_seq]
                         .items()):
                        # make a copy of the record
                        record_copy = r2_transposon_record.copy()
                        # add additional restriction enzyme info
                        record_copy.update({
                            'restriction_enzyme':
                            restriction_enzyme,
                            'count': count})
                        r2_transposon_summary.append(record_copy)

        return r1_primer_summary, r2_transposon_summary

    def write(self,
              raw: bool = False,
              component_dict: dict = None,
              output_dirpath: str = ".",
              filename: str = "barcode_qc",
              suffix: str = "") -> None:
        """Write a pickle and/or a comma-delimited file summarizing the
        barcode QC metrics.

        Args:
            raw (bool, optional): If True, pickles the object.
                Defaults to False.
            component_dict (dict, optional): A dictionary containing keys
                for 'tf', 'r1_primers', and 'r2_transposons', and their
                respective lists of values. If provided, writes summaries
                for each component. Defaults to None.
            output_dirpath (str, optional): The output directory path where
                the files will be saved. Defaults to the current directory.
            filename (str, optional): The base filename for the output files.
                Defaults to "barcode_qc".
            suffix (str, optional): A suffix to be appended to the base
                filename. Defaults to an empty string.
        """
        # check that the output_dirpath is a valid directory
        if not os.path.join(output_dirpath):
            raise ValueError("output_dirpath must be a valid directory")
        # if raw is true, then pickle the object
        if raw:
            pickle_path = os.path.join(
                output_dirpath, filename + '_' + suffix + ".pickle")
            logger.info("pickling barcode_qc object to %s{pick_path}")
            with open(pickle_path, "wb") as pickle_file:
                pickle.dump(self, pickle_file)

        # if component_dict is passed
        if component_dict:
            # input checks
            if not isinstance(component_dict, dict):
                raise TypeError("component_dict must be a dictionary")
            if not {'tf', 'r1_primer', 'r2_transposon'} == \
                    set(list(component_dict.keys())):
                raise ValueError("component_dict must be a dictionary "
                                 "where the keys are 'tf', 'r1_primers', "
                                 "'r2_transposons' and the values are "
                                 "lists of the same length. The index of "
                                 "each list corresponds to the same "
                                 "transcription factor.")
            for k, v in component_dict.items():
                if not isinstance(v, list):
                    raise TypeError("component_dict values must be lists")
            if len({len(x) for x in component_dict.values()}) != 1:
                raise ValueError("component_dict values must be lists of "
                                 "the same length")
            # extract summaries from the metrics
            r1_primer_summary, r2_transposon_summary = \
                self._summarize_by_tf(component_dict)

            # write r1_primer_summary to file
            append_suffix = '_' + suffix if suffix else ''
            r1_primer_basename = \
                filename + "_r1_primer_summary" + append_suffix + ".csv"
            r1_primer_summary_path = os.path.join(
                output_dirpath, r1_primer_basename)
            r1_primer_summary_df = pd.DataFrame(r1_primer_summary)
            logger.info("writing r1_primer_summary "
                        "to %s{r1_primer_summary_path}")
            r1_primer_summary_df.to_csv(r1_primer_summary_path, index=False)

            # write r2_transposon summary to file
            r2_transposon_summary_basename = \
                filename + "_r2_transposon_summary" + append_suffix + ".csv"
            r2_transposon_summary_path = os.path.join(
                output_dirpath, r2_transposon_summary_basename)
            r2_transposon_summary_df = pd.DataFrame(r2_transposon_summary)
            logger.info("writing r2_transposon_summary "
                        "to %s{r2_transposon_summary_path}")
            r2_transposon_summary_df.to_csv(
                r2_transposon_summary_path, index=False)