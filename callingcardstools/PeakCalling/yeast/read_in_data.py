import os

import pandas as pd
import pyranges as pr


def qbed_df_to_pyranges(df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert a qbed dataframe to a pyranges dataframe.

    :param df: The qbed dataframe.
    :type df: pd.DataFrame
    :return: The pyranges dataframe.
    :rtype: pd.DataFrame

    :Example:

    >>> import pandas as pd
    >>> df = pd.DataFrame({'chr': ['chr1', 'chr1'],
    ...                    'start': [1, 2],
    ...                    'end': [2, 3],
    ...                    'strand': ['+', '-'],
    ...                    'depth': [1, 2]})
    >>> pyranges_df = qbed_df_to_pyranges(df)
    >>> list(pyranges_df.columns) == ['Chromosome', 'Start', 'End', 'Strand', 'depth']
    True
    """
    # check input
    if not isinstance(df, pd.DataFrame):
        raise ValueError("df must be a DataFrame")
    if not all(col in df.columns for col in ["chr", "start", "end", "strand", "depth"]):
        raise ValueError(
            "df must have columns `chr`, `start`, `end`, `strand`, " "and `depth`"
        )

    return pr.PyRanges(
        df.rename(
            {"chr": "Chromosome", "start": "Start", "end": "End", "strand": "Strand"},
            axis=1,
        )
    )


def read_in_chrmap(chrmap_data_path: str, required_col_set: set) -> pd.DataFrame:
    """
    Read in the chrmap data from the `chrmap_data_path` and check that the
    `required_col_set` is a subset of the columns in the chrmap dataframe.

    :param chrmap_data_path: Path to the chrmap file.
    :type chrmap_data_path: str
    :param required_col_set: The required columns in the chrmap dataframe.
    :type required_col_set: set
    :return: The chrmap dataframe.
    :rtype: pd.DataFrame

    :raises ValueError: If the `chrmap_data_path` does not exist or
        is not a file; if the `required_col_set` is not a subset of the
        columns in the chrmap dataframe.

    :Example:

    >>> import pandas as pd
    >>> import os
    >>> import tempfile
    >>> tmp_chrmap = tempfile.NamedTemporaryFile(suffix='.csv').name
    >>> with open(tmp_chrmap, 'w') as f:
    ...    _ = f.write('experiment_orig_chr_convention,'
    ...                'promoter_orig_chr_convention,'
    ...                'background_orig_chr_convention,'
    ...                'unified_chr_convention\\n')
    ...    _ = f.write('chr1,chr1,chr1,chrI\\n')
    >>> chrmap_df = read_in_chrmap(tmp_chrmap,
    ...                            {'experiment_orig_chr_convention',
    ...                             'promoter_orig_chr_convention',
    ...                             'background_orig_chr_convention',
    ...                             'unified_chr_convention'})
    >>> list(chrmap_df.columns) == ['experiment_orig_chr_convention',
    ...                             'promoter_orig_chr_convention',
    ...                             'background_orig_chr_convention',
    ...                             'unified_chr_convention']
    True
    """
    if not os.path.exists(chrmap_data_path):
        raise ValueError(f"chrmap_data_path does " f"not exist: {chrmap_data_path}")
    required_col_set.add("type")
    # read in the chr map
    chrmap_df = pd.read_csv(chrmap_data_path)

    # raise an error if the nrow of the chrmap is 0
    if chrmap_df.shape[0] == 0:
        raise ValueError("The chrmap file is empty")
    # raise an error if the ncol of the chrmap is less than 2
    if chrmap_df.shape[1] < 1:
        raise ValueError(
            "The chrmap file must have at least two column, ",
            "one for the original chromosome name and one for ",
            "the new chromosome name",
        )

    # raise an error if the experiment_orig_chr_convention,
    # promoter_orig_chr_convention,
    # background_orig_chr_convention and unified_chr_convention are not in
    # the chrmap_df columns
    missing_chr_cols = required_col_set.difference(chrmap_df.columns)
    if len(missing_chr_cols) > 0:
        raise ValueError(
            f"The following chromosome columns are missing "
            f"from the chrmap file: {missing_chr_cols}"
        )

    # cast all columns to str
    chrmap_df = chrmap_df.astype("str")

    return chrmap_df


def relabel_chr_column(
    data_df: pd.DataFrame,
    chrmap_df: pd.DataFrame,
    curr_chr_name_convention: str,
    new_chr_name_convention: str,
    genomic_only: bool = True,
) -> pd.DataFrame:
    """
    Given a `data_df` with column `chr`, a `curr_chr_name_convention` and
    a `new_chr_name_convention`, that are both columns of `chrmap_df`, join
    the `chrmap` to the `data_df` based on the `curr_chr_name_convention` and
    swap the values in the `chr` column to the `new_chr_name_convention`.
    relabel the `new_chr_name_convention` to `chr` and return the dataframe
    with columns in the same order as the input dataframe. If `geonmic_only` is
    set to `True`, only those chromosomes with type == 'genomic' are returned.

    :param df: The dataframe to relabel.
    :type df: pd.DataFrame
    :param curr_chr_name_convention: The current chromosome name convention.
    :type curr_chr_name_convention: str
    :param new_chr_name_convention: The new chromosome name convention.
    :type new_chr_name_convention: str
    :param genomic_only: Whether to return only records with type == "genomic".
    :type genomic_only: bool

    :return: The relabeled dataframe.
    :rtype: pd.DataFrame

    :raises ValueError: If the `curr_chr_name_convention` or
        `new_chr_name_convention` are not columns in `chrmap_df`.

    :Example:

    >>> import pandas as pd
    >>> data_df = pd.DataFrame({'chr': ['chr1', 'chr2', 'chr3'],
    ...                         'start': [1, 2, 3],})
    >>> chrmap_df = pd.DataFrame({'curr_chr_name_convention':
    ...                            ['chr1', 'chr2', 'chr3'],
    ...                           'new_chr_name_convention':
    ...                            ['chrI', 'chrII', 'chrIII']})
    >>> relabeled_df = relabel_chr_column(data_df, chrmap_df,
    ...                                   'curr_chr_name_convention',
    ...                                   'new_chr_name_convention')
    >>> list(relabeled_df.columns) == ['chr', 'start']
    True
    >>> list(relabeled_df['chr']) == ['chrI', 'chrII', 'chrIII']
    True
    """
    # check input
    if "chr_curr" in chrmap_df.columns:
        raise ValueError(
            "chr_curr cannot be a column in chrmap_df for the "
            "purposes of relabelling. rename that column in "
            "chrmap_df and resubmit"
        )
    if curr_chr_name_convention not in chrmap_df.columns:
        raise ValueError("curr_chr_name_convention " "must be a column in chrmap_df")
    if new_chr_name_convention not in chrmap_df.columns:
        raise ValueError("new_chr_name_convention " "must be a column in chrmap_df")

    # rename the current chr column to chr_curr to avoid any errors in
    # joining, if the old/new format is called 'chr'
    data_df = data_df.rename(columns={"chr": "chr_curr"})

    if curr_chr_name_convention != new_chr_name_convention:
        # join a subset of chrmap_df -- only the columns we need -- to data_df
        data_df = data_df.merge(
            chrmap_df[[curr_chr_name_convention, new_chr_name_convention, "type"]],
            left_on="chr_curr",
            right_on=curr_chr_name_convention,
        )
        # swap values in chr column
        data_df["chr_curr"] = data_df[new_chr_name_convention]
    else:
        data_df = data_df.merge(
            chrmap_df[[curr_chr_name_convention, "type"]],
            left_on="chr_curr",
            right_on=curr_chr_name_convention,
        )

    if genomic_only:
        data_df = data_df.query("type=='genomic'")

    return data_df.drop(
        columns=[new_chr_name_convention, curr_chr_name_convention, "type"]
    ).rename(columns={"chr_curr": "chr"})


def read_in_experiment_data(
    experiment_data_path: str,
    curr_chr_name_convention: pd.DataFrame,
    new_chr_name_convention: pd.DataFrame,
    chrmap_df: str,
    deduplicate: bool = True,
    **kwargs,
) -> pd.DataFrame:
    """
    Read in experiment (hops) data from a qbed file. The qbed file may be
    plain text or gzipped and may or may not have column headers. If the
    column headers are present, they must be in the following order:
    `chr`, `start`, `end`, `strand`, `depth`. If the column headers are
    not present, the columns must be in same order and number. Datatypes
    are checked but will not be coerced -- errors are raised if they do not
    match the expected datatypes. the `chr` column is relabeled from the
    `curr_chr_name_convention` to the `new_chr_name_convention` using the
    `chrmap_df`.

    Additional keyword arguments:
        - genomic_only (bool): Whether to return only records with type == "genomic".
            See `relabel_chr_column` for more information. Defaults to True.

    :param experiment_data_path: Path to the qbed file, plain text or gzipped,
        with or without column headers
    :type experiment_data_path: str
    :param curr_chr_name_convention: The current chromosome name convention.
    :type curr_chr_name_convention: str
    :param new_chr_name_convention: The new chromosome name convention.
    :type new_chr_name_convention: str
    :param chrmap_df: The chrmap dataframe.
    :type chrmap_df: pd.DataFrame
    :param deduplicate: Whether to deduplicate the experiment data based on
        `chr`, `start`, `end` such that if an insertion occurs at the same
        coordinate but on opposite strands, only one record is retained.
    :type deduplicate: bool

    :return: The experiment data as a dataframe with the `chr` column
        refactored to the `new_chr_name_convention`
    :rtype: pd.DataFrame

    :raises ValueError: If the `experiment_data_path` does not exist or
        is not a file; if the column headers exist but do not match expectation
        or if the datatypes do not match expectation.

    :Example:

    >>> import pandas as pd
    >>> import os
    >>> import tempfile
    >>> tmp_qbed = tempfile.NamedTemporaryFile(suffix='.qbed').name
    >>> with open(tmp_qbed, 'w') as f:
    ...    _ = f.write('chr\\tstart\\tend\\tstrand\\tdepth\\n')
    ...    _ = f.write('chr1\\t1\\t2\\t+\\t1\\n')
    >>> # create a temporary chrmap file
    >>> tmp_chrmap = tempfile.NamedTemporaryFile(suffix='.csv').name
    >>> chrmap_df = pd.DataFrame({'curr_chr_name_convention':
    ...                            ['chr1', 'chr2', 'chr3'],
    ...                           'new_chr_name_convention':
    ...                            ['chrI', 'chrII', 'chrIII']})
    >>> # read in the data
    >>> experiment_df, experiment_total_hops = read_in_experiment_data(
    ...   tmp_qbed,
    ...   'curr_chr_name_convention',
    ...   'new_chr_name_convention',
    ...    chrmap_df)
    >>> list(experiment_df.columns) == ['chr', 'start', 'end', 'depth',
    ...                                 'strand']
    True
    >>> experiment_total_hops
    1
    """
    # check input
    if not os.path.exists(experiment_data_path):
        raise ValueError("experiment_data_path must exist")
    if not os.path.isfile(experiment_data_path):
        raise ValueError("experiment_data_path must be a file")

    # check if data is gzipped
    gzipped = str(experiment_data_path).endswith(".gz")
    # check if data has column headers
    header = pd.read_csv(
        experiment_data_path, sep="\t", compression="gzip" if gzipped else None, nrows=0
    )
    if header.columns.tolist() != ["chr", "start", "end", "depth", "strand"]:
        header = None
    else:
        header = 0
    # read in data
    try:
        experiment_df = pd.read_csv(
            experiment_data_path,
            sep="\t",
            header=header,
            names=["chr", "start", "end", "depth", "strand"],
            dtype={"chr": str, "start": int, "end": int, "depth": int, "strand": str},
            compression="gzip" if gzipped else None,
        )
    except ValueError as e:
        raise ValueError(
            "experiment_data_path must be a qbed file "
            "with columns `chr`, `start`, `end`, `strand`, "
            "and `depth`"
        ) from e

    # if the file is empty, raise an error. This needs to be handled in an outer
    # function in order to exit gracefully if there is no actual data to do
    # calculations on
    if experiment_df.shape[0] == 0:
        raise ValueError("The experiment file is empty -- no data to process")

    # relabel chr column
    experiment_df = relabel_chr_column(
        experiment_df,
        chrmap_df,
        curr_chr_name_convention,
        new_chr_name_convention,
        kwargs.get("genomic_only", True),
    )

    if deduplicate:
        # deduplicate based on `chr`, `start`, `end` such that if an insertion occurs
        # at the same coordinate but on opposite strands, only one record is retained
        experiment_df.drop_duplicates(subset=["chr", "start", "end"], inplace=True)
        # set strand to '*'
        experiment_df["strand"] = "*"

    return qbed_df_to_pyranges(experiment_df), len(experiment_df)


def read_in_promoter_data(
    promoter_data_path: str,
    curr_chr_name_convention: pd.DataFrame,
    new_chr_name_convention: pd.DataFrame,
    chrmap_df: str,
) -> pd.DataFrame:
    """
    Read in the promoter data. The promoter data should be a tsv with extension
    `.bed` or `.bed.gz` and should have the following columns:
    `chr` `start`   `end`  `name` `score` `strand`. The `chr` column is
    refactored from the `curr_chr_name_convention` to the
    `new_chr_name_convention` using the `chrmap_df`.

    :param promoter_data_path: Path to the promoter bed file (plus colnames)
    :type promoter_data_path: str
    :param curr_chr_name_convention: The current chromosome name convention.
    :type curr_chr_name_convention: str
    :param new_chr_name_convention: The new chromosome name convention.
    :type new_chr_name_convention: str
    :param chrmap_df: The chrmap dataframe.
    :type chrmap_df: pd.DataFrame
    :return: The promoter data as a dataframe with the `chr` column
        refactored to the `new_chr_name_convention`
    :rtype: pd.DataFrame

    :raises ValueError: If the `promoter_data_path` does not exist or
        is not a file; if the column headers exist but do not match expectation
        or if the datatypes do not match expectation.

    :Example:

    >>> import pandas as pd
    >>> import os
    >>> import tempfile
    >>> tmp_bed = tempfile.NamedTemporaryFile(suffix='.bed').name
    >>> with open(tmp_bed, 'w') as f:
    ...    _ = f.write('chr\\tstart\\tend\\tname\\tscore\\tstrand\\n')
    ...    _ = f.write('chr1\\t1\\t2\\ttest\\t1\\t+\\n')
    >>> chrmap_df = pd.DataFrame({'curr_chr_name_convention':
    ...                            ['chr1', 'chr2', 'chr3'],
    ...                           'new_chr_name_convention':
    ...                            ['chrI', 'chrII', 'chrIII']})
    >>> promoter_df = read_in_promoter_data(
    ... tmp_bed,
    ... 'curr_chr_name_convention',
    ... 'new_chr_name_convention',
    ... chrmap_df)
    >>> list(promoter_df.columns) == ['chr', 'start', 'end', 'name',
    ...                               'score', 'strand']
    True
    >>> len(promoter_df) == 1
    True
    """
    # check input
    if not os.path.exists(promoter_data_path):
        raise ValueError("promoter_data_path must exist")
    if not os.path.isfile(promoter_data_path):
        raise ValueError("promoter_data_path must be a file")

    # check if data is gzipped
    gzipped = str(promoter_data_path).endswith(".gz")
    # check if data has column headers
    header = pd.read_csv(
        promoter_data_path, sep="\t", compression="gzip" if gzipped else None, nrows=0
    )
    if header.columns.tolist() != ["chr", "start", "end", "name", "score", "strand"]:
        header = None
    else:
        header = 0
    # read in data
    try:
        promoter_df = pd.read_csv(
            promoter_data_path,
            sep="\t",
            header=header,
            names=["chr", "start", "end", "name", "score", "strand"],
            dtype={
                "chr": str,
                "start": int,
                "end": int,
                "name": str,
                "score": float,
                "strand": str,
            },
            compression="gzip" if gzipped else None,
        )
    except ValueError as e:
        raise ValueError(
            "promoter_data_path must be a bed file "
            "with columns `chr`, `start`, `end`, `name`, "
            "`score`, and `strand`"
        ) from e

    # if the file is empty, raise an error.
    if promoter_df.shape[0] == 0:
        raise ValueError("The promoter file is empty -- no data to process")

    # relabel chr column
    chr_relabeled_promoter_df = relabel_chr_column(
        promoter_df, chrmap_df, curr_chr_name_convention, new_chr_name_convention
    )

    return chr_relabeled_promoter_df


def read_in_background_data(
    background_data_path: str,
    curr_chr_name_convention: pd.DataFrame,
    new_chr_name_convention: pd.DataFrame,
    chrmap_df: str,
    **kwargs,
) -> pd.DataFrame:
    """
    Read in background (hops) data from a qbed file. The qbed file may be
    plain text or gzipped and may or may not have column headers. If the
    column headers are present, they must be in the following order:
    `chr`, `start`, `end`, `strand`, `depth`. If the column headers are
    not present, the columns must be in same order and number. Datatypes
    are checked but will not be coerced -- errors are raised if they do not
    match the expected datatypes. the `chr` column is relabeled from the
    `curr_chr_name_convention` to the `new_chr_name_convention` using the
    `chrmap_df`. NOTE: unlike the experiment df, there is no option to deduplicate
    as the background file is expected to be the combination of multiple experiments
    at this point.

    Additional keyword arguments:
        - genomic_only (bool): Whether to return only records with type == "genomic".
            See `relabel_chr_column` for more information. Defaults to True.

    :param background_data_path: Path to the background data qbed file, plain
        text or gzipped, with or without column headers
    :type background_data_path: str
    :param curr_chr_name_convention: The current chromosome name convention
    :type curr_chr_name_convention: str
    :param new_chr_name_convention: The new chromosome name convention
    :type new_chr_name_convention: str
    :param chrmap_df: The chrmap dataframe
    :type chrmap_df: pd.DataFrame
    :return: The background data.
    :rtype: pd.DataFrame

    :raises ValueError: If the `background_data_path` does not exist or
        is not a file; if the column headers exist but do not match expectation
        or if the datatypes do not match expectation.

    :Example:

    >>> import pandas as pd
    >>> import os
    >>> import tempfile
    >>> tmp_qbed = tempfile.NamedTemporaryFile(suffix='.qbed').name
    >>> with open(tmp_qbed, 'w') as f:
    ...    _ = f.write('chr\\tstart\\tend\\tstrand\\tdepth\\n')
    ...    _ = f.write('chr1\\t1\\t2\\t+\\t1\\n')
    >>> # create a temporary chrmap file
    >>> chrmap_df = pd.DataFrame({'curr_chr_name_convention':
    ...                            ['chr1', 'chr2', 'chr3'],
    ...                           'new_chr_name_convention':
    ...                            ['chrI', 'chrII', 'chrIII']})
    >>> # read in the data
    >>> background_df, background_total_hops = read_in_background_data(
    ...   tmp_qbed,
    ...   'curr_chr_name_convention',
    ...   'new_chr_name_convention',
    ...    chrmap_df)
    >>> list(background_df.columns) == ['chr', 'start', 'end', 'depth',
    ...                                  'strand']
    True
    >>> background_total_hops
    1
    """
    # check input
    if not os.path.exists(background_data_path):
        raise ValueError("background_data_path must exist")
    if not os.path.isfile(background_data_path):
        raise ValueError("background_data_path must be a file")

    # check if data is gzipped
    gzipped = str(background_data_path).endswith(".gz")
    # check if data has column headers
    header = pd.read_csv(background_data_path, sep="\t", nrows=0)
    if header.columns.tolist() != ["chr", "start", "end", "depth", "strand"]:
        header = None
    else:
        header = 0

    # read in data
    try:
        background_df = pd.read_csv(
            background_data_path,
            sep="\t",
            header=header,
            names=["chr", "start", "end", "depth", "strand"],
            dtype={
                "chr": str,
                "start": int,
                "end": int,
                "depth": "int64",
                "strand": str,
            },
            compression="gzip" if gzipped else None,
        )
    except ValueError as e:
        raise ValueError(
            "background_data_path must be a qbed file "
            "with columns `chr`, `start`, `end`, `depth`, "
            "and `strand`"
        ) from e

    # if the file is empty, raise an error. Background data should never be empty
    if background_df.shape[0] == 0:
        raise ValueError("The background file is empty -- no data to process")

    # relabel chr column
    background_df = relabel_chr_column(
        background_df,
        chrmap_df,
        curr_chr_name_convention,
        new_chr_name_convention,
        kwargs.get("genomic_only", True),
    )

    return qbed_df_to_pyranges(background_df), len(background_df)
