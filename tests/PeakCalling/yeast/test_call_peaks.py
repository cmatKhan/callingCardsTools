import argparse
import os

import pandas as pd
import pyranges as pr

from callingcardstools.PeakCalling.yeast.call_peaks import count_hops
from callingcardstools.PeakCalling.yeast.call_peaks import \
    main as call_peaks_main


def test_count_hops():
    # Setup
    promoter_df = pd.DataFrame(
        {
            "chr": ["chr1", "chr1"],
            "start": [100, 300],
            "end": [200, 400],
            "strand": ["+", "-"],
            "name": ["prom1", "prom2"],
        }
    )

    promoter_pr = pr.PyRanges(
        promoter_df.copy().rename(
            columns={
                "chr": "Chromosome",
                "start": "Start",
                "end": "End",
                "strand": "Strand",
            }
        )
    )

    # the first and last overlap the two promoter regions respectively on the same
    # strand, the second overlaps the first promoter region on the opposite strand
    experiment_qbed_df = pd.DataFrame(
        {
            "chr": ["chr1", "chr1", "chr1"],
            "start": [150, 152, 300],
            "end": [151, 153, 301],
            "strand": ["+", "-", "-"],
        }
    )

    experiment_pr = pr.PyRanges(
        experiment_qbed_df.rename(
            columns={
                "chr": "Chromosome",
                "start": "Start",
                "end": "End",
                "strand": "Strand",
            }
        )
    )

    # the first promoter does not overlap any background regions, the second overlaps
    # the the second promoter region
    background_qbed_df = pd.DataFrame(
        {
            "chr": ["chr1", "chr1"],
            "start": [200, 400],
            "end": [201, 401],
            "strand": ["-", "-"],
        }
    )
    background_pr = pr.PyRanges(
        background_qbed_df.rename(
            columns={
                "chr": "Chromosome",
                "start": "Start",
                "end": "End",
                "strand": "Strand",
            }
        )
    )

    expected_result_stranded = pd.DataFrame(
        {
            "name": ["prom1", "prom2"],
            "chr": ["chr1", "chr1"],
            "start": [100, 300],
            "end": [200, 400],
            "strand": ["+", "-"],
            "experiment_hops": [1, 1],
            "background_hops": [0, 1],
        }
    )
    experiment_result_stranded = count_hops(
        promoter_pr, experiment_pr, "experiment_hops", strandedness="same"
    )

    background_result_stranded = count_hops(
        promoter_pr, background_pr, "background_hops", strandedness="same"
    )

    result_stranded = (
        promoter_df.set_index("name")
        .join(
            [
                experiment_result_stranded.set_index("name", drop=True),
                background_result_stranded.set_index("name", drop=True),
            ],
            how="left",
        )
        .reset_index()
        .fillna(0)
    )

    assert (
        result_stranded.columns
        == [
            "name",
            "chr",
            "start",
            "end",
            "strand",
            "experiment_hops",
            "background_hops",
        ]
    ).all()

    pd.testing.assert_frame_equal(
        result_stranded, expected_result_stranded, check_dtype=False
    )

    expected_result_unstranded = pd.DataFrame(
        {
            "name": ["prom1", "prom2"],
            "chr": ["chr1", "chr1"],
            "start": [100, 300],
            "end": [200, 400],
            "strand": ["+", "-"],
            "experiment_hops": [2, 1],
            "background_hops": [1, 1],
        }
    )
    experiment_hops_unstranded = count_hops(
        promoter_pr, experiment_pr, "experiment_hops", strandedness=False
    )
    background_hops_unstranded = count_hops(
        promoter_pr, background_pr, "background_hops", strandedness=False
    )

    result_unstranded = (
        promoter_df.set_index("name", drop=True)
        .join(
            [
                experiment_hops_unstranded.set_index("name", drop=True),
                background_hops_unstranded.set_index("name", drop=True),
            ],
            how="left",
        )
        .reset_index()
        .fillna(0)
    )

    pd.testing.assert_frame_equal(
        result_unstranded, expected_result_unstranded, check_dtype=False
    )


def test_with_data(tmpdir):
    test_data_directory = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
        "test_data/yeast/Analysis",
    )

    assert os.path.isdir(test_data_directory) is True

    output_path = os.path.join(tmpdir, "test_call_peaks_output.csv")

    args = argparse.Namespace(
        experiment_data_path=os.path.join(test_data_directory, "hap5_expr17.qbed.gz"),
        experiment_orig_chr_convention="id",
        promoter_data_path=os.path.join(test_data_directory, "yiming_promoters.bed.gz"),
        promoter_orig_chr_convention="id",
        background_data_path=os.path.join(
            test_data_directory, "adh1_background.qbed.gz"
        ),
        background_orig_chr_convention="ucsc",
        chrmap_data_path=os.path.join(test_data_directory, "chrmap.csv.gz"),
        unified_chr_convention="ucsc",
        output_path=output_path,
        deduplicate_experiment=True,
        pseudocount=0.2,
        compress_output=False,
    )

    call_peaks_main(args)

    assert os.path.exists(output_path) is True

    output_df = pd.read_csv(output_path)
    experiment_df = pd.read_csv(args.experiment_data_path, sep="\t")
    background_df = pd.read_csv(args.background_data_path, sep="\t")

    assert (
        output_df.columns
        == [
            "name",
            "chr",
            "start",
            "end",
            "strand",
            "experiment_hops",
            "background_hops",
            "background_total_hops",
            "experiment_total_hops",
            "callingcards_enrichment",
            "poisson_pval",
            "hypergeometric_pval",
        ]
    ).all()

    # check that the deduplication worked as expected
    assert (output_df["experiment_total_hops"] != experiment_df.shape[0]).all()
    assert (
        output_df["experiment_total_hops"]
        == experiment_df.drop_duplicates(subset=["chr", "start", "end"]).shape[0]
    ).all()

    # filter experiment_df to only include rows where chr==2 and the start is
    # between 36350 and 37050 inclusive
    experiment_df_subset = experiment_df[
        (experiment_df["chr"] == 2)
        & (experiment_df["start"] >= 36350)
        & (experiment_df["start"] <= 37050)
    ]
    # do the same with the background_df
    background_df_subset = background_df[
        (background_df["chr"] == "chrII")
        & (background_df["start"] >= 36350)
        & (background_df["start"] <= 37050)
    ]

    # filter output_df to chr==chrII, start== 36350 and end == 37050
    # and check that the output_df_subset['experiment_hops'] is equal to the deduplicated
    # number of rows in experiment_df, and not equal to un-deduplicated number of rows
    output_df_subset = output_df[
        (output_df["chr"] == "chrII")
        & (output_df["start"] == 36350)
        & (output_df["end"] == 37050)
    ]

    assert (output_df["background_total_hops"] == background_df.shape[0]).all()
    assert (output_df_subset["background_hops"] == background_df_subset.shape[0]).all()
    assert (
        output_df_subset["experiment_hops"]
        == experiment_df_subset.drop_duplicates(subset=["chr", "start", "end"]).shape[0]
    ).all()
    assert (output_df_subset["experiment_hops"] != experiment_df_subset.shape[0]).all()


def test_combine_replicates(tmpdir):
    test_data_directory = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
        "test_data/yeast/Analysis",
    )

    assert os.path.isdir(test_data_directory) is True

    cc_rep1 = os.path.join(
        test_data_directory, "callingcards_replicates/ccexperiment_292_hap5_chrI.csv.gz"
    )
    assert os.path.exists(cc_rep1) is True
    cc_rep2 = os.path.join(
        test_data_directory, "callingcards_replicates/ccexperiment_297_hap5_chrI.csv.gz"
    )
    assert os.path.exists(cc_rep2) is True

    output_path1 = os.path.join(tmpdir, "peaks_rep1.csv")

    args = argparse.Namespace(
        experiment_data_path=cc_rep1,
        experiment_orig_chr_convention="ucsc",
        promoter_data_path=os.path.join(test_data_directory, "yiming_promoters.bed.gz"),
        promoter_orig_chr_convention="id",
        background_data_path=os.path.join(
            test_data_directory, "adh1_background.qbed.gz"
        ),
        background_orig_chr_convention="ucsc",
        chrmap_data_path=os.path.join(test_data_directory, "chrmap.csv.gz"),
        unified_chr_convention="ucsc",
        output_path=output_path1,
        deduplicate_experiment=True,
        pseudocount=0.2,
        compress_output=False,
    )

    call_peaks_main(args)

    assert os.path.exists(output_path1) is True

    output_path2 = os.path.join(tmpdir, "peaks_rep2.csv")

    args = argparse.Namespace(
        experiment_data_path=cc_rep2,
        experiment_orig_chr_convention="ucsc",
        promoter_data_path=os.path.join(test_data_directory, "yiming_promoters.bed.gz"),
        promoter_orig_chr_convention="id",
        background_data_path=os.path.join(
            test_data_directory, "adh1_background.qbed.gz"
        ),
        background_orig_chr_convention="ucsc",
        chrmap_data_path=os.path.join(test_data_directory, "chrmap.csv.gz"),
        unified_chr_convention="ucsc",
        output_path=output_path2,
        deduplicate_experiment=True,
        pseudocount=0.2,
        compress_output=False,
    )

    call_peaks_main(args)

    assert os.path.exists(output_path2) is True

    # note: ensure that the combined peaks are deduplicated based on the chr, start, and
    # end columns
    df1 = pd.read_csv(cc_rep1, sep="\t")
    df1 = df1.drop_duplicates(subset=["chr", "start", "end"])
    df2 = pd.read_csv(cc_rep2, sep="\t")
    df2 = df2.drop_duplicates(subset=["chr", "start", "end"])

    combined_df = pd.concat([df1, df2])

    combined_qbed_path = os.path.join(tmpdir, "combined_qbed.qbed")
    combined_df.to_csv(combined_qbed_path, sep="\t", index=False)

    output_path3 = os.path.join(tmpdir, "peaks_combined.csv")

    args = argparse.Namespace(
        experiment_data_path=combined_qbed_path,
        experiment_orig_chr_convention="ucsc",
        promoter_data_path=os.path.join(test_data_directory, "yiming_promoters.bed.gz"),
        promoter_orig_chr_convention="id",
        background_data_path=os.path.join(
            test_data_directory, "adh1_background.qbed.gz"
        ),
        background_orig_chr_convention="ucsc",
        chrmap_data_path=os.path.join(test_data_directory, "chrmap.csv.gz"),
        unified_chr_convention="ucsc",
        output_path=output_path3,
        deduplicate_experiment=False,
        pseudocount=0.2,
        compress_output=False,
    )

    call_peaks_main(args)

    assert os.path.exists(output_path3) is True

    # read in output_path1, output_path2 and output_path3 as dataframes.
    # join on `chr`, `start`, `end`, `score`, `strand`. the output_path3
    # `experiment_hops` and `background_hops` and `experiment_total_hops` should
    # be the same as the sum of the `experiment_hops`, `background_hops` and
    # `experiment_total_hops` columns in output_path1 and output_path2.
    # `background_total_hops` should be the same for all three sources.
    df1 = pd.read_csv(output_path1)
    df2 = pd.read_csv(output_path2)
    df3 = pd.read_csv(output_path3)

    combined_df = df1.merge(
        df2, on=["name", "chr", "start", "end", "strand"], suffixes=("_1", "_2")
    )
    combined_df = combined_df.merge(df3, on=["name", "chr", "start", "end", "strand"])

    assert (
        combined_df["experiment_hops_1"] + combined_df["experiment_hops_2"]
        == combined_df["experiment_hops"]
    ).all()
    assert (
        combined_df["experiment_total_hops_1"] + combined_df["experiment_total_hops_2"]
        == combined_df["experiment_total_hops"]
    ).all()
    assert (combined_df["background_hops_1"] == combined_df["background_hops"]).all()
    assert (combined_df["background_hops_2"] == combined_df["background_hops"]).all()
    assert (
        combined_df["background_total_hops_1"] == combined_df["background_total_hops"]
    ).all()
    assert (
        combined_df["background_total_hops_2"] == combined_df["background_total_hops"]
    ).all()
