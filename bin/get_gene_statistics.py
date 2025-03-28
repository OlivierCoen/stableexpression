#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import polars as pl
from pathlib import Path
from dataclasses import dataclass, field
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# nb of top stable genes to select and to display at the end
DEFAULT_NB_TOP_STABLE_GENES = 1000
# we want to select samples that show a particularly low nb of genes
MIN_RATIO_GENE_COUNT_TO_MEAN = 0.75  # experimentally chosen
WEIGHT_RATIO_NB_NULLS = 1


# outfile names
TOP_STABLE_GENE_SUMMARY_OUTFILENAME = "top_stable_genes_summary.csv"
ALL_GENES_RESULT_OUTFILENAME = "stats_all_genes.csv"
ALL_COUNTS_FILTERED_PARQUET_OUTFILENAME = "all_counts_filtered.parquet"
TOP_STABLE_GENES_COUNTS_OUTFILENAME = "top_stable_genes_transposed_counts_filtered.csv"

# column names
RANK_COLNAME = "Rank"
ORIGINAL_GENE_ID_COLNAME = "original_gene_id"
ORIGINAL_GENE_IDS_COLNAME = "original_gene_ids"
ENSEMBL_GENE_ID_COLNAME = "ensembl_gene_id"
M_MEASURE_COLNAME = "m_measure"
GENE_NAME_COLNAME = "name"
GENE_DESCRIPTION_COLNAME = "description"
VARIATION_COEFFICIENT_COLNAME = "variation_coefficient"
STANDARD_DEVIATION_COLNAME = "standard_deviation"
MEAN_COLNAME = "mean"
EXPRESSION_LEVEL_QUANTILE_INTERVAL_COLNAME = "expression_level_quantile_interval"
EXPRESSION_LEVEL_STATUS_COLNAME = "expression_level_status"
GENE_COUNT_COLNAME = "count"
SAMPLE_COLNAME = "sample"
NB_NULLS_COLNAME = "total_nb_nulls"
NB_NULLS_VALID_SAMPLES_COLNAME = "nb_nulls_valid_samples"
NB_ZEROS_COLNAME = "nb_zeros"
STABILITY_SCORE_COLNAME = "stability_score"
KS_TEST_COLNAME = "kolmogorov_smirnov_to_uniform_dist_pvalue"

STATISTICS_COLS = [
    RANK_COLNAME,
    ENSEMBL_GENE_ID_COLNAME,
    STABILITY_SCORE_COLNAME,
    M_MEASURE_COLNAME,
    STANDARD_DEVIATION_COLNAME,
    VARIATION_COEFFICIENT_COLNAME,
    MEAN_COLNAME,
    EXPRESSION_LEVEL_STATUS_COLNAME,
    NB_NULLS_COLNAME,
    NB_NULLS_VALID_SAMPLES_COLNAME,
    GENE_NAME_COLNAME,
    GENE_DESCRIPTION_COLNAME,
    ORIGINAL_GENE_IDS_COLNAME,
]

ALL_GENES_STATS_COLS = [
    ENSEMBL_GENE_ID_COLNAME,
    STABILITY_SCORE_COLNAME,
    M_MEASURE_COLNAME,
    MEAN_COLNAME,
    STANDARD_DEVIATION_COLNAME,
    VARIATION_COEFFICIENT_COLNAME,
]

# quantile intervals
NB_QUANTILES = 100

NB_TOP_GENES_TO_SHOW_IN_LOG_COUNTS = 100


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Get statistics from count data for each gene"
    )
    parser.add_argument(
        "--counts", type=Path, dest="count_file", required=True, help="Count file"
    )
    parser.add_argument(
        "--metadata",
        type=str,
        dest="metadata_files",
        required=True,
        help="Metadata file",
    )
    parser.add_argument(
        "--mappings", type=str, dest="mapping_files", required=True, help="Mapping file"
    )
    parser.add_argument(
        "--m-measures",
        type=str,
        dest="m_measure_file",
        required=True,
        help="M-measure file",
    )
    parser.add_argument(
        "--nb-top-stable-genes",
        type=int,
        dest="nb_top_stable_genes",
        required=True,
        help="Number of top stable genes to show",
    )
    parser.add_argument(
        "--ks-stats",
        type=Path,
        dest="ks_stats_file",
        required=True,
        help="KS stats file",
    )
    parser.add_argument(
        "--ks-pvalue-threshold",
        type=str,
        dest="ks_pvalue_threshold",
        required=True,
        help="KS p-value threshold",
    )
    return parser.parse_args()


def is_valid_lf(lf: pl.LazyFrame, file: Path) -> bool:
    """Check if a LazyFrame is valid.

    A LazyFrame is considered valid if it contains at least one row.
    """
    try:
        return not lf.limit(1).collect().is_empty()
    except FileNotFoundError:
        # strangely enough we get this error for some files existing but empty
        logger.error(f"Could not find file {str(file)}")
        return False
    except pl.exceptions.NoDataError as err:
        logger.error(f"File {str(file)} is empty: {err}")
        return False


def get_valid_lazy_lfs(files: list[Path]) -> list[pl.LazyFrame]:
    """Get a list of valid LazyFrames from a list of files.

    A LazyFrame is considered valid if it contains at least one row.
    """
    lf_dict = {file: pl.scan_csv(file) for file in files}
    return [lf for file, lf in lf_dict.items() if is_valid_lf(lf, file)]


def cast_cols_to_string(lf: pl.LazyFrame) -> pl.LazyFrame:
    return lf.select(
        [pl.col(column).cast(pl.String) for column in lf.collect_schema().names()]
    )


def concat_cast_to_string_and_drop_duplicates(files: list[Path]) -> pl.LazyFrame:
    """Concatenate LazyFrames, cast all columns to String, and drop duplicates.

    The first step is to concatenate the LazyFrames. Then, the dataframe is cast
    to String to ensure that all columns have the same data type. Finally, duplicate
    rows are dropped.
    """
    lfs = get_valid_lazy_lfs(files)
    lfs = [cast_cols_to_string(lf) for lf in lfs]
    concat_lf = pl.concat(lfs)
    # dropping duplicates
    # casting all columns to String
    return concat_lf.unique()


def get_count_columns(lf: pl.LazyFrame) -> list[str]:
    """Get all column names except the ENSEMBL_GENE_ID_COLNAME column.

    The ENSEMBL_GENE_ID_COLNAME column contains only gene IDs.
    """
    return lf.select(pl.exclude(ENSEMBL_GENE_ID_COLNAME)).collect_schema().names()


def cast_count_columns_to_float32(lf: pl.LazyFrame) -> pl.LazyFrame:
    return lf.select(
        [pl.col(ENSEMBL_GENE_ID_COLNAME)]
        + [pl.col(column).cast(pl.Float32) for column in get_count_columns(lf)]
    )


def get_counts(
    file: Path, ks_stats_file: Path, ks_pvalue_threshold: str
) -> pl.LazyFrame:
    # sorting dataframe (necessary to get consistent output)
    count_lf = pl.scan_parquet(file).sort(ENSEMBL_GENE_ID_COLNAME, descending=False)
    ks_stats_df = pl.read_csv(
        ks_stats_file, has_header=False, new_columns=[SAMPLE_COLNAME, KS_TEST_COLNAME]
    )

    # parsing threshold
    try:
        ks_pvalue_threshold = float(ks_pvalue_threshold)
    except ValueError:
        raise ValueError(
            f"KS p-value threshold {ks_pvalue_threshold} could not be cast to float"
        )

    # logging number of samples excluded from analysis
    not_valid_samples = ks_stats_df.filter(
        ks_stats_df[KS_TEST_COLNAME] <= ks_pvalue_threshold
    )[SAMPLE_COLNAME].to_list()
    logger.warning(
        f"Excluded {len(not_valid_samples)} samples showing a KS p-value below {ks_pvalue_threshold}"
    )

    # getting samples for which the Kolmogorov-Smirnov test pvalue is above the threshold
    valid_samples = ks_stats_df.filter(
        ks_stats_df[KS_TEST_COLNAME] > ks_pvalue_threshold
    )[SAMPLE_COLNAME].to_list()
    # filtering the count dataframe to keep only the valid samples
    return count_lf.select([ENSEMBL_GENE_ID_COLNAME] + valid_samples)


def get_metadata(metadata_files: list[Path]) -> pl.LazyFrame:
    """Retrieve and concatenate metadata from a list of metadata files."""
    return concat_cast_to_string_and_drop_duplicates(metadata_files)


def get_mappings(mapping_files: list[Path]) -> pl.LazyFrame:
    concat_lf = concat_cast_to_string_and_drop_duplicates(mapping_files)
    # group by new gene IDs and gets the lis
    """Group by new gene IDs, get the list of distinct original gene IDs and convert to a string representation."""
    # t of distinct original gene IDs for each group
    # convert the list column to a string representation
    # separate the original gene IDs with a semicolon
    return concat_lf.group_by(ENSEMBL_GENE_ID_COLNAME).agg(
        pl.col(ORIGINAL_GENE_ID_COLNAME)
        .unique()
        .sort()
        .str.join(";")
        .alias(ORIGINAL_GENE_IDS_COLNAME)
    )


def add_m_measures(stat_lf: pl.LazyFrame, m_measure_file: Path) -> pl.LazyFrame:
    if m_measure_file != "none" and Path(m_measure_file).exists():
        stat_lf = stat_lf.join(
            pl.scan_csv(m_measure_file), on=ENSEMBL_GENE_ID_COLNAME, how="left"
        ).sort(M_MEASURE_COLNAME, descending=False)
    return stat_lf


def merge_data(
    stat_lf: pl.LazyFrame, metadata_lf: pl.LazyFrame, mapping_lf: pl.LazyFrame
) -> pl.LazyFrame:
    """Merge the statistics dataframe with the metadata dataframe and the mapping dataframe."""
    # we need to ensure that the index of stat_lf are strings
    return stat_lf.join(metadata_lf, on=ENSEMBL_GENE_ID_COLNAME, how="left").join(
        mapping_lf, on=ENSEMBL_GENE_ID_COLNAME, how="left"
    )


def sort_dataframe(lf: pl.LazyFrame) -> pl.LazyFrame:
    if M_MEASURE_COLNAME in lf.collect_schema().names():
        lf = lf.sort(M_MEASURE_COLNAME, descending=False, nulls_last=True)
    else:
        lf = lf.sort(STABILITY_SCORE_COLNAME, descending=False, nulls_last=True)
    return (
        lf.with_row_index(name="index")
        .with_columns((pl.col("index") + 1).alias("Rank"))
        .drop("index")
    )


def get_status(quantile_interval: int) -> str:
    """Return the expression level status of the gene given its quantile interval."""
    if NB_QUANTILES - 5 <= quantile_interval:
        return "Very high expression"
    elif NB_QUANTILES - 10 <= quantile_interval < NB_QUANTILES - 5:
        return "High expression"
    elif 4 < quantile_interval <= 9:
        return "Low expression"
    elif quantile_interval <= 4:
        return "Very low expression"
    else:
        return "Medium range"


def get_top_stable_gene_summary(
    stat_lf: pl.LazyFrame, nb_top_stable_genes: int
) -> pl.LazyFrame:
    """
    Extract the most stable genes from the statistics dataframe.
    """
    logger.info("Getting most stable genes per quantile interval")
    mapping_dict = {
        quantile_interval: get_status(quantile_interval)
        for quantile_interval in range(NB_QUANTILES)
    }
    lf = stat_lf.head(nb_top_stable_genes).with_columns(
        pl.col(EXPRESSION_LEVEL_QUANTILE_INTERVAL_COLNAME)
        .replace_strict(mapping_dict)
        .alias(EXPRESSION_LEVEL_STATUS_COLNAME)
    )
    return lf.select(
        [column for column in STATISTICS_COLS if column in lf.collect_schema().names()]
    )


def format_all_genes_statistics(stat_lf: pl.LazyFrame) -> pl.LazyFrame:
    """
    Format the dataframe containing statistics for all genes by selecting the right columns
    and sorting the dataframe by gene ID.
    """
    return stat_lf.select(
        [
            column
            for column in ALL_GENES_STATS_COLS
            if column in stat_lf.collect_schema().names()
        ]
    ).sort(STABILITY_SCORE_COLNAME, descending=False)


def get_top_stable_genes_counts(
    log_count_lf: pl.LazyFrame, top_stable_genes_summary_lf: pl.LazyFrame
) -> pl.DataFrame:
    # getting list of top stable genes in the order
    sorted_stable_genes = (
        top_stable_genes_summary_lf.head(NB_TOP_GENES_TO_SHOW_IN_LOG_COUNTS)
        .select(ENSEMBL_GENE_ID_COLNAME)
        .collect()
        .to_series()
        .to_list()
    )
    mapping_dict = {item: index for index, item in enumerate(sorted_stable_genes)}

    # extracting log counts of top stable genes
    sorted_transposed_counts_df = (
        log_count_lf.filter(pl.col(ENSEMBL_GENE_ID_COLNAME).is_in(sorted_stable_genes))
        .with_columns(
            pl.col(ENSEMBL_GENE_ID_COLNAME)
            .replace_strict(mapping_dict)
            .alias("sort_order")
        )
        .sort("sort_order", descending=False)
        .drop(["sort_order", ENSEMBL_GENE_ID_COLNAME])
    ).collect()

    return sorted_transposed_counts_df.transpose(column_names=sorted_stable_genes)


def export_data(
    top_stable_genes_summary_lf: pl.LazyFrame,
    formated_stat_lf: pl.LazyFrame,
    all_counts_lf: pl.LazyFrame,
    top_stable_genes_counts_df: pl.DataFrame,
):
    """Export gene expression data to CSV files."""
    logger.info(
        f"Exporting statistics of the top stable genes to: {TOP_STABLE_GENE_SUMMARY_OUTFILENAME}"
    )
    top_stable_genes_summary_lf.collect().write_csv(TOP_STABLE_GENE_SUMMARY_OUTFILENAME)

    logger.info(
        f"Exporting statistics for all genes to: {ALL_GENES_RESULT_OUTFILENAME}"
    )
    formated_stat_lf.collect().write_csv(ALL_GENES_RESULT_OUTFILENAME)

    logger.info(f"Exporting all counts to: {ALL_COUNTS_FILTERED_PARQUET_OUTFILENAME}")
    all_counts_lf.collect().write_parquet(ALL_COUNTS_FILTERED_PARQUET_OUTFILENAME)

    logger.info(
        f"Exporting counts of the top stable genes to: {TOP_STABLE_GENES_COUNTS_OUTFILENAME}"
    )
    top_stable_genes_counts_df.write_csv(TOP_STABLE_GENES_COUNTS_OUTFILENAME)

    logger.info("Done")


#####################################################
#####################################################
# CLASSES
#####################################################
#####################################################


@dataclass
class StabilityScorer:
    count_lf: pl.LazyFrame

    gene_count_per_sample_df: pl.DataFrame = field(init=False)
    stat_lf: pl.LazyFrame = field(init=False)
    count_columns: list[str] = field(init=False)
    samples_with_low_gene_count: list[str] = field(init=False)

    def __post_init__(self):
        self.count_columns = get_count_columns(self.count_lf)
        self.gene_count_per_sample_df = self.get_gene_counts_per_sample()
        self.samples_with_low_gene_count = self.get_samples_with_low_gene_count()

    def get_valid_counts(self) -> pl.LazyFrame:
        return self.count_lf.select(pl.exclude(ENSEMBL_GENE_ID_COLNAME))

    def get_gene_counts_per_sample(self) -> pl.DataFrame:
        """
        Get the number of non-null values per sample.
        :return:
        A polars dataframe containing 2 columns:
            - sample: name of the sample
            - nb_not_nulls: number of non-null values
        """
        return (
            self.count_lf.select(pl.exclude(ENSEMBL_GENE_ID_COLNAME))
            .count()
            .collect()
            .transpose(
                include_header=True, header_name="sample", column_names=["count"]
            )
        )

    def get_samples_with_low_gene_count(self) -> list[str]:
        mean_gene_count = self.gene_count_per_sample_df[GENE_COUNT_COLNAME].mean()
        return (
            self.gene_count_per_sample_df.filter(
                (pl.col(GENE_COUNT_COLNAME) / mean_gene_count)
                < MIN_RATIO_GENE_COUNT_TO_MEAN
            )
            .select(SAMPLE_COLNAME)
            .to_series()
            .to_list()
        )

    def get_main_statistics(self) -> pl.LazyFrame:
        """
        Compute count descriptive statistics for each gene in the count dataframe.
        """
        logger.info("Getting descriptive statistics")
        # computing main stats
        augmented_count_lf = self.count_lf.with_columns(
            mean=pl.concat_list(self.count_columns).list.drop_nulls().list.mean(),
            std=pl.concat_list(self.count_columns).list.drop_nulls().list.std(),
        )
        return augmented_count_lf.select(
            pl.col(ENSEMBL_GENE_ID_COLNAME),
            pl.col("mean").alias(MEAN_COLNAME),
            pl.col("std").alias(STANDARD_DEVIATION_COLNAME),
            (pl.col("std") / pl.col("mean")).alias(VARIATION_COEFFICIENT_COLNAME),
        )

    def compute_nb_null_values(self):
        # the samples showing a low gene count will not be taken into account for the zero count penalty
        cols_to_exclude = [ENSEMBL_GENE_ID_COLNAME] + self.samples_with_low_gene_count
        total_nb_nulls = (
            self.count_lf.select(pl.exclude(ENSEMBL_GENE_ID_COLNAME).is_null())
            .collect()
            .sum_horizontal()
        )
        nb_nulls_valid_samples = (
            self.count_lf.select(pl.exclude(cols_to_exclude).is_null())
            .collect()
            .sum_horizontal()
        )
        self.stat_lf = self.stat_lf.with_columns(
            total_nb_nulls.alias(NB_NULLS_COLNAME),
            nb_nulls_valid_samples.alias(NB_NULLS_VALID_SAMPLES_COLNAME),
        )

    def get_quantile_intervals(self):
        """
        Compute the quantile intervals for the mean expression levels of each gene in the dataframe.

        The function assigns to each gene a quantile interval of its mean cpm compared to all genes.
        """
        logger.info("Getting cpm quantiles")
        self.stat_lf = self.stat_lf.with_columns(
            (pl.col(MEAN_COLNAME).rank() / pl.col(MEAN_COLNAME).count() * NB_QUANTILES)
            .floor()
            .cast(pl.Int8)
            # we want the only value = NB_QUANTILES to be NB_QUANTILES - 1
            # because the last quantile interval is [NB_QUANTILES - 1, NB_QUANTILES]
            .replace({NB_QUANTILES: NB_QUANTILES - 1})
            .alias(EXPRESSION_LEVEL_QUANTILE_INTERVAL_COLNAME)
        )

    def compute_stability_score(self):
        logger.info("Computing stability score")
        nb_valid_samples = self.gene_count_per_sample_df.select(pl.len()).item() - len(
            self.samples_with_low_gene_count
        )
        ratio_nb_nulls = (
            self.stat_lf.select(
                pl.col(NB_NULLS_VALID_SAMPLES_COLNAME) / nb_valid_samples
            )
            .collect()
            .to_series()
        )
        expr = (
            pl.col(STANDARD_DEVIATION_COLNAME) + ratio_nb_nulls * WEIGHT_RATIO_NB_NULLS
        )
        self.stat_lf = self.stat_lf.with_columns(expr.alias(STABILITY_SCORE_COLNAME))

    def compute_statistics_and_score(self) -> pl.LazyFrame:
        logger.info("Computing statistics and stability score")
        # getting expression statistics
        self.stat_lf = self.get_main_statistics()
        # adding column for nb of null values for each gene
        self.compute_nb_null_values()
        # computing stability score
        self.compute_stability_score()
        # getting quantile intervals
        self.get_quantile_intervals()

        return self.stat_lf


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()
    metadata_files = [Path(file) for file in args.metadata_files.split(" ")]
    mapping_files = [Path(file) for file in args.mapping_files.split(" ")]

    # putting all counts into a single dataframe
    count_lf = get_counts(args.count_file, args.ks_stats_file, args.ks_pvalue_threshold)

    # getting metadata and mappings
    metadata_lf = get_metadata(metadata_files)
    mapping_lf = get_mappings(mapping_files)

    # computing statistics (mean, standard deviation, coefficient of variation, quantiles)
    stability_scorer = StabilityScorer(count_lf)
    stat_lf = stability_scorer.compute_statistics_and_score()

    # adding other statistics (example: m-measure)
    stat_lf = add_m_measures(stat_lf, args.m_measure_file)

    # add gene name, description and original gene IDs
    stat_lf = merge_data(stat_lf, metadata_lf, mapping_lf)

    # sort genes according to the metrics present in the dataframe
    stat_lf = sort_dataframe(stat_lf)

    # getting the most stable genes
    # we don't want to exceed 1000 (for multiqc)
    nb_top_stable_genes = min(args.nb_top_stable_genes, DEFAULT_NB_TOP_STABLE_GENES)
    top_stable_genes_summary_lf = get_top_stable_gene_summary(
        stat_lf, nb_top_stable_genes
    )

    formated_stat_lf = format_all_genes_statistics(stat_lf)

    # reducing dataframe size (it is only used for plotting by MultiQC)
    count_lf = cast_count_columns_to_float32(count_lf)

    top_stable_genes_counts_df = get_top_stable_genes_counts(
        count_lf, top_stable_genes_summary_lf
    )

    # exporting computed data
    export_data(
        top_stable_genes_summary_lf,
        formated_stat_lf,
        count_lf,
        top_stable_genes_counts_df,
    )


if __name__ == "__main__":
    main()
