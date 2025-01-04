#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import polars as pl
from pathlib import Path
import logging
from functools import reduce

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

ORIGINAL_GENE_ID_COLNAME = "original_gene_id"
ENSEMBL_GENE_ID_COLNAME = "ensembl_gene_id"
GENE_VARIATION_OUTFILENAME = "gene_variations.csv"
COUNT_OUTFILENAME = "all_normalised_counts.csv"

GENE_VARIATION_COLNAME = "variation_coefficient"
STANDARD_DEVIATION_COLNAME = "standard_deviation"
METHOD_NAMES = {"std": STANDARD_DEVIATION_COLNAME, "cv": GENE_VARIATION_COLNAME}

# Get variation coefficient from count data for each gene
# The variation coefficient is the ratio of the standard deviation to the mean
# We want genes that are neither expressed too much nor too little
# Metadata (name and description) are used to annotate the results
# Likewise, mappings (original gene ids) are used to better associate gene ids with their original gene ids
# Usage:
# python get_variation_coefficients.py --counts <count_file> --metadata <metadata_file> --mapping <mapping_file>


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Get variation from count data for each gene"
    )
    parser.add_argument(
        "--counts", type=str, dest="count_files", required=True, help="Count file"
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
        "--method",
        type=str,
        dest="gene_variation_method",
        required=True,
        choices=["std", "cv"],
        help="Method for computation of the gene variation",
    )
    return parser.parse_args()


def is_valid_df(df: pl.LazyFrame, file: Path) -> bool:
    try:
        return not df.limit(1).collect().is_empty()
    except (
        FileNotFoundError
    ):  # strangely enough we get this error for files existing but empty
        logger.error(f"Could not find file {str(file)}")
        return False
    except pl.exceptions.NoDataError as err:
        logger.error(f"File {str(file)} is empty: {err}")
        return False


def get_valid_lazy_dfs(files: list[Path]) -> list[pl.LazyFrame]:
    df_dict = {file: pl.scan_csv(file) for file in files}
    return [df for file, df in df_dict.items() if is_valid_df(df, file)]


def join_dfs(df1: pl.LazyFrame, df2: pl.LazyFrame):
    return df1.join(df2, on=ENSEMBL_GENE_ID_COLNAME, how="full", coalesce=True)


def concat_cast_to_string_and_drop_duplicates(files: list[Path]) -> pl.LazyFrame:
    dfs = get_valid_lazy_dfs(files)
    concat_df = pl.concat(dfs)
    # dropping duplicates
    # casting all columns to String
    return concat_df.unique().with_columns(
        [
            pl.col(column).cast(pl.String)
            for column in concat_df.collect_schema().names()
        ]
    )


def get_count_columns(df: pl.LazyFrame) -> list[str]:
    return [
        col for col in df.collect_schema().names() if col != ENSEMBL_GENE_ID_COLNAME
    ]


def get_counts(files: list[Path]) -> pl.LazyFrame:
    # lazy loading
    dfs = get_valid_lazy_dfs(files)
    # joining all count files
    merged_df = reduce(join_dfs, dfs)
    # casting count columns to Float64
    # casting gene id column to String
    count_columns = get_count_columns(merged_df)
    return merged_df.with_columns(
        [pl.col(column).cast(pl.Float64) for column in count_columns]
    ).with_columns([pl.col(ENSEMBL_GENE_ID_COLNAME).cast(pl.String)])


def get_metadata(metadata_files: list[Path]) -> pl.LazyFrame:
    return concat_cast_to_string_and_drop_duplicates(metadata_files)


def get_mappings(mapping_files: list[Path]) -> pl.LazyFrame:
    concat_df = concat_cast_to_string_and_drop_duplicates(mapping_files)
    # group by new gene IDs and gets the list of distinct original gene IDs for each group
    # convert the list column to a string representation
    # separate the original gene IDs with a semicolon
    return concat_df.group_by(ENSEMBL_GENE_ID_COLNAME).agg(
        pl.col(ORIGINAL_GENE_ID_COLNAME)
        .unique()
        .sort()
        .str.join(";")
        .alias("original_gene_ids")
    )


def get_standard_deviation(count_df: pl.LazyFrame) -> pl.LazyFrame:
    """
    Compute variation for each gene in the count dataframe.
    """
    logger.info("Getting standard deviations")
    count_columns = get_count_columns(count_df)
    return count_df.with_columns(
        std=pl.concat_list(count_columns).list.std(),
    ).select(
        pl.col(ENSEMBL_GENE_ID_COLNAME),
        pl.col("std").alias(STANDARD_DEVIATION_COLNAME),
    )


def get_coefficient_of_variation(count_df: pl.LazyFrame) -> pl.LazyFrame:
    """
    Compute variation for each gene in the count dataframe.
    """
    logger.info("Getting coefficients of variation")
    count_columns = get_count_columns(count_df)
    return count_df.with_columns(
        mean=pl.concat_list(count_columns).list.mean(),
        std=pl.concat_list(count_columns).list.std(),
    ).select(
        pl.col(ENSEMBL_GENE_ID_COLNAME),
        (pl.col("std") / pl.col("mean")).alias(GENE_VARIATION_COLNAME),
    )


def get_average_log2(count_df: pl.LazyFrame) -> pl.LazyFrame:
    # the dataframe has already been filtered to exclude rows where mean is 0
    # adds 1 to avoid log(0) and to stabilize variance
    logger.info("Getting average log2 counts")
    count_columns = get_count_columns(count_df)
    transformed_cols = [pl.col(ENSEMBL_GENE_ID_COLNAME)] + [
        (pl.col(column) + 1).log(2).alias(column) for column in count_columns
    ]
    return (
        count_df.select(transformed_cols)
        .with_columns(average_log2_count=pl.concat_list(count_columns).list.mean())
        .select(pl.col(ENSEMBL_GENE_ID_COLNAME), pl.col("average_log2_count"))
    )


def aggregate_expression_values(
    count_df: pl.LazyFrame, gene_variation_method: str
) -> pl.LazyFrame:
    # handling NA values (genes that are not found in all datasets)
    # replace NaN values with 0
    count_df = count_df.fill_null(0)

    # getting column containing the gene variations
    if gene_variation_method == GENE_VARIATION_COLNAME:
        gene_variation = get_coefficient_of_variation(count_df)
    else:
        gene_variation = get_standard_deviation(count_df)

    # get the average log cpm value for each gene
    # to get an idea of the overall expression level of each gene
    average_log_cpm = get_average_log2(count_df)

    return gene_variation.join(average_log_cpm, on=ENSEMBL_GENE_ID_COLNAME, how="left")


def merge_data(
    cv_df: pl.LazyFrame, metadata_df: pl.LazyFrame, mapping_df: pl.LazyFrame
) -> pl.LazyFrame:
    # we need to ensure that the index of cv_df are strings
    return (
        cv_df.join(metadata_df, on=ENSEMBL_GENE_ID_COLNAME, how="left")
        .join(mapping_df, on=ENSEMBL_GENE_ID_COLNAME, how="left")
        .unique()  # just in case
    )


def export_data(cv_df: pl.LazyFrame, count_df: pl.LazyFrame):
    logger.info(f"Exporting normalised counts to: {COUNT_OUTFILENAME}")
    count_df.collect().write_csv(COUNT_OUTFILENAME)

    logger.info(f"Exporting gene variations to: {GENE_VARIATION_OUTFILENAME}")
    cv_df.collect().write_csv(GENE_VARIATION_OUTFILENAME)


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()
    count_files = [Path(file) for file in args.count_files.split(" ")]
    metadata_files = [Path(file) for file in args.metadata_files.split(" ")]
    mapping_files = [Path(file) for file in args.mapping_files.split(" ")]
    gene_variation_method = METHOD_NAMES[args.gene_variation_method]

    count_df = get_counts(count_files)
    cv_df = aggregate_expression_values(count_df, gene_variation_method)

    metadata_df = get_metadata(metadata_files)
    mapping_df = get_mappings(mapping_files)

    cv_df = merge_data(cv_df, metadata_df, mapping_df)
    cv_df = cv_df.sort(gene_variation_method, descending=False)
    export_data(cv_df, count_df)


if __name__ == "__main__":
    main()