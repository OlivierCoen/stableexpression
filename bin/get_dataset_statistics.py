#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
from pathlib import Path
from scipy import stats
import pandas as pd
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

QUANT_NORM_SUFFIX = ".quant_norm.parquet"
DATASET_STATISTICS_SUFFIX = ".dataset_stats.csv"

ENSEMBL_GENE_ID_COLNAME = "ensembl_gene_id"


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute general statistics from count data for each sample"
    )
    parser.add_argument(
        "--counts", type=Path, dest="count_file", required=True, help="Count file"
    )
    return parser.parse_args()


def compute_kolmogorov_smirnov_test_to_uniform_distribution(count_df: pd.DataFrame):
    """Compute Kolmogorov-Smirnov test to uniform distribution."""
    ks_tests = pd.Series(index=count_df.columns)
    for col in count_df.columns:
        ks = stats.ks_1samp(count_df[col], stats.uniform.cdf, nan_policy="omit")
        ks_tests[col] = ks.pvalue
    return ks_tests


def compute_dataset_statistics(count_df: pd.DataFrame):
    dataset_stats_df = count_df.describe()
    dataset_stats_df.loc["skewness"] = count_df.skew()
    # for each sample, test distance to uniform distribution
    ks_tests = compute_kolmogorov_smirnov_test_to_uniform_distribution(count_df)
    dataset_stats_df.loc["ks_uniform"] = ks_tests
    return dataset_stats_df


def export_count_data(dataset_stats_df: pd.DataFrame, count_file: Path):
    """Export dataset statistics to CSV files."""
    outfilename = count_file.name.replace(QUANT_NORM_SUFFIX, DATASET_STATISTICS_SUFFIX)
    logger.info(f"Exporting dataset statistics counts to: {outfilename}")
    dataset_stats_df.to_csv(outfilename)


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()
    count_file = args.count_file

    logger.info(f"Computing dataset statistics for {count_file.name}")
    count_df = pd.read_parquet(count_file)
    count_df.set_index(ENSEMBL_GENE_ID_COLNAME, inplace=True)

    dataset_stats_df = compute_dataset_statistics(count_df)

    export_count_data(dataset_stats_df, count_file)


if __name__ == "__main__":
    main()
