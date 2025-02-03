#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
from pathlib import Path
import pandas as pd
from sklearn.preprocessing import QuantileTransformer
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

QUANT_NORM_SUFFIX = ".quant_norm.parquet"

ENSEMBL_GENE_ID_COLNAME = "ensembl_gene_id"
N_QUANTILES = 1000
OUTPUT_DISTRIBUTION = "uniform"


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute mean from count data for each gene"
    )
    parser.add_argument(
        "--counts", type=Path, dest="count_file", required=True, help="Count file"
    )
    return parser.parse_args()


def quantile_normalize(data: pd.DataFrame):
    """
    Quantile normalize a data matrix based on a target distribution.
    """
    transformer = QuantileTransformer(
        n_quantiles=N_QUANTILES, output_distribution=OUTPUT_DISTRIBUTION
    )

    normalised_data = pd.DataFrame(index=data.index, columns=data.columns)
    for col in data.columns:
        normalised_data[col] = transformer.fit_transform(data[col].to_frame())

    return normalised_data


def export_count_data(quantile_normalized_counts: pd.DataFrame, count_file: Path):
    """Export gene expression data to CSV files."""
    outfilename = count_file.name.replace(".csv", QUANT_NORM_SUFFIX)
    logger.info(f"Exporting quantile normalised counts to: {outfilename}")
    quantile_normalized_counts.reset_index().to_parquet(outfilename)


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()
    count_file = args.count_file

    logger.info(f"Quantile normalising {count_file.name}")
    count_df = pd.read_csv(count_file, index_col=0)
    count_df.index.name = ENSEMBL_GENE_ID_COLNAME

    quantile_normalized_counts = quantile_normalize(count_df)

    export_count_data(quantile_normalized_counts, count_file)


if __name__ == "__main__":
    main()
