#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path

import pandas as pd
import polars as pl
import matplotlib.pyplot as plt
from genorm import m_measure

bin_script_file = Path(__file__).parents[1] / "bin"
sys.path.append(str(bin_script_file))
from normfinder import NormFinder

"""
The input is a file containing CT results, which should be in CSV format and have the following columns:
    sample,group,gene1,gene2,gene3,gene4,...
    S1,G1,16,57,18,39,...
    S2,G1,16,64,18,89,...
    S3,G1,17,75,19,23,...
    S4,G2,17,61,18,74,...
    S5,G2,17,59,18,33,...
    S6,G2,17,26,18,94,...
"""

def parse_args():
    parser = argparse.ArgumentParser(description="Analyse qPCR results")
    parser.add_argument(
        "--ct", type=Path, dest="ct_file", required=True, help="File containing CT results"
    )
    return parser.parse_args()


def get_genorm_measures(counts: pd.DataFrame) -> pd.Series:
    prepared_counts = counts.drop(columns=["group"]).set_index("sample")
    m_measures_series = m_measure(prepared_counts)
    m_measures_series = m_measures_series.rename("m_measure").sort_values()
    m_measures_series.index.name = "gene_id"
    return m_measures_series


def main():
    args = parse_args()
    ct_file = args.ct_file
    counts = pd.read_csv(ct_file)
    
    # set column types
    counts["sample"] = counts["sample"].astype(str)
    counts["group"] = counts["group"].astype(str)
    count_columns = [col for col in counts.columns if col not in ["sample", "group"]]
    for col in count_columns:
        counts[col] = counts[col].astype(float)

    ##########################################################
    # GENORM
    ##########################################################
    
    m_measures = get_genorm_measures(counts)
    print(f"M values:\n{m_measures}")

    ##########################################################
    # NORMFINDER
    ##########################################################

    samples = counts["sample"].to_list()
    count_df = pl.from_pandas(counts[count_columns]).transpose(
        include_header=True, 
        header_name="gene_id",
        column_names=samples,
    )

    design_df = counts[["sample", "group"]]
    design_df = pl.from_pandas(design_df).select(
        pl.lit("batch").alias("batch"),
        pl.col("group").alias("condition"),
        pl.col("sample")
    )

    # the NormFinder object expects a polars lazy DataFrame and a polars design DataFrame
    nfd = NormFinder(count_df.lazy(), design_df)
    stabilities = nfd.compute_stability_scoring()
    stabilities = stabilities.sort(by="normfinder_stability_value")
    print(f"Stabilities:\n{stabilities}")
    
    ##########################################################
    # EXPORTING
    ##########################################################

    m_measures.to_csv("m_measures.csv")
    stabilities.write_csv("stabilities.csv")
    

if __name__ == "__main__":
    main()