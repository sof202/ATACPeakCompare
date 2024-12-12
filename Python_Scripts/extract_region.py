import numpy as np
import pandas as pd
from IO import BedGraph, BedBase


def subset_bedgraph(bedgraph: BedGraph,
                    chromosome: str,
                    start: int,
                    end: int) -> BedGraph:
    """Subset a BedGraph to a given region.

    Returns:
        A BedGraph covering the region selected
    """
    bedgraph = bedgraph.get()
    # Ranges in bedgraph do not overlap and so these should always yield a
    # single index each.
    start_index = bedgraph.loc[
        (bedgraph["CHR"] == chromosome) &
        (bedgraph["START"] <= start) &
        (bedgraph["END"] >= start)
    ].index[0]
    end_index = bedgraph.loc[
        (bedgraph["CHR"] == chromosome) &
        (bedgraph["START"] <= end) &
        (bedgraph["END"] >= end)
    ].index[0]
    bedgraph = bedgraph.iloc[start_index:end_index+1]
    bedgraph = BedGraph(
        CHR=bedgraph["CHR"],
        START=bedgraph["START"],
        END=bedgraph["END"],
        SCORE=bedgraph["SCORE"]
    )
    return bedgraph


def convert_to_bedbase(bedgraph: BedGraph,
                       chromosome: str,
                       start: int,
                       end: int) -> BedBase:
    bases = pd.Series(list(range(start, end+1)))
    chromosome_column = pd.Series([chromosome] * len(bases))
    bins = list(bedgraph.get("START")) + [bedgraph.get("END").iloc[-1]]
    score = pd.cut(
        bases,
        bins=bins,
        labels=bedgraph.get("SCORE"),
        include_lowest=True,
        ordered=False,
        right=False
    )
    score = np.nan_to_num(score, nan=np.nanmin(score))
    bedbase = BedBase(
        CHR=chromosome_column,
        BASE=bases,
        SCORE=score
    )
    return bedbase


def extract_bedbase_region(bedgraph: BedGraph,
                           chromosome: str,
                           start: int,
                           end: int) -> BedBase:
    """Extract a region of a bedgraph data frame and convert it into bedbase
    format

    Args
        bedgraph (BedGraph): Bedgraph to extract from.
        chromosome (str): Chromosome to extract.
        start (int): Start of region.
        end (int): End of region.

    Returns:
        A BedBase covering the region selected
    """
    bedgraph = subset_bedgraph(bedgraph, chromosome, start, end)
    bedbase = convert_to_bedbase(bedgraph, chromosome, start, end)
    return bedbase
