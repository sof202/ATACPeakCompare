import pandas as pd


def subset_bedgraph(bedgraph: pd.DataFrame,
                    chromosome: str,
                    start: int,
                    end: int) -> pd.DataFrame:
    """Subset a bedgraph to a given region.

    Returns:
        A pandas.DataFrame covering the region selected in bedgraph format
    """
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
    return bedgraph


def convert_to_bedbase(bedgraph: pd.DataFrame,
                       chromosome: str,
                       start: int,
                       end: int) -> pd.DataFrame:
    """Convert a bedgraph dataframe into bedbase format

    Returns:
        A pandas.DataFrame covering the region selected in bedbase format
    """
    bases = list(range(start, end+1))
    chromosome_column = [chromosome] * len(bases)
    bedbase = pd.DataFrame(
        zip(chromosome_column, bases),
        columns=["CHR", "BASE"]
    )
    bins = list(bedgraph["START"]) + [bedgraph["END"].iloc[-1]]
    bedbase["SCORE"] = pd.cut(
        bedbase["BASE"],
        bins=bins,
        labels=bedgraph["SCORE"],
        include_lowest=True,
        right=False
    )
    return bedbase


def extract_bedbase_region(bedgraph: pd.DataFrame,
                           chromosome: str,
                           start: int,
                           end: int) -> pd.DataFrame:
    """Extract a region of a bedgraph data frame and convert it into bedbase
    format

    Args
        bedgraph (pandas.DataFrame): Bedgraph to extract from.
        chromosome (str): Chromosome to extract.
        start (int): Start of region.
        end (int): End of region.

    Returns:
        A pandas.DataFrame covering the region selected in bedbase format
    """
    bedgraph = subset_bedgraph(bedgraph, chromosome, start, end)
    bedbase = convert_to_bedbase(bedgraph, chromosome, start, end)
    return bedbase
