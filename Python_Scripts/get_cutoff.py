import pandas as pd
import argparse


def get_cutoff(filepath: str, average_peak_length: int) -> float:
    cutoffs = pd.read_csv(filepath, sep="\t", header=0)
    cutoff_to_use = cutoffs.loc[cutoffs['avelpeak']
                                >= average_peak_length].iloc[0].score
    return cutoff_to_use


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("cutoffs_file_path", type=str)
    parser.add_argument("average_peak_length", type=int)
    args = parser.parse_args()

    cutoff = get_cutoff(args.cutoffs_file_path, args.average_peak_length)
    print(cutoff)
