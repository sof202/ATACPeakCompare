import pandas as pd
from typing import Optional, Union


class IncompatabilityError(Exception):
    """Exception for when positions in two files don't align"""
    pass


class GenomicData:
    """
    Base class for genomic data representations (e.g., BedBase, BedGraph).
    """

    def __init__(self, df: pd.DataFrame):
        """
        Initializes a GenomicData object.

        Args:
            df (pd.DataFrame): The underlying pandas DataFrame.
        """
        self.df = df

    def get(self, column_name: str = None) -> Union[pd.DataFrame, pd.Series]:
        """
        Returns underlying pandas DataFrame
        """
        if column_name is None:
            return self.df
        else:
            return self.df[column_name]

    def write_file(self, file_path: str) -> None:
        """Writes a pandas DataFrame to a file in tab-separated format.

        Args:
            data (pd.DataFrame): The DataFrame to write to the file.
            file_path (str): The path to the output file.
        """
        try:
            with open(file_path, 'w') as file:
                self.df.to_csv(file, sep="\t", header=False, index=False)
        except (FileNotFoundError, IOError):
            print(f"Data could not be written to {file_path}.")
        except IsADirectoryError:
            print(f"{file_path} is a directory.")
        except PermissionError:
            print(f"Permission denied for {file_path}")
        except OSError as e:
            print(f"OS error occurred: {e}")


class BedGraph(GenomicData):
    """
    Represents a BedGraph DataFrame with specific columns: CHR, START, END,
    SCORE.
    """

    def __init__(self, CHR, START, END, SCORE):
        """
        Initializes a BedBase object.

        Args:
            CHR (pd.Series): Series representing the chromosome column.
            START (pd.Series): Series representing start of region.
            END (pd.Series): Series representing end of region.
            SCORE (pd.Series): Series representing the score column. Score can
                be anything, such as number of reads or p-value.
        """
        self.df = pd.DataFrame({
            "CHR": CHR,
            "START": START,
            "END": END,
            "SCORE": SCORE
        })
        self.df = self.df.astype({
            "CHR": 'object',
            "START": 'int64',
            "END": 'int64',
            "SCORE": 'float64'
        })

    def has_same_positions(self, comparison: 'BedGraph') -> bool:
        has_same_chromosome = self.df["CHR"].equals(comparison.df["CHR"])
        has_same_start = self.df["START"].equals(comparison.df["START"])
        has_same_end = self.df["END"].equals(comparison.df["END"])
        if has_same_chromosome and has_same_start and has_same_end:
            return True
        return False

    @classmethod
    def read_from_file(cls, file_path: str) -> Optional["BedGraph"]:
        """Reads a bedgraph file into a pandas DataFrame.

        Args:
            file_path (str): The path to the bedgraph file.

        Returns:
            Optional[pd.DataFrame]: The DataFrame containing the bedgraph data,
            or None if an error occurred.
        """
        try:
            with open(file_path, 'r') as file:
                first_line = file.readline()

            # Some bedgraph files start with a meta data line
            if first_line.startswith("track"):
                bedgraph = pd.read_table(file_path, sep="\t", skiprows=1)
            else:
                bedgraph = pd.read_table(file_path, sep="\t")

            if bedgraph.shape[1] != 4:
                raise ValueError(f"Bedgraph file at {file_path}"
                                 " does not have exactly 4 columns.")
            bedgraph.columns = ["CHR", "START", "END", "SCORE"]
            return cls(
                bedgraph["CHR"],
                bedgraph["START"],
                bedgraph["END"],
                bedgraph["SCORE"]
            )
        except (FileNotFoundError, IOError):
            print(f"{file_path} does not exist or could not be read.")
            return None
        except IsADirectoryError:
            print(f"{file_path} is a directory.")
            return None
        except PermissionError:
            print(f"Permission denied for {file_path}")
            return None
        except OSError as e:
            print(f"OS error occurred: {e}")
            return None


class BedBaseCI(GenomicData):
    """
    Represents a BedBase with confidence intervals with specific columns:
    CHR, BASE, LOWER_SCORE, UPPER_SCORE.
    """

    def __init__(self, CHR, BASE, LOWER_SCORE, UPPER_SCORE):
        """
        Initializes a BedBase object.

        Args:
            CHR (pd.Series): Series representing the chromosome column.
            BASE (pd.Series): Series representing the base column.
            LOWER_SCORE (pd.Series): Series representing the score column for
                lower bound of confidence interval. Score can be anything,
                such as number of reads or p-value.
            UPPER_SCORE (pd.Series): Series representing the score column for
                upper bound of confidence interval.
        """
        self.df = pd.DataFrame({
            "CHR": CHR,
            "BASE": BASE,
            "LOWER_SCORE": LOWER_SCORE,
            "UPPER_SCORE": UPPER_SCORE
        })
        self.df = self.df.astype({
            "CHR": 'object',
            "BASE": 'int64',
            "LOWER_SCORE": 'float64',
            "UPPER_SCORE": 'float64'
        })

    def has_same_positions(self, comparison: 'BedBaseCI') -> bool:
        has_same_chromosome = self.df["CHR"].equals(comparison.df["CHR"])
        has_same_bases = self.df["BASE"].equals(comparison.df["BASE"])
        if has_same_chromosome and has_same_bases:
            return True
        return False


class BedBase(GenomicData):
    """
    Represents a BedBase DataFrame with specific columns: CHR, BASE, SCORE.
    """

    def __init__(self, CHR, BASE, SCORE):
        """
        Initializes a BedBase object.

        Args:
            CHR (pd.Series): Series representing the chromosome column.
            BASE (pd.Series): Series representing the base column.
            SCORE (pd.Series): Series representing the score column. Score can
                be anything, such as number of reads or p-value.
        """
        self.df = pd.DataFrame({
            "CHR": CHR,
            "BASE": BASE,
            "SCORE": SCORE
        })
        self.df = self.df.astype({
            "CHR": 'object',
            "BASE": 'int64',
            "SCORE": 'float64'
        })

    def has_same_positions(self, comparison: 'BedBase') -> bool:
        has_same_chromosome = self.df["CHR"].equals(comparison.df["CHR"])
        has_same_bases = self.df["BASE"].equals(comparison.df["BASE"])
        if has_same_chromosome and has_same_bases:
            return True
        return False


class Bed(GenomicData):
    """
    Represents a Bed DataFrame with specific columns: CHR, START, END
    """

    def __init__(self, CHR, START, END):
        """
        Initializes a Bed object.

        Args:
            CHR (pd.Series): Series representing the chromosome column.
            START (pd.Series): Series representing start of region.
            END (pd.Series): Series representing end of region.
        """
        self.df = pd.DataFrame({
            "CHR": CHR,
            "START": START,
            "END": END
        })
        self.df = self.df.astype({
            "CHR": 'object',
            "START": 'int64',
            "END": 'int64'
        })

    def has_same_positions(self, comparison: 'Bed') -> bool:
        has_same_chromosome = self.df["CHR"].equals(comparison.df["CHR"])
        has_same_start = self.df["START"].equals(comparison.df["START"])
        has_same_end = self.df["END"].equals(comparison.df["END"])
        if has_same_chromosome and has_same_start and has_same_end:
            return True
        return False

    @classmethod
    def read_from_file(cls, file_path: str) -> Optional["Bed"]:
        """Reads a BED3+7 file from MACS into a pandas DataFrame.

        Args:
            file_path (str): The path to the bedbase file.

        Returns:
            Optional[pd.DataFrame]: The DataFrame containing the bed data,
            or None if an error occurred.
        """
        try:
            bed = pd.read_table(file_path, sep="\t", skiprows=1)
            if bed.shape[1] < 3:
                raise ValueError(f"Bed file at {file_path}"
                                 " does not have enough columns.")
            # Remaining columns that might exist are useless
            bed = bed.iloc[:, 0:3]
            bed.columns = ["CHR", "START", "END"]
            return cls(
                bed["CHR"],
                bed["START"],
                bed["END"]
            )
        except (FileNotFoundError, IOError):
            print(f"{file_path} does not exist or could not be read.")
            return None
        except IsADirectoryError:
            print(f"{file_path} is a directory.")
            return None
        except PermissionError:
            print(f"Permission denied for {file_path}")
            return None
        except OSError as e:
            print(f"OS error occurred: {e}")
            return None
