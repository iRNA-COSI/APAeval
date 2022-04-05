"""
Convert .csv to .tsv
"""
import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter

def parse_arguments():
    parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        "--csv",
        dest="in_csv",
        required=True,
        help="Path to the input CSV file.",
    )
    parser.add_argument(
        "--tsv",
        dest="out_tsv",
        required=True,
        help="Path for the output TSV file.",
    )
    return parser


def main():
    try:
        csv = pd.read_csv(options.in_csv)
        csv.to_csv(options.out_tsv, sep='\t', index=False)
        
    except Exception as e:
        raise e


if __name__ == "__main__":

    try:
        # parse the command-line arguments
        options = parse_arguments().parse_args()
        main()
        
    except Exception as e:
        raise e
