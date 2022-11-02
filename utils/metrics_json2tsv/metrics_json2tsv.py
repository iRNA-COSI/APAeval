from argparse import ArgumentParser, RawTextHelpFormatter
import json
import logging

import pandas as pd

logging.basicConfig(
            level=logging.DEBUG, 
            format='%(asctime)s %(levelname)s:%(message)s', datefmt='%Y-%m-%d %H:%M:%S ')
LOGGER = logging.getLogger(__name__)


def parse_arguments():
    '''
    Parser of command-line arguments.
    '''
    parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        "-f", 
        "--file-list", 
        help="list of OEB assessment_datasets.json filepaths.", required=True, nargs="+"
    )
    parser.add_argument(
        "-o",
        "--output",
        help="path to output samples table",
        required=True
    )
    return parser


def main():
    # input parameters
    file_list = options.file_list
    outfile = options.output
    metrics_df = pd.DataFrame()

    # for each json file in input list
    for f in file_list:
        LOGGER.debug(f"Filename: {f}")

        with open(f, mode='r', encoding="utf-8") as f:
            try:
                jmetrics = json.load(f)
            except Exception as exc:
                LOGGER.error(f"Can't read json file, please check format.")
                raise exc

        # Make sure we only consider assessment objects
        jmetrics = [j for j in jmetrics if j["type"] == "assessment"]
        
        # Convert json to df
        tmetrics = pd.json_normalize(jmetrics)
        
        # Add to the DataFrame for all tools
        metrics_df = pd.concat([metrics_df, tmetrics], ignore_index=True)

    LOGGER.debug(f"metrics_df: {metrics_df}")

    # Split metric_id to get window size and return type
    # First get metric name in new col
    metrics_df[['metric','temp']] = metrics_df['metrics.metric_id'].str.split(':', 1, expand=True)

    # Now swap return type and window (neccessary to get NA if return type is missing)
    metrics_df['temp'] = metrics_df['temp'].apply(lambda x: str_swap(x,':') if x is not None else "nt:all" )

    # type and window in new col
    metrics_df[['window_size', 'site_set']] = metrics_df['temp'].str.split(':', 1, expand=True)
    # drop obsolete cols
    metrics_df.drop(['temp','metrics.metric_id'],axis=1, inplace=True)
    # Get rid of "nt" in window size
    metrics_df['window_size'] = metrics_df['window_size'].str.strip("nt")


    # And write metrics table to file
    with open(outfile, mode='w', encoding='utf-8') as o:
        metrics_df.to_csv(o,sep="\t",header=True)

    return


def str_swap(x, sep):
    l = x.split(sep)
    x = sep.join(l[::-1])
    return x


if __name__ == '__main__':
    # parse the command-line arguments
    options = parse_arguments().parse_args()
    
    # execute the body of the script
    LOGGER.info("Starting script")

    main()

    LOGGER.info("Finished script successfully.")
