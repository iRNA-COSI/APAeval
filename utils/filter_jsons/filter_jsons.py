from argparse import ArgumentParser, RawTextHelpFormatter
import json
import logging

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
        help="list of OEB .json (assessments, validation or consolidation) filepaths.", required=True, nargs="+"
    )
    parser.add_argument(
        "-o",
        "--out_prefix",
        help="prefix for output files",
        required=True
    )
    parser.add_argument(
        "-m",
        "--b_metrics",
        help="list of metrics to be removed",
        nargs='?',
        default = ""
    )
    parser.add_argument(
        "-c",
        "--b_challenges",
        help="list of challenges to be removed",
        nargs='?',
        default= ""
    )
    return parser


def main():
    # input parameters
    file_list = options.file_list
    out_prefix = options.out_prefix
    bad_metrics = options.b_metrics.strip().split(" ")
    bad_challenges = options.b_challenges.strip().split(" ")

    LOGGER.info(f"badlisted metric elements: {bad_metrics}")
    LOGGER.info(f"badlisted challenge elements: {bad_challenges}")

    # for each json file in input list
    for f in file_list:
        LOGGER.debug(f"Filename: {f}")
        new_jobjects = []
        with open(f, mode='r', encoding="utf-8") as infile:
            try:
                jobjects = json.load(infile)
            except Exception as exc:
                LOGGER.error(f"Can't read json file, please check format.")
                raise exc

        # for each object in json file
        for j in jobjects:
            LOGGER.debug(f"Object: {j}")
            # object is from the manifest
            if not 'type' in j:
                bad_c = []
                if bad_challenges[0]:
                    bad_c = [True for c in bad_challenges if c in j['id']]
                if not any(bad_c):
                    new_jobjects.append(j)

            # object is participant object
            elif j['type'] == 'participant':
                # From a list of challenges only keep those that don't contain a bad element
                challenges = []
                bad_c = []
                for c in j['challenge_id']:
                    if bad_challenges[0]:
                        bad_c = [True for b in bad_challenges if b in c]
                    if not any(bad_c):
                        challenges.append(c)
                j['challenge_id'] = challenges
                new_jobjects.append(j)

            # object is assessment or aggregation object
            elif (j['type'] == 'assessment') or (j['type'] == 'aggregation'):
                # Check for bad metric elements
                bad_m = []
                if bad_metrics[0]:
                    bad_m = [True for b in bad_metrics if b in j['_id']]
                # Check for bad challenge elements
                bad_c = []
                if bad_challenges[0]:
                    bad_c = [True for c in bad_challenges if c in j['_id']]
                # Decide whether to keep object
                if (not any(bad_m)) and (not any(bad_c)):
                    new_jobjects.append(j)

            
        # And write to file
        with open(out_prefix + f, mode='w', encoding='utf-8') as outfile:
            jdata = json.dumps(new_jobjects, sort_keys=True, indent=4, separators=(',', ': '))
            outfile.write(jdata)

    return





if __name__ == '__main__':
    # parse the command-line arguments
    options = parse_arguments().parse_args()
    
    # execute the body of the script
    LOGGER.info("Starting script")

    main()

    LOGGER.info("Finished script successfully.")