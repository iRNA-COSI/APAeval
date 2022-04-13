
import os
import json
import logging
import datetime
from argparse import ArgumentParser, RawTextHelpFormatter


def parse_arguments():
    '''
    Parser of command-line arguments.
    '''
    parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        "-a", 
        "--assessment_data", 
        help="path to participant assessment.json", required=True
    )
    parser.add_argument(
        "-b", 
        "--benchmark_dir",
        help="dir where the data for the benchmark are stored", 
        required=True
    )
    parser.add_argument(
        "-o",
        "--output",
        help="output directory where the manifest and output JSON files will be written",
        required=True
    )
    parser.add_argument(
        "-d",
        "--event_date",
        help="passes in the event_date defined in nextflow.config",
        required=True
    )
    parser.add_argument(
        "--offline",
        help="offline mode; if set to 1, existing benchmarking datasets will be read from local dir instead of OEB DB",
        default=0,
        required=False,
        type=int
    )
    return parser


def main():

    # input parameters
    benchmark_dir = options.benchmark_dir
    assessment_data = options.assessment_data
    output_dir = options.output
    event_date = options.event_date
    offline = options.offline

    # Assuring the output directory for participant does exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    ########################################################
    # 1. Get assesment (=metrics) for current participant
    ########################################################

    # Load assessment file into dict of challenges with dicts of metric IDs
    community_id, participant_id, challenges = get_metrics_per_challenge(assessment_data)
    # Get a list of challenge ids
    challenge_ids = list(challenges.keys())

    logging.debug(f"Participant {participant_id}, challenges {challenges}, challenge_ids {challenge_ids}")

    ########################################################
    # 2. Handle aggregation file(s)
    # Loop over challenge IDs, for each challenge:
    # a) download results of previous benchmark from OEB
    # b) or get it from local data directory
    # c) or start fresh with the provided template
    # d) Get list of participants
    # e) add participant's metrics to aggregation
    # f) write aggregation file
    # g) store info for summary file (manifest) 
    ########################################################
    
    # Store info for summary file
    manifest = []

    for challenge_id in challenge_ids:
        # 2.a) TO DO ###########################
        # If not offline:
        # Load aggregation file from OEB, convert, and write to benchmark_dir
        # Does this have to be performed on every run of the summary workflow? Only makes sense if we're synchronized with the OEB DB?
        # see manage_assessment_data.py for functions to query OEB, and to convert the result to aggregation file format
        ####################################


        # 2.b) If aggregation file in benchmark_dir
        # Load aggregation file from there
        aggregation_old = os.path.join(benchmark_dir,challenge_id + ".json")
        # Or if necessary template from here
        aggregation_template= os.path.join(os.path.dirname(os.path.realpath(__file__)), "aggregation_template_Q.json")
        
        try:
            # First try to open an existing aggregation file
            try:
                with open(aggregation_old, mode='r', encoding="utf-8") as f:
                    aggregation = json.load(f)
            # 2.c) If there's none, open the aggregation template instead        
            except FileNotFoundError as e:
                logging.warning(f"Couldn't find an existing aggregation file for challenge {challenge_id}. Creating a new file from {aggregation_template}!")

                with open(aggregation_template, mode='r', encoding="utf-8") as t:
                    aggregation = json.load(t)
                # We still need to fill the template with ID and challenge ID, data will be appended later
                for item in aggregation:
                    item["_id"] = f"{community_id}:{event_date}_{challenge_id}_aggregation"
                    item["challenge_ids"] = [challenge_id]
               
        # if something else than the file missing went wrong
        except Exception:
            raise
            
        logging.debug(f"aggregation on load: {aggregation}")

        # 2.d) Get a list of participants from aggregation (needed for manifest)
        participants = []
        
        # Already recorded participants
        # Dirty: we're only looking at the first aggregation object, assuming the participants are the same for all objects (=plots)
        for item in aggregation[0]["datalink"]["inline_data"]["challenge_participants"]:
            participants.append(item["participant_id"])
        # Current participant
        participants.append(participant_id)

        logging.debug(f"participants: {participants}")

        # 2.e) Add the current participant's metrics to the aggregation for the current challenge_id
        new_aggregation = add_to_aggregation(aggregation, participant_id, challenges[challenge_id])

        logging.debug(f"aggregation after update: {aggregation}")

        # 2.f) Write aggregation file to benchmark_dir (and TO DO local results dir)
        aggregation_file = os.path.join(benchmark_dir, challenge_id + ".json")
        with open(aggregation_file, mode='w', encoding="utf-8") as f:
            json.dump(new_aggregation, f, sort_keys=True, indent=4, separators=(',', ': '))
        
        # 2.g) store manifest object for current challenge
        mani_obj = {
            "id" : challenge_id,
            "participants": participants,
            'timestamp': datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0).isoformat()
        }

        manifest.append(mani_obj)


    # After we have updated all aggregation files for all challenges, save the summary manifest
    with open(os.path.join(output_dir, "Manifest.json"), mode='w', encoding="utf-8") as f:
        json.dump(manifest, f, sort_keys=True, indent=4, separators=(',', ': '))


# Function definitions
def get_metrics_per_challenge(assessment_data):
    '''
    From the assessment file collect all challenges and all metrics per challenge
    Returns:
    participant ID
    dict of challenges, per challenge dict of metric IDs

    challenges = {
        "Qa": {
            "correlation": 1,
            "percent_matched": 1
        },
        "Qb": {
            "correlation": 0.9,
            "percent_matched": 0.5
        }
    }

    QUESTION: should we rather load the total content of the json files, instead of picking only what we need. Because what we need might change later on and then we'd have to make bigger changes to the code

    '''
    # read assessment file
    with open(assessment_data, mode='r', encoding="utf-8") as f:
        assessments = json.load(f)
    
    # initialize
    challenges = {}
    participant_id = assessments[0]["participant_id"]
    community_id = assessments[0]["community_id"]

    # go through all assessment objects
    for item in assessments:

        # make sure we're dealing with assessment objects
        assert_object_type(item, "assessment")
        # get metric
        metric_id = item["metrics"]["metric_id"]
        metric_val = item["metrics"]["value"]
        # add metric to dict of its challenge
        challenges.setdefault(item['challenge_id'], {})[metric_id] = metric_val

        # make sure all objects belong to same participant
        if item["participant_id"] != participant_id:
            raise ValueError(f"Something went wrong, not all metrics in the assessment file belong to the same participant.")
    
    return community_id, participant_id, challenges


def assert_object_type(json_obj, curr_type):
    '''
    Check OEB json object type
    Input:
    json object
    desired object type (e.g. "assessment", "aggregation", "participant")
    Returns:
    None; raises TypeError if object is not of given type
    '''
    type_field = json_obj["type"]
    if type_field != curr_type:
        raise TypeError(f"json object is of type {type_field}, should be {curr_type}")



def add_to_aggregation(aggregation, participant_id, challenge):
    '''
    Add the metrics for the current challenge to the challenge's aggregation file. Aggregation file can have more than one aggregation object, one per plot type.
    '''
    for item in aggregation:
        assert_object_type(item, "aggregation")

        # get the current visualization object
        plot = item["datalink"]["inline_data"]["visualization"]

        # Depending on the type of plot we'll need to create different participant objects
        participant = {}
        participant["participant_id"] = participant_id
        if plot["type"] == "2D-plot":
            try:
                participant["metric_x"] = challenge[plot["x_axis"]]
                participant["metric_y"] = challenge[plot["y_axis"]]
            except Exception as e:
                logging.exception(str(e))
                raise e
        
        elif plot["type"] == "bar-plot":
            try:
                participant["metric_value"] = challenge[plot["metric"]]
            except KeyError as e:
                logging.exception(f"The assessment file does not contain data for metric {plot['metric']}.")
                raise e

        # Append current participant to the list of participant objects
        item["datalink"]["inline_data"]["challenge_participants"].append(participant)
    
    return aggregation




if __name__ == '__main__':
    try:
        # parse the command-line arguments
        options = parse_arguments().parse_args()

        # set up logging during the execution
        logging.basicConfig(level=logging.DEBUG, 
                    format='%(asctime)s %(levelname)s:%(message)s', datefmt='%Y-%m-%d %H:%M:%S ')
        # execute the body of the script
        logging.info("Starting script")
        main()
        logging.info("Finished script successfully.")
    
    except Exception as e:
        logging.exception(str(e))
        raise e