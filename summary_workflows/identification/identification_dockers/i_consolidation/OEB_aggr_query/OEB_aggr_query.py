#!/usr/bin/env python3

from __future__ import division
import requests
import json
import os
import logging
import sys
from argparse import ArgumentParser

DEFAULT_OEB_API = "https://dev-openebench.bsc.es/api/scientific/graphql"
DEFAULT_bench_event_id = "OEBE0070000001"


def main(args):

    # input parameters
    data_dir = args.benchmark_data
    challenge_name = args.challenge

    if not args.event_id:
        event_id = DEFAULT_bench_event_id
    else:
        event_id = args.event_id

    logging.info(f"Passed event_id: {args.event_id}")
    logging.info(f"Set event_id: {event_id}")
    logging.info("Querying OEB database...")

    # Load data from OEB DB
    response = query_OEB_DB(event_id, challenge_name)

    # For debugging:
    with open(os.path.join(data_dir, "response.json"), mode='w', encoding="utf-8") as f:
            json.dump(response, f, sort_keys=True, indent=4, separators=(',', ': '))

    # from the loaded data, create an aggregation file and write it to data_dir
    getOEBAggregations(response, data_dir)

## get existing aggregation datasets from OEB DB
def query_OEB_DB(bench_event_id, challenge_name):
    logging.info(f"Querying OEB for aggregation objects of event {bench_event_id} and challenge {challenge_name}...")

    json_query = {
        'query': 
            """query AggregationQuery($bench_event_id: String, $challenge_name: String) {
            getChallenges(challengeFilters: {benchmarking_event_id: $bench_event_id, acronym: $challenge_name}){
                acronym
                datasets(datasetFilters: {type: "aggregation"}) {
                    orig_id
                    datalink {
                        inline_data
                    }
                }
            }
            }""",
        'variables': {
            'bench_event_id': bench_event_id,
            'challenge_name': challenge_name
        }
    }
    try:
        url = DEFAULT_OEB_API
        # get challenges and input datasets for provided benchmarking event
        r = requests.post(url=url, json=json_query, headers={'Content-Type': 'application/json'})
        response = r.json()
        data = response.get("data")
        if data is None:
            logging.error("For {} got response error from graphql query: {}".format(bench_event_id, r.text))
            sys.exit(6)
        if len(data["getChallenges"]) == 0:
            logging.warning("No challenges associated to benchmarking event " + bench_event_id +
                          " in OEB. Please contact OpenEBench support for information about how to open a new challenge")
            return
        else:
            return data.get('getChallenges')
    except Exception as e:

        logging.exception(e)
        
# function to populate assess_dir with existing aggregations
def getOEBAggregations(response, output_dir):

    # The response contains data for all plots of one challenge
    challenge = response[0]
    challenge_name = challenge["acronym"]
    challenge_aggregation_objects = []

    # For each plot
    for item in challenge["datasets"]:
        aggregation_data = json.loads(item["datalink"]["inline_data"])

        # replace key "tool_id" for "participant_id" (for the visualization)
        for i in aggregation_data["challenge_participants"]:
            i["participant_id"] = i.pop("tool_id")
        
        # aggregation object compatible with our workflow
        aggregation_object = {
            "_id": item["orig_id"],
            "challenge_ids": [
                 challenge_name
            ],
            "datalink": {
                "inline_data": aggregation_data
                },
            "type": "aggregation"
        }

        challenge_aggregation_objects.append(aggregation_object)


    with open(os.path.join(output_dir, challenge_name+ ".json"), mode='w', encoding="utf-8") as f:
        json.dump(challenge_aggregation_objects, f, sort_keys=True, indent=4, separators=(',', ': '))



if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument("-b", "--benchmark_data", help="dir where the data for the benchmark are stored", required=True)
    parser.add_argument("-e", "--event_id", help="OEB benchmarking event ID for which data shall be downloaded", required=False)
    parser.add_argument("-c", "--challenge", help="challenge name for which data shall be downloaded", required=True)
    args = parser.parse_args()

    # set up logging during the execution
    logging.basicConfig(level=logging.DEBUG, 
                format='%(asctime)s %(levelname)s:%(message)s', datefmt='%Y-%m-%d %H:%M:%S ')
    main(args)