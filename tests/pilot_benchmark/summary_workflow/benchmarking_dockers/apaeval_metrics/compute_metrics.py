from __future__ import division
import io
import os
import pandas
import math
import json
from argparse import ArgumentParser
from JSON_templates import JSON_templates


def main(args):

    # input parameters
    input_participant = args.participant_data
    gold_standards_dir = args.metrics_ref
    cancer_types = args.cancer_types
    participant = args.participant_name
    community = args.community_name
    out_path = args.output

    # Assuring the output path does exist
    if not os.path.exists(os.path.dirname(out_path)):
        try:
            os.makedirs(os.path.dirname(out_path))
            with open(out_path, mode="a"):
                pass
        except OSError as exc:
            print("OS error: {0}".format(exc) + "\nCould not create output path: " + out_path)

    compute_metrics(input_participant,  gold_standards_dir, cancer_types, participant, community, out_path)



def compute_metrics(input_participant,  gold_standards_dir, cancer_types, participant, community, out_path):

    # get participant dataset
    participant_data = pandas.read_csv(input_participant, sep='\t',
                                       comment="#", header=0)
    
    # filter data by q-value
    if participant == "MutSig2CV":

        filtered_data = participant_data.loc[participant_data['qvalue'] <= 0.1]

        predicted_genes = filtered_data.iloc[:, 0].values

    elif participant == "ActiveDriver":

        filtered_data = participant_data.loc[participant_data['qvalue'] <= 0.0001]

        predicted_genes = filtered_data.iloc[:, 0].values

    elif participant == "MuSiC":

        filtered_data = participant_data.loc[participant_data['pvalue'] <= math.exp(-8)]
        filtered_data = filtered_data[filtered_data['info'] == "FILTER=PASS"]

        predicted_genes = filtered_data.iloc[:, 0].values

    else:

        filtered_data = participant_data.loc[participant_data['qvalue'] <= 0.05]

        predicted_genes = filtered_data.iloc[:, 0].values
    
    # define array that will hold the full set of assessment datasets
    ALL_ASSESSMENTS = []

    for cancer in cancer_types:
        # get metrics dataset
        metrics_data = pandas.read_csv(os.path.join(gold_standards_dir, cancer + ".txt"),
                                       comment="#", header=None)
        gold_standard = metrics_data.iloc[:, 0].values


        # TRUE POSITIVE RATE
        overlapping_genes = set(predicted_genes).intersection(gold_standard)
        TPR = len(overlapping_genes) / len(gold_standard)

        # ACCURACY/ PRECISION
        if len(predicted_genes) == 0:
            acc = 0
        else:
            acc = len(overlapping_genes) / len(predicted_genes)

        assessment_data = {'toolname': participant, 'x': TPR, 'y': acc, 'e': 0, 'cancer_type': cancer}

        # get json assessment file for both metrics
        data_id_1 = community + ":" + cancer + "_TPR_" + participant + "_A"
        std_error= 0
        assessment_TPR = JSON_templates.write_assessment_dataset(data_id_1, community, cancer, participant, "TPR", TPR, std_error)

        data_id_2 = community + ":" + cancer + "_precision_" + participant + "_A"
        assessment_precision = JSON_templates.write_assessment_dataset(data_id_2, community, cancer, participant, "precision", acc, std_error)

        # push the two assessment datasets to the main dataset array
        ALL_ASSESSMENTS.extend([assessment_TPR, assessment_precision])

    # once all assessments have been added, print to json file
    with io.open(out_path,
                 mode='w', encoding="utf-8") as f:
        jdata = json.dumps(ALL_ASSESSMENTS, sort_keys=True, indent=4, separators=(',', ': '))
        f.write(unicode(jdata,"utf-8"))


if __name__ == '__main__':
    
    parser = ArgumentParser()
    parser.add_argument("-i", "--participant_data", help="list of cancer genes prediction", required=True)
    parser.add_argument("-c", "--cancer_types", nargs='+', help="list of types of cancer selected by the user, separated by spaces", required=True)
    parser.add_argument("-m", "--metrics_ref", help="dir that contains metrics reference datasets for all cancer types", required=True)
    parser.add_argument("-p", "--participant_name", help="name of the tool used for prediction", required=True)
    parser.add_argument("-com", "--community_name", help="name/id of benchmarking community", required=True)
    parser.add_argument("-o", "--output", help="output path where assessment JSON files will be written", required=True)
    
    args = parser.parse_args()

    
    main(args)



