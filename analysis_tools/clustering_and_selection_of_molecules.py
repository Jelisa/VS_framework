"""
This script will cluster the molecules according to their CHEMBL/ZINC id and pick the best one according
to a given scoring function. To do this the program will require of:
- activities_simulation_file : This one is obtained from the script compound_sim_activity_parser.py and it contains
                            the relationship between the simulation ID and the ZINC/CHEMBL ids
- sfs_order: one or multiple files containing the SFs values that the user wants to use to pick the
            best molecule for each compound.
The program will generate one output .csv file for each sfs_order file given, this file will contain one sim ID per
line, and this ID will correspond to the best molecule according to the sf used to obtain the order
"""

import csv
import re
import sys
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-activities_simulation_file", required=True)
parser.add_argument("-sfs_file", required=True, nargs='+')
parser.add_argument("-output_prefix", default="selected_according_to")
args = parser.parse_args()

input_files = {}
out_filenames = {}

# Extract the tags for the filenames
for filename in args.sfs_file:
    if re.match(r".*mm[-_]*gbsa[_-]*.*", filename, re.IGNORECASE):
        scoring_function = "mm_gbsa"
    elif re.match(".*glide.*", filename, re.IGNORECASE):
        if re.match(".*glideSP.*", filename, re.IGNORECASE):
            scoring_function = "glideSP"
        elif re.match(".*glideXP.*", filename, re.IGNORECASE):
            scoring_function = "glideXP"
        else:
            scoring_function = "glide"
    elif re.match(".*vina.*", filename, re.IGNORECASE):
        scoring_function = "vina"
    elif re.match(".*Xscore.*", filename, re.IGNORECASE):
        scoring_function = "xscore"
    elif re.match(".*dsx.*", filename, re.IGNORECASE):
        scoring_function = "dsx"
    elif re.match(".*RF-Score.*", filename, re.IGNORECASE):
        scoring_function = "RF-Score"
    elif re.match(".*nn.*", filename, re.IGNORECASE):
        if re.match(".*average.*", filename, re.IGNORECASE):
            scoring_function = "nn_score_average"
        elif re.match(".*best.*", filename, re.IGNORECASE):
            scoring_function = "nn_score_best"
        else:
            scoring_function = raw_input("Which name do you want to use for this file {}?\n".format(filename))
    elif re.match(".*rdock.*", filename, re.IGNORECASE):
        scoring_function = "rdock"
    elif re.match(".*PELE.*", filename, re.IGNORECASE):
        scoring_function = "PELE"
    else:
        scoring_function = raw_input("Which name do you want to use for this file {}?\n".format(filename))
    # input_tags.append(scoring_function)
    input_files[scoring_function] = filename
    out_filenames[scoring_function] = "{0}_{1}.csv".format(args.output_prefix, scoring_function)

print "\n".join(["{} : {}".format(filename, tag) for tag, filename in input_files.iteritems()])
negative_answers = ['no', 'n', 'not', 'nop', 'nope', "N"]
positive_answers = ['yes', 'y', 'yep', 'si', "Y", "S", "s"]
check = None
while check not in negative_answers + positive_answers:
    check = raw_input("Is the above list of filename-tag correct?[Y/n]")
    if check == "":
        check = "y"
if check in negative_answers:
    # TODO: write the list of files and get the tags manually...
    sys.exit()

with open(args.activities_simulation_file) as infile:
    csv_parser = csv.DictReader(infile, delimiter=",")
    sim_id_name_dictio = {line['sim_id']: line['name'] for line in csv_parser}
    # sim_id_name_dictio = {}
    # for line in csv_parser:
    #     try:
    #         sim_id_name_dictio[line['sim_id']].append(line['sim_id'])
    #     except KeyError:
    #         sim_id_name_dictio[line['name']] = [line['sim_id']]

output_dictionary = {}
for tag, filename in input_files.iteritems():
    output_ids = {}
    tmp_dict = {}
    with open(filename) as infile:
        csv_parser = csv.reader(infile, delimiter=",")
        first_line = True
        for line in csv_parser:
            if first_line:
                first_line =False
                continue
            # a = re.search(".*_(\d+)_*", line[0])
            # print line[0]
            # if a is None:
            sim_id = line[0].strip()
            # else:
            #     sim_id = a.group(1)
            # print sim_id
            # print sim_id_name_dictio.keys()
            compound_name = sim_id_name_dictio[sim_id]
            output_ids[sim_id] = line[0].strip()
            try:
                tmp_dict[compound_name]
            except KeyError:
                tmp_dict[compound_name] = {}
            tmp_dict[compound_name][line[0]] = float(line[1].strip())
    with open(out_filenames[tag], 'w') as outfile:
        sf_name = "{0}_value".format(tag)
        csv_writer = csv.DictWriter(outfile, fieldnames=["ID", sf_name], delimiter=",")
        csv_writer.writeheader()
        for compound, values in tmp_dict.iteritems():
            best_molecule_id_value = sorted(values.items(), key=lambda x: x[1])[0]
            csv_writer.writerow({"ID": best_molecule_id_value[0], sf_name: best_molecule_id_value[1]})

