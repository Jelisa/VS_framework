"""
A script to parse the report files and compute the mean average of a value from the PELE report file.
"""

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import numpy as np
import os
import csv
import sys
from re import search

def parse_pele_report_file(report_filename):
    """
    This function reads a report file from PELE and parses it.
    :param report_filename: a string containing the name of the file to read
    :return: report_array: a numpy array containing all the values from the report and
    report_header: a list containing the header of the file
    """
    report_values = []
    first_line = True
    report_header = ""
    with open(report_filename, 'r') as report_file:
        first_line = True
        for line in report_file.readlines():
            line = line.strip()
            if not line:
                continue
            split_line = line.split("    ")
            if first_line:
                first_line = False
                report_header = split_line
            else:
                values = []
                for v in split_line:
                    try:
                        values.append(int(v))
                    except ValueError:
                        try:
                            values.append(float(v))
                        except ValueError:
                            values.append(v)
                report_values.append(values)
    report_array = np.asarray(report_values)

    return report_array, report_header

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("-input_files", nargs='+', required=True)
parser.add_argument("-output_file", default="mean_PELE_BE.csv")
parser.add_argument("-element_to_extract", default="Binding Energy")
parser.add_argument("-initial_model", default=1, type=int)
args = parser.parse_args()

sim_value_dictio = {}

for filename in args.input_files:
    if not os.path.isfile(filename):
        print "The file {0} doesn't exist".format(filename)
        continue
    values, header = parse_pele_report_file(filename)
    if header == "":
        print "The file {0} is empty".format(filename)
        continue
    try:
        b_e_index = header.index(args.element_to_extract)
    except ValueError:
        print "The file {0} doesn't have the Binding energy".format(filename)
        continue
    pele_mean_energy = values[args.initial_model:, b_e_index].mean()
    name = filename.split(os.sep)[-1]
    pattern = search(r"(.*)_report", name)
    if pattern:
        name = pattern.group(1)
    sim_value_dictio[name] = pele_mean_energy

with open(args.output_file, 'w') as outfile:
    field_names = ["ID", "PELE_mean_BE"]
    csv_writter = csv.DictWriter(outfile, fieldnames=field_names)
    csv_writter.writeheader()
    for k, v in sorted(sim_value_dictio.iteritems()):
        csv_writter.writerow({"ID": k, "PELE_mean_BE": v})
