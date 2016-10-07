"""
A script to read all the output from pele simulations, should be previously generated, and detect which ones have
failed and generate a report.
"""

from argparse import ArgumentParser
import os
import re

import parameters_help

parser = ArgumentParser()
parser.add_argument("-input", required=True, nargs="+", help=parameters_help.input_desc)
parser.add_argument("-report", default="sims_report.txt", help=parameters_help.report_desc)
args = parser.parse_args()

total_systems = 0
systems_ok = 0
failing_systems = []
ok_systems = []
incomplete_systems = 0
different_errors = {'segmentation_fault': [0, []], 'missing_nma_ter_mark_missing': [0, []], 'mismatch': [0, []],
                    'configuration_file_problem': [0, []], 'missing_conf_file': [0, []], "missing_atom": [0, []],
                    'others': [0, []], 'missing_template': [0, []], "incomplete_simulation": [0, []]}

for filename in args.input:
    if not os.path.isfile(filename):
        print "The filename {} isn't a file, it'll be skipped."
        continue
    total_systems += 1
    keep_lines = False
    error_text = []
    maximum_lines = 0
    others = []
    with open(filename, 'r') as infile:
        output_text = infile.read()
    if re.search("Macro ended successfully", output_text):
        found_error = True
        systems_ok += 1
    elif re.search('/!', output_text):
        failing_systems.append(filename)
        if re.search(r"Segmentation fault", output_text):
            found_error = True
            try:
                different_errors["segmentation_fault"][0] += 1
                different_errors["segmentation_fault"][1].append(filename)
            except KeyError:
                different_errors["segmentation_fault"] = [1, [filename]]
        elif re.search("configuration file doesn't exist", output_text):
            found_error = True
            try:
                different_errors["missing_conf_file"][0] += 1
                different_errors["missing_conf_file"][1].append(filename)
            except KeyError:
                different_errors["missing_conf_file"] = [1, [filename]]
        elif re.search(r"Template atom", output_text):
            error_pattern = re.search(r"Template atom '(\w+)' in template '(\w+)' not found in residue '(...:\d{,3}:\w)'", output_text)
            atom_name = error_pattern.group(1)
            template_name = error_pattern.group(2)
            residue_name = error_pattern.group(3)
            try:
                different_errors["missing_atom"][0] += 1
                different_errors["missing_atom"][1].append([filename, atom_name, residue_name])
            except KeyError:
                different_errors["missing_atom"][0] = [1, [[filename, atom_name, residue_name]]]
        elif re.search(r"Atom name not found in template", output_text):
            error_pattern = re.search(r"Atom name not found in template '(\w+)'; atom is '(\s*\w+\s*)' in residue '(...:\d{,3}:\w)'", output_text)
            template_name = error_pattern.group(1)
            atom_name = error_pattern.group(2)
            residue_name = error_pattern.group(3)
            if template_name == "NMA" and re.search(r"HH3[123]", atom_name):
                try:
                    different_errors["missing_nma_ter_mark_missing"][0] += 1
                    different_errors["missing_nma_ter_mark_missing"][1].append([filename, residue_name])
                except KeyError:
                    different_errors["missing_nma_ter_mark_missing"][0] = [1, [[filename, residue_name]]]
            else:
                try:
                    different_errors["mismatch"][0] += 1
                    different_errors["mismatch"][1].append([filename, atom_name, residue_name])
                except KeyError:
                    different_errors["mismatch"][0] = [1, [[filename, atom_name, residue_name]]]
        elif re.search(r'Error getting structural template', output_text):
            error_pattern = re.search("\*(\w+\s*)\*", output_text)
            missing_template = error_pattern.group(1)
            try:
                different_errors["missing_template"][0] += 1
                different_errors["missing_template"][1].append([filename, missing_template])
            except KeyError:
                different_errors["missing_template"][0] = [1, [[filename, missing_template]]]
        elif re.search(r'validation:FAIL|Errors:', output_text):
            different_errors["configuration_file_problem"][0] += 1
            different_errors["configuration_file_problem"][1].append([filename])
    else:
        incomplete_systems += 1
        try:
            different_errors["incomplete_simulation"][0] += 1
            different_errors["incomplete_simulation"][1].append(filename)
        except KeyError:
            different_errors["incomplete_simulation"][0] = 1
            different_errors["incomplete_simulation"][1] = [filename]

failed_systems = total_systems - systems_ok - incomplete_systems
header = "{} files have been analyzed.\n * {} have finished correctly\n * {} have failed\n * {} are incomplete\n"
report_text = header.format(total_systems, systems_ok, failed_systems, incomplete_systems)
report_text += "The summary of failing systems is:\n"
failing_systems_classification = []
failing_systems_specifics = []
for keyword, values in different_errors.iteritems():
    if keyword == "mismatch":
        failing_systems_classification.append(" * {} systems failed due to {} between the template and the residue, "
                                              "this can be either more atoms than they should, bad naming or a ter"
                                              "mark missing after an standard aa.".format(values[0], keyword))
    else:
        failing_systems_classification.append(" * {} systems failed due to {}".format(values[0], keyword))
    if values[0] != 0:
        failing_systems_specifics.append("The following systems have failed due to {} error:".format(
            keyword.replace("_", " ")))
        for value in sorted(values[1]):
            if keyword in ['segmentation_fault', 'missing_conf_file', "incomplete_simulation"]:
                failing_systems_specifics.append(" - {}".format(value))
            elif keyword == "missing_nma_ter":
                failing_systems_specifics.append(" - {0[0]} : missing TER mark after residue {0[1}".format(value))
            elif keyword == "mismatch":
                failing_systems_specifics.append(" - {0[0]} : the problematic atom is {0[1]} in residue "
                                                      "{0[2]}".format(value))
            elif keyword == "missing_atom":
                failing_systems_specifics.append(" - {0[0} : is missing the atom {0[1]} in the residue "
                                                      "{0[2]}".format(value))
            elif keyword == 'configuration_file_problem':
                failing_systems_specifics.append(" - {} : the configuration file format is incorrect"
                                                      ".".format(value))
            elif keyword == 'missing_template':
                failing_systems_specifics.append("-  {0[0]} : is missing the template {0[1]} .".format(value))
report_text += "\n".join(failing_systems_classification)
report_text += "\nFailing systems specifications:\n"
report_text += "\n".join(failing_systems_specifics)

with open(args.report, 'w') as outfile:
    outfile.write(report_text)
