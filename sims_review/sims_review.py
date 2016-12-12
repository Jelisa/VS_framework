"""
A script to read all the output from pele simulations, should be previously generated, and detect which ones have
failed and generate a report.
"""

from argparse import ArgumentParser
import os
import re

import parameters_help

parser = ArgumentParser()
parser.add_argument("-input", required=True, nargs="+", help=parameters_help.sims_review_input_desc)
parser.add_argument("-report", default="sims_report.txt", help=parameters_help.sims_review_report_desc)
args = parser.parse_args()

total_systems = 0
systems_ok = 0
failing_systems = []
ok_systems = []
incomplete_systems = 0
different_errors = {'segmentation_fault': [0, []], 'missing_nma_ter_mark_missing': [0, []], 'mismatch': [0, []],
                    'configuration_file_problem': [0, []], 'missing_conf_file': [0, []], "missing_atom": [0, []],
                    'others': [0, []], 'missing_template': [0, []], "incomplete_simulation": [0, []],
                    "cannot_write_report": [0, []], "missing_obc_parameters": [0, []],
                    "rotamers_template_missmatch": [0, []]}

for filename in args.input:
    if not os.path.isfile(filename):
        print "The filename {0} isn't a file, it'll be skipped.".format(filename)
        continue
    total_systems += 1
    keep_lines = False
    error_text = []
    maximum_lines = 0
    others = []
    with open(filename, 'r') as infile:
        output_text = infile.read()
    if re.search("Macro ended successfully", output_text, re.IGNORECASE):
        found_error = True
        systems_ok += 1
    elif re.search('/!', output_text):
        failing_systems.append(filename)
        if re.search(r"Segmentation fault", output_text, re.IGNORECASE):
            found_error = True
            different_errors["segmentation_fault"][0] += 1
            different_errors["segmentation_fault"][1].append(filename)
        elif re.search("configuration file doesn't exist", output_text, re.IGNORECASE):
            found_error = True
            different_errors["missing_conf_file"][0] += 1
            different_errors["missing_conf_file"][1].append(filename)
        elif re.search(r"Template atom", output_text):
            error_pattern = re.search(
                r"Template atom '(\w+)' in template '(\w+)' not found in residue '(...:\d{,3}:\w)'", output_text,
                re.IGNORECASE)
            if error_pattern:
                atom_name = error_pattern.group(1)
                template_name = error_pattern.group(2)
                residue_name = error_pattern.group(3)
                different_errors["missing_atom"][0] += 1
                different_errors["missing_atom"][1].append([filename, atom_name, residue_name])
            error_pattern = re.search(r"Rotamer library refers to atoms '(\s*\w+\s*)' and '(\s*\w+\s*)' in link"
                                      r" '(\w+)(\d+)\.(\w)', but the atom to apply the rotamer could not be found",
                                      output_text)
            if error_pattern:
                atom_name1 = error_pattern.group(1)
                atom_name2 = error_pattern.group(2)
                residue_name = error_pattern.group(3)
                residue_number = error_pattern.group(4)
                residue_chain = error_pattern.group(5)
                different_errors["rotamers_template_missmatch"][0] += 1
                different_errors["rotamers_template_missmatch"][1].append([filename, atom_name1, atom_name2,
                                                                           residue_name, residue_number, residue_chain])
        elif re.search(r"Atom name not found in template", output_text):
            error_pattern = re.search(
                r"Atom name not found in template '(\w+)'; atom is '(\s*\w+\s*)' in residue '(...:\d{,3}:\w)'",
                output_text, re.IGNORECASE)
            template_name = error_pattern.group(1)
            atom_name = error_pattern.group(2)
            residue_name = error_pattern.group(3)
            if template_name == "NMA" and re.search(r"HH3[123]", atom_name):
                different_errors["missing_nma_ter_mark_missing"][0] += 1
                different_errors["missing_nma_ter_mark_missing"][1].append([filename, residue_name])
            else:
                different_errors["mismatch"][0] += 1
                different_errors["mismatch"][1].append([filename, atom_name, residue_name])
        elif re.search(r"ImplementationObcAlphaSasaUpdater::checkAtomSolventId:", output_text):
            error_pattern = re.search(r"Error atom parameters NOT found in template for (\w{,3})([a-z])_(\w+)",
                                      output_text, re.IGNORECASE)
            residue_name = error_pattern.group(1)
            chain_name = error_pattern.group(2)
            atom_name = error_pattern.group(3)
            different_errors["missing_obc_parameters"][0] += 1
            different_errors["missing_obc_parameters"][1].append([filename, atom_name, residue_name, chain_name])
        elif re.search(r'Error getting structural template', output_text):
            error_pattern = re.search("\*(\w+\s*)\*", output_text)
            missing_template = error_pattern.group(1)
            different_errors["missing_template"][0] += 1
            different_errors["missing_template"][1].append([filename, missing_template])
        elif re.search(r'validation:FAIL|Errors:', output_text):
            different_errors["configuration_file_problem"][0] += 1
            different_errors["configuration_file_problem"][1].append([filename])
        elif re.search("Could't create PeleReport file.*(\(\w+.\w+)\)", output_text):
            different_errors["cannot_write_report"][0] += 1
            different_errors["cannot_write_report"][1].append([filename])
        else:
            different_errors["others"][0] += 1
            different_errors["others"][1].append(filename)
    elif re.search(r"Running macro\.\.\.", output_text):
        incomplete_systems += 1
        try:
            different_errors["incomplete_simulation"][0] += 1
            different_errors["incomplete_simulation"][1].append(filename)
        except KeyError:
            different_errors["incomplete_simulation"][0] = 1
            different_errors["incomplete_simulation"][1] = [filename]
    else:
        different_errors["others"][0] += 1
        different_errors["others"][1].append(filename)
failed_systems = total_systems - systems_ok - incomplete_systems
header = "{0} files have been analyzed.\n * {1} have finished correctly\n * {2} have failed\n * {3} are incomplete\n"
report_text = header.format(total_systems, systems_ok, failed_systems, incomplete_systems)
report_text += "The summary of failing systems is:\n"
failing_systems_classification = []
failing_systems_specifics = []
for keyword, values in different_errors.iteritems():
    if keyword == "mismatch":
        failing_systems_classification.append(" * {0} systems failed due to {1} between the template and the residue, "
                                              "this can be either more atoms than they should, bad naming or a ter"
                                              "mark missing after an standard aa.".format(values[0], keyword))
    elif keyword == "cannot_write_report":
        failing_systems_classification.append(" * {0} systems failed due to {1} , meaning that the program is not able"
                                              " to write the report file, probably the output folder is missing, "
                                              "or you don't have the permissions to write on "
                                              "the folder.".format(values[0], keyword))
    elif keyword == "others":
        failing_systems_classification.append(" * {0} systems failed due to {1} errors, probably segmentation "
                                              "fault".format(values[0], keyword))
    else:
        failing_systems_classification.append(" * {0} systems failed due to {1}".format(values[0], keyword))
    if values[0] != 0:
        failing_systems_specifics.append("The following systems have failed due to {0} error:".format(
            keyword.replace("_", " ")))
        for value in sorted(values[1]):
            if keyword in ['segmentation_fault', 'missing_conf_file', "incomplete_simulation", "others",
                           "cannot_write_report"]:
                failing_systems_specifics.append(" + {0}".format(value))
            elif keyword == "missing_nma_ter":
                failing_systems_specifics.append(" - {0[0]} : missing TER mark after residue {0[1}".format(value))
            elif keyword == "mismatch":
                failing_systems_specifics.append(" - {0[0]} : the problematic atom is {0[1]} in residue "
                                                 "{0[2]}".format(value))
            elif keyword == "missing_atom":
                failing_systems_specifics.append(" - {0[0} : is missing the atom {0[1]} in the residue "
                                                 "{0[2]}".format(value))
            elif keyword == 'configuration_file_problem':
                failing_systems_specifics.append(" - {0} : the configuration file format is incorrect"
                                                 ".".format(value))
            elif keyword == 'missing_template':
                failing_systems_specifics.append("-  {0[0]} : is missing the template {0[1]} .".format(value))
            elif keyword == "missing_obc_parameters":
                failing_systems_specifics.append("- {0[0]} : is missing the obc parameters for the atom {0[1]} in "
                                                 "residue {0[2]} in chain {0[3]}".format(value))
            elif keyword == "rotamers_template_missmatch" :
                failing_systems_specifics.append(" - {0[0]} : has a problem in the atoms '{0[1]}' , '{0[2]}' in the "
                                                 "residue {0[3]}:{0[4]}:{0[5]} (resname:resnum:chain)".format(value))

report_text += "\n".join(failing_systems_classification)
report_text += "\nFailing systems specifications:\n"
report_text += "\n".join(failing_systems_specifics) + "\n"
with open(args.report, 'w') as outfile:
    outfile.write(report_text)
