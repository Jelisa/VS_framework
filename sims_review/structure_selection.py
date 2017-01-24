"""
author: jelisa iglesias
e-mail: jelisa.iglesias@gmail.com
"""

from argparse import ArgumentParser
import os
# import sys
import numpy as np
import re
import datetime
import logging

import prody

import parameters_help


def extract_minimum_energy_models(report_filename, number_skipped_systems, deltag_limit):
    # This block will analyze the report file
    report_header = []
    report_values = []
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
    if not report_values:
        number_skipped_systems += 1
        logging.error(" The report file was empty but for the header.")
        return False
    report_values_array = np.asarray(report_values)
    # model_number_index = report_header.index("AcceptedSteps")
    binding_energy_index = report_header.index("Binding Energy")
    minimum_energy_value = report_values_array[:, binding_energy_index].min()
    minimum_energy_model = report_values_array[:, binding_energy_index].argmin()
    # The following instruction returns the models number to use as the best ones including the one with minimum energy.
    # print 'b,', minimum_energy_value, minimum_energy_value + args.deltag
    # print np.where(report_values_array[:, binding_energy_index] <= minimum_energy_value +
    #                                              args.deltag)[0]
    systems2use = [v for v in np.where(report_values_array[:, binding_energy_index] <= minimum_energy_value +
                                       deltag_limit)[0] if v != minimum_energy_model]

    if len(systems2use) == 0:
        logging.info(" There's no other pose within 5 Kcal of the best pose.")
    logging.info(" The models to use are {}".format(", ".join(map(str, systems2use + [minimum_energy_model]))))
    return systems2use, minimum_energy_model, report_values_array, binding_energy_index, minimum_energy_value


parser = ArgumentParser()
parser.add_argument("-input", required=True, nargs="+", help=parameters_help.structure_selection_input_desc)
parser.add_argument("-deltag", default=5, type=int, help=parameters_help.structure_selection_deltag_desc)
parser.add_argument("-ligand_chain", default="Z", help=parameters_help.structure_selection_ligand_chain_desc)
parser.add_argument("-single_structure", default=False, action="store_true",
                    help=parameters_help.structure_selection_single_structure)
parser.add_argument("-initial_energies", action="store_true", help=parameters_help.initial_energies)
parser.add_argument("-output_folder", default="./")
parser.add_argument("-only_statistics", default=False, action="store_true",
                    help=parameters_help.structure_selection_only_statistics)
args = parser.parse_args()

logging.basicConfig(filename="structures_extraction.log", format="%(message)s", level=logging.INFO, filemode="w")
logging.info("{} : Program starting".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M")))

if args.output_folder != "./":
    output_folder = args.output_folder
else:
    output_folder = "./selected_models"
if not os.path.isdir(output_folder) and not args.only_statistics:
    os.mkdir(output_folder)
if output_folder[-1] != os.sep:
    output_folder += os.sep

skipped_systems = 0
statistics = {}
energies_values = {}
initial_energies_values = {}
for folder in args.input:
    # Adding the OS separator to the folder path
    if folder[-1] != os.sep:
        folder += os.sep
    logging.info("Working with the folder: {0}".format(folder))
    if not os.path.isdir(folder):
        logging.error(" The folder {0} doesn't exist, thus it will be ignored.".format(folder))
        skipped_systems += 1
        continue
    # Checking if the folder to look for exists if not ignore the system.
    pele_report_filename = ""
    trajectory_filename = ""
    for filename in os.listdir(folder):
        if "report" in filename:
            pele_report_filename = folder + filename
        elif "trajectory" in filename:
            trajectory_filename = folder + filename
    if not pele_report_filename or not trajectory_filename:
        logging.error(" Report or trajectory file missing the system will be discontinued.")
        skipped_systems += 1
        continue
    # obtain the simulation id from the folder path.
    pattern = re.search(r'{0}*(\w+_\d+){0}.*'.format(os.sep), folder)
    if pattern is None:
        if folder.split('/')[-2] != "output":
            folder_name = folder.split('/')[-2]
        else:
            folder_name = folder.split('/')[-3]
    else:
        folder_name = pattern.group(1)

    try:
        models_2_use, pele_minimum_energy_model, pele_energies_values_array,\
            pele_binding_energy_index, pele_minimum_energy_value = extract_minimum_energy_models(pele_report_filename,
                                                                                                 skipped_systems,
                                                                                                 args.deltag)
    except ValueError:
        print "error"
        continue
    initial_energies_values[folder_name] = pele_energies_values_array[0][pele_binding_energy_index]
    trajectory = prody.parsePDB(trajectory_filename)
    if args.single_structure:
        models_2_use = [pele_minimum_energy_model]
    else:
        # The following block loads the trajectory and computes the rmsd between the different models to use,
        # using as reference the one with minimum energy.
        ligand = trajectory.select("(within 6 of chain {0}) and not hydrogen".format(args.ligand_chain))
        # print minimum_energy_model
        reference_coordinates = ligand.getCoordsets(pele_minimum_energy_model)
        if models_2_use:
            minimim_rmsd_value = 0.25
            maximum_rmsd_value = 2.5
            for index, x in enumerate(models_2_use):
                rmsd = prody.calcRMSD(reference_coordinates, ligand.getCoordsets(x))
                # print 'c', minimum_energy_model, x, rmsd
                if rmsd < minimim_rmsd_value or rmsd > maximum_rmsd_value:
                    models_2_use.pop(index)
        models_2_use.insert(0, pele_minimum_energy_model)
    pattern = re.search(r"([\w-]*)trajectory", trajectory_filename.split(os.sep)[-1])
    system_name = pattern.group(1)
    if system_name[-1] != '_' and system_name[-1] != "-":
        system_name += "_"
    current_output_folder = output_folder + folder_name[:-1]
    general_output_name = current_output_folder + system_name
    binding_energy_text = []
    for index in models_2_use:
        general_output_name = current_output_folder + "_" + str(index + 1) + "/"
        # print general_output_name
        if not os.path.isdir(general_output_name) and not args.only_statistics:
            os.mkdir(general_output_name)
        output_filename = general_output_name + system_name + str(index + 1)
        if not args.only_statistics:
            prody.writePDB(output_filename, trajectory, csets=index)
            pele_bind_energy_filename = general_output_name + "pele_binding_energy.csv"
            with open(pele_bind_energy_filename, 'w') as outfile:
                # print binding_energy_text
                outfile.write("\n".join(binding_energy_text))
        binding_energy_text.append("{0},{1:.3f}".format(output_filename,
                                                        pele_energies_values_array[index, pele_binding_energy_index]))

    statistics[general_output_name] = len(models_2_use)
    # print 'aa', general_output_name
    energies_values[output_filename] = pele_minimum_energy_value
if args.only_statistics:
    statistics_filename = "general_statistics.csv"
    energies_filename = "minimum_PELE_binding_energies.csv"
    initial_energies_filename = "initial_PELE_binding_energies.csv"
else:
    statistics_filename = output_folder + "general_statistics.csv"
    energies_filename = output_folder + "minimum_PELE_binding_energies.csv"
    initial_energies_filename = output_folder + "initial_PELE_binding_energies.csv"
with open(statistics_filename, 'w') as outfile:
    outfile.write("\n".join(["{0},{1}".format(key[:-1], value) for key, value in statistics.iteritems()]))
# print energies_values
with open(energies_filename, 'w') as outfile:
    outfile.write("ID,pele_score\n")
    outfile.write("\n".join(["{0},{1}".format(key, value) for key, value in energies_values.iteritems()]))
if args.initial_energies:
    with open(initial_energies_filename, 'w') as outfile:
        outfile.write("ID,pele_score\n")
        outfile.write("\n".join(["{0},{1}".format(key, value)
                                 for key, value in initial_energies_values.iteritems()]))
logging.info("{} : Program finished correctly.".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M")))
logging.shutdown()
