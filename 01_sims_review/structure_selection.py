"""
author: jelisa iglesia
e-mail: jelisa.iglesias@gmail.com
"""

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import os
import glob
import pandas as pd
import sys
import numpy as np
import re
import datetime
import logging
import prody
import parameters_help


def raise_parsing_output_error(directory, type_of_simulation):
    raise IOError("The directory {0} doesn't contain the right information"
                  "for a PELE simulation done with the protocol: {1}".format(directory,
                                                                             type_of_simulation))


def parse_adaptive_sampling_reports(directory):
    """
    This function reads all the report files generated during the simulation
    and returns a pandas dataframe with all the data
    :param directory: a string containing the name of the output folder of the simulation
    :return: pandas dataframe with all the values from the reports in it.
    """
    files2analyze = glob.glob(directory + '*/*report*')
    files2analyze.sort()
    if not files2analyze:
        return False
    reports = pd.DataFrame()
    for name in files2analyze:
        next_report = pd.read_csv(name, sep="    ", engine='python')
        epoch = name.split(os.sep)[-2]
        processor = name.split(os.sep)[-1].split('_')[-1]
        try:
            next_report['trajectory'] = glob.glob("{}{}/*trajectory_{}.pdb".format(directory, epoch, processor))[0]
        except IndexError:
            print "{}{}/*trajectory_{}.pdb".format(directory, epoch, processor)
            raise_parsing_output_error(directory, 'adaptive')
            return False
        reports = reports.append(next_report, ignore_index=True)
    reports.rename(columns={'numberOfAcceptedPeleSteps': 'accepted_steps'}, inplace=True)
    return reports


def parse_mpi_simulation_reports(directory):
    """
    This function reads all the files generated by an mpi PELE simulation and
    returns a pandas dataframe with all the data.
    :param directory: a string containing the output folder of the simulation
    :return:
    """
    files2analyze = glob.glob(directory + '*report*')
    files2analyze.sort()
    reports = pd.DataFrame()
    for name in files2analyze:
        next_report = pd.read_csv(name, sep="    ", engine='python')
        processor = name.split(os.sep)[-1].split('_')[-1]
        try:
            next_report['trajectory'] = glob.glob("{}*trajectory_{}.pdb".format(directory, processor))[0]
        except IndexError:
            # raise_parsing_output_error(directory, 'mpi')
            return False
        reports = reports.append(next_report, ignore_index=True)
    reports.rename(columns={'numberOfAcceptedPeleSteps': 'accepted_steps'}, inplace=True)
    return reports


def parse_single_core_simulation_report(directory):
    """
    This function reads the output files contained in a single processor PELE simulation
    and returns the report values and the name of the trajectory file in a pandas dataframe.
    :param directory: a string containing the path to the output folder of a PELE
    single core simulation.
    :return: Pandas DataFrame object containing the reports values and the name
    of their corresponding  trajectory files
    """
    report_file = glob.glob(directory + '*report*')
    # if len(report_file) > 1:
    #     raise_parsing_output_error(directory, 'single_core (too many report files)')
    # elif len(report_file) == 0:
    #     raise_parsing_output_error(directory, 'single_core (no report file present)')
    if len(report_file) != 1:
        return False
    reports = pd.read_csv(report_file[0], sep="    ", engine="python")
    try:
        reports['trajectory'] = glob.glob("{}*trajectory*.pdb".format(directory))[0]
    except IndexError:
        return False
    reports.rename(columns={'AcceptedSteps': 'accepted_steps'}, inplace=True)
    return reports


def write_selected_structure(new_directory, new_structure_filename, selected_structure_text):
    """
    This function creates a new directory and a file containing the selected_structure
    :param new_directory: a string containing the new folder to create
    :param new_structure_filename: a string containing the name of the new file to be created.
    :param selected_structure_text: a string containing the new selected structure text.
    :return:
    """
    try:
        os.mkdir(new_directory)
    except OSError as e:
        if e[0] == 17:
            # THis means the file already exists
            pass
        else:
            raise OSError(new_directory)
    with open(new_structure_filename, 'w') as outfile:
        outfile.write(selected_structure_text)


def main(args):
    """
    The main function of the program that takes care of minor checks and calls to the functions.
    :param args: the arguments from the command line.
    """
    logging.basicConfig(filename=args.log_file, format="%(message)s",
                        level=logging.INFO, filemode="w")
    logging.info("{} : Program starting".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M")))

    # Check for the existence of the output folder.
    # If it doesn't exist create it.
    if args.output_folder != "./":
        output_folder = args.output_folder
    else:
        output_folder = "./selected_models"
    if not os.path.isdir(output_folder) and not args.only_statistics:
        try:
            os.mkdir(output_folder)
        except OSError as e:
            # Since the new folder should be written into an existing path
            # the program will fail if the path to the folder doesn't exist.
            if e[0] == 2:
                logging.critical("The path to your output folder doesn't exist."
                                 "\n{0}".format(output_folder))
                logging.info("{0} : Abrupte termination of the program"
                             ".".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M")))
                logging.shutdown()
                raise OSError("The path to your output folder doesn't exist."
                              "\n{0}".format(output_folder))
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
        # Checking the existence of the folder, if it doesn't exist ignore the system.
        if not os.path.isdir(folder):
            logging.error("ERROR: The folder {0} doesn't exist, thus it will be ignored.".format(folder))
            skipped_systems += 1
            continue
        # Obtain the simulation name (or id) from the folder path.
        pattern = re.search(r'{0}*(\w+_\d+){0}.*'.format(os.sep), folder)
        if pattern is None:
            if "output" not in folder.split('/')[-2]:
                system_general_name = folder.split('/')[-2]
            else:
                system_general_name = folder.split('/')[-3]
        else:
            system_general_name = pattern.group(1)
        # Set the outputs names
        current_output_folder = output_folder + system_general_name
        # Parse the reports files.
        if args.simulation_type == "mpi":
            reports_values = parse_mpi_simulation_reports(folder)
        elif args.simulation_type == "adaptive":
            print 'h'
            reports_values = parse_adaptive_sampling_reports(folder)
            print reports_values
        else:
            reports_values = parse_single_core_simulation_report(folder)
        if reports_values is False:
            logging.error("WARNING: The report or trajectory file are missing, the system will be discontinued.")
            skipped_systems += 1
            continue

        if args.type_of_selection == "rmsd_clustering":
            raise LookupError("The selection using the rmsd clustering isn't implemented yet.")
        elif args.type_of_selection == "energy_clustering":
            raise LookupError("The selection using the energt clustering isn't implemented yet.")
        else:
            if args.energy_type == "binding_energy":
                # print reports_values
                model_2_use = reports_values.iloc[reports_values['Binding Energy'].argmin()]
            else:
                model_2_use = reports_values.iloc[reports_values['currentEnergy'].argmin()]
            trajectory_model = model_2_use['accepted_steps'] + 1
            trajectory_filename = model_2_use['trajectory']
            with open(trajectory_filename) as infile:
                trajectory_text = infile.read()
            pattern = re.search(r'MODEL +{}\n(.*?)\nENDMDL'.format(trajectory_model),
                                trajectory_text, re.DOTALL)
            if pattern is None:
                logging.warning("WARNING: There's some problem with the trajectory. The selected "
                                "model {0} cannot be found in {1}".format(trajectory_model,
                                                                          trajectory_filename))
                skipped_systems += 1
                continue
            else:
                selected_model_text = pattern.group(1)
                current_output_folder += "_{0}/".format(trajectory_model)
                output_filename = current_output_folder + \
                                  "{0}_{1}.pdb".format(system_general_name, trajectory_model)
                write_selected_structure(current_output_folder, output_filename, selected_model_text)

    logging.info("{} : Program finished correctly.\n"
                 "With {} warnings".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), skipped_systems))
    logging.shutdown()


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter,
                            description=parameters_help.structure_Selection_program_desc)
    parser.add_argument("-input", required=True, nargs="+",
                        help=parameters_help.structure_Selection_input_desc)
    parser.add_argument("-output_folder", default="", required=True,
                        help=parameters_help.structure_Selection_output_folder_desc)
    parser.add_argument("-ligand_chain", default="Z",
                        help=parameters_help.structure_selection_ligand_chain_desc)
    parser.add_argument("-only_statistics", default=False, action="store_true",
                        help=parameters_help.structure_selection_only_statistics)
    parser.add_argument("-initial_energies", action="store_true",
                        help=parameters_help.sims_review_initial_energies)
    parser.add_argument("-simulation_type", default="single_core",
                        choices=["single_core", "adaptive", "mpi"],
                        help=parameters_help.structure_Selection_simulation_type_desc)
    parser.add_argument("-log_file", default="structures_extraction.log",
                        help=parameters_help.log_file_desc)

    subparsers = parser.add_subparsers(dest='type_of_selection',
                                       help=parameters_help.structure_Selection_selection_criteria)
    parser_minimum_energy = subparsers.add_parser('minimum_energy',
                                                  description=parameters_help.structure_Selection_min_ener_desc)
    parser_minimum_energy.add_argument('-energy_type', default="binding_energy",
                                       const='binding_energy', nargs="?",
                                       choices=["binding_energy", "pele_energy"],
                                       help=parameters_help.ss_min_energy_type_desc)
    parser_rmsd_clust = subparsers.add_parser('rmsd_clustering',
                                              description=parameters_help.structure_Selection_rmsd_clust_desc)
    parser_energy_clust = subparsers.add_parser('energy_clustering',
                                                description=parameters_help.structure_Selection_ener_clust_desc)
    parser_energy_clust.add_argument("-deltag", default=5, type=int,
                                     help=parameters_help.ss_energy_c_deltag_desc)
    # parser.add_argument("-single_structure", action="store_true",
    #                     help=new_parameters_help.structure_selection_single_structure)

    arguments = parser.parse_args()
    main(arguments)
