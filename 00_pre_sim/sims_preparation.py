#!usr/bin/python

# TODO: Change the call to something that grabs the outputs to process them

""""
This program prepares the complexes and files to launch PELE simulations.
It requires the independent programs:
    - mutations_program.py script
    - PlopRotTemp.py
    - pdbconvert utility from Schrodinger
It has the python dependencies:
    - numpy
    - prody (from mutations_program.py)
author: Jelisa Iglesias
mail: jelisa.iglesias@bsc.es
"""

import datetime
import fnmatch
import logging
import os
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from re import search, match, findall, finditer, IGNORECASE
import re
import shutil
from string import Template
from subprocess import check_output, CalledProcessError, check_call, STDOUT
import subprocess
from glob import glob
from tempfile import mkdtemp
import numpy as np
# from modules import module

import enviroment_parameters
import external_software_calls as esc
import help_descriptions
import constant_values as const_values


def createsymboliclink(original, link):
    """
    This function tries to create the folder and catches the errors.
    :param original: a string containing the path of the original folder.
    :param link: a string containing the path  where the link should be generated.
    :return: void
    """
    try:
        os.symlink(original, link)
    except OSError as e:
        if e[1] == "File exists":
            logging.info(" - The link {} already exists. So it won't be created again.".format(link))
        else:
            logging.error(" - The link {} couldn't be created.".format(link))
    else:
        logging.info("   - Link created correctly")


""" NEW FUNCTIONS FROM HERE DOWNWARDS, EVERYTHING ABOVE HAS TO BE RECHECKED """


def check_output_folder_and_create_it(output_directory, interactive):
    """
    This function checks whether the output directory has an appropiate format (ends with the OS separator)
    and if it exists or not, if it doesn't exist it will ask the user if he/she wants to create it and gives
    the possibility to change it, as long as the program is in interactive mode.
    :param output_directory: a string containing the path where the output should be created.
    :param interactive: a boolean indicating whether the program should ask questions to the user or not. Default:True
    :return: a string containing the path for the output with the OS separator as the last character.
    """
    # check if the path from the option subfolder_path, ends in the system separator otherwise it adds it,
    # by default the program creates the subfolders in the directory where the program is executed.
    if output_directory != os.sep:
        output_directory = output_directory + os.sep
    if not os.path.isdir(output_directory):
        if interactive:
            logging.critical("ERROR: The path to create the subfolders doesn't "
                         "exist, it should be created previously by the user."
                         " The program will be terminated.")
            logging.shutdown()
            raise IOError("In the not interactive mode the user should have "
                          "created previously the general folder where the new "
                          "folders should be created.\nTerminating the program")
        else:
            print "The path to create the subfolders doesn't exist"
            answer = raw_input("Do you want to create the folder {}? Y/n\n".format(output_directory))
            if search("n.*", answer, IGNORECASE):
                answer = raw_input("Do you want to create and use another path? y/[N]?\n")
                if search("n", answer, IGNORECASE) or answer == "":
                    logging.critical("ERROR: The path to create the subfolders doesn't exist\n"
                                 "The user has chosen to terminate the program")
                    sys.exit("The program will be terminated.")
                else:
                    new_subfolders_path = raw_input(
                        "Please introduce the path where the folders should be generated:\n")
                    print "Creating the folder {}".format(new_subfolders_path)
                    os.mkdir(new_subfolders_path)
                    current_directory = os.getcwd()
                    output_directory = current_directory + new_subfolders_path
                    if not os.path.isdir(output_directory):
                        output_directory = new_subfolders_path
                    logging.info("The folder {} has been created and the new folders will be generated in it.")
            elif search("y.*", answer, IGNORECASE) or answer == "":
                os.mkdir(output_directory)
                print "Creating the folder {}".format(output_directory)
                logging.info("The folder {} has been created and the new folders will be generated in it.")
    return output_directory


def schrodinger_mae2pdb_converter(file_without_extension, terminate_if_fail):
    """
    This function receives a .mae file and converts it to .pdb format using the
    schrodinger pdb converter utility.
    :param file_without_extension: a string containing the name of the .mae file without extension.
    :return: a string containing the name of the file in .pdb format, and a boolean indicating if
    there have been errors or not
    :param terminate_if_fail: a boolean indicating whether the program should fail or not
    in case of error when executing the command.
    """
    maefile = "{0}.mae".format(file_without_extension)
    pdbfile = "{0}.pdb".format(file_without_extension)
    command2call = esc.schrodinger_mae2pdb_convert_command.format(maefile, pdbfile)
    try:
        check_call(command2call.split())
    except CalledProcessError as error_message:
        if terminate_if_fail:
            logging.critical("ERROR: The pdbconverter from schrodinger couldn't be executed.")
            logging.info("The error is in the command ' {} '".format(command2call))
            logging.shutdown()
            raise OSError("Program terminated due to the error:\n{}".format(error_message))
        else:
            logging.error(" - ERROR: The pdbconverter from schrodinger couldn't be executed.")
            error = True
    else:
        error = False
    return pdbfile, error


def check_mutations_program_output(command_call, output_directory, terminate_if_fail, verbose=True):
    """
    This function uses the subprocess module to call the external mutations program using the command
    provided and checks if it executes correctly, and returns the new file name. If the external
    program fails it will generate a file with the output of the failing command.
    :param command_call: a string containing a call to the external program mutations_program.py
    :param output_directory: a string containing the path where the output should be generated.
    :return: a string containing the name of the new file generated and a boolean indicating whether
    the mutations_program.py has failed or not.
    :param terminate_if_fail: a boolean indicating whether the program should fail or not
    in case of error when executing the command.
    """

    call_output = subprocess.check_output(command_call.split(), stderr=subprocess.STDOUT)
    generated_file = ""
    errors = False
    # We want to use the file generated by the mutations program since it will have the correct format for PELE.
    correct_output = re.search(r"Writing the structure to (.*.pdb)", call_output, re.IGNORECASE)
    if correct_output is None:
        error_filename = output_directory + "receptor_preparation_error_log.txt"
        logging.error(
            " - ERROR: The mutations program hasn't been able to process the file.\n"
            " - More information about this error in the file {}".format(error_filename))
        errors = True
        with open(error_filename, 'w') as error_file:
            error_file.write(call_output)
        if terminate_if_fail:
            logging.info("The error is in the command ' {} '\nThe program will be terminated.".format(command_call))
            logging.shutdown()
            raise OSError("Program terminated due to an error executing the command: {}".format(command_call))
    else:
        generated_file = correct_output.group(1)
        if verbose:
            logging.info(" - The file ready for pele it's called: {0}".format(generated_file))
    return generated_file, errors


def receptor_check(receptor_filename, output_directory, logging_file, log):
    """
    This function checks the format of the receptor file only looking at its extension, it should be .mae or .pdb,
    if it's a .mae file it is converted to .pdb using the schrodinger pdbconverter utility, once in .pdb format it's
    copied into the output directory, in case that a file with the same name exists in the output directory the
    program will add the _1 suffix at the end of the name but before the extension. Finally it will run the copied
    file through the mutations_program.py so it has the correct format for the PELE simulations.
    :param logging_file:
    :param receptor_filename: a string containing the relative path to the file containing the receptor in
                             PDB or MAE formats.
    :param output_directory: a string containing the output directory where the receptor should be copied.
    :param log: the logger object that registers a log file throughout the simulation.
    :return: a string containing the receptor information in pdb format, and another string
    containing the format of the original receptor file.
    """
    # check the format of the receptor file if it exists, if it isn't valid the program is terminated.
    format_ok = findall(const_values.accepted_formats, receptor_filename)
    if not format_ok or len(format_ok) > 1:
        log.critical("ERROR: The format of the receptor file isn't valid.\n\Terminating the program.")
        raise TypeError("ERROR: The format of the receptor file isn't valid.\n\Terminating the program.")
    else:
        log.info("The receptor file specified is: {}\n".format(receptor_filename))
        receptor_format = format_ok[0]
    log.info("INFO: The input files will be treated as single ligands and the receptor "
             "file will be used to form all the complexes\n")

    if receptor_format == '.mae':
        log.info("The receptor file will be converted to .pdb format using schrodinger's converter."
                 " This new file will be used as the receptor.")
        receptor_name = search(r"(.*)\.mae", receptor_filename)
        if receptor_name is None:
            raise Exception("Something went really wrong the receptor name is {} and it doesn't match the pattern"
                            "'(.*)\.mae' when it should.")
        receptor_filename, __ = schrodinger_mae2pdb_converter(receptor_name.group(1), True)
    # This option should be chosen only if the user has specified a receptor file.
    # If a receptor is given copy it to the new folder.
    if os.sep in receptor_filename:
        receptor_copy = output_directory + receptor_filename.split(os.sep)[-1]
    else:
        receptor_copy = output_directory + receptor_filename
    if os.path.isfile(receptor_copy):
        """If in the output folder exists a file with the same name as the receptor
        the receptor name will be changed adding '_1' to the name before the extension."""
        log.info(" - INFO: The receptor file it's already in the new subfolder,"
                 " the receptor will be renamed")
        pattern4new_receptor = search("(.*)(\.mae|\.pdb)", receptor_copy)
        if pattern4new_receptor is None:
            # This shouldn't happen ever since a similar pattern has been checked previously...
            log.critical("There's a problem in the code, tell the developer code L281")
            log.shutdown()
            sys.exit("You got somewhere where it should be impossible to get. Code L281")
        else:
            new_receptor_name = "{0}_1{1}".format(pattern4new_receptor.group(1), pattern4new_receptor.group(2))
            receptor_copy = new_receptor_name
            log.info("INFO: The receptor file in the folder {} is named {}".format(output_directory,
                                                                                   new_receptor_name))
    shutil.copyfile(receptor_filename, receptor_copy)

    log.info(" - Calling the 'mutations_program.py' to format the receptor file.")
    mutations_program_call = esc.mutations_program_command_receptor.format(receptor_copy)
    processed_receptor_filename, errors = check_mutations_program_output(mutations_program_call, output_directory,
                                                                         True)
    if not processed_receptor_filename:
        log.critical("")
        log.shutdown()
        sys.exit("There's been a problem while processing the receptor file with the mutations program.\n"
                 "Check the log file {0} for more information about the problem.".format(logging_file))
    text = ""
    with open(processed_receptor_filename, 'r') as receptor_file:
        for line in receptor_file:
            if line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("TER"):
                text += line
    return text, receptor_format


def input_format_check(tentative_input, logging):
    """
    This function iterates over the input files and checks that they have the appropriate format (it just checks the
    extension) eliminating the files without the right extension from the list of input files.
    :param tentative_input: a list containing the names of the files to use as input in string format.
    :param logging: the logger object that registers a log file throughout the simulation.
    :return: a list containing the files with the right extension if none has it the program is terminated.
    """
    # check if the input files are in the right format, if any of them isn't it will be skipped
    correct_files = []
    number_of_skipped_files = 0
    heading = "The following files don't have the right extension (the program will skip them):"
    # This loop iterates over the filenames and discards those that don't have the right extension.
    for filename in tentative_input:
        format_ok = findall(const_values.accepted_formats, filename)
        if not format_ok or len(format_ok) > 1:
            number_of_skipped_files += 1
            if number_of_skipped_files == 0:
                line = heading + " {}. {}".format(number_of_skipped_files, filename)
            else:
                line = " {}. {}".format(number_of_skipped_files, filename)
            logging.info(line)
        else:
            correct_files.append(filename)
    if number_of_skipped_files == 0:
        logging.info("All the input files have the right format.")
    logging.info("The input files are:\n{}".format("\n".join([filename for filename in correct_files])))
    return correct_files


def parse_template_text(text, keywords_pattern=const_values.templates_keywords_pattern):
    """
    This function searches for the specific pattern for what the program and the templates module consider
    keywords in the input text and then creates a python template object from the input.
    This keywords are the ones the program will need to fulfill when it uses the template object.
    :param keywords_pattern:
    :param text: a string containing the text that needs to be converted into a template.
    :return: a list containing strings with the keywords found in the text and a template object.
    """
    template_keywords = set(pattern.group(1) for pattern in
                            finditer(keywords_pattern, text))
    template = Template(text)
    return template_keywords, template


def parse_template_file(filename, template_type):
    """
    This function reads a file and converts it into a python Template object,
    at the same time that extracts the needed keywords for the template, and if
    the template is for a configuration file for PELE it also extracts the
    solvent name.
    :param filename: a string containing the path to the file to read.
    :param template_type: a string containing the type of template (PELE, adaptive)
    :return: a python Template object, a list containing the names of the keywords
    for the template, a string containing the kind of solvent for the simulation.
    """
    with open(filename, 'r') as template_file:
        template_or_text = "".join(template_file.readlines())
    template_keywords, template = parse_template_text(template_or_text)
    if template_type == 'pele':
        solvent = re.search(r'"solventType" : "(.*)",*', template_or_text)
        if solvent is None:
            logging.critical("ERROR: The template doesn't specify the solvent type PELE won't work.")
            logging.info("Terminating the program")
            logging.shutdown()
            sys.exit("Terminating the program.\nThe configuration file doesn't specify the solvent type.")
        else:
            solvent = solvent.group(1).lower()
        output_folder = re.search(r'"trajectoryPath" : "(.*)"', template_or_text)
        if output_folder is None:
            logging.critical("ERROR: The template doesn't specify the trajectory file PELE won't work.")
            logging.info("Terminating the program")
            logging.shutdown()
            sys.exit("Terminating the program.\nThe configuration file doesn't specify the trajectory file.")
        else:
            filename = output_folder.group(1)
            out_path = os.sep.join(filename.split(os.sep)[:-1])
    else:
        solvent = ""
        out_path = ""
    return template, template_keywords, solvent, out_path


def check_ligand(lig_filename, ligand_chain, output_directory, terminate_if_fail):
    """
    A function to read the ligand filename, check and fix the ligand so it
    is in the chain specified as the ligand chain and has a name different than UNK.
    :param output_directory:
    :param terminate_if_fail:
    :param lig_filename: a string containing the path to the ligand file in PDB
                        format to be read.
    :param ligand_chain: a string containing one letter to be used as the chain name
                        in the PDB file.
    :return: peptidic_ligand: a boolean indicating if the ligand is peptidic or not
             lig_name: and a string containing the ligand name.
             errors: a boolean indicating if there's been an error during the checks or not.
    """
    changes = False
    lig_name = ""
    errors = False
    with open(lig_filename, 'r') as lig_file:
        peptidic_ligand = False
        new_lig = ""
        for l in lig_file:
            if l.startswith("ATOM"):
                peptidic_ligand = True
            elif l.startswith("HETATM"):
                if not lig_name:
                    lig_name = l[17:20]
                if l[21] != ligand_chain:
                    changes = True
                    l = l[:21] + ligand_chain + l[22:]
                if search(r'<.>', l[17:20]) or search(r'UNK', l[17:20]):
                    changes = True
                    l = l[:17] + "LIG" + l[20:]
                    lig_name = 'LIG'
                if peptidic_ligand:
                    message = "WARNING: This ligand has peptidic and non peptidic residues. It will be skipped.\n" \
                              "Remember that the peptide-like ligands should be considered as one single ligand, " \
                              "with only one residue."
                    logging.info(message)
                    print "Check the system {} its ligand is peptide-like and it has several residues.".format(
                        lig_filename)
                    errors = True
                    break
            new_lig += l
    if not errors:
        if changes:
            with open(lig_filename, 'w') as lig_file:
                lig_file.write(new_lig)
        command_call = esc.mutations_program_command_lig_or_comp.format(lig_filename,
                                                                        ligand_chain) + " -opdb {}".format(lig_filename)
        output_filename, errors = check_mutations_program_output(command_call, output_directory,
                                                                 terminate_if_fail, False)
        if output_filename != lig_filename:
            logging.error("The mutations_program.py isn't working correctly.")
            errors = True

    return peptidic_ligand, lig_name, errors


def createfolder(folder2create_name):
    """
    This function creates the folder with all intermediates and catches the errors.
    :param folder2create_name:a string containing the folder to create.
    :return: a boolean indicating whether the creation has been possible or not
    """
    error = False
    try:
        os.makedirs(folder2create_name)
    except OSError as e:
        if e[1] == "File exists":
            logging.info(" - The folder {0} already exists. So it won't be created again.".format(folder2create_name))
        else:
            logging.error(" - The folder {0} couldn't be created.".format(folder2create_name))
            error = True
    else:
        logging.info("   - Folder {0} created correctly".format(folder2create_name))
    return error


def input_preprocess(original_file, original_output_directory, rewrite, terminate_if_errors_dict):
    """
    This function uses the name without extension of the original file as an ID to create a new subfolder
    in the output folder, then copies the original file into this new folder without rewriting the file
    if it exists unless explicitly said so.
    If the file is in mae format it's converted to pdb format using the pdb converter utility
    from Schrodinger.
    converts it to pdb format
    :param original_file: a string containing the file name of the original input
    :param original_output_directory: a string containing the main folder where all the output should be generated.
    :param rewrite: a boolean indicating whether the files should be rewritten or not in case they already exist.
    :param terminate_if_errors_dict: a dictionary containing the name of external programs as keys and booleans as
    values that are used to indicate whether or not the main program should terminate if the external program fails.
    :return: system_id: a string containing a substring of the original_file to be used as ID.
             new_copy: another string containing the path to the new file generated.
             new_output: a string containing the path where the output for this system should be written.
             error: a boolean indicating if there's been errors during the function.
    """
    error = False
    filename = original_file.split(os.sep)[-1]
    # The next instruction removes the extension of the file while keeping all possible . in the name.
    system_id = ".".join(filename.split(".")[:-1])
    new_output = original_output_directory + system_id + os.sep
    # Create the subfolder for the input file
    s_error = createfolder(new_output)
    if s_error:
        new_copy = ""
        error = True
    else:
        new_copy = new_output + filename
        if rewrite or not os.path.isfile(new_copy):
            shutil.copyfile(original_file, new_copy)
        # If the file is in .mae format it'll be converted to .pdb format to be able to do the checks
        # and modifications needed in a simpler way.
        is_mae_file = search("(.*)\.mae", new_copy)
        if is_mae_file:
            #  The input file has .mae format so it should be transformed to .pdb format
            new_copy, s_error = schrodinger_mae2pdb_converter(is_mae_file.group(1),
                                                              terminate_if_errors_dict['mae2pdb_convert'])
            if s_error:
                error = True
            else:
                terminate_if_errors_dict['mae2pdb_convert'] = False
    return system_id, new_copy, new_output, error


def create_complex_file(general_text, filename_lig, output_path, syst_id):
    """
    This function merges the protein information with the ligand information
    generating the complex file with the protein-ligand complex information in it.
    :param output_path:
    :param syst_id:
    :param general_text: a string containing all the protein information
    :param filename_lig: a string containing the name of the ligand file.
    :return: the filename of the complex, in case of failure it'll returns None
    """
    logging.info(" - Generating the complex file.")
    complex_text = general_text
    with open(filename_lig, 'r') as ligand_file:
        for line in ligand_file:
            if line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("TER"):
                complex_text += line
    if complex_text:
        new_file_name = output_path + syst_id + "_complex.pdb"
        with open(new_file_name, 'w') as complex_file:
            complex_file.write(complex_text)
    else:
        # If the code gets here is due to the fact that the receptor and the ligand files are empty so
        # it shouldn't continue.
        logging.critical("There's a problem in the code, tell the developer code L446")
        logging.shutdown()
        print "BU!! NOT WORKING!"  # TODO:Change this message to something more serious
        sys.exit("See what's happening..LINE 277 ")
    return new_file_name


def split_complex(filename, ligand_chain, output_path, syst_id):
    """
    This function reads a .pdb file and extracts whatever forms the chain
    specified by the ligand_chain argument, it should be the one where the ligand
    is present and it should not contain anything else.
    :param ligand_chain: a string containing the name of the chain to look for
    :param filename: a string containing the path to the file to read.
    :param output_path: a string containing the path where the new files should be generated.
    :param syst_id: a string containing the name of the system.
    :return: new_ligand_filename: the name of the file containing the extracted ligand
             receptor: a string containing the text for the part of the complex that isn't the ligand
             error: a boolean indicating whether the subroutine has failed or not.
             #lig_name: a string containing the ligand name
    """
    error = False
    with open(filename, 'r') as complex_text:
        ligand = ""
        receptor = ""
        # lig_name = ""
        for line in complex_text:
            if "ATOM" in line or "HETATM" in line:
                if line[21] in [ligand_chain.upper(), ligand_chain.lower()]:
                    ligand += line
                    # if "HETATM" in line and not lig_name:
                    #     lig_name = line[17:20]
                else:
                    receptor += line
            elif line == "TER\n":
                if ligand:
                    ligand += line
                else:
                    receptor += line
    if ligand:
        logging.info(" - Extracting the ligand from the complex.")
        new_ligand_filename = output_path + syst_id + "_ligand.pdb"
        with open(new_ligand_filename, 'w') as new_ligand_file:
            new_ligand_file.write(ligand)
            # In this case we don't really need the receptor in a separated file since we already have the complex.
            # receptor_filename = new_general_subfolder + new_folder_name + "_receptor.pdb"
            # with open(receptor_filename, 'w') as receptor_file:
            #     receptor_file.write(receptor)
    else:
        # If the program hasn't found the ligand chain it won't generate the files
        logging.info(" - There's nothing in chain {0} for the file {1}.".format(ligand_chain, filename))
        logging.warning(" -WARNING: The procedure for the folder {0} will be interrupted from now on.".format(
            output_path))
        error = True
        new_ligand_filename = ""
    return new_ligand_filename, receptor, error # , lig_name


def process_complex(complex_filename, ligand_chain, output_path, syst_id, terminate_if_fail):
    """

    :param complex_filename: a string containing the path to the complex file in PDB format.
    :param ligand_chain: a string containing the name of the chain where the ligand should be.
    :param output_path: a string containing the path where the new files should be generated.
    :param syst_id: a string containing the name of the system.
    :param terminate_if_fail: a boolean indicating whether the program should be terminated or not
                            in case of an error in the subroutine.
    :return:ligand_filename: a string containing the file with the ligand in PDB format
            ligand_name: a string contaiing the name of the ligand molecule.
            receptor_text: a string containing the information of the receptor in PDB format.
            no_need_for_template; a boolean indicating whether or not the template of the ligand should be generated.
            error: a string indicating the type of error ocurred if any.
    """
    error = ""
    ligand_name = ""
    no_need_for_template = False
    ligand_filename, receptor_text, s_warning = split_complex(complex_filename,
                                                              ligand_chain,
                                                              output_path, syst_id)
    if s_warning:
        error = "warning"
    else:
        no_need_for_template, ligand_name, s_error = check_ligand(ligand_filename, ligand_chain,
                                                                  output_path, terminate_if_fail)
        if s_error:
            error = "error"
    return ligand_filename, ligand_name, receptor_text,  no_need_for_template, error


def create_template_and_rotamerlib(initial_pdb, template_folder, rotamer_library_folder, terminate_if_errors_dict):
    """
    This function calls the PlopRotTemp program to generate the template and
    rotamer libraries for PELE simulations, and places them into their
    respective subdirectories of DataLocal
    :param initial_pdb: a string containing the relative path to the file
                        containing the ligand in PDB format.
    :param template_folder: a string containing the relative path to the folder
                            where the template file should be placed.
    :param rotamer_library_folder: a string containing the relative path to the folder
                                  where the rotamer library file should be placed.
    :param out: a boolean indicating whether the program should be terminated or not
                in case of failure of the PlopRotTemp script. Default=False
    # :param e_counter: an integer counting the number of error occurred during the
    #                 program execution.
    :return: a string containing the template filename (including the relative path), and the
            updated e_counter.
    """
    error = False
    template_filename = ""
    ploprottemp_wd = mkdtemp(prefix="plop_exec_")
    ploprottemp_wd += os.sep
    ligand_filename_mae = ploprottemp_wd + initial_pdb.split(os.sep)[-1].split(".pdb")[0] + ".mae"
    # Convert the pdb file to .mae format
    command2call2 = esc.schrodinger_pdb2mae_convert_command.format(initial_pdb, ligand_filename_mae)
    try:
        check_call(command2call2.split())
    except CalledProcessError as error_message:
        error = True
        if terminate_if_errors_dict['pdb2mae_convert']:
            print "The pdbconverter of maestro failed when being execute. The program will be terminated."
            logging.critical("ERROR: The pdbconverter from schrodinger couldn't be executed.")
            logging.info("The error is in the command ' {} '".format(command2call2))
            logging.shutdown()
            raise IOError("Program terminated due to the error:\n{}".format(error_message))
        else:
            logging.warning("The templates won't be generated for this system, due to the fail od the pdbconverter"
                            "from Schrodinger to convert the ligand to a .mae file.")
    else:
        terminate_if_errors_dict['pdb2mae_convert'] = False
    if error:
        return template_filename, error
    # Let's generate the template and the rotamer library.
    command2call2 = esc.ploprottemp_command.format(ligand_filename_mae)
    try:
        ploprottemp_output = check_output(command2call2.split(), stderr=STDOUT, cwd=ploprottemp_wd)
    except CalledProcessError as error_message:
        if terminate_if_errors_dict['ploprottemp']:
            print "The PlopRotTemp program couldn't be execute. The program will be terminated."
            logging.critical("ERROR: The PlopRotTemp program couldn't be executed.")
            logging.info("The error is in the command ' {} '".format(command2call2))
            logging.shutdown()
            raise IOError("Program terminated due to the error:\n{}".format(error_message))
        else:
            logging.error(" - ERROR: The command ' {} ' couldn't be executed. Skipping system".format(command2call2))
            error = True
    else:
        terminate_if_errors_dict['ploprottemp'] = False
        if "rotamer library has been successfully created in" not in ploprottemp_output:
            error = True
            logging.error(" - The PlopRotTemp.py hasn't finished correctly, so the template and/or the rotamer "
                          "library may be missing.")
            logging.info(" - This system will be discontinued.")
    if not error:
        for name4file in os.listdir(ploprottemp_wd):
            if fnmatch.fnmatch(name4file, "???z"):
                template_filename = template_folder + name4file.split(os.sep)[-1]
                initial_name_temp = ploprottemp_wd + name4file.split(os.sep)[-1]
                shutil.move(initial_name_temp, template_filename)
            elif fnmatch.fnmatch(name4file, "*.rot.assign"):
                initial_name_rotlib = ploprottemp_wd + name4file
                with open(initial_name_rotlib, 'r') as filein:
                    # This checks whether the rotamer library is empty or not, sometimes
                    # ploprottemp generates libraries with just the header.
                    number_of_lines = len(filein.readlines())
                    if number_of_lines <= 1:
                        logging.info("WARNING: The rotamer library generated by PlopRotTemp is empty, removing it.")
                        os.remove(initial_name_rotlib)
                    else:
                        destination_file = rotamer_library_folder + name4file.split(os.sep)[-1]

                        logging.info("The library of rotamers has been correctly generated and is: {}".format(
                            destination_file))
                        shutil.move(initial_name_rotlib, destination_file)
    shutil.rmtree(ploprottemp_wd)
    return template_filename, error


def check_previous_next(prev_element, next_element, current_element, n, reset):
    """
    This function is used to check if the previous element and the following one are ok and if the
    current element should be constrained or not.
    :param prev_element: a string containing the information of the previous element
    :param next_element: a string containing the information of the following element
    :param current_element: a string containing the information of the current element
    :param n: an integer that acts as a counter
    :param reset: an integer specifying the threshold where the counter should be reset.
    :return: n: the updated counter
             atom_id: a string containing the id for PELE of the current element in case n%reset = 0,
                     otherwise it'll be empty.
    """
    # if args.debug:
    #     print type(n), n, type(reset), reset
    chain_id = current_element[21]
    resnum = current_element[22:26].strip()
    atom_name = current_element[12:16].replace(" ", "_")
    if prev_element == "TER":
        atom_id = '"{}:{}:{}"'.format(chain_id, resnum, atom_name)
        n = 0
    elif next_element == "TER":
        atom_id = '"{}:{}:{}"'.format(chain_id, resnum, atom_name)
        n = 0
    elif (n % reset) == 0:
        try:
            atom_id = '"{}:{}:{}"'.format(chain_id, resnum, atom_name)
        except IndexError:
            print current_element
            raise IOError("Bad terminations")
        n = 0
    else:
        n += 1
        atom_id = ""
    return n, atom_id


def obtain_constraints_from_pdb(pdb_text, every, constraint, atoms2restrain):
    """
    This function reads a string containing information about the receptor in PDB format and process it
    to generate the constrains lines for the PELE configuration file to run simulations.
    :param pdb_text:a string containing the receptor information in PDB format.
    :param every: an integer specifiying every how many atoms should the program put a resntrain.
    :param constraint: an integer specifying the stregnth of the constraint
    :param atoms2restrain: a list containing the names of the atoms to restrain.
    :return: constraints_text: a string containing all the contrains to write in the configuration file.
    """
    constraints_text = ""
    counter = 0
    pdb_text = [lin for lin in pdb_text.split('\n') if lin.startswith("TER") or lin[12:16] in atoms2restrain and
                (lin.startswith("ATOM") or lin.startswith("HETATM"))]
    constrain_line = '\n{0:6}{{ "type" : "constrainAtomToPosition", "springConstant":{1}, "equilibriumDistance": 0.0,' \
                     ' "constrainThisAtom": {2} }}'
    if pdb_text[0].strip() == "TER":
        logging.error(" - WARNING: Review the receptor: the initial residue from the "
                      "first chain is preceded by a TER mark.")
    for index, element in enumerate(pdb_text):
        if element.strip() == "TER":
            continue
        else:
            chain_id = element[21]
            resnum = element[22:26].strip()
            atom_name = element[12:16].replace(" ", "_")
            if index == 0:
                if pdb_text[index + 1].strip() == "TER":
                    logging.error(" - WARNING: Review the receptor: the initial residue from the "
                                  "first chain is followed by a TER mark.")
                counter, atom_id = check_previous_next("", pdb_text[index + 1], element, counter, every)
                constraints_text += constrain_line.format(" ", constraint, atom_id)
                counter += 1
                continue
            if 0 < index < len(pdb_text) - 1:
                counter, atom_id = check_previous_next(pdb_text[index - 1], pdb_text[index + 1],
                                                       element, counter, every)
            else:
                counter, atom_id = check_previous_next(pdb_text[index - 1], "", element, counter, every)
            if counter == 0:
                counter += 1
            if atom_id:
                constraints_text += ',{0}'.format(constrain_line.format(" ", constraint, atom_id))
            else:
                continue
    if not constraints_text:
        logging.error(" - ERROR: The constraints haven't been created for this system. It'll be discontinued.")
    else:
        constraints_text += '\n'
    return constraints_text


def compute_center_of_mass(lig_filename):
    """
    This function computes the center of mass of the ligand provided in the parameters.
    :param lig_filename: a string containing the path to the ligand file in PDB format.
    :return: a string containing the coordinates of the ligand center of mass formatted for a
            PELE configuration file.
    """
    with open(lig_filename, 'r') as filein:
        molecule_atoms_text = [l for l in filein if l.startswith("ATOM") or l.startswith("HETATM")]
    masses_sum = 0.0
    coordinates_sumatory = np.zeros(3)
    for atom_line in molecule_atoms_text:
        resname = atom_line[12:16].lower()
        if resname in const_values.atomic_table_weights.keys():
            atom_mass = const_values.atomic_table_weights[resname]
        else:
            resname = resname.strip()
            if resname in const_values.atomic_table_weights.keys():
                atom_mass = const_values.atomic_table_weights[resname]
            else:
                if resname[:2] in const_values.atomic_table_weights.keys():
                    atom_mass = const_values.atomic_table_weights[resname[:2]]
                elif resname[0] in const_values.atomic_table_weights.keys():
                    atom_mass = const_values.atomic_table_weights[resname[0]]
                else:
                    print "The molecule has an atom not present in the atomic_weights table."
                    print "The atom is: {}".format(resname)
                    return False
        coordinates = np.asarray([float(x) for x in atom_line[30:54].split()])
        coordinates_sumatory += coordinates * atom_mass
        masses_sum += atom_mass
    center_of_mass = coordinates_sumatory / masses_sum

    return "{0[0]:.4f}, {0[1]:.4f}, {0[2]:.4f}".format(center_of_mass)


def generate_pele_conf_values(template_keywords, complex_filename, receptor_text, ligand_complete_path, gen_param,
                              syst_id, hbonds_dictio={}, adaptive_boolean=False):
    """
    This function generates the values needed for each of the keywords in the PELE configuration template file.
    :param hbonds_dictio: a dictionary containing the atoms to use with the hbond keyword
    :param template_keywords: a list containing the keywords present in the PELE configuration template file.
    :param complex_filename: a string with the path to the complex file in PDB format to be used as
                            input for the simulation.
    :param receptor_text: a string containing the receptor information in PDB format.
    :param ligand_complete_path: a string containing the path to the ligand file (in PDB format) to be used
                                for the simulation.
    :param gen_param: an ArgumentParser object containing all the options for the program.
    :param syst_id: a string containing the name of the system.
    :param adaptive_boolean: a boolean indicating whether or not the type of PELE simulation should be
                            an adaptive sampling protocol or not.
    :return:keywords_values: a dictionary containing the values for each of the keys in the PELE
                            configuration template file.
            error: a boolean indicating whether or not there has been errors in the subroutine.
    """
    keywords_values = {}
    error = False
    for keyword in template_keywords:
        if adaptive_boolean and keyword in enviroment_parameters.adaptive_sampling_keywords:
            keywords_values[keyword] = "${0}".format(keyword)
        elif search(".*constraints", keyword, IGNORECASE):
            # obtain constraints from a pdb file
            constraints = obtain_constraints_from_pdb(receptor_text, gen_param.fix_every_x_atoms,
                                                      gen_param.constraint_strength, gen_param.atoms2constraint)
            if not constraints:
                print "Couldn't generate the constraints for file {}".format(receptor_text)
                logging.error("  ERROR: Couldn't generate the constraints for file {}".format(receptor_text))
                logging.info(" This system will be discontinued since constraints couldn't be generated.")
                error = True
                break
            keywords_values[keyword] = constraints
        elif search("center_*of_*mass", keyword, IGNORECASE):
            mass_center = compute_center_of_mass(ligand_complete_path)
            if not mass_center:
                logging.warning("  WARNING: This system will be discontinued. Couldn't compute the center of mass.")
                logging.info("INFO: to solve this problem talk with he developer and provide the ligand file.")
                error = True
                break
            keywords_values[keyword] = mass_center
        elif search("input.*", keyword, IGNORECASE):
            keywords_values[keyword] = complex_filename.split('/')[-1]
        elif search("license.*", keyword, IGNORECASE):
            keywords_values[keyword] = gen_param.pele_license
        elif search("system.*", keyword, IGNORECASE):
            keywords_values[keyword] = syst_id
        elif search("hbond2", keyword, IGNORECASE):
            keywords_values[keyword] = hbonds_dictio[syst_id][1]
        elif search("hbond1", keyword, IGNORECASE):
            keywords_values[keyword] = hbonds_dictio[syst_id][0]
        elif search("COMPLEXES", keyword):
            raise IOError("You're using a keyword from the adaptive simulations without using the adaptive parameters.")
        else:
            print "The keyword '{}' isn't a valid one".format(keyword)
            logging.warning("  WARNING: This system will be discontinued. The keyword {} isn't a valid one".format(
                keyword))
            logging.info("INFO: Check your keyword if it should be valid talk with he developer and provide "
                         "the template file.")
            error = True
            break
            # Every time a new keyword appears a new search should be added.
    return keywords_values, error


def generate_adaptive_conf_values(template_keywords, complex_filename, lig_name, new_pele_conf_file_filename,
                                  pele_license_path, syst_id):
    """
    This function generates the values for each of the keywords in the adaptive configuration file template.
    :param template_keywords: a list containing the keywords in the adaptive configuration file template.
    :param complex_filename: a string with the path to the complex file in PDB format to be used as
                            input for the simulation.
    :param lig_name: a string containing the name of the ligand.
    :param new_pele_conf_file_filename: a string containing the path to the PELE configuration file to use.
    :param pele_license_path: a string containing the path to the license file for PELE.
    :param syst_id: a string containing the genral name from the system.
    :return: adaptive_conf_key_val; a dictionary containing the  values for each of the keys in the adaptive_sampling.py
                            configuration template file.
    """
    adaptive_conf_key_val = {}
    for k in template_keywords:
        if search("input.*", k, IGNORECASE):
            adaptive_conf_key_val[k] = complex_filename.split('/')[-1]
        elif search("ligand_name", k, IGNORECASE):
            adaptive_conf_key_val[k] = lig_name
        elif search("pele_*conf_*file", k, IGNORECASE):
            adaptive_conf_key_val[k] = new_pele_conf_file_filename.split(os.sep)[-1]
        elif search("license*", k, IGNORECASE):
            adaptive_conf_key_val[k] = pele_license_path
        elif search("system.*", k, IGNORECASE):
            adaptive_conf_key_val[k] = syst_id
    return adaptive_conf_key_val


def generate_obc_parameters(lig_template, output_folder, terminate_if_fail):
    logging.info(" - Checking for the existence of the OBC parameters.")
    error = False
    original_obc_param = enviroment_parameters.pele_folders + os.sep + "Data/OBC/solventParamsHCTOBC.txt"
    copy_obc_folder = output_folder + const_values.datalocal_obc_param
    copy_obc_param = copy_obc_folder + "solventParamsHCTOBC.txt"
    if os.path.isdir(copy_obc_folder) and os.path.isfile(copy_obc_param):
        # Checks if the folder and file already exist otherwise the program creates them.
        logging.info(" - The parameters already exist.")
    else:
        logging.info(" - Generating the solvent parameters.")
        command2call = esc.obc_param_command.format(lig_template)
        # print command2call
        try:
            __ = check_output(command2call.split(), stderr=STDOUT)
        except CalledProcessError as error_message:
            if terminate_if_fail:
                print "The program failed to execute the OBC script. The program will be terminated."
                logging.critical("ERROR: The script to generate the OBC parameters couldn't be executed.")
                logging.info("The error is in the command ' {} '".format(command2call))
                logging.shutdown()
                sys.exit("Program terminated due to the error:\n{}".format(error_message))
            else:
                error = True
        else:
            obc_command_output_file = lig_template + "_OBCParams.txt"
            try:
                os.mkdir(copy_obc_folder)
            except OSError as e:
                if search(r"exists", e[1]):
                    pass
                else:
                    print "There's been a problem creating the OBC folder, review it.", e
                    error = True
            if not error:
                shutil.copyfile(original_obc_param, copy_obc_param)
                if os.path.isfile(obc_command_output_file):
                    with open(obc_command_output_file, 'r') as obc_template_file:
                        obc_text = "".join(obc_template_file.readlines()) + "\n"
                    with open(copy_obc_param, 'a') as obc_param_file:
                        obc_param_file.write(obc_text)
                    os.remove(obc_command_output_file)
    return error


def update_counter(error_bool, counter):
    if error_bool:
        counter += 1
    return error_bool, counter


def parse_hbond_file (filename):
    hbond_dictio = {}
    with open(filename) as infile:
        for line in infile.readlines():
            line = line.strip().split()
            if not line:
                continue
            if len(line) == 1:
                raise IndexError("Each line in the hbonds file should have at least two elements.")
            hbond_dictio[line[0]] = line[1:]
    return hbond_dictio


def main(args, log):
    command_first_execution = {"ploprottemp": True, 'mae2pdb_convert': True, 'pdb2mae_convert': True,
                               'mutations_program': True, 'obc_param_gen': True}

    # Set up a counter for non-critical errors and warnings to do a simple summary in the log file.
    total_errors_counter = 0
    total_warnings_counter = 0
    general_output_folder = check_output_folder_and_create_it(args.subfolders_path, args.not_interactive)
    if args.receptor:
        receptor_text, ori_format = receptor_check(args.receptor, general_output_folder, args.log_file, log)
        command_first_execution['mutations_program'] = False
        if ori_format == '.mae':
            command_first_execution['mae2pdb_convert'] = False
        create_complex = True
    else:
        create_complex = False
    # check if the input files are in the right format, if any of them isn't it will be skipped
    input_files = input_format_check(args.input_files, log)
    if not input_files:
        log.critical("\nERROR: None of the files specified in the input has the right format."
                     "\nTerminating the program.")
        log.info("The input files are {}".format("\n".join([filename for filename in args.input_files])))
        log.shutdown()
        raise IOError("ERROR: No valid input file.\nTerminating the program.")

    # Check the PELE configuration file template to look for the parameters to generate.
    if match("none", args.conf_template, IGNORECASE):
        log.info("The program won't generate any configuration file. Because the template has been specified as none.")
        print "No configuration file will be generated."
        generate_configuration_file_template = False
        keywords_in_the_pele_conf_template = {}
        solvent_type = ""
    else:
        log.info(" - Reading the configuration file Template.")
        generate_configuration_file_template = True
        pele_conf_template, keywords_in_the_pele_conf_template, \
        solvent_type, pele_sim_output_path = parse_template_file(args.conf_template, 'pele')
        if args.hb_from_file and not ("hbond1" not in keywords_in_the_pele_conf_template or
            "hbond2" not in keywords_in_the_pele_conf_template):
            raise IOError("When the option hb_from_file is used the template should have either the"
                          "keywork 'hbond1' or 'hbond2'.")
        elif args.hb_from_file:
            hbonds_dictionary = parse_hbond_file(args.hb_from_file)
        elif ("hbond1" not in keywords_in_the_pele_conf_template or
              "hbond2" not in keywords_in_the_pele_conf_template) and not args.hb_from_file:
            raise IOError("The pele template has the hbond keywords but none has been provided.")
        else:
            hbonds_dictionary = {}
    # Check the existence of configuration file template for the adaptive and generate a template and
    #  look for the parameters to generate.
    if args.adaptive_sampling:
        log.info(" - The program will prepare the input for the adaptiveSampling.py script.")
        if os.path.isfile(args.adaptive_sampling):
            log.info(" - Reading the configuration template for the adaptive sampling.")
            adaptive_conf_template, keywords_in_the_adaptive_conf_template, \
            __, pele_sim_output_path = parse_template_file(args.adaptive_sampling, 'adaptive')
            adaptive_simulation = True
        else:
            log.critical(
                "The file you have specified as the template for the adaptive control file doesn't exist.\n"
                "Terminating the program.")
            raise IOError("Terminating the program.\nThe template file for the adaptive "
                          "'{0}' doesn't exist.".format(args.adaptive_sampling))
    else:
        adaptive_simulation = False
    # This variable will store the name of the configuration files created.
    # This huge for loop takes care of all the preparation for each system. From creating the folders
    # to generate the configuration file needed for PELE, it also
    for filename in input_files:
        log.info("Working with file: {}".format(filename))
        input_id, input_filename, current_output_path, error = input_preprocess(filename, general_output_folder,
                                                                                args.rewrite, command_first_execution)
        if error:
            total_errors_counter += 1
            logging.error(" - This system will be discontinued.")
            continue
        # Check if a complex in PDB format already exists in the current output folder.
        if args.receptor:
            possible_complex_filename = current_output_path + const_values.complexes_pre_defined_name
            existing_complexes_files = glob(possible_complex_filename)
        else:
            possible_complex_filename = current_output_path + input_id + const_values.complexes_pre_created_name
            if args.debug:
                print possible_complex_filename
            existing_complexes_files = glob(possible_complex_filename)
        if args.debug:
            print '{}, {}'.format(0, existing_complexes_files)
        if args.debug:
            print '{}, {}, {}'.format(1, existing_complexes_files, args.rewrite)
        if not existing_complexes_files or args.rewrite:
            if args.debug:
                print '{}, {}'.format(2, create_complex)
            if create_complex:
                no_need_for_template, ligand_name, s_error = check_ligand(input_filename, args.ligand_chain,
                                                                          current_output_path,
                                                                          command_first_execution['mutations_program'])
                if command_first_execution['mutations_program']:
                    command_first_execution['mutations_program'] = False
                is_error, total_errors_counter = update_counter(s_error, total_errors_counter)
                if is_error:
                    continue
                complex_filename = create_complex_file(receptor_text, input_filename, current_output_path, input_id)
                ligand_filename = input_filename
            else:
                if args.debug:
                    print '2'
                command_to_call = esc.mutations_program_command_lig_or_comp.format(input_filename, args.ligand_chain)
                complex_filename, s_error = check_mutations_program_output(command_to_call,
                                                                           current_output_path,
                                                                           command_first_execution['mutations_program'])
                if command_first_execution['mutations_program']:
                    command_first_execution['mutations_program'] = False
                if s_error:
                    total_errors_counter += 1
                    continue
                needed_variables = process_complex(complex_filename, args.ligand_chain, current_output_path,
                                                   input_id, command_first_execution['mutations_program'])
                ligand_filename, ligand_name, receptor_text, no_need_for_template, s_error = needed_variables
                if s_error:
                    if s_error == "error":
                        total_errors_counter += 1
                    else:
                        total_warnings_counter += 1
                    continue
            if args.no_templates:
                no_need_for_template = True
        else:
            complex_filename = existing_complexes_files[0]
            needed_variables = process_complex(complex_filename, args.ligand_chain, current_output_path, input_id,
                                               command_first_execution['mutations_program'])
            # if args.debug:
            #     print needed_variables
            ligand_filename, ligand_name, receptor_text, no_need_for_template, s_error = needed_variables
            if s_error:
                if s_error == "error":
                    total_errors_counter += 1
                else:
                    total_warnings_counter += 1
                continue
            if args.no_templates:
                no_need_for_template = True
        # Lets create the template, rotamer library and the obc parameters if needed.
        if not no_need_for_template:
            current_pele_template_folder = current_output_path + const_values.datalocal_templates
            s_error = createfolder(current_pele_template_folder)
            if s_error:
                total_errors_counter += 1
                continue
            current_rotamerlibs_folder = current_output_path + const_values.datalocal_rotamerlib
            s_error = createfolder(current_rotamerlibs_folder)
            if s_error:
                total_errors_counter += 1
                continue
            # print command_first_execution
            lig_template, s_error = create_template_and_rotamerlib(ligand_filename, current_pele_template_folder,
                                                                   current_rotamerlibs_folder,
                                                                   command_first_execution)
            if s_error:
                total_errors_counter += 1
                continue
            if solvent_type == "obc":
                s_error = generate_obc_parameters(lig_template, current_output_path,
                                                  command_first_execution['obc_param_gen'])
                if command_first_execution['obc_param_gen']:
                    command_first_execution['obc_param_gen'] = False
                if s_error:
                    total_errors_counter +=1
                    continue
        # This block creates the configuration file into the folder.
        if generate_configuration_file_template:
            if args.debug:
                print '{}: {}'.format(3, receptor_text)
            pele_conf_key_val, s_errors = generate_pele_conf_values(keywords_in_the_pele_conf_template,
                                                                    complex_filename, receptor_text, ligand_filename,
                                                                    args, input_id, hbonds_dictionary,
                                                                    adaptive_simulation)
            if s_error:
                total_errors_counter += 1
                continue
            try:
                new_pel_conf_file_text = pele_conf_template.substitute(pele_conf_key_val)
            except (KeyError, ValueError) as e:
                log.critical("ERROR: the configuration file template isn't able to obtain all the needed values.")
                log.critical("The program has been terminated.")
                log.shutdown()
                print "Programing error. Talk with the developer. Code L549"
                if args.debug:
                    print keywords_in_the_pele_conf_template
                    print pele_conf_key_val.keys()
                print 0, e
                raise IOError("Error when creating the configuration file, missing keywords.")
            else:
                if args.adaptive_sampling:
                    new_pele_conf_file_filename = "{0}{1}_pele" \
                                                  "_adaptive_sampling_{2}.conf".format(current_output_path,
                                                                                   input_id,
                                                                                   args.conf_file_suffix)
                else:
                    new_pele_conf_file_filename = "{0}{1}{2}.conf".format(current_output_path, input_id,
                                                                          args.conf_file_suffix)
                    if args.debug:
                        print "{}: {}".format(4, new_pele_conf_file_filename)
                with open(new_pele_conf_file_filename, 'w') as new_conf_file:
                    new_conf_file.write(new_pel_conf_file_text)
            if pele_sim_output_path:
                pele_output = current_output_path + pele_sim_output_path
                s_error = False
                if os.path.isdir(pele_output):
                    logging.warning(" - WARNING: The output folder in the template already exists.")
                    total_warnings_counter += 1
                else:
                    s_error = createfolder(pele_output)
                    if s_error:
                        raise IOError("You don't have permissions to create the output folder {}".format(pele_output))
        if args.adaptive_sampling:
            adaptive_conf_key_val = generate_adaptive_conf_values(keywords_in_the_adaptive_conf_template,
                                                                  complex_filename, ligand_name,
                                                                  new_pele_conf_file_filename,
                                                                  args.pele_license, input_id)
            try:
                adaptive_conf_file_text = adaptive_conf_template.substitute(adaptive_conf_key_val)
            except (KeyError, ValueError) as e:
                log.critical("ERROR: the configuration file template isn't able to obtain all the needed values.")
                log.critical("The program has been terminated.")
                log.shutdown()
                print "Programing error. Talk with the developer. Code L903"
                if args.debug:
                    print keywords_in_the_adaptive_conf_template
                    print adaptive_conf_key_val.keys()
                print e
                sys.exit("Error when creating the configuration file, missing keywords.")
            else:
                if args.conf_file_suffix:
                    new_adaptive_conf_file_filename = "{0}{1}_{3}_adaptive_sampling_{2}.conf".format(current_output_path,
                                                                                                     input_id,
                                                                                                     args.conf_file_suffix,
                                                                                                     solvent_type)

                else:
                    new_adaptive_conf_file_filename = "{0}{1}_{2}_adaptive_sampling.conf".format(current_output_path,
                                                                                                 input_id,
                                                                                                 solvent_type)
                with open(new_adaptive_conf_file_filename, 'w') as new_conf_file:
                    new_conf_file.write(adaptive_conf_file_text)
            adaptive_output = search(r'"outputPath" *: *"(.*)",', adaptive_conf_file_text)
            if adaptive_output is None:
                log.critical("ERROR: the template for the adaptive configuration file doesn't present"
                             "and outputPath option.\nThe program won't work.")
                log.critical("The program has been terminated.")
                log.shutdown()
                sys.exit("The template for the adaptive configuration file doesn't present"
                         "and outputPath option.\nThe program won't work.")
            else:
                adaptive_output_path = "{0}{1}".format(current_output_path, adaptive_output.group(1))
            if not os.path.isdir(adaptive_output_path):
                s_error = createfolder(adaptive_output_path)


    # Now we have to create a file to submit the jobs. The file will depend on where they want to be run, so
    # we should use a template, also there are some requirements that change depending on BSC/AZ
    # In AZ the nodes have 12 physical cores so all the jobs should be  multiple of 24 or 12
    # with open(args.sub_template, 'r') as sub_template_file:
    #     sub_template_text = "".join(sub_template_file.readlines())

    print "Finished Correctly."
    log.info("{} : Program finished normally with {} warnings and {} non-critical errors.".format(
        datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), total_warnings_counter, total_errors_counter))
    log.shutdown()


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-input_files", "-input", required=True, nargs="+",
                        help=help_descriptions.input_files_desc)
    parser.add_argument("-receptor", default=enviroment_parameters.receptor_file_default,
                        help=help_descriptions.receptor_desc)
    parser.add_argument("-subfolders_path", default=enviroment_parameters.subfolders_path_default,
                        help=help_descriptions.subfolder_path_desc)
    parser.add_argument("-conf_template", "-conf_file_template",
                        default=enviroment_parameters.conformational_template,
                        help=help_descriptions.conf_template_desc)
    parser.add_argument("-conf_file_suffix", default=enviroment_parameters.conf_file_suffix,
                        help=help_descriptions.conf_file_suffix_desc)
    parser.add_argument("-adaptive_sampling", default=enviroment_parameters.adaptive_sampling_def,
                        help=help_descriptions.adaptive_sampling_desc)
    parser.add_argument("-ligand_chain", default=enviroment_parameters.ligand_chain,
                        help=help_descriptions.ligand_chain_desc)
    parser.add_argument("-no_templates", action="store_true", help=help_descriptions.no_templates)
    parser.add_argument("-rewrite", action="store_true", help=help_descriptions.rewrite_desc)
    parser.add_argument("-not_interactive", action="store_true", help=help_descriptions.not_interactive_desc)
    parser.add_argument("-log_file", default=enviroment_parameters.log_file_default,
                        help=help_descriptions.log_file_desc)
    parser.add_argument("-debug", action="store_true", help=help_descriptions.debug_desc)
    pele_group = parser.add_argument_group("PELE options", "The options to use for the PELE configuration file")
    pele_group.add_argument("-pele_folders", default=enviroment_parameters.pele_folders,
                            help=help_descriptions.pele_folder_desc)
    pele_group.add_argument("-pele_license", default=enviroment_parameters.pele_license_path,
                            help=help_descriptions.pele_license_desc)
    pele_group.add_argument("-every", "--fix_every_x_atoms", default=enviroment_parameters.default_every,
                            type=int, help=help_descriptions.every_desc)
    pele_group.add_argument("-constraint", "--constraint_strength", default=enviroment_parameters.constraints_strength,
                            help=help_descriptions.constraint_desc)
    pele_group.add_argument("-atoms2constraint", default=enviroment_parameters.atoms2apply_constraints_default,
                            nargs='+', help=help_descriptions.atoms2constraint_desc)
    pele_group.add_argument("-hb_from_file", default=enviroment_parameters.hb_from_file_def,
                            help=help_descriptions.hb_from_file_desc)
    external_soft_group = parser.add_argument_group("External Software path", "The options to specify the path "
                                                                              "for the external software.")
    external_soft_group.add_argument("-schrodinger_path",
                                     default=enviroment_parameters.schrodinger_path,
                                     help=help_descriptions.schrodinger_path_desc)
    external_soft_group.add_argument("-plop_path", default=enviroment_parameters.plop_path,
                                     help=help_descriptions.plop_path_desc)
    external_soft_group.add_argument("-mutations_program_path", default=enviroment_parameters.mutations_program_path,
                                     help=help_descriptions.mutations_program_path_desc)
    external_soft_group.add_argument("-obc_param_generator", default=enviroment_parameters.obc_param_generator_path,
                                     help=help_descriptions.obc_param_generator_desc)
    arguments = parser.parse_args()
    # The program will always generate a log file.
    logging.basicConfig(filename=arguments.log_file, format="%(message)s", level=logging.INFO, filemode="w")
    logging.info("{} : Program starting".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M")))

    main(arguments, logging)
# The first step to do is to check if the arguments are in the right format:
# the input files and the receptor should be either in .pdb or .mae formats,
# to check for this a pattern will be used.

# These are the commands that will be executed during the program, they change depending
# on where the program is being executed AZ or BSC
# Commands AZ
# This line is mandatory only for AZ, as the import
# The first thing is to load the schrodinger module
# module("load", "schrodinger/2016.01")
# Then we load the commands from the constant values depending on where we are.
# schrodinger_mae2pdb_convert_command = "pdbconvert -imae {} -opdb {}"
# schrodinger_pdb2mae_convert_command = "pdbconvert  -ipdb {} -omae {}"
# ploprottemp = "{}/utilities/python {}".format(args.schrodinger_path, args.plop_path)
# ploprottemp_command = ploprottemp + " {} -mae_charges=no -mtor=5 -g=30 -clean=yes"
# mutations_program_command = "python " + args.mutations_program_path + \
#                             " -ipdb {} -make_unique Z -gaps_ter "
# Commands BSC
