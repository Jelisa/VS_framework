"""
This is the script that takes care of preparing and launching the SF computation.
This script will take as input a list of files and will copy them into a folder called scoring_functions_values.
 (or one given by the user) and will generate one folder per structure with the same name as the input file without the
  format extension.
Then it will prepare all the files needed to compute the SFs chosen by the user, by default it will compute:
 - Glide
 - vina
 - xscore
 - dsx
 - rf-score
 - mmgbsa (energy descriptors)
 - binana descriptors (structural descriptors)
After computing all the scoring functions it will parse the different outputs from all the sf, and will generate a file
 with all the values. If the user specifies the energies values they'll be added to the final file.
Due to how glide works if this scoring functions has to be computed the program will generate the following folders:
 - a folder called glide inside the scoring_functions_values folder where the glide computation will be done
 - a structures folder inside the glide folder, to put all the symbolic links to the structures.

author: jelisa iglesias
email: jelisa.iglesias@gmail.com
"""

import argparse
import re
import os
import sys
from string import Template
from subprocess32 import check_output, STDOUT, CalledProcessError, check_call, TimeoutExpired
import logging
import datetime

import external_software_paths
# import modules  # This is for AZ.


def extract_ligand(pdb_filename, general_name, ligand_chain, executing_folder):
    """
    This function creates 4 different files in .pdb format
    :param pdb_filename:
    :param general_name:
    :param ligand_chain:
    :return:
    :param executing_folder:
    """
    # receptor_with_waters_text = ""
    receptor_text = ""
    waters_text = ""
    ligand_text = ""
    ligand_filename = general_name + "_ligand.pdb"
    with open(pdb_filename, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                if line[21] == ligand_chain:
                    ligand_text += line
                else:
                    if line[17:20] == "HOH":
                        waters_text += line
                        # receptor_with_waters_text += line
                    else:
                        receptor_text += line
    if ligand_text == "":
        print "Something went wrong when extracting the ligand."
        return False
    elif receptor_text == "":
        print "Something went wrong when extracting the receptor."
        return False
    else:
        logging.info(" - Ligand and receptor extracted correctly.")
    with open(executing_folder + ligand_filename, 'w') as ligand_file:
        ligand_file.write(ligand_text)
    if waters_text:
        waters_filename = general_name + "_waters.pdb"
        with open(executing_folder + waters_filename, 'w') as waters_file:
            waters_file.write(waters_text)
        receptor_with_waters_filename = general_name + "_protein.pdb"
        with open(executing_folder + receptor_with_waters_filename, 'w') as receptor_with_waters_file:
            receptor_with_waters_file.write(receptor_text + waters_text)
        receptor_no_waters_filename = general_name + "_receptor_no_waters.pdb"
        with open(executing_folder + receptor_no_waters_filename, 'w') as receptor_file:
            receptor_file.write(receptor_text)
    else:
        receptor_no_waters_filename = general_name + "_protein.pdb"
        with open(executing_folder + receptor_no_waters_filename, 'w') as receptor_file:
            receptor_file.write(receptor_text)
        receptor_with_waters_filename = ""
        waters_filename = ""
    return receptor_no_waters_filename, receptor_with_waters_filename, waters_filename, ligand_filename


def command_execution_failed(command, error_message):
    print "The program failed to execute the command:\n{0}\nWith the following error:\n{1}\n" \
          "The program will be terminated.".format(command, error_message)
    logging.critical("CRITICAL: Couldn't execute the command:\n{}".format(command))
    logging.info("Terminating program.")
    logging.shutdown()
    sys.exit("Terminating the program.")


def structures_conversion_pdb_to_format(initial_filename, output_format, conversor_command, warnings_count,
                                        executing_folder, keyword, rewrite):
    if args.debug:
        print "wd:", executing_folder
    general_name = initial_filename.split(".pdb")[0]
    output_name = general_name + output_format
    if os.path.isfile(executing_folder + output_name):
        if not rewrite:
            logging.info(" - The file {0} already exists it won't be rewritten".format(output_name))
            return output_name, warnings_count
        else:
            logging.info(" - The file {0} already exists it will be rewritten".format(output_name))
        logging.info(" - The file already exists")
    command2execute = conversor_command.format(initial_filename, output_name)
    try:
        if args.debug:
            print command2execute
        conversor_output = check_output(command2execute.split(), cwd=executing_folder, stderr=STDOUT)
    except CalledProcessError as e:
        if first_time_execution[keyword]:
            command_execution_failed(command2execute, e)
        else:
            logging.warning(" - WARNING: Command {0} failed.".format(command2execute))
            warnings_count += 1
    else:
        first_time_execution[keyword] = False
        if output_format in [".mol2", ".sdf"]:
            if "0 molecules converted" in conversor_output:
                obabel_error = "Obabel hasn't converted any molecule from file {0}".format(initial_filename)
                logging.warning("Warning: Obabel hasn't converted any molecule from file {0}.".format(initial_filename))
                warnings_count += 1

    return output_name, warnings_count


def get_geometric_center_and_dimensions(filename, extra_space, executing_folder):
    coordinates_found = False
    min_x = None
    min_y = None
    min_z = None
    max_x = None
    max_y = None
    max_z = None
    with open(executing_folder + filename, 'r') as pdb_infile:
        for line in pdb_infile:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    sys.exit("Couldn't read the coordinates from the file {0}".format(filename))
                if min_x is None:
                    coordinates_found = True
                    min_x = x
                    max_x = x
                    min_y = y
                    max_y = y
                    min_z = z
                    max_z = z
                else:
                    if x < min_x:
                        min_x = x
                    elif x > max_x:
                        max_x = x
                    if y < min_y:
                        min_y = y
                    elif y > max_y:
                        max_y = y
                    if z < min_z:
                        min_z = z
                    elif z > max_z:
                        max_z = z
    if not coordinates_found:
        sys.exit("The file {0} doesn't have any coordinates".format(filename))

    # THe center is medium point in each axe.
    center = ((min_x + max_x) / 2,
              (min_y + max_y) / 2,
              (min_z + max_z) / 2,)
    dimensions = (max_x - min_x + extra_space,
                  max_y - min_y + extra_space,
                  max_z - min_z + extra_space)

    return center, dimensions


def vina_execution(ligand_pdbqt, receptor_pdbqt, vina_box_distance_to_ligand, working_folder, template_keywords,
                   template):
    """
    This function prepares the configuration file for computing Vina scoring function and executes it.
    :param ligand_pdbqt: The .pdbqt containing the ligand, this will be used as the ligand
                        for the scoring for computation and to compute the grid center and dimensions.
    :param receptor_pdbqt: The .pdbqt file containing the
    :param vina_box_distance_to_ligand: The size of the box (from the ligand to the limit) for vina calculations.
    :param working_folder: The folder where the config.txt file for vina should be created and where
                            Vina should be computed.
    :param template_keywords: The keywords present in the config.txt template file for vina computations.
    :param template: A Template object containing the config.txt template string for vina computations.
    :return:
    """
    if args.debug:
        print "executing vina"
    center, dimensions = get_geometric_center_and_dimensions(ligand_pdbqt, vina_box_distance_to_ligand, working_folder)
    vina_keywords_values = {}
    for keyword in template_keywords:
        if re.search(r"receptor", keyword):
            vina_keywords_values[keyword] = receptor_pdbqt
        elif re.search(r"ligand", keyword):
            vina_keywords_values[keyword] = ligand_pdbqt
        elif re.search(r"center.*", keyword):
            if "x" in keyword:
                vina_keywords_values[keyword] = "{0:.3f}".format(center[0])
            elif "y" in keyword:
                vina_keywords_values[keyword] = "{0:.3f}".format(center[1])
            elif "z" in keyword:
                vina_keywords_values[keyword] = "{0:.3f}".format(center[2])
            else:
                print "{0} is an invalid keyword for the vina template."
                print "If it should be accepted talk with the developer."
                sys.exit("Invalid vina template")
        elif re.search(r"size.*", keyword):
            if "x" in keyword:
                vina_keywords_values[keyword] = "{0:.3f}".format(dimensions[0])
            elif "y" in keyword:
                vina_keywords_values[keyword] = "{0:.3f}".format(dimensions[1])
            elif "z" in keyword:
                vina_keywords_values[keyword] = "{0:.3f}".format(dimensions[2])
            else:
                print "{0} is an invalid keyword for the vina template."
                print "If it should be accepted talk with the developer."
                sys.exit("Invalid vina template")
        else:
            print "{0} is an invalid keyword for the vina template."
            print "If it should be accepted talk with the developer."
            sys.exit("Invalid vina template")
    vina_configuration_file = "vina_config.txt"
    if template:
        config_text = template.substitute(vina_keywords_values)
        with open(working_folder + vina_configuration_file, 'w') as fileout:
            fileout.write(config_text)
    else:
        print "There's no vina template defined, and it should. Check it."
        logging.critical("CRITICAL: There's no vina template and it should exist.")
        logging.info("Terminating program.")
        logging.shutdown()
        sys.exit("Terminating the program.")
    if args.debug:
        print "wd: ", working_folder
    vina_command = "{0}vina --config {1} --out vina_out.txt " \
                   "--log vina_log.txt".format(external_software_paths.vina_executable_path, vina_configuration_file)
    try:
        vina_output = check_output(vina_command.split(), cwd=working_folder)
    except CalledProcessError as e:
        command_execution_failed(vina_command, e)
    else:
        if "Error" in vina_output:
            command_execution_failed(vina_command, vina_output)


def dsx_execution(receptor_filename, ligand_filename, execution_folder, waters_filename=""):
    """
    This function computes the dsx scoring funtion.
    :param receptor_filename: The .pdb file containing the protein
    :param ligand_filename: The .mol2 file containing the ligand molecule
    :param execution_folder: The folder where the scoring function has to be computed and the results writen
    :param waters_filename: The .mol2 file containing the waters present in the system if there are any.
    :return:
    """
    if args.debug:
        print "executing dsx"
    if waters_filename:
        dsx_command = "{0}linux64/dsx_linux_64.lnx -P {1} " \
                      "-L {2} -W {3} -D {0}pdb_pot_0511/".format(external_software_paths.dsx_path, receptor_filename,
                                                                 ligand_filename, waters_filename)
    else:
        dsx_command = "{0}linux64/dsx_linux_64.lnx -P {1} " \
                      "-L {2}  -D {0}pdb_pot_0511/".format(external_software_paths.dsx_path, receptor_filename,
                                                           ligand_filename)
    try:
        if args.debug:
            print "a", execution_folder
        check_output(dsx_command.split(), cwd=execution_folder, stderr=STDOUT)
    except CalledProcessError as e:
        command_execution_failed(dsx_command, e)


def xscore_execution(receptor_filename, ligand_filename, execution_folder):
    """
    This function computes the xscore scoring function, which has to be correctly installed previously.
    :param receptor_filename: the name of the pdb file containing the protein and waters (if present)
    :param ligand_filename:  the name of the mol2 file containing the ligand molecule
    :param execution_folder: the folder where the scoring function should be computed.
    :return:
    """
    if args.debug:
        print "executing xscore"
        print "receptor: '{0}'".format(receptor_filename)
        print "ligand: '{0}'".format(ligand_filename)
        print "wd: ", execution_folder
    xscore_command = "xscore -score {0} {1}".format(receptor_filename, ligand_filename)
    try:
        xscore_output = check_output(xscore_command.split(), cwd=execution_folder, stderr=STDOUT)
    except (CalledProcessError, OSError) as e:
        if first_time_execution["xscore"]:
            print execution_folder
            command_execution_failed(xscore_command, e)
        else:
            logging.info(" - Error while computing Xscore.")
    else:
        if "Error:" in xscore_output:
            if first_time_execution["xscore"]:
                command_execution_failed(xscore_command, xscore_output)
            else:
                logging.info(" - Error while computing Xscore.")
        else:
            output_filename = execution_folder + "xscore_out.txt"
            with open(output_filename, "w") as outfile:
                outfile.write(xscore_output)
        first_time_execution["xscore"] = False



def binana_execution(receptor_filename, ligand_filename, execution_folder):
    if args.debug:
        print "Executing binana"
    binana_command = "python {0}binana_1_2_0.py -receptor {1} -ligand {2}".format(external_software_paths.binana_path,
                                                                           receptor_filename, ligand_filename)
    if args.debug:
        print "comand:", binana_command
        print "wd:", execution_folder
    try:
        binana_output = check_output(binana_command.split(), cwd=execution_folder, stderr=STDOUT)
    except CalledProcessError as e:
        command_execution_failed(binana_command, e)
    else:
        binana_output_filename = execution_folder + "binana_output.txt"
        with open(binana_output_filename, 'w') as binana_outfile:
            binana_outfile.write(binana_output)


def mmgbsa_execution(complex_filename, execution_folder, host, number_cpus, lig_chain):
    if args.debug:
        print "Executing mmgbsa"
    mmgbsa_command = "{0}prime_mmgbsa;{1};-csv;yes;-ligand;'chain.name  {4} ';" \
                     "-HOST;{2}:{3};-WAIT".format(external_software_paths.schrodinger_path,
                                                   complex_filename, host, number_cpus, lig_chain)
    if args.debug:
        print mmgbsa_command.split(';')
    try:
        check_call(mmgbsa_command.split(';'), cwd=execution_folder, timeout=7200)
    except CalledProcessError as e:
        if first_time_execution['mmgbsa']:
            command_execution_failed(mmgbsa_command.replace(";", " "), e)
        else:
            logging.info(" - Error while computing mmgbsa.")
    except TimeoutExpired as e:
        logging.info(" - mmgbsa execution takes more than 2 hours for this system, it won't be computed.")


parser = argparse.ArgumentParser()
parser.add_argument("-input_files", nargs="+", required=True)
parser.add_argument("-folders_path", default="./scoring_functions_values")
#parser.add_argument("-output_general_name", default="all")
parser.add_argument("-scoring_functions", nargs="+",
                    default=["glide", "vina", "xscore", "dsx", "mmgbsa", "binana", "rf_score"])
parser.add_argument("-schrodinger_host", default="Calculon_slurm")
parser.add_argument("-schrodinger_cpus", default="1")
parser.add_argument("-vina_box_distance_to_ligand", default=20, type=int)
parser.add_argument("-experimental_deltag", default="")
parser.add_argument("-rf_score_output_file", default=-2)
parser.add_argument("-debug", default=False, action="store_true")
parser.add_argument("-ligand_chain", default="Z")
parser.add_argument("-log_file", default="scoring_function_computations_log.txt")
parser.add_argument("-rewrite", action="store_true")
args = parser.parse_args()

logging.basicConfig(filename=args.log_file, format="%(message)s", level=logging.INFO, filemode="w")
logging.info("{} : Program starting".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M")))
warnings_counter = 0
old_warnings = 0
errors_counter = 0
# Some constants definitions and environmental commands needed to use other applications. They'll
# change depending on the environment where the program is launched (BSC, AZ, etc.)
# AZ
# modules.module("load", "schrodinger/2016.01")
pdb_convert_path = external_software_paths.schrodinger_path
schrodinger_pdb2mae_convert_command = pdb_convert_path + "utilities/pdbconvert -ipdb {} -omae {}"
# modules.module("load", "openbabel/python2.7-2.3.1")
obabel_command = "obabel {} -O {}"
prepare_ligand_path = "{0}bin/pythonsh {1}prepare_ligand4.py ".format(external_software_paths.vina_tools,
                                                                      external_software_paths.vina_utilities)
prepare_receptor_path = "{0}bin/pythonsh {1}prepare_receptor4.py ".format(external_software_paths.vina_tools,
                                                                          external_software_paths.vina_utilities)
prepare_ligand_vina_command = prepare_ligand_path + "-l {0} -o {1}"
prepare_receptor_vina_command = prepare_receptor_path + "-U nphs -U lps -U nonstdres -r {0} -o {1}"
rf_descriptors_extractor_path = "{0}RF-descriptors_extraction ".format(
    external_software_paths.rf_score_descriptors_extractor_path)
rf_descriptors_commands = rf_descriptors_extractor_path + "--input_from_file {0} --output_file {1}"




# Check for the folder separator in the new_folders path, add it at the end if it isn't.

compute_autodock_vina = False
compute_xscore = False
compute_dsx = False
compute_rf_score = False
compute_mmgbsa = False
compute_binana = False
compute_obabel_desc = False
glides2compute = []
# Parse the scoring functions to compute and check that all of them are supported.
first_time_execution = {}
for sf in args.scoring_functions:
    if re.search(r"(?:autodock)*_*-*vina", sf, re.IGNORECASE):
        compute_autodock_vina = True
        first_time_execution['vina'] = True
        continue
    elif re.search(r"xscore", sf, re.IGNORECASE):
        compute_xscore = True
        first_time_execution['xscore'] = True
    elif re.search(r"dsx", sf, re.IGNORECASE):
        compute_dsx = True
    elif re.search(r"rf[_-]*score", sf, re.IGNORECASE):
        compute_rf_score = True
        first_time_execution['rf_score'] = True
    elif re.search(r"mm[_-]*gbsa", sf, re.IGNORECASE):
        compute_mmgbsa = True
        first_time_execution['mmgbsa'] = True
    elif re.search(r"binana", sf, re.IGNORECASE):
        compute_binana = True
    elif re.search(r"glide.*", sf, re.IGNORECASE):
        a = re.search(r"glide-*_*(.P)*|glide-*_*(HTS)", sf, re.IGNORECASE)
        if a.group().lower() == "glide":
            glides2compute.append("glide_SP")
            first_time_execution['glide_SP'] = True
        else:
            glides2compute.append("glide_"+a.group(1).upper())
            first_time_execution["glide_"+a.group(1).upper()] = True
    elif re.search(r"obabel[_-]*desc(?:riptors)*", sf, re.IGNORECASE):
        compute_obabel_desc = True
first_time_execution["obabel"] = True
first_time_execution["pdbqt"] = True
first_time_execution["mae"] = True
# If glide is one of the scoring functions to compute create the folders needed for glide computation
# into the general subfolder, and write the input files for glide.
if glides2compute:
    if args.folders_path[-1] != os.sep:
        new_folders_path = args.folders_path + os.sep
    else:
        new_folders_path = args.folders_path
    # Check whether the subfolders path exists or not. If it doesn't exist create it
    if not os.path.exists(new_folders_path):
        os.mkdir(new_folders_path)
    # all the glide versions can be computed on the same folder. They need slightly different input files though.
    glide_subfolder = new_folders_path + "glide_calculations" + os.sep
    try:
        os.mkdir(glide_subfolder)
    except OSError as e:
        if "File exists" not in e[1]:
            sys.exit(e)
    if args.debug:
        print 'bla'
    glide_structures_subfolder = glide_subfolder + "structures" + os.sep
    if args.debug:
        print glide_structures_subfolder
    try:
        os.mkdir(glide_structures_subfolder)
    except OSError as e:
        if "File exists" not in e[1]:
            sys.exit(e)
    glide_template_filename = external_software_paths.templates_directory + "glide_mininplace_template.in"
    with open(glide_template_filename, 'r') as infile:
        glide_template_text = infile.read()
    glides_template = Template(glide_template_text)
    glide_input_files = []
    for glide_model in glides2compute:
        glide_precision = glide_model.split('_')[-1]
        glide_input_text = glides_template.substitute(precision=glide_precision)
        glide_working_directory = glide_subfolder + "glide{}".format(glide_precision) + os.sep
        try:
            os.mkdir(glide_working_directory)
        except OSError as e:
            if "File exists" not in e[1]:
                sys.exit(e)
        glide_input_filename = glide_working_directory + "glide_{}_mininplace.in".format(glide_precision)
        glide_input_files.append(glide_input_filename)
        with open(glide_input_filename, 'w') as outfile:
            outfile.write(glide_input_text)
else:
    glide_subfolder = ""
keywords_in_the_vina_template = {}
vina_template = ""
if compute_autodock_vina:
    vina_template_filename = external_software_paths.templates_directory + "vina_score_template.txt"
    with open(vina_template_filename, 'r') as infile:
        vina_template_text = infile.read()
    keywords_in_the_vina_template = set(pattern.group(1) for pattern in re.finditer("\$\{*(\w*_*w*)\}?",
                                                                                    vina_template_text))
    vina_template = Template(vina_template_text)

deltag_values = {}
if compute_rf_score:
    logging.info("Remember that for RF-Score you still have to do the computation, this just prepares the input.")
    rf_score_name_deltag_dictio = {}
    if args.experimental_deltag:
        with open(args.experimental_deltag, 'r') as infile:
            first_line = True
            for line in infile:
                line = line.strip()
                if "," in line:
                    line = line.split(",")
                elif ";" in line:
                    line = line.split(";")
                else:
                    line = line.split()
                if first_line:
                    first_line = False
                    try:
                        int(line[1])
                    except ValueError:
                        sim_index = line.index("sim_id")
                        activity_index = line.index("activity")
                        continue
                    else:
                        sim_index = 1
                        activity_index = 2
                deltag_values[line[sim_index].strip()] = line[activity_index].strip()
    else:
        logging.info("INFO: Note that no energies file has been given, since RF-Score needs energy values as "
                     "input the program will assign a dummy energy of 0 to all the systems")

minimum_filename = ""
for filename in args.input_files:
    logging.info("Working with the file: {0}".format(filename))
    if ".pdb" not in filename:
        logging.error(" - ERROR: This program expects to work with .pdb files since that is the output format of PELE.")
        logging.error("- INFO: Skipping file {0}".format(filename))
        errors_counter += 1
    basic_filename = filename.split(os.sep)[-1]
    minimum_filename = basic_filename.split('.pdb')[0]
    pattern = re.search(r'([a-z0-9]+_\d+)_complex_processed', minimum_filename, re.IGNORECASE)
    if pattern:
        minimum_filename = pattern.group(1)
    # Since the structures have been extracted from the trajectory previously and the folder structure has been already
    # created we just need to extract all the files into the same folder as the input file.

    working_folder = os.sep.join(filename.split(os.sep)[:-1]) + os.sep
    if working_folder == os.sep:
        working_folder = '.' + os.sep
    if glide_subfolder:
        # noinspection PyUnboundLocalVariable
        new_link = glide_structures_subfolder + basic_filename
        try:
            os.symlink(os.path.abspath(filename), new_link)
        except OSError as e:
            if "File exists" not in e[1]:
                sys.exit(e)
        else:
            logging.info(" - linking the structure {0}".format(new_link))
    logging.info(" - Separating ligand and protein")
    try:
        protein_filename_pdb, protein_with_waters_pdb, \
        waters_filename_pdb, ligand_filename_pdb = extract_ligand(filename, minimum_filename, args.ligand_chain, working_folder)
    except TypeError:
        logging.error(" - This system will be discontinued, due to a problem in the "
                      "separation of the ligand and protein.")
        errors_counter += 1
        continue
    if waters_filename_pdb == "":
        waters_mol2_filename = ""
    else:
        if ".pdb" not in waters_filename_pdb:
            logging.critical("CRITICAL: There's been a programming error talk with the developer.")
            logging.info("ERROR CODE: CSF202, give this to the developer.")
            logging.shutdown()
            sys.exit("PROGRAMMING ERROR: CSF202")
        waters_mol2_filename = ""
        if compute_dsx:
            logging.info(" - Converting waters to .mol2")
            waters_mol2_filename, warnings_counter = structures_conversion_pdb_to_format(waters_filename_pdb, ".mol2",
                                                                                         obabel_command,
                                                                                         warnings_counter,
                                                                                         working_folder, "mol2",
                                                                                         args.rewrite)
        if warnings_counter > old_warnings:
            old_warnings = warnings_counter
    if ".pdb" not in ligand_filename_pdb:
        logging.critical("CRITICAL: There's been a programming error talk with the developer.")
        logging.info("ERROR CODE: CSF215, give this to the developer.")
        logging.shutdown()
        sys.exit("PROGRAMMING ERROR: CSF215")
    if ".pdb" not in protein_filename_pdb:
        logging.critical("CRITICAL: There's been a programming error talk with the developer.")
        logging.info("ERROR CODE: CSF220, give this to the developer.")
        logging.shutdown()
        sys.exit("PROGRAMMING ERROR: CSF220")
    ligand_general_filename = ligand_filename_pdb.split(".pdb")[0]
    if compute_dsx and compute_xscore:
        if args.debug:
            print 1
        logging.info(" - Converting ligand to .mol2")
        ligand_filename_mol2, warnings_counter = structures_conversion_pdb_to_format(ligand_filename_pdb, ".mol2",
                                                                                     obabel_command, warnings_counter,
                                                                                     working_folder, "obabel",
                                                                                     args.rewrite)
        if warnings_counter > old_warnings:
            old_warnings = warnings_counter
        else:
            logging.info(" - Executing DSX")
            dsx_execution(protein_filename_pdb, ligand_filename_mol2, working_folder, waters_mol2_filename)
            logging.info(" - Executing Xscore")
            if waters_filename_pdb:
                xscore_execution(protein_with_waters_pdb, ligand_filename_mol2, working_folder)
            else:
                xscore_execution(protein_filename_pdb, ligand_filename_mol2, working_folder)
    elif compute_dsx:
        if args.debug:
            print 2
        logging.info(" - Converting ligand to .mol2")
        ligand_filename_mol2, warnings_counter = structures_conversion_pdb_to_format(ligand_filename_pdb, ".mol2",
                                                                                     obabel_command, warnings_counter,
                                                                                     working_folder, "obabel",
                                                                                     args.rewrite)
        if warnings_counter > old_warnings:
            old_warnings = warnings_counter
        else:
            logging.info(" - Executing DSX")
            dsx_execution(protein_filename_pdb, ligand_filename_mol2, working_folder, waters_mol2_filename)
    elif compute_xscore:
        if args.debug:
            print 3
        ligand_filename_mol2, warnings_counter = structures_conversion_pdb_to_format(ligand_filename_pdb, ".mol2",
                                                                                     obabel_command, warnings_counter,
                                                                                     working_folder, "obabel",
                                                                                     args.rewrite)
        if warnings_counter > old_warnings:
            old_warnings = warnings_counter
        else:
            logging.info(" - Executing Xscore")
            if waters_filename_pdb:
                xscore_execution(protein_with_waters_pdb, ligand_filename_mol2, working_folder)
            else:
                xscore_execution(protein_filename_pdb, ligand_filename_mol2, working_folder)
    if compute_rf_score:
        if args.debug:
            print 4
        logging.info(" - Converting ligand to .sdf")
        ligand_filename_sdf, warnings_counter = structures_conversion_pdb_to_format(ligand_filename_pdb, ".sdf",
                                                                                    obabel_command, warnings_counter,
                                                                                    working_folder, "obabel",
                                                                                    args.rewrite)
        if args.debug:
            print 9
        # The real computation of RF-score is done manually.
        logging.info(" - Preparing RF deltaG")
        if deltag_values:
            # noinspection PyUnboundLocalVariable
            if minimum_filename in deltag_values.keys():
                if args.debug:
                    print 10
                rf_score_name_deltag_dictio[minimum_filename] = deltag_values[minimum_filename]
            else:
                if args.debug:
                    print 11
                pattern = re.search(r"\w+_(\d+)_", minimum_filename)
                if pattern is None:
                    logging.warning(" - WARNING: The name couldn't be found in the deltaG file, "
                                    "a dummy deltaG of 0 will be used for the system.")
                    rf_score_name_deltag_dictio[minimum_filename] = 0
                else:
                    try:
                        rf_score_name_deltag_dictio[minimum_filename] = deltag_values[pattern.group(1)]
                    except KeyError:
                       logging.warning(" - WARNING: The name couldn't be found in the deltaG file, "
                                    "a dummy deltaG of 0 will be used for the system.")
                       rf_score_name_deltag_dictio[minimum_filename] = 0

        else:
            if args.debug:
                print 12
            rf_score_name_deltag_dictio[minimum_filename] = 0
    if compute_autodock_vina and compute_binana:
        if args.debug:
            print 5
        logging.info(" - Converting ligand to .pdbqt")
        ligand_filename_pdbqt, warnings_counter = structures_conversion_pdb_to_format(ligand_filename_pdb, ".pdbqt",
                                                                                      prepare_ligand_vina_command,
                                                                                      warnings_counter, working_folder,
                                                                                      "pdbqt", args.rewrite)
        if warnings_counter > old_warnings:
            old_warnings = warnings_counter
        else:
            if args.debug:
                print 'receptor_prep'
            logging.info(" - Converting receptor to .pdbqt")
            receptor_filename_pdbqt, warnings_counter = structures_conversion_pdb_to_format(protein_filename_pdb,
                                                                                            ".pdbqt",
                                                                                            prepare_receptor_vina_command,
                                                                                            warnings_counter,
                                                                                            working_folder, "pdbqt",
                                                                                            args.rewrite)
            if warnings_counter > old_warnings:
                old_warnings = warnings_counter
            else:
                logging.info(" - Executing Vina")
                vina_execution(ligand_filename_pdbqt, receptor_filename_pdbqt, args.vina_box_distance_to_ligand,
                               working_folder, keywords_in_the_vina_template, vina_template)
                logging.info(" - Executing Binana")
                if args.debug:
                    print receptor_filename_pdbqt
                    print ligand_filename_pdbqt
                    print working_folder
                binana_execution(receptor_filename_pdbqt, ligand_filename_pdbqt, working_folder)
    elif compute_autodock_vina:
        if args.debug:
            print 6
        logging.info(" - Converting ligand to .pdbqt")
        ligand_filename_pdbqt, warnings_counter = structures_conversion_pdb_to_format(ligand_filename_pdb, ".pdbqt",
                                                                                      prepare_ligand_vina_command,
                                                                                      warnings_counter, working_folder,
                                                                                      "pdbqt", args.rewrite)
        if warnings_counter > old_warnings:
            old_warnings = warnings_counter
        else:
            logging.info(" - Converting receptor to .pdbqt")
            receptor_filename_pdbqt, warnings_counter = structures_conversion_pdb_to_format(protein_filename_pdb,
                                                                                            ".pdbqt",
                                                                                            prepare_receptor_vina_command,
                                                                                            warnings_counter,
                                                                                            working_folder, "pdbqt",
                                                                                            args.rewrite)
            if warnings_counter > old_warnings:
                old_warnings = warnings_counter
            else:
                logging.info(" - Executing Vina")
                vina_execution(ligand_filename_pdbqt, receptor_filename_pdbqt, args.vina_box_distance_to_ligand,
                               working_folder, keywords_in_the_vina_template, vina_template)
    elif compute_binana:
        if args.debug:
            print 7
        logging.info(" - Converting ligand to .pdbqt")
        ligand_filename_pdbqt, warnings_counter = structures_conversion_pdb_to_format(ligand_filename_pdb, ".pdbqt",
                                                                                      prepare_ligand_vina_command,
                                                                                      warnings_counter, working_folder,
                                                                                      "pdbqt", args.rewrite)
        if warnings_counter > old_warnings:
            old_warnings = warnings_counter
        else:
            logging.info(" - Converting receptor to .pdbqt")
            receptor_filename_pdbqt, warnings_counter = structures_conversion_pdb_to_format(protein_filename_pdb,
                                                                                            ".pdbqt",
                                                                                            prepare_receptor_vina_command,
                                                                                            warnings_counter,
                                                                                            working_folder, "pdbqt",
                                                                                            args.rewrite)
            if warnings_counter > old_warnings:
                old_warnings = warnings_counter
            else:
                logging.info(" - Executing Binana")
                binana_execution(receptor_filename_pdbqt, ligand_filename_pdbqt, working_folder)

    if compute_mmgbsa:
        if args.debug:
            print 8
        logging.info(" - Converting complex to .mae")
        complex_filename_mae, warnings_counter = structures_conversion_pdb_to_format(basic_filename, ".mae",
                                                                                     schrodinger_pdb2mae_convert_command,
                                                                                     warnings_counter, working_folder,
                                                                                     "mae", args.rewrite)
        if warnings_counter > old_warnings:
            old_warnings = warnings_counter
        else:
            logging.info(" - Executing mmgbsa")
            mmgbsa_execution(complex_filename_mae, working_folder, args.schrodinger_host, args.schrodinger_cpus,
                             args.ligand_chain)

if errors_counter == len(args.input_files):
    logging.critical("None of the systems could be processed.\nTerminating the program.")
    logging.shutdown()
    sys.exit("Couldn't process any of the inputs.")

if not minimum_filename:
    logging.critical("CRITICAL: There's no filename how did it get to here??")
    logging.shutdown()
if glides2compute:
    if args.debug:
        print 13
    logging.info(" - Executing glide")
    # noinspection PyUnboundLocalVariable
    for input_file in glide_input_files:
        infile = input_file.split(os.sep)[-1]
        glide_working_directory = os.sep.join(input_file.split(os.sep)[:-1])
        glide_command = "{0}run xglide.py {1} -HOST {2}:{3}".format(external_software_paths.schrodinger_path,
                                                                    infile, args.schrodinger_host,
                                                                    args.schrodinger_cpus)
        if args.debug:
            print glide_working_directory
        try:
            check_call(glide_command.split(), cwd=glide_working_directory)
        except CalledProcessError as e:
            command_execution_failed(glide_command, e)
if compute_rf_score:
    if args.debug:
        print 14
    if args.rf_score_output_file == -2:
        rf_output_folder = "/".join(minimum_filename.split("/")[:-1])
    else:
        rf_output_folder = args.rf_score_output_file
    rf_output_filename = rf_output_folder + "rf_descriptors_input_file.txt"
    rf_descriptors_output_filename = rf_output_folder + "rf_descriptors_output_file.csv"
    logging.info("Writing the initial file for RF-Score: {0}".format(rf_output_filename))
    with open(rf_output_filename, 'w') as outfile:
        outfile.write("\n".join(["{0},{1}".format(key, value) for key, value in
                                 rf_score_name_deltag_dictio.iteritems()]))
    rf_command = rf_descriptors_commands.format(rf_output_filename, rf_descriptors_output_filename)
    try:
        rf_output = check_output(rf_command.split())
    except (CalledProcessError, OSError )as e:
        command_execution_failed(rf_command, e)
    else:
        if "Filename" not in rf_output:
            command_execution_failed(rf_command, e)
    os.remove(rf_output_filename)
    # logging.info("Remember yo have to execute RF yourself!")
    # print "Remember yo have to execute RF yourself!"
logging.info("{} : Program finished normally with {} warnings and {} non-critical errors.".format(
    datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), warnings_counter, errors_counter))
logging.shutdown()
