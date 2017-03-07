"""
A script to extract the scoring functions (SFs) values from the output generated from each SF. and generate
the following series of files:
- A file containing the SF value and the experimental energy, for each SF.
- A file containing the values of all the SFs and their respective experimental energies, in this
case only those systems for which all the SFs have been correctly computed will be taken into account.
- A file containing all the descriptors and the experimental energies.
- A file containing only the structural descriptors and their experimental energies.
- A file containing only the energy descriptors and their experimental energies.
The experimental energies will be added only if a file containing the experimental energies is provided.
If present the SFs xscore and nn-score will be converted to deltaG values since they produce pkd values

DEVELOPERS NOTES:
- The program some times wirtes out things like CODE:LXXX where XXX is a number, that number
corresponds to the line of code (when it was written) basically it's an easy way to find the cause of error,
since they're unique codes.
Jelisa Iglesias 24/11/2016
mail: jelisa.iglesias@gmail.com
"""

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import os
import logging
import datetime
import re
from math import log
import sys
import csv
import pybel as py

import extract_sf_help


def binana_extraction(filename, ensemble_flag):
    system_id = filename.split('/')[-1]
    filein = open(filename, 'r')
    file_lines = filein.readlines()
    filein.close()
    values_to_extract = {"2.5_contacts": 0, "4.0_contacts": 0, "hydrophobic": 0, "pi-cation": 0,
                         "salt_bridge": 0, "pi-t": 0, "pi-pi": 0, "sidechain_alpha": 0,
                         "sidechain_beta": 0, "sidechain_other": 0, "backbone_alpha": 0,
                         "backbone_beta": 0, "backbone_other": 0, "H_B_receptor": 0, "H_B_ligand": 0}
    read_table = False
    skip = 0
    key = None
    # flexible_dictionary_absolute_values = {"sidechain_alpha": 0.0, "sidechain_beta": 0.0, "sidechain_other": 0.0,
    #                                        "backbone_alpha": 0.0, "backbone_beta": 0.0, "backbone_other": 0.0}
    h_bonds_dictionary = {"receptor": 0, "ligand": 0}
    total_count = 0
    if ensemble_flag:
        name_pattern = r"[_{0}]*(\w+_\d+)_.*ligand".format(os.sep)
    else:
        name_pattern = r"[_{0}]*([a-z0-9]+_\d+)_.*ligand".format(os.sep)
    ligand_file_pattern = r"\s*ligand\s*\|\s+(\w+)\."
    for line in file_lines:
        if line.strip() == "":
            if read_table:
                read_table = False
                key = None
                count, total_count = 0, 0
            continue
        pattern = re.search(ligand_file_pattern, line, re.IGNORECASE)
        if pattern:
            name = pattern.group(1)
            pattern = re.search(name_pattern, name, re.IGNORECASE)
            if pattern is None:
                logging.warning(" # WARNING: The file {0} doesn't match the naming convention.".format(filename))
                return False
            system_id = pattern.group(1)
            continue
        if skip > 0:
            skip -= 1
            continue
        if "Atom-type" in line:
            if "2.5" in line:
                key = "2.5_contacts"
            elif "4.0" in line:
                key = "4.0_contacts"
            else:
                print "Where are you?? How did you get HERE??"
                continue
        elif " flexibility:" in line:
            key = "flexibility"
        elif "Hydrogen bonds" in line:
            key = "h_bonds"
        elif "Hydrophobic contacts" in line:
            key = "hydrophobic"
        elif "pi-pi stacking interactions:" in line:
            # print 'here', line
            key = "pi-pi"
        elif "T-stacking (" in line:
            key = "pi-t"
        elif "Cation-pi" in line:
            key = "pi-cation"
        elif "Salt Bridges:" in line:
            key = "salt_bridge"
        elif line.strip() == "Raw data:":
            # print "here", read_table
            if read_table:
                read_table = False
                key = None
                total_count = 0
            continue
        if not read_table and key is not None:
            read_table = True
            skip = 2
            continue
        elif read_table and key is None:
            print "WTH?? How did you get here... Let's recheck things!"
            continue
        if read_table:
            line = line.strip("REMARK").strip()
            if key in ["2.5_contacts", "4.0_contacts", "hydrophobic", "pi-cation"]:
                at_tp_1, at_tp_2, count = [element.strip() for element in line.split('|')]
                try:
                    # noinspection PyUnboundLocalVariable
                    total_count += int(count)
                except NameError:
                    print "WTH?"
                    total_count = int(count)
                values_to_extract[key] = str(total_count)
            elif key == "flexibility":
                # print line
                sd_bb, structure_3d, count = [element.lower().strip() for element in line.split('|')]
                flex_key = "_".join([sd_bb, structure_3d])
                values_to_extract[flex_key] = count
            elif key == "h_bonds":
                hbond_dictio_key, sd_bb, structure_3d, count = [element.lower().strip() for element in line.split('|')]
                h_bonds_dictionary[hbond_dictio_key] += int(count)
            elif key in ["salt_bridge", "pi-t", "pi-pi"]:
                sc_struct, count = [element.strip() for element in line.split('|')]
                total_count += int(count)
                values_to_extract[key] = str(total_count)
                # print key,  str(total_count)
    values_to_extract['H_B_ligand'] = str(h_bonds_dictionary['ligand'])
    values_to_extract['H_B_ligand'] = str(h_bonds_dictionary['receptor'])
    # values_to_extract["flexibility"] = flexible_dictionary_absolute_values
    # for key, value in values_to_extract.iteritems():
    #     values_to_extract[key] = str(key)

    # values_to_extract["h_bonds"] = h_bonds_dictionary
    return system_id, values_to_extract


def mmgbsa_extraction(filename, ensemble_flag):
    # Old version that also extracts the total energy and the values discarded
    # for the consensus as descriptors.
    # keywords2look4 = ['r_psp_MMGBSA_dG_Bind', 'r_psp_MMGBSA_dG_Bind_Coulomb', 'r_psp_MMGBSA_dG_Bind_Covalent',
    #                   'r_psp_MMGBSA_dG_Bind_Hbond', 'r_psp_MMGBSA_dG_Bind_Lipo', 'r_psp_MMGBSA_dG_Bind_Packing',
    #                   'r_psp_MMGBSA_dG_Bind_SelfCont', 'r_psp_MMGBSA_dG_Bind_Solv_GB', 'r_psp_MMGBSA_dG_Bind_vdW',
    #                   'r_psp_Lig_Strain_Energy']
    keywords2look4 = ['r_psp_MMGBSA_dG_Bind', 'r_psp_MMGBSA_dG_Bind_Coulomb', 'r_psp_MMGBSA_dG_Bind_Covalent',
                      'r_psp_MMGBSA_dG_Bind_Hbond', 'r_psp_MMGBSA_dG_Bind_Lipo', 'r_psp_MMGBSA_dG_Bind_Packing',
                      'r_psp_MMGBSA_dG_Bind_Solv_GB', 'r_psp_MMGBSA_dG_Bind_vdW', 'r_psp_Lig_Strain_Energy']
    if ensemble_flag:
        name_pattern = r"[_{0}]*(\w+_\d+)_.*".format(os.sep)
    else:
        name_pattern = r"[_{0}]*(\w+_\d+)_.*".format(os.sep)
    complete_name = filename.split(os.sep)[-1].split("-out.csv")[0]
    pattern = re.search(name_pattern, complete_name)
    if pattern:
        system_id = pattern.group(1)
    else:
        system_id = complete_name
    all_descriptors = {}
    with open(filename) as infile:
        file_text = infile.readlines()
        if len(file_text) == 1:
            if file_text[0].strip() == "":
                logging.error(" ERROR: The file is empty!")
            else:
                logging.error(" ERROR: This file has only the header but no values.")
                return False
        elif len(file_text) != 2:
            logging.error(" ERROR: This file has more than two lines, it shouldn't happen...")
            return False
        else:
            skip_this_file = [line for line in file_text if len(line.split(",")) <= 1]
            if skip_this_file:
                logging.error(" ERROR:This file isn't separated with comas, it isn't valid!")
                return False
            all_descriptors = {key: "{0:.5f}".format(float(value)) for key, value in zip(file_text[0].split(','),
                                                                                         file_text[1].split(','))
                               if value != "" and key != "title"}
    # return keywords2look4, all_descriptors
    return system_id, {key[6:]: all_descriptors[key] for key in keywords2look4[1:]}, all_descriptors[keywords2look4[0]]


def dsx_extraction(filename, ensemble_flag):
    """
    A function to extract the DSX value from an output file of this scoring function were only on ligand was treated.
    :param ensemble_flag:
    :param filename: a string containing the name of the file to parse.
    :return: Two strings on contains the id of the system and the other the scoring function value
    """
    with open(filename, 'r') as infile:
        score_txt = infile.read()
    score_pattern = r"\d+\s+\|\s+\w+.*\s+\|.*\|\s+(.*)\s+\|.*\|.*\|.*\|.*\|"
    pattern = re.search(score_pattern, score_txt, re.IGNORECASE)
    if pattern is None:
        # print 'c'
        logging.warning(" # WARNING: The file {0} doesn't have the needed values check it.".format(filename))
        return False
    score = pattern.group(1).strip()
    # THis pattern will have to change to adapt to different names. another option for files
    # where the filename has a path is commented and provided
    # ligand_filename_pattern = r"ligands file\s+:\s+.*{0}+(\w+)\.(\w+)".format(os.sep)
    ligand_filename_pattern = r"ligands file\s+:\s+(\w+)\.(\w+)"
    pattern = re.search(ligand_filename_pattern, score_txt, re.IGNORECASE)
    if pattern is None:
        logging.warning(" # WARNING: The file {0} doesn't have the needed values check the file, if it's OK "
                        "tell the developer to review the patterns. CODE:L187".format(filename))
        # print 'a'
        return False
    name = pattern.group(1)
    if ensemble_flag:
        pattern2find_in_name = r"[_{0}]*(\w+_\d+)_.*l{,1}i{,1}g{,1}a{,1}n{,1}d{,1}"
    else:
        pattern2find_in_name = r"[_{0}]*([a-z0-9]+_\d+)_\d+.*l{,1}i{,1}g{,1}a{,1}n{,1}d{,1}"
        # If there's any problem with the naming just change this patter. But
    # it's quite generic for the DSX score, since we use the suffix _ligand for all the process.
    pattern = re.search(pattern2find_in_name, name)
    if pattern is None:
        logging.warning(" - WARNING: The pattern {0} doesn't match the file name {1}, "
                        "check both, CODE:L193".format(pattern2find_in_name, name))
        logging.info(" - The program will use the whole name as key for this scoring functions")
        # print 'b'
        system_id = name
    else:
        system_id = pattern.group(1)

    return system_id, score


def xscore_extraction(filename, rt_factor, ensemble_flag):
    if ensemble_flag:
        id_pattern = r".* ligand from '[_{0}]*(\w+_\d+).*_ligand.mol2'".format(os.sep)
    else:
        id_pattern = r".* ligand from '[_{0}]*([a-z0-9]+_\d+).*_ligand.mol2'".format(os.sep)
    scores_patterns = r"(H[M,P,S]SCORE) .* = (-*\d.\d*)|(average) .* = (-*\d.\d*)"
    with open(filename, "r") as infile:
        text = infile.read()
    pattern = re.search(id_pattern, text, re.IGNORECASE)
    if pattern is None:
        print 1
        logging.warning(" # WARNING: Something went really wrong when computing Xscore for this system, It'll be skip.")
        return False
    system_id = pattern.group(1)
    values = {}
    for finding in re.finditer(scores_patterns, text, re.IGNORECASE):
        if "average" in finding.group():
            values[finding.group(3)] = "{0:0.5f}".format(rt_factor * (log(10 ** (-float(finding.group(4))))))
        else:
            values[finding.group(1)] = "{0:0.5f}".format(rt_factor * (log(10 ** (-float(finding.group(2))))))
    if not values:
        print 2
        logging.warning(" # WARNING: Something went wrong when computing Xscore for this system, It'll be skip.")
        return False
    return system_id, values


def vina_extraction(filename, ensemble_flag):
    pattern2look4 = r"Affinity:\s*(-\d+\.\d+)"
    system_id = filename.split(os.sep)[-2]
    if not ensemble_flag:
        name_pattern = r"([a-z0-9]+_\d+).*"
        pattern = re.search(name_pattern, system_id, re.IGNORECASE)
        if pattern:
            system_id = pattern.group(1)
    with open(filename, 'r') as infile:
        score_txt = infile.read()
    pattern = re.search(pattern2look4, score_txt)
    if pattern is None:
        logging.warning(" # WARNING: Something went wrong when computing vina for this system, It'll be skip.")
        return False
    score = pattern.group(1)
    return system_id, score


def parse_csv_file(filename, ensemble_flag):
    """
    A function meant to parse csv files where we only want two columns,
     the first as ID and the second as value.
    :param ensemble_flag:
    :param filename: the name of the file to open
    :return: a dictionary containing the first column as keys and the second one as values.
    """
    with open(filename) as infile:
        csv_parser = csv.reader(infile, delimiter=",")
        csv_parser.next()  # With this line we get to skip the header
        dictio2return = {}
        for line in csv_parser:
            if ensemble_flag:
                name_pattern = r"(\w+_\d+).*"
            else:
                name_pattern = r"([a-z0-9]+_\d+).*"
            pattern = re.search(name_pattern, line[0])
            if pattern:
                system_id = pattern.group(1)
            else:
                system_id = line[0].split(os.sep)[-1]
            dictio2return[system_id] = line[1].strip()
    return dictio2return


def process_values(dictionary2analyze, keyword2dictionary, possible_name):
    try:
        dictionary2analyze[keyword2dictionary]
    except KeyError:
        print "dictionary:", dictionary2analyze
        print "dict_keyword:", keyword2dictionary
        print 'scoring function/descriptor analyzed:', possible_name
        sys.exit("OK, this shouldn't ever happen. There has been an error with the data you've provided. Review it.")
        # return False
    else:
        elements = dictionary2analyze[keyword2dictionary]
    # print elements
    keys2return = ""
    values2return = ""
    if type(elements) is dict:
        # print 'a'
        for key, value in sorted(elements.iteritems()):
            if type(value) is dict:
                values2return = ",".join([str(element) for element in value.values()])
                keys2return = ",".join([str(element) for element in value.keys()])
            else:
                values2return += "," + str(value)
                keys2return += "," + str(key)
    else:
        # print 'a2'
        values2return = "{0},".format(elements)
        keys2return = str(possible_name)
    return keys2return, values2return


def merge_data(data2use, type_of_data, systems2use, energies_dictio, output_name):
    not_concatenate_values = True
    tmp_list_all = []
    header_all = "ID"
    unknowns_counter = 0
    idx_pattern = r"\w+_(\d+)_"
    for keyword, dictionary in sorted(data2use.iteritems()):
        current_values = []
        text = ""
        for idx, system in enumerate(sorted(systems2use)):
            try:
                # print system, keyword, dictionary
                keys2use, values2use = process_values(dictionary, system, keyword)
                if not keys2use:
                    return False
            except TypeError:
                print 'dictio', dictionary.keys()[:10], dictionary.keys()[-10:]
                print 'key', system
                print 'sf/desc', keyword
                return False
            if not_concatenate_values:
                if values2use[0] == ",":
                    tmp_list_all.append("{0}{1}".format(system, values2use))
                else:
                    tmp_list_all.append("{0},{1}".format(system, values2use))
            else:
                if tmp_list_all[idx][-1] == "," and values2use[0] == ",":
                    tmp_list_all[idx] += values2use[1:]
                else:
                    tmp_list_all[idx] += values2use
            if energies_dictio:
                try:
                    energies_dictio[system]
                except KeyError:
                    pattern = re.search(idx_pattern, system)
                    if pattern is None:
                        # print energies_dictio, system
                        logging.critical("Missing the energy for the system {0}".format(system))
                        logging.info("If there's an energy file the keywords should match all the keywords in the "
                                     "scoring and descriptors files.")
                        logging.info("Terminating the program")
                        sys.exit("If there's an energy file the keywords should match all the keywords in the "
                                 "scoring and descriptors files.\nProgram terminated")
                    else:
                        try:
                            energies_dictio[pattern.group(1)]
                        except KeyError:
                            # print energies_dictio, system
                            logging.critical("Missing the energy for the system {0}".format(system))
                            logging.info("After trying the index approach it still fails."
                                         "If there's an energy file the keywords should match all the keywords in the "
                                         "scoring and descriptors files.")
                            logging.info("Terminating the program")
                            sys.exit("If there's an energy file the keywords should match all the keywords in the "
                                     "scoring and descriptors files.\nProgram terminated")
                        else:
                            # print pattern.group(1), energies_dictio[pattern.group(1)]
                            energy_value = energies_dictio[pattern.group(1)]
                            # print pattern.group(1), energies_dictio
                else:
                    energy_value = energies_dictio[system]
                if values2use[0] == ",":
                    current_values.append("{0}{1},{2}".format(system, values2use, energy_value))
                else:
                    current_values.append("{0},{1},{2}".format(system, values2use, energy_value))
                if not text:
                    if keys2use[0] == ",":
                        text = "ID{0},Exp_energy\n".format(keys2use)
                    else:
                        text = "ID,{0},Exp_energy\n".format(keys2use)
            else:
                if values2use[0] == ',':
                    current_values.append("{0}{1}".format(system, values2use))
                else:
                    current_values.append("{0},{1}".format(system, values2use))
                if not text:
                    if keys2use[0] == ",":
                        text = "ID{0}\n".format(keys2use)
                    else:
                        text = "ID,{0}\n".format(keys2use)
        if not text or not current_values:
            logging.error("ERROR: There's no values for the keyword {0}. Just HOW?".format(keyword))
            continue
        if keys2use[0] != ",":
            header_all += "," + keys2use
        else:
            header_all += keys2use
        not_concatenate_values = False
        text += "\n".join(current_values)
        if type_of_data == "descriptors":
            if keyword == "structural":
                output_filename = output_name + "_structural_descriptors_common_systems.csv"
            elif keyword == "mmgbsa":
                output_filename = output_name + "_energy_descriptors_common_systems.csv"
            else:
                # print keyword
                unknowns_counter += 1
                output_filename = output_name + "_unknown_descriptor_{0:02d}.csv".format(unknowns_counter)
                logging.error("ERROR: The descriptors {0} doesn't agree with any of the predetermined descriptors it"
                              "will be written to the file {1}".format(keyword, output_filename))
            with open(output_filename, 'w') as outfile:
                outfile.write(text)
    if energies_dictio:
        all_descriptors_energies = []
        for system in sorted(systems2use):
            try:
                energies_dictio[system]
            except KeyError:
                pattern = re.search(idx_pattern, system)
                if pattern:
                    try:
                        energies_dictio[pattern.group(1)]
                    except KeyError:
                        pass
                    else:
                        all_descriptors_energies.append(energies_dictio[pattern.group(1)])
            else:
                all_descriptors_energies.append(energies_dictio[system])

        tmp_list_all = ["{0},{1}".format(line, energy_value) for line, energy_value in
                        zip(tmp_list_all, all_descriptors_energies)]
        header_all += ",Exp_energy"
    merged_info_text = header_all + "\n" + "\n".join(tmp_list_all)
    if type_of_data == "descriptors":
        output_filename = output_name + "_all_descriptors_merged.csv"
    else:
        output_filename = output_name + "_all_sfs_merged.csv"
    with open(output_filename, 'w') as outfile:
        outfile.write(merged_info_text)
    return True


parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("-dsx_files", default=False, nargs="+", help=extract_sf_help.dsx_files)
parser.add_argument("-xscore_files", default=False, nargs="+", help=extract_sf_help.xscore_files)
parser.add_argument("-binana_files", default=False, nargs="+", help=extract_sf_help.binana_files)
parser.add_argument("-mmgbsa_files", default=False, nargs="+", help=extract_sf_help.mmgbsa_files)
parser.add_argument("-rf_score_file", default=False, help=extract_sf_help.rf_score_file)
parser.add_argument("-vina_files", default=False, nargs="+", help=extract_sf_help.vina_files)
parser.add_argument("-xglide_file", default=False, help=extract_sf_help.xglide_file)
parser.add_argument("-glide_ranking_csv_file", default=False, nargs=2,
                    help=extract_sf_help.glide_ranking_csv_file)
parser.add_argument("-pele_file", default=False, help=extract_sf_help.pele_file)
parser.add_argument("-mmgbsa_as_sf", action="store_true", help=extract_sf_help.mmgbsa_as_sf)
parser.add_argument("-rotable_bonds_files", default=False, nargs="+", help=extract_sf_help.rotable_bonds_files)
parser.add_argument("-energies_file", default=False, help=extract_sf_help.energies_file)
parser.add_argument("-obabel_desc", default=False, nargs='+', help=extract_sf_help.obabel_desc)
parser.add_argument("-convert", "-convert2deltaG", nargs="+", default=["xscore", "nn_score"],
                    help=extract_sf_help.convert)
parser.add_argument("-conversion_temperature", "-temperature", type=int, default=300,
                    help=extract_sf_help.conversion_temperature)
parser.add_argument("-conversion_r_value", "-R", default=0.002, help=extract_sf_help.conversion_r_value)
parser.add_argument("-output_general_name", required=True, help=extract_sf_help.output_general_name)
parser.add_argument("-ensemble", action="store_true", help=extract_sf_help.ensemble_data)
parser.add_argument("-log_file", default="sf_extraction_log.txt", help=extract_sf_help.log_file)
args = parser.parse_args()

# with this line I check if any of the scoring functions has been entered.
if not {key: value for key, value in vars(args).iteritems() if "file" in key and 'log' not in key and value}:
    parser.error("No action requested. At least one of the scoring functions or descriptors should be provided.")

logging.basicConfig(filename="sf_extraction_log.txt", format="%(message)s", level=logging.INFO, filemode="w")
logging.info("{} : Program starting".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M")))

warnings_counter = 0
errors_counter = 0

energy_dictionary = {}
if args.energies_file:
    logging.info("Energy file to use: {0}".format(args.energies_file))
    if not os.path.isfile(args.energies_file):
        logging.error("ERROR: The energy file specified doesn't exist.")
    else:
        with open(args.energies_file) as infile:
            for l in infile:
                l = l.strip()
                if not l:
                    continue
                if len(l.split(",")) == 1:
                    l = l.split(";")
                else:
                    l = l.split(",")
                if len(l) == 2 or len(l) > 3:
                    energy_dictionary[l[0]] = l[1]
                elif len(l) == 3:
                    energy_dictionary[l[1]] = l[2]
                else:
                    logging.error("Check the energies files it should contain at least the id and the energy "
                                  "in the two first columns and the separator should be  ',' or ';'")
                    errors_counter += 1

descriptors = {}
scoring_functions = {}
if args.binana_files:
    descriptors["structural"] = {}
    logging.info("Extracting binana values.")
    for filename in args.binana_files:
        logging.info(" - Working with file: {0}".format(filename))
        if not os.path.isfile(filename):
            logging.warning("## WARNING: The file {0} doesn't exist. It'll be skipped.".format(filename))
            warnings_counter += 1
        try:
            pdb_code, binana_dictio = binana_extraction(filename, args.ensemble)
        except TypeError:
            continue
        else:
            descriptors["structural"][pdb_code] = binana_dictio

if args.mmgbsa_files:
    descriptors["mmgbsa"] = {}
    if args.mmgbsa_as_sf:
        scoring_functions["mmgbsa"] = {}
    logging.info("Extracting mmgbsa values.")
    for filename in args.mmgbsa_files:
        if not os.path.isfile(filename):
            logging.warning(" # WARNING: The file {0} doesn't exist. It'll be skipped.".format(filename))
        if filename[-4:] != ".csv":
            logging.error(" # ERROR: The program expect .csv files separated by comas.")
            logging.info(" - The following file won't be parsed:\n{0}".format(filename))
            continue
        logging.info(" - Working with file: {0}".format(filename))
        if not os.path.isfile(filename):
            print "## WARNING: The file {0} doesn't exist. It'll be skipped.".format(filename)
            warnings_counter += 1
        try:
            pdb_code, mmgbsa_values, mmgbsa_score = mmgbsa_extraction(filename, args.ensemble)
        except TypeError:
            pass
        else:
            descriptors["mmgbsa"][pdb_code] = mmgbsa_values
            if args.mmgbsa_as_sf:
                scoring_functions["mmgbsa"][pdb_code] = mmgbsa_score

if args.obabel_desc:
    logging.info("Computing the pybel descriptors.")
    desired_descriptors = ["MW", "TPSA", "logP"]
    try:
        descriptors["structural"]
    except KeyError:
        descriptors["structural"] = {}
    for filename in args.obabel_desc:
        if '.pdb' not in filename:
            logging.warning("WARNING: The file {0} doesn't have the right extension."
                            " It'll be skipped.".format(filename))
        ligand_molecule = py.readfile('pdb', filename).next()
        ligand_descriptors = ligand_molecule.calcdesc()
        if args.ensemble:
            name_pattern = r"[_{0}]*(\w+_\d+)_.*ligand".format(os.sep)
        else:
            name_pattern = r"[_{0}]*([a-z0-9]+_\d+)_.*ligand".format(os.sep)
        pattern = re.search(name_pattern, filename, re.IGNORECASE)
        if pattern is None:
            logging.warning(" # WARNING: The file {0} doesn't match the naming convention.".format(filename))
            continue
        else:
            pdb_code = pattern.group(1)
        for name in desired_descriptors:
            try:
                descriptors["structural"][pdb_code]
            except:
                pass
            else:
                descriptors["structural"][pdb_code][name] = ligand_descriptors[name]


if args.rotable_bonds_files:
    try:
        descriptors["structural"]
    except KeyError:
        descriptors["structural"] = {}
    logging.info("Extracting rotable bonds.")
    for filename in args.rotable_bonds_files:
        pdb_code = filename.split(os.sep)[-1].split("_ligand.pdbqt")[0]
        with open(filename, 'r') as infile:
            text = infile.read()
        rotable_bonds_pattern = r"REMARK\s+(\d+)\s+active torsions:"
        pattern = re.search(rotable_bonds_pattern, text)
        if pattern is None:
            logging.error("The ligand has no rotable bonds...")
            continue
        else:
            try:
                descriptors["structural"][pdb_code]
            except KeyError:
                pass
            else:
                descriptors["structural"][pdb_code]["rotable_bonds"] = pattern.group(1)

if args.rf_score_file:
    logging.info("Extracting rf_score values")
    scoring_functions["rf_score"] = {}
    if not os.path.isfile(args.rf_score_file):
        logging.warning(" # WARNING: The file {0} doesn't exist. It'll be skipped.".format(args.rf_score_file))
    else:
        with open(args.rf_score_file, 'r') as infile:
            text = infile.read()
            rf_values_in_output_pattern = r"\"(\w+)(?:_ligand)\",(-*\d+\.\d+)"
            first = True
            if "rf_score" in args.convert:
                scoring_functions["rf_score"] = {pattern.group(1): "{0:0.5f}".format(args.conversion_r_value *
                                                                                     args.conversion_temperature *
                                                                                     (log(10 ** (
                                                                                         -float(pattern.group(2))))))
                                                 for pattern in re.finditer(rf_values_in_output_pattern, text,
                                                                            re.IGNORECASE)}
            else:
                scoring_functions["rf_score"] = {pattern.group(1): pattern.group(2)
                                                 for pattern in re.finditer(rf_values_in_output_pattern, text,
                                                                            re.IGNORECASE)}

if args.dsx_files:
    logging.info("Extracting DSX score")
    scoring_functions["dsx"] = {}
    for filename in args.dsx_files:
        if not os.path.isfile(filename):
            logging.warning(" # WARNING: The file {0} doesn't exist. It'll be skipped.".format(filename))
            continue
        logging.info(" - Working with file {0}".format(filename))
        try:
            pdb_code, dsx_values = dsx_extraction(filename, args.ensemble)
        except TypeError:
            pass
        else:
            scoring_functions["dsx"][pdb_code] = dsx_values

if args.xscore_files:
    # print 'here'
    scoring_functions["xscore"] = {}
    logging.info("Extracting xscore values")
    for filename in args.xscore_files:
        if not os.path.isfile(filename):
            logging.warning(" # WARNING: The file {0} doesn't exist. It'll be skipped.".format(filename))
            continue
        logging.info(" - Working with file {0}".format(filename))
        try:
            pdb_code, xscore_values = xscore_extraction(filename, args.conversion_r_value * args.conversion_temperature,
                                                        args.ensemble)
        except TypeError:
            pass
        else:
            scoring_functions["xscore"][pdb_code] = xscore_values

if args.vina_files:
    logging.info("Extracting vina score")
    scoring_functions["vina"] = {}
    for filename in args.vina_files:
        if not os.path.isfile(filename):
            logging.warning(" # WARNING: The file {0} doesn't exist. It'll be skipped.".format(filename))
            continue
        logging.info(" - Working with file {0}".format(filename))
        try:
            pdb_code, vina_score = vina_extraction(filename, args.ensemble)
        except TypeError:
            pass
        else:
            scoring_functions["vina"][pdb_code] = vina_score

if args.xglide_file:
    logging.info("Extracting glide score")
    if not os.path.isfile(args.xglide_file):
        logging.warning(" # WARNING: The file {0} doesn't exist. It'll be skipped.".format(args.xglide_file))
    else:
        with open(args.xglide_file, 'r') as infile:
            text = infile.read()
        if args.ensemble:
            xglide_log_pattern = r"(\w+_\d+).*\s+\w+\s+\(.*\)\s+\d+\.\d+\s+(-*\d+\.\d+)"
        else:
            xglide_log_pattern = r"([a-z0-9]+_\d+).*\s+\w+\s+\(.*\)\s+\d+\.\d+\s+(-*\d+\.\d+)"
        scoring_functions["glide"] = {x.group(1): x.group(2) for x in re.finditer(xglide_log_pattern, text)}
        if not scoring_functions["glide"]:
            logging.error(" # ERROR: Couldn't find any line with the pattern used in this script to extract "
                          "the scores from the log file of xglide script.")

if args.pele_file:
    logging.info("Extracting pele score")
    if not os.path.isfile(args.pele_file):
        logging.warning(" # WARNING: The file {0} doesn't exist. It'll be skipped.".format(args.pele_file))
    else:
        scoring_functions["pele"] = parse_csv_file(args.pele_file, args.ensemble)

if args.glide_ranking_csv_file:
    if not (os.path.isfile(args.glide_ranking_csv_file[0]) or os.path.isfile(args.glide_ranking_csv_file[1])):
        print "At least one of the two elements given to the option glide_ranking_csv_file should be a file."
        logging.ERROR("No file has been provides for this option so this SF won't be extracted.")
        errors_counter += 1
    else:
        if os.path.isfile(args.glide_ranking_csv_file[0]):
            filename, prefix = args.glide_ranking_csv_file
        else:
            prefix, filename = args.glide_ranking_csv_file[1]
        if prefix[-1] != "_":
            prefix += "_"
        glide_ranking = {}
        with open(filename) as csvfile:
            reader = csv.DictReader(csvfile, delimiter=",")
            # for index, line in enumerate(reader, 1):
            for sim_num, line in enumerate(reader, 1):
                if not ("Title" in line or "docking score" in line):
                    logging.ERROR("The file {0} doesn't have either the field Title or the docking score field"
                                  "which are mandatory, so this SF won't be extracted.")
                    errors_counter += 1
                    break
                else:
                    sim_id = prefix + str(sim_num)
                    if "receptor" in line['Title']:
                        continue
                    glide_ranking[sim_id] = line['docking score']
        scoring_functions["glide"] = glide_ranking

logging.info("Finding the systems with all the scoring functions")
systems_with_all_sfs_descriptors = set()
all_systems = set()
for dictionary in scoring_functions.values() + descriptors.values():
    if not systems_with_all_sfs_descriptors:
        systems_with_all_sfs_descriptors = set(dictionary.keys())
        all_systems = set(dictionary.keys())
    else:
        systems_with_all_sfs_descriptors = systems_with_all_sfs_descriptors.intersection(set(dictionary.keys()))
        all_systems = all_systems.union(set(dictionary.keys()))

logging.info("{0:d} systems out of {1:d} have all the scoring functions and descriptors correctly "
             "computed".format(len(systems_with_all_sfs_descriptors), len(all_systems)))
if systems_with_all_sfs_descriptors:
    if descriptors:
        logging.info("Writing the merged file for the descriptors.")
        merge_data(descriptors, "descriptors", systems_with_all_sfs_descriptors, energy_dictionary,
                   args.output_general_name)
    if scoring_functions:
        logging.info("Writing the merged file for the scoring functions.")
        merge_data(scoring_functions, "sfs", systems_with_all_sfs_descriptors, energy_dictionary,
                   args.output_general_name)

unknown_descriptors_counter = 0
for keyword, dictionary in sorted(descriptors.iteritems()):
    current_values = []
    header = ""
    logging.info("Processing the values from the {0} descriptors.".format(keyword))
    for system, elements in sorted(dictionary.iteritems()):
        keys2use, values2use = process_values(dictionary, system, keyword)
        if energy_dictionary:
            # print 'here'
            try:
                current_values.append("{0},{1},{2}".format(system, values2use, energy_dictionary[system]))
            except KeyError as e:
                idx_pattern = r"\w+_(\d+)_"
                pattern = re.search(idx_pattern, system)
                if pattern is None:
                    logging.critical("Missing the energy for the system {0}".format(e[0]))
                    logging.info("If there's an energy file the keywords should match all the keywords in the "
                                 "scoring and descriptors files.")
                    logging.info("Terminating the program")
                    sys.exit("If there's an energy file the keywords should match all the keywords in the "
                             "scoring and descriptors files.\nProgram terminated")
                else:
                    try:
                        energy = energy_dictionary[pattern.group(1)]
                    except KeyError:
                        logging.critical("Missing the energy for the system {0}".format(e[0]))
                        logging.info("If there's an energy file the keywords should match all the keywords in the "
                                     "scoring and descriptors files.")
                        logging.info("Terminating the program")
                        sys.exit("If there's an energy file the keywords should match all the keywords in the "
                                 "scoring and descriptors files.\nProgram terminated")
                    else:
                        if values2use[-1] == ",":
                            current_values.append("{0},{1}{2}".format(system, values2use, energy))
                        if values2use[0] == ",":
                            current_values.append("{0}{1},{2}".format(system, values2use, energy))
                        else:
                            current_values.append("{0},{1},{2}".format(system, values2use, energy))
                    if not header:
                        header = "ID{0},Exp_energy\n".format(keys2use)
            else:
                if not header:
                    header = "ID{0},Exp_energy\n".format(keys2use)
        else:
            if values2use[-1] == ",":
                current_values.append("{0},{1}".format(system, values2use[:-1]))
            elif values2use[0] == ",":
                current_values.append("{0}{1}".format(system, values2use))
            else:
                current_values.append("{0},{1}".format(system, values2use))
                # current_values = ["{0},{1}".format(system, ",".join(map(str, dictio.values())))
                #                   for system, dictio in dictionary.iteritems()]
        if not header:
            # print 'c', keys2use
            header = "{0}{1}\n".format("ID", keys2use)
    text = header + "\n".join(current_values)
    if keyword == "structural":
        output_filename = args.output_general_name + "_structural_descriptors_all_systems.csv"
    elif keyword == "mmgbsa":
        output_filename = args.output_general_name + "_energy_descriptors_all_systems.csv"
    else:
        unknown_descriptors_counter += 1
        output_filename = args.output_general_name + \
                          "_unknown_descriptor_{0:02d}_all_systems.csv".format(unknown_descriptors_counter)
        logging.warning("WARNING: The descriptors {0} doesn't agree with any of the predetermined descriptors it"
                        "will be written to the file {1}".format(keyword, output_filename))
    with open(output_filename, 'w') as outfile:
        outfile.write(text)

for keyword, dictionary in sorted(scoring_functions.iteritems()):
    header = None
    current_values = []
    logging.info("Processing the values from the scoring function: {0} .".format(keyword))
    for system in sorted(dictionary.keys()):
        keys2use, values2use = process_values(dictionary, system, keyword)
        if header is None:
            if keys2use[0] == ",":
                header = "ID" + keys2use
            else:
                header = "ID," + keys2use
        if energy_dictionary:
            try:
                current_values.append("{0},{1}{2}".format(system, values2use, energy_dictionary[system]))
            except KeyError as e:
                idx_pattern = r"\w+_(\d+)_"
                pattern = re.search(idx_pattern, system)
                if pattern is None:
                    logging.critical("Missing the energy for the system {0}".format(e[0]))
                    logging.info("If there's an energy file the keywords should match all the keywords in the "
                                 "scoring and descriptors files.")
                    logging.info("Terminating the program")
                    sys.exit("If there's an energy file the keywords should match all the keywords in the "
                             "scoring and descriptors files.\nProgram terminated")
                else:
                    try:
                        energy = energy_dictionary[pattern.group(1)]
                        # print 'a', current_values, energy
                    except KeyError:
                        logging.critical("Missing the energy for the system {0}".format(e[0]))
                        logging.info("If there's an energy file the keywords should match all the keywords in the "
                                     "scoring and descriptors files.")
                        logging.info("Terminating the program")
                        sys.exit("If there's an energy file the keywords should match all the keywords in the "
                                 "scoring and descriptors files.\nProgram terminated")
                    else:
                        if values2use[-1] == ",":
                            current_values.append("{0},{1}{2}".format(system, values2use, energy))
                        elif values2use[0] == ',':
                            current_values.append("{0}{1},{2}".format(system, values2use, energy))
                        else:
                            current_values.append("{0},{1},{2}".format(system, values2use, energy))
        else:
            if values2use[-1] == ",":
                current_values.append("{0},{1}".format(system, values2use[:-1]))
            elif values2use[0] == ",":
                current_values.append("{0}{1}".format(system, values2use))
            else:
                current_values.append("{0},{1}".format(system, values2use))
    if energy_dictionary:
        if header[-1] == ",":
            header += "Exp_energy"
        else:
            header += ",Exp_energy"
    text = header + "\n" + "\n".join(current_values)
    output_filename = args.output_general_name + "_" + keyword + "_all_systems.csv"
    with open(output_filename, 'w') as outfile:
        outfile.write(text)

logging.info("{} : Program finished correctly.".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M")))
logging.shutdown()
