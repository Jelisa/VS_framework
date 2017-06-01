sims_review_input_desc = "A list of files containing the output from PELE simulations to analyze."
sims_review_report_desc = "The complete path (including name) of the file where the report will be written. "
structure_Selection_criteria = "The one with minimum binding energy and the ones that have a bionding energy within " \
                               "5 Kcal of the minimum one, and then those that have an rmsd between 0.25 and 2 to the " \
                               "structure with the minimum energy, computed using the heavy atoms of the ligand and " \
                               "those of the protein within 6 A."
structure_selection_input_desc = "A list of folders where to look for the PELE output files to analyze"
structure_selection_deltag_desc = "The value in PELE Energy units (Kcal) to use as threshold to select multiple " \
                                  "poses from PELE simulations as the best poses."
structure_selection_ligand_chain_desc = "The name of the chain where the ligand can be found.\n" \
                                        "Remember that the ligand should be the only molecule in this chain."
structure_selection_single_structure = "This option acts as a flag, if present the program will select one single " \
                                       "structure otherwise it will select multiple structures using the following " \
                                       "criteria {}".format(structure_Selection_criteria)
structure_selection_only_statistics = "This option acts as a flag, if present the program will only compute the " \
                                      "statistics with the selected criteria."
initial_energies = "When this option is present the script will write the file initial_PELE_binding_energies.csv " \
                   "that contains the PELE binding energy for the initial complex given to PELE."
max_steps = "The x first steps of the PELE simulation to study. I.E.: To use only the first 200 steps of a 1000 " \
            "steps simulation"
log_file = "The name for the log file to write."
# structure_selection_output_desc = "The"