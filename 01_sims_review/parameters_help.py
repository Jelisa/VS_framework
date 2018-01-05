structure_Selection_program_desc = "This programs reviews the output from any kind of PELE simulations " \
                                   "and selects structures from the trajectory based on serveral " \
                                   "criteria, such as minimum binding energy, distances, RMSD or " \
                                   "energy clustering, etc. And creates one folder for each selected " \
                                   "structure containing a pdb file with the structure. " \
                                   "It also generates a file containing how many structures have " \
                                   "been selected for each folder in the input, a file containing" \
                                   "the PELE energy (binding or normal) of the selected structure" \
                                   "and if the option is chosen it will also extract the initial" \
                                   "PELE energy (binding or normal)."
structure_Selection_input_desc = "A list containing the path to the output folders " \
                                 "from any kind of PELE simulations to analyze. \n" \
                                 "I.E.: /path/to/pele/simulation/system_1/output_obc/"
structure_Selection_output_folder_desc = "The directory where the program will " \
                                         "create the folders containing the " \
                                         "selected structures (one for each " \
                                         "structure selected). It's mandatory"
structure_selection_ligand_chain_desc = "The name of the chain where the ligand" \
                                        " can be found.\nRemember that the ligand" \
                                        " should be the only molecule in this chain."
log_file_desc = "The name for the log file."

structure_Selection_selection_criteria = "The different criteria to select the structure(s)."
structure_Selection_simulation_type_desc = "This option specifies which kind of PELE simulation" \
                                           "has generated the output folder given in the input:" \
                                           "single_core refers to a PELE simulation run with " \
                                           "only one procesor (traditional PELE), mpi refers" \
                                           "to the mpi version of PELE and adaptive to the" \
                                           "adaptive script to launch PELE."
structure_selection_only_statistics = "This option acts as a flag, if present the program " \
                                      "will only compute the number of structures selected " \
                                      "and their energies according to the selection criteria." \
                                      "It won't extract the structures."
structure_Selection_min_ener_desc = "When this option is chosen the program will select " \
                                    "the structure with the minimum (binding or pele) " \
                                    "energy over all the simulation."
ss_min_energy_type_desc = "This option allows the user to choose which type of PELE " \
                          "energy to use: pele_energy is the PELE energy and " \
                          "binding_energy is the B.E. computed by PELE."
structure_Selection_rmsd_clust_desc = "To be developed"
structure_Selection_ener_clust_desc = "To be developed"
ss_energy_c_deltag_desc = "The value in PELE Energy units (Kcal) to use as threshold to do " \
                          "the energy clustering."
# structure_selection_single_structure = "This option acts as a flag, if present the program will select one single " \
#                                        "structure otherwise it will select multiple structures using the following " \
#                                        "criteria {}\nDEFAULT: False".format(structure_Selection_criteria)
# structure_Selection_criteria = "The one with minimum binding energy and the ones that have a bionding energy within " \
#                                "5 Kcal of the minimum one, and then those that have an rmsd between 0.25 and 2 to the " \
#                                "structure with the minimum energy, computed using the heavy atoms of the ligand and " \
#                                "those of the protein within 6 A."
sims_review_initial_energies = "When this option is present the script will write the file initial_PELE_binding_energies.csv " \
                   "that contains the PELE binding energy for the initial complex given to PELE."
sims_review_input_desc = "A list of files containing the output from PELE simulations to analyze."
sims_review_report_desc = "The complete path (including name) of the file where the report will be written. "
# structure_selection_output_desc = "The"