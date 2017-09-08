input_files_desc = "The files containing either a complex formed by a ligand and a protein or just one ligand. " \
                   "The file can be either in .mae or in .pdb format."
receptor_desc = "This option specifies the file containing the protein to use to form the complexes (one per file in " \
                "input_file) in .mae or .pdb format. If it isn't present the program will asume that the input files " \
                "are complexes with the ligand in the chain Z and try to extract them." \
                "The docking process should have been done previously, and it should be ready to form " \
                "complexes with the files specified in input_files, just by simply merging the files."
subfolder_path_desc = "This option allows the user to specify the path to create the subfolders. By default the " \
                      "subfolder will be generated in the same directory where the program is executed."
log_file_desc = "This option specifies the file where to write the log if the program."
pele_folder_desc = "This option specifies the path to look for the Data and Documents folders needed to launch PELE " \
                   "simulations."
debug_desc = "This option is mostly for the developers, when present the program prints more messages in every step, " \
             "that may help to find errors."
template_desc = "The configuration file used as template, to generate all the configurations files needed."
not_interactive_desc = "If this option is present the program won't make any requests to the user, " \
                       "and will assume certain things (see the README.cmd file)."
schrodinger_path_desc = "The complete path to schrodinger software, without final separator."
plop_path_desc = "The complete path to the PlopRotTemp.py script, including the script name."
pele_license_desc = "The complete path to use in the configuration file for PELE indicating where the" \
                    " license is present."
every_desc = "If in the template file the constraints keyword is specified, this option will " \
             "indicate how many occurrences of the atoms to constraint are between constraints"
constraint_desc = "This option specifies the constraint value to use when extracting constraints from a pdb" \
                  " file for the PELE configuration file."
atoms2constraint_desc = "This option is a list of the atom names (blank spaces included) to apply contraints when" \
                        " generating the configuration file for PELE."
mutations_program_path_desc = "The complete path to the mutations program script to prepare pdb files for PELE."
obc_param_generator_desc = "The complete path to the solventOBCParamsGenerator script to generate the OBC parameters" \
                           "for the ligand for PELE."
ligand_chain_desc = "The chain where the ligand is present. This chain should contain only the ligand and unless it " \
                    "is a peptidic ligand it should be formed by one single residue."
no_templates = "This option acts as a flag when present the program won't generate the template nor the ligand " \
               "rotamer library."
adaptive_sampling_desc = "If used, it should specify the path to the templatized control file to be used " \
                         "to create the control files for the adaptiveSampling.py script."
