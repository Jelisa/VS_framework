input_files_desc = "The files containing either a complex formed by a ligand and a protein or just one ligand. " \
                   "The file can be either in .mae or in .pdb format."

receptor_desc = "This option specifies the file containing the protein to use to form the complexes (one per file in " \
                "input_file) in .mae or .pdb format. If it isn't present the program will asume that the input files " \
                "are complexes with the ligand in the chain Z and try to extract them." \
                "The docking process should have been done previously, and it should be ready to form " \
                "complexes with the files specified in input_files, just by simply merging the files." \

subfolder_path_desc = "This option allows the user to specify the path to create the subfolders. By default the " \
                      "subfolder will be generated in the same directory where the program is executed."

log_file_desc = "This option specifies the file where to write the log if the program. DEFAULT: log.txt"

pele_folder_desc = "This option specifies the path to look for the Data and Documents folders needed to launch PELE " \
                   "simulations."
debug_desc = "This option is mostly for the developers, when present the program prints more messages in every step, " \
             "that may help to find errors."
template_desc = "The configuration file used as template, to generate all the configurations files needed." \
                " DEFAULT: ./template.conf"
not_interactive_desc = "If this option is present the program won't make any requests to the user, " \
                       "and will assume certain things (see the README.cmd file). DEFAULT: False"
schrodinger_path_desc = "The complete path to schrodinger software, without final separator:" \
                        "DEFAULT: /home/jelisa/programs/Schrodinger_Suites_2016-1_Advanced_Linux-x86_64"
plop_path_desc = "The complete path to the PlopRotTemp.py script, including the script name.\n" \
                 "DEFAULT: /home/jelisa/PhD_work/VS/working_plop/PlopRotTermp.py"
pele_license_desc = "The complete path to use in the configuration file for PELE indicating where the" \
                    " license is present.\nDEFAULT: /gpfs/projects/bsc72/PELE++/license"
every_desc = "If in the template file the constraints keyword is specified, this option will " \
             "indicate how many occurrencies of the atoms to constraint are between constraints\n" \
             "DEFAULT: 1"
constraint_desc = "This option specifies the constraint value to use when extracting constraints from a pdb file for the PELE" \
             "configuration file.\nDEFAULT: 0.2"
atoms2constraint_desc = "This option is a list of the atom names (blank spaces included) to apply contraints when" \
                        " generating the configuration file for PELE.\nDEFAULT: ' CA '"
mutations_program_path_desc = "The complete path to the mutations program script to prepare pdb files for PELE.\n" \
                              "DEFAULT: /home/kqtw353/stay/work/mutations_program/mutations_program.py"
obc_param_generator_desc = "The complete path to the solventOBCParamsGenerator script to generate the OBC parameters" \
                           "for the ligand for PELE.\nDEFAULT: " \
                           "/home/kqtw353/stay/work/OBCparamsgen/solventOBCParamsGenerator.py"
