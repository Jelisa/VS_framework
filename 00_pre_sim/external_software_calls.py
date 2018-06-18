"""
A module to host all the functions that call external programs.
"""

import enviroment_parameters

schrodinger_converter = "{}/utilities/pdbconvert ".format(enviroment_parameters.schrodinger_path)
schrodinger_mae2pdb_convert_command = schrodinger_converter + " -imae {} -opdb {}"
schrodinger_pdb2mae_convert_command = schrodinger_converter + " -ipdb {} -omae {}"
ploprottemp = "{}/utilities/python {}".format(enviroment_parameters.schrodinger_path,
                                              enviroment_parameters.plop_path)
ploprottemp_command = ploprottemp + " {}  --mtor 5 --gridres 30 --clean"
mutations_program_command_lig_or_comp = "python " + enviroment_parameters.mutations_program_path + \
                                        " -ipdb {} -make_unique {}"

mutations_program_command_receptor = "python " + enviroment_parameters.mutations_program_path + \
                                     " -ipdb {} "
obc_param_generator = "python {}".format(enviroment_parameters.obc_param_generator_path)
obc_param_command = obc_param_generator + " {}"

obabel_convert_to_pdb = "obabel -i{} {} -O {}"
obabel_convert_to_pdb_and_gen_3d = "obabel -i{} {} -O {} -h --gen3D"
