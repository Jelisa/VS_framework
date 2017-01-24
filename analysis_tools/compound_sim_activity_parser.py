import csv
from argparse import ArgumentParser
import sys

parser = ArgumentParser()
parser.add_argument("-glide_docking_rank", required=True)
parser.add_argument("-actives_file", required=True)
parser.add_argument("-inactives_file", default=False)
parser.add_argument("-output_file", required=True)
args = parser.parse_args()

# simulations_ids = []
molecule_names = []
with open(args.glide_docking_rank) as csvfile:
    reader = csv.DictReader(csvfile)
    # for index, line in enumerate(reader, 1):
    for line in reader:
        try:
            line["Title"]
        except KeyError:
            print "FATAL: The glide file doesn't have the 'Title' file which is mandatory for this program."
            print "Interrupting the program"
        else:
            if "receptor" in line['Title']:
                continue
            molecule_names.append(line['Title'])
            # simulations_ids.append(str(index))

name_index = 2  # Value hardcoded for the DUD-E dataset activities files .ism format
# In the .ism file the values are in the form of type_of_experiment = experimental_value , thus I'm using the
# = as the marker to extract these values.
name_experiment_type_dictio = {}
name_experiment_value_dictio = {}
name_activity_dictio = {}
with open(args.actives_file) as csvfile:
    reader = csv.reader(csvfile, delimiter=" ")
    for x, line in enumerate(reader):
        try:
            name = line[name_index]
            equal_index = line.index('=')
            experiment_type = line[equal_index - 1]
            experimental_value = line[equal_index + 1]
        except IndexError:
            print "There are problems with the columns in the file {0} in line number {1}".format(args.actives_file,
                                                                                                  x)
            print "ignoring this line."
        else:
            if name in name_experiment_value_dictio.keys() or name in name_experiment_type_dictio.keys():
                print "There are duplicates in the file, did you know it? Talk to the developer to tell her what to do."
                sys.exit()
            else:
                name_experiment_type_dictio[name] = experiment_type
                name_experiment_value_dictio[name] = experimental_value
                name_activity_dictio[name] = '1'

if args.inactives_file:
    reader = csv.reader(csvfile, delimiter=" ")
    for x, line in enumerate(reader):
        try:
            name = line[name_index]
            equal_index = line.index('>')
            experiment_type = line[equal_index - 1]
            experimental_value = line[equal_index + 1]
        except IndexError:
            print "There are problems with the columns in the file {0} in line number {1}".format(args.actives_file,
                                                                                                  x)
            print "ignoring this line."
        else:
            if name in name_experiment_value_dictio.keys() or name in name_experiment_type_dictio.keys():
                print "There are duplicates in the file, did you know it? Talk to the developer to tell her what to do."
                sys.exit()
            else:
                name_experiment_type_dictio[name] = experiment_type
                name_experiment_value_dictio[name] = experimental_value
                name_activity_dictio[name] = '0'

unset_molecules = [molecule for molecule in molecule_names if molecule not in name_activity_dictio.keys()]
for molecule in unset_molecules:
    if molecule in name_experiment_type_dictio.keys():
        continue
    # print molecule,
    name_experiment_type_dictio[molecule] = "UNK"
    name_experiment_value_dictio[molecule] = "UNK"
    name_activity_dictio[molecule] = '0'

output_text = ["name,sim_id,activity,exp_val,experiment"]
for index, molecule in enumerate(molecule_names, 1):
    output_line = ",".join([molecule, str(index), name_activity_dictio[molecule],
                            name_experiment_value_dictio[molecule], name_experiment_type_dictio[molecule]])
    output_text.append(output_line)

with open(args.output_file, 'w') as outfile:
    outfile.write("\n".join(output_text))

n_actives = len([name for name, value in name_activity_dictio.iteritems() if value != '0'])
print "There are {0} actives out of {1} molecules.".format(n_actives, len(name_activity_dictio.keys()))
