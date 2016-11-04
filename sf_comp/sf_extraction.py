"""
author: Jelisa Iglesias
email: jelisa.iglesias@gmail.com
"""

import argparse
from re import search, IGNORECASE

parser = argparse.ArgumentParser()
parser.add_argument("-input_files", nargs="+", required=True)
parser.add_argument("-folders_path", default="./scoring_functions_values")
parser.add_argument("-output_general_name", default="all")
parser.add_argument("-scoring_functions", nargs="+",
                    default=["glide", "vina", "xscore", "dsx", "mmgbsa", "binana", "rf_score"])
args = parser.parse_args()

scoring_functions_values = {name.lower(): {} for name in args.scoring_functions}

for filename in args.input_files:
    if search(r"glide.*\.log", filename):
        glide_keys = [x for x in scoring_functions_values.keys() if search(r"glide", x)]
        if len(glide_keys) > 1:
            pattern = search(r"(glide)[-_]*([SX]+p)", filename, IGNORECASE)
            if pattern is None:
                print "UPS! The filename {} doesn't match the naming rules when computing more than one glide score."
                continue
            else:
                for x in glide_keys:
                    if pattern.group(2) in x:
                        key = x
                        break
                else:
                    print "The file {} doesn't agree with any of the glide keywords check it!".format(filename)
                    continue
        else:
            key = glide_keys[0]

        with open(filename, 'r') as xglide_log_file:
            start_reading = 0
            gscore_index = None
            receptor_index = None
            for line in xglide_log_file:
                if search(r"gscore", line, IGNORECASE):
                    start_reading += 1
                    # line = line.strip()
                    # for i, element in enumerate(line.split()):
                    #     if search(r"receptor", element, IGNORECASE):
                    #         receptor_index = i
                    #     elif search(r"gscore", element, IGNORECASE):
                    #         gscore_index = i
                    # continue
                elif start_reading == 1:
                    start_reading += 1
                    continue
                elif start_reading == 2:
                    # if receptor_index is None or gscore_index is None:
                    #     print "Something went wrong the patterns couldn't be found."
                    #     print "The program couldn't process the file {}".format(filename)
                    #     break
                    # else:
                    line = line.strip()
                    if line == "":
                        break
                    pattern = search(r'(\w+)\s+(\w+)\s+\(.*\)\s+\d+.*\d+\s+(-*\d+.*\d+)', line)
                    receptor = pattern.group(1)
                    score = float(pattern.group(3))
                    scoring_functions_values[key][receptor] = score
                else:
                    continue

for key, values in scoring_functions_values.iteritems():
    if not values:
        continue
    fileout_name = key + "_values.csv"
    with open(fileout_name, 'w') as outfile:
        outfile.write("\n".join(["{0},{1}".format(receptor, value) for receptor, value in values.iteritems()]))
