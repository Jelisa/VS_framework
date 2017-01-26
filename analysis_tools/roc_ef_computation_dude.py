import csv
from argparse import ArgumentParser
from sklearn.metrics import auc
import sys
import numpy as np
import pylab as pl
import re
import os

import roc_ef_help


def obtain_energies_from_activity_file(ids2pick, activity_dictionary):
    activities = []
    for identifier in ids2pick:
        pattern = re.search(r'\w+_(\d+)_', identifier)
        if pattern is None:
            pattern = re.search(r'\w+_(\d+)_*', identifier)
            if pattern is None:
                number = identifier
            else:
                number = pattern.group(1)
        else:
            number = pattern.group(1)
        activities.append(activity_dictionary[number])
    return np.asarray(activities)

# np.asarray([simulation_activity_dictio[sim_id] for sim_id in systems_ids1])


def select_systems(sims2select, sfs2use, sfs_systems):
    indexes2pick = [index for index, sfs_system in enumerate(sfs_systems) if sfs_system in sims2select]
    sims_present = [sim for sim in sims2select if sim in sfs_systems]
    sfs2return = np.asarray([sfs2use[x] for x in indexes2pick])
    return sfs2return, sims_present


def read_scoring_functions_file(filename, keys=None):
    sfs_array = []
    ids = []
    with open(filename, 'r') as infile:
        for line in infile:
            if keys is None:
                keys = line.strip().split(',')[1:]
            else:
                lin = line.strip().split(',')[1:]
                current_id = line.strip().split(',')[0]
                sfs_array.append([float(x) for x in lin])
                a = re.search(r"(\w+_\d+)_", current_id)
                if a:
                    ids.append(a.group(1))
                else:
                    ids.append(current_id)
    sfs_array = np.asarray(sfs_array)
    return keys, sfs_array, ids


def boolean_matrix_check(matrix):
    tmp_mat = np.asarray([x for x in matrix if x == 0 or x == 1])
    if tmp_mat.shape == matrix.shape:
        is_boolean = True
    else:
        is_boolean = False
    return is_boolean


def plot_all_sfs(keys, values_matrix, delta_g_values, subfix, ef_thresholds):
    roc_values = {}
    ef_dictionary = {}
    # Plot the image with all the SFs.
    title = "{} {} ".format(args.general_title, "All SFS")
    fig, ax = pl.subplots()
    color = pl.cm.rainbow(np.linspace(0, 1, len(keys)))
    # TODO: check if the energies matrix is boolean or not.
    if boolean_matrix_check(delta_g_values):
        actives = delta_g_values[delta_g_values.nonzero()].shape[0]
    else:
        tmp = delta_g_values[delta_g_values.nonzero()]
        actives = tmp[tmp < 0].shape[0]
    inactives = delta_g_values.shape[0] - actives
    # print keys
    for score, c in zip(keys, color):
        print score,
        sf2use2 = keys.index(score)
        # print sf2use2
        values2study = [[values_matrix[i, sf2use2], delta_g_values[i]] for i in range(len(delta_g_values))]
        # print values2study2
        experimental_value_prediction_ordered = [x[1] for x in sorted(values2study)]
        tpr, fpr, auc_val = compute_roc(experimental_value_prediction_ordered, actives, inactives)
        roc_values[score] = [fpr, tpr, auc_val]
        ax.plot(fpr, tpr, label='{0} ROC curve (area = {1:0.2f})'.format(score, auc_val), c=c)
        print "auc:", auc_val
        efs = compute_enrichment_factor(experimental_value_prediction_ordered, ef_thresholds, actives)
        print score, "EF:", efs
        ef_dictionary[score] = efs
    ax.plot([0, 1], [0, 1], color='navy', linestyle='--', label="Random selection")
    ax.legend(loc="lower right")
    ax.set_ylim([0, 1])
    ax.set_xlim([0, 1])
    ax.set_title(title)
    ax.xaxis.set_label_text("FPR")
    ax.yaxis.set_label_text("TPR")
    pl.show()
    fig.savefig(args.output_prefix + "_{}_all_sfs.png".format(subfix), bbox_inches='tight')
    return roc_values, ef_dictionary


def compute_roc(ordered_list, actives, inactives):
    """
    Function to compute the ROC curve
    :param ordered_list: a list containing the activity in a boolean manner: 0 (inactive) or 1 (active)
    :param actives: the number of actives compounds
    :param inactives: the number of inactives compounds
    :return: returns two lists containing: true positives ratio, false positives ratio, each and the AUC value.
    """
    # deltag_non_zero_matrix2 = actives[actives.nonzero()]
    # actives = deltag_non_zero_matrix2.shape[0]
    # inactives2 = actives.shape[0] - actives
    # print actives2, inactives2
    tpr, fpr = [], []
    tp = 0.0
    fp = 0.0
    for value in ordered_list:
        if value:
            tp += 1
        else:
            fp += 1
        tpr.append(tp / float(actives))
        # print deltag_non_zero_matrix2
        fpr.append(fp / float(inactives))
    # print tp, fp
    auc_value = auc(fpr, tpr)
    return tpr, fpr, auc_value


def compute_enrichment_factor(ordered_list, thresholds_2_use, actives):
    enrichment_factor = {}
    enrichment_factor['max'] = 1.0 / (float(actives) / len(ordered_list))
    for n in thresholds_2_use:
        true_positives = [x for x in ordered_list[:n] if x <= 0]
        enrichment_factor_val = (float(len(true_positives)) / n) / (float(actives) / len(ordered_list))
        enrichment_factor[n] = enrichment_factor_val
    return enrichment_factor


def write_roc_files(general_dictionary, filename_prefix):
    for score, values in general_dictionary.iteritems():
        filename = "{0}_roc_curve_values_{1}.csv".format(filename_prefix, score)
        fpr, tpr, auc = values
        with open(filename, 'w') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter=',')
            csvwriter.writerow(["FPR", "TPR", "AUC"])
            csvwriter.writerow([fpr[0], tpr[0], auc])
            csvwriter.writerows([[x, y] for x, y in zip(fpr[1:], tpr[1:])])


def write_ef_files(general_dictionary, filename_prefix):
    # for score, values in general_dictionary.iteritems():
    filename = "{0}_ef_values.csv".format(filename_prefix)
    with open(filename, 'w') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter=',')
        # The following line access the first dictionary inside general_dictionary to obtain its keys
        # and thus obtain the right order of the keys.
        csvwriter.writerow(['score'] + sorted(general_dictionary[general_dictionary.keys()[0]]))
        csvwriter.writerows([sc] + [ef_dict[x] for x in sorted(ef_dict.keys())]
                            for sc, ef_dict in general_dictionary.iteritems())


parser = ArgumentParser()
parser.add_argument("-all_sfs_file", required=True, help=roc_ef_help.all_sfs_file_help)
parser.add_argument("-previous_scoring_functions_file", help=roc_ef_help.pele_best_energy_file_help)
parser.add_argument("-general_title", default="", help=roc_ef_help.general_title_help)
parser.add_argument("-output_prefix", default="roc_ef_comp", help=roc_ef_help.output_prefix_help)
parser.add_argument("-activities_file", required=True, help=roc_ef_help.activities_file_help)
parser.add_argument("-ef_thresholds", default=[10, 20, 50, 100, 250, 500], type=int, nargs='+')
parser.add_argument("-output_folder", default=".")
parser.add_argument("-systems2use", default=False)
args = parser.parse_args()

# Piece of code to set the images size.
# Get current size
fig_size = pl.rcParams["figure.figsize"]

# Prints: [8.0, 6.0]
print "Current size:", fig_size

# Set figure width to 12 and height to 9
fig_size[0] = 13
fig_size[1] = 10
pl.rcParams["figure.figsize"] = fig_size
print "Current size:", fig_size

roc_values0 = {}

if args.output_folder != "'.":
    if args.output_folder[-1] != os.sep:
        args.output_folder += os.sep
    if not os.path.isdir(args.output_folder):
        os.mkdir(args.output_folder)
    args.output_prefix = args.output_folder + args.output_prefix

simulation_activity_dictio = {}
simulation_name_dictio = {}
with open(args.activities_file) as csvfile:
    csvparser = csv.DictReader(csvfile)
    for line in csvparser:
        molecule_name = line['name']
        sim_id = line['sim_id']
        activity = int(line['activity'])
        simulation_activity_dictio[sim_id] = activity
        simulation_name_dictio[sim_id] = molecule_name

ids2select = []
if args.systems2use:
    with open(args.systems2use) as infile:
        csv_parser = csv.DictReader(infile, delimiter=',')
        ids2select = [line['ID'] for line in csv_parser]


# Defines a matrix from the _all_sfs_merged.csv with all the sfs values and a to lists
# that are used later for the rest of the script as the dictionary with the systems to use
# keys Will be used as the dictionary to substract the different SFs.
# all_sf_values1 Will become a matrix containing the values for all the SFs.
# systems_ids1 Will be used as the dictionary containing the systems to use.
# print all_sf_filename
print "Working with the file: {0}".format(args.all_sfs_file)
keys1, all_sf_values1, systems_ids1 = read_scoring_functions_file(args.all_sfs_file)
if ids2select:
    all_sf_values1, systems_ids1 = select_systems(ids2select, all_sf_values1, systems_ids1)
if "Exp_energy" in keys1:
    delta_g_values1 = all_sf_values1[:, keys1.index("Exp_energy")]
else:
    delta_g_values1 = obtain_energies_from_activity_file(systems_ids1, simulation_activity_dictio)

roc_values1, ef_dicitonaries1 = plot_all_sfs(keys1, all_sf_values1, delta_g_values1, "if", args.ef_thresholds)

write_roc_files(roc_values1, args.output_prefix + "_if")
write_ef_files(ef_dicitonaries1, args.output_prefix + "_if")

# print 0, keys1
# print len(keys1), all_sf_values1.shape
if args.previous_scoring_functions_file:
    print "Working with the file: {0}".format(args.all_sfs_file)
    keys0, all_sf_values0, systems_ids0 = read_scoring_functions_file(args.all_sfs_file)
    if ids2select:
        all_sf_values0, systems_ids0 = select_systems(ids2select, all_sf_values0, systems_ids1)
    if "Exp_energy" in keys1:
        delta_g_values0 = all_sf_values0[:, keys0.index("Exp_energy")]
    else:
        delta_g_values0 = obtain_energies_from_activity_file(systems_ids0, simulation_activity_dictio)
    roc_values0, ef_dicitonaries0 = plot_all_sfs(keys0, all_sf_values0, delta_g_values0, 'init', args.ef_thresholds)
    write_roc_files(roc_values0, args.output_prefix + "_init")
    write_ef_files(ef_dicitonaries0, args.output_prefix + "_init")
    common_scores = set.intersection(set(keys1), set(keys0))
    for score in common_scores:
        title = "{0} (SF:{1}) ".format(args.general_title, score)
        fig_tmp, ax_tmp = pl.subplots()
        print score,
        # print roc_values0
        fpr0, tpr0, auc_val0 = roc_values0[score]
        fpr1, tpr1, auc_val1 = roc_values1[score]
        ax_tmp.plot(fpr0, tpr0, label='initial {0} ROC curve (area = {1:0.2f})'.format(score, auc_val0), c='b')
        ax_tmp.plot(fpr1, tpr1, label='i.f. {0} ROC curve (area = {1:0.2f})'.format(score, auc_val1), c='r')
        print "delta auc: {0:6.4}".format(auc_val0 - auc_val1)
        ax_tmp.plot([0, 1], [0, 1], color='navy', linestyle='--', label="Random selection")
        ax_tmp.legend(loc="lower right")
        ax_tmp.set_ylim([0, 1])
        ax_tmp.set_xlim([0, 1])
        ax_tmp.set_title(title)
        ax_tmp.xaxis.set_label_text("FPR")
        ax_tmp.yaxis.set_label_text("TPR")
        pl.show()
        fig_tmp.savefig(args.output_prefix + "_{}_comparison.png".format(score), bbox_inches='tight')
    delta_ef = {score: {th: ef_dicitonaries1[score][th] - ef_dicitonaries0[score][th] for th in args.ef_thresholds}
                for score in common_scores}
    write_ef_files(delta_ef, args.output_prefix + "_delta")
