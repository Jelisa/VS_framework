import csv
import pandas as pd
from argparse import ArgumentParser
from sklearn.metrics import auc
import sys
import numpy as np
import pylab as pl
import re
import os

import roc_ef_help


def boolean_matrix_check(matrix):
    tmp_mat = np.asarray([x for x in matrix if x == 0 or x == 1])
    if tmp_mat.shape == matrix.shape:
        is_boolean = True
    else:
        is_boolean = False
    return is_boolean


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
    # print ordered_list
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


def compute_enrichment_factor(ordered_list, thresholds_2_use, actives, boolean_values):
    enrichment_factor = {}
    # enrichment_factor['max'] = 1.0 / (float(actives) / len(ordered_list))
    for n in thresholds_2_use:
        if boolean_values:
            true_positives = [x for x in ordered_list[:n] if x == 1]
        else:
            true_positives = [x for x in ordered_list[:n] if x <= 0]
        enrichment_factor_val = (float(len(true_positives)) / n) / (float(actives) / len(ordered_list))
        enrichment_factor[n] = enrichment_factor_val
    return enrichment_factor


def plot_all_sfs(df, columns2study, subfix, ef_thresholds, gen_title):
    roc_values = {}
    ef_dictionary = {}
    # Plot the image with all the SFs.
    title = "{0} {1} {2}".format(gen_title, "All SFS", subfix)
    fig, ax = pl.subplots()
    color = pl.cm.rainbow(np.linspace(0, 1, len(columns2study)))
    number_of_m = df.shape[0]
    number_of_actives = df['activity'].nonzero()[0].shape[0]
    number_of_inactives = number_of_m - number_of_actives
    if df['activity'].isin([0, 1]).shape == df['activity'].shape:
        boolean_matrix = True
    else:
        boolean_matrix = False
    df.sort_values(["activity"], inplace=True, ascending=False)
    ef_dictionary["Exp_energy"] = compute_enrichment_factor(df['activity'], ef_thresholds,
                                                            number_of_actives, boolean_matrix)
    for score, c in zip(columns2study, color):
        df.sort_values([score], inplace=True)
        compound_names = list(df['name'])
        sim_id = list(df['ID'])
        activities = list(df['activity'])
        tpr, fpr, auc_val = compute_roc(df['activity'].values, number_of_actives, number_of_inactives)
        roc_values[score] = [sim_id, activities, compound_names, fpr, tpr, auc_val]
        efs = compute_enrichment_factor(df['activity'], ef_thresholds, number_of_actives, boolean_matrix)
        ef_dictionary[score] = efs
        ax.plot(fpr, tpr, label='{0} ROC curve (area = {1:0.2f})'.format(score, auc_val), c=c)
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


def write_ef_files(general_dictionary, filename_prefix):
    filename = "{0}_ef_values.csv".format(filename_prefix)
    fieldnames = general_dictionary['Exp_energy'].keys()
    fieldnames.sort()
    fieldnames.insert(0, 'Score')
    ordered_sfs = ['Exp_energy'] + sorted([x for x in general_dictionary.keys() if x != "Exp_energy"])
    with open(filename, 'w') as outfile:
        csvwriter = csv.DictWriter(outfile, fieldnames=fieldnames)
        csvwriter.writeheader()
        for x in ordered_sfs:
            tmp = {fieldnames[0]: x}
            tmp.update(general_dictionary[x])
            csvwriter.writerow(tmp)


def write_roc_files(general_dictionary, filename_prefix):
    for score, values in general_dictionary.iteritems():
        filename = "{0}_roc_curve_values_{1}.csv".format(filename_prefix, score)
        sim_ids, activity, names, fpr, tpr, auc = values
        with open(filename, 'w') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter=',')
            csvwriter.writerow(["sim_ID", "Activity", "Name", "FPR", "TPR", "AUC"])
            csvwriter.writerow([sim_ids[0], activity[0], names[0], fpr[0], tpr[0], auc])
            csvwriter.writerows([[x, y, p, z, q] for x, y, p, z, q in zip(sim_ids[1:], activity[1:], names[1:],
                                                                          fpr[1:], tpr[1:])])


parser = ArgumentParser()
parser.add_argument("-all_sfs_file", required=True, help=roc_ef_help.all_sfs_file_help)
parser.add_argument("-previous_scoring_functions_file", help=roc_ef_help.pele_best_energy_file_help)
parser.add_argument("-general_title", default="", help=roc_ef_help.general_title_help)
parser.add_argument("-output_prefix", default="roc_ef_comp", help=roc_ef_help.output_prefix_help)
parser.add_argument("-activities_file", required=True, help=roc_ef_help.activities_file_help)
parser.add_argument("-ef_thresholds", default=[10, 20, 50, 100, 250, 500], type=int, nargs='+')
parser.add_argument("-output_folder", default=".")
parser.add_argument("-systems2use1", default=False)
parser.add_argument("-systems2use0", default=False)
parser.add_argument('-debug')
args = parser.parse_args()

fig_size = pl.rcParams["figure.figsize"]

# Set figure width to 12 and height to 9
fig_size[0] = 13
fig_size[1] = 10
pl.rcParams["figure.figsize"] = fig_size
pl.ion()

if args.output_folder != "'.":
    if args.output_folder[-1] != os.sep:
        args.output_folder += os.sep
    if not os.path.isdir(args.output_folder):
        os.mkdir(args.output_folder)
    args.output_prefix = args.output_folder + args.output_prefix

name_sim_activity_df = pd.read_csv(args.activities_file)

if_sfs_df = pd.read_csv(args.all_sfs_file)
ids2select_df1 = None
ids2select_df0 = None
compounds = None
if args.systems2use1:
    # print 0
    ids2select_df1 = pd.read_csv(args.systems2use1)
    # print ids2select_df1
    if_sfs_df = if_sfs_df.loc[if_sfs_df['ID'].isin(ids2select_df1['ID'])]
    ids2select_df0 = ids2select_df1
    compounds = name_sim_activity_df['name'].loc[name_sim_activity_df['sim_id'].isin(ids2select_df1['ID'])]
if args.systems2use0:
    # print 1
    ids2select_df0 = pd.read_csv(args.systems2use0)
    if len(ids2select_df1) > len(ids2select_df0):
        compounds = name_sim_activity_df.loc[name_sim_activity_df['sim_id'].isin(ids2select_df0['ID'])]
# print 2
# print if_sfs_df[:5]['HMScore']
if "Exp_energy" in if_sfs_df.columns:
    if args.debug:
        print 'a'
    if_df = if_sfs_df.merge(name_sim_activity_df[['sim_id', 'name']], left_on="ID", right_on="sim_id")
    if_df.rename(columns={"Exp_energy": "activity"}, inplace=True)
else:
    if args.debug:
        print 'b'
    if_df = if_sfs_df.merge(name_sim_activity_df[['sim_id', 'name', "activity"]], left_on="ID", right_on="sim_id")
if_df.drop(['sim_id'], axis=1, inplace=True)
# print if_df[:5]['HMScore']
if compounds is not None:
    if_df = if_df.loc[if_df['name'].isin(compounds)]
    if len(if_df) < len(compounds):
        compounds = if_df['name']

scores1 = [score for score in if_df.columns if score not in ["ID", "activity", "name"]]
if args.debug:
    print if_df.shape, if_sfs_df.shape
roc_values1, ef_dictionary1 = plot_all_sfs(if_df, scores1, 'if', args.ef_thresholds, args.general_title)
write_roc_files(roc_values1, args.output_prefix + "_if")
# print ef_dictionary1.keys()
write_ef_files(ef_dictionary1, args.output_prefix + "_if")
if args.debug:
    print 4, roc_values1['HMScore'][0][:5]
    print 5, roc_values1['HMScore'][3][:5]
if args.previous_scoring_functions_file:
    ini_sfs_df = pd.read_csv(args.previous_scoring_functions_file)
    if ids2select_df0 is not None:
        ini_sfs_df = ini_sfs_df.loc[ini_sfs_df['ID'].isin(ids2select_df0['ID'])]
    if "Exp_energy" in ini_sfs_df.columns:
        ini_df = ini_sfs_df.merge(name_sim_activity_df[['sim_id', 'name']], left_on="ID", right_on="sim_id")
        ini_df.rename(columns={"Exp_energy": "activity"}, inplace=True)
    else:
        ini_df = ini_sfs_df.merge(name_sim_activity_df[['sim_id', 'name', "activity"]], left_on="ID", right_on="sim_id")
    ini_df.drop(['sim_id'], axis=1, inplace=True)
    if compounds is not None:
        ini_df = ini_df.loc[ini_df['name'].isin(compounds)]
    if args.debug:
        print ini_df.shape
    scores0 = [score for score in ini_df.columns if score not in ["ID", "activity", "name"]]
    roc_values0, ef_dictionary0 = plot_all_sfs(ini_df, scores0, 'init', args.ef_thresholds, args.general_title)
    if args.debug:
        print 9, roc_values0['HMScore'][3][:5]
    write_roc_files(roc_values0, args.output_prefix + "_init")
    write_ef_files(ef_dictionary0, args.output_prefix + "_init")
    common_scores = set.intersection(set(scores1), set(scores0))
    if args.debug:
        print 3, common_scores
    for score in sorted(common_scores):
        title = "{0}  (SF:{1} comparison) ".format(args.general_title, score)
        fig_tmp, ax_tmp = pl.subplots()
        sim_ids1, activities1, names1, fpr1, tpr1, auc_val1 = roc_values1[score]
        if args.debug:
            print 6, sim_ids1[:5]
            print 7, fpr1[:5]
        sim_ids0, activities0, names0, fpr0, tpr0, auc_val0 = roc_values0[score]
        if args.debug:
            print 8, fpr0[:5]
        ax_tmp.plot(fpr0, tpr0, label='initial {0} ROC curve (area = {1:0.2f})'.format(score, auc_val0), c='b')
        ax_tmp.plot(fpr1, tpr1, label='i.f. {0} ROC curve (area = {1:0.2f})'.format(score, auc_val1), c='r')
        print "{1} delta auc: {0:6.4}".format( auc_val1 - auc_val0, score)
        ax_tmp.plot([0, 1], [0, 1], color='navy', linestyle='--', label="Random selection")
        ax_tmp.legend(loc="lower right")
        ax_tmp.set_ylim([0, 1])
        ax_tmp.set_xlim([0, 1])
        ax_tmp.set_title(title)
        ax_tmp.xaxis.set_label_text("FPR")
        ax_tmp.yaxis.set_label_text("TPR")
        pl.show()
        fig_tmp.savefig(args.output_prefix + "_{}_comparison.png".format(score), bbox_inches='tight')
    common_scores.update({"Exp_energy"})
    delta_ef = {score: {th: ef_dictionary1[score][th] - ef_dictionary0[score][th] for th in args.ef_thresholds}
                for score in common_scores}
    # print delta_ef.keys()
    write_ef_files(delta_ef, args.output_prefix + "_delta")
