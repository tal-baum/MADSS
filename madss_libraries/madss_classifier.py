"""
MADSS v1.1, Updated December 7, 2015
Citation:
Lorberbaum, T., Nasir, M., Keiser, M., Vilar, S., Hripcsak, G. and Tatonetti, N. (2015),
Systems Pharmacology Augments Drug Safety Surveillance. Clinical Pharmacology & Therapeutics,
97: 151-158. doi: 10.1002/cpt.2
http://onlinelibrary.wiley.com/doi/10.1002/cpt.2/abstract

Copyright (C) 2014-2015, Tatonetti Lab
Tal Lorberbaum <tal.lorberbaum@columbia.edu>
Nicholas P. Tatonetti <nick.tatonetti@columbia.edu>
All rights reserved.

This script is part of MADSS which is released under a CC BY-NC-SA 4.0 license.
For full license details see LICENSE.txt or go to:
http://creativecommons.org/licenses/by-nc-sa/4.0/

------------------------------------------------------------------------
Script to train random forest classifier using connectivity scores as
features. Exports classifier probabilities and generates ROC plot.

"""

import os
import sys
import csv
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn import ensemble
from collections import defaultdict

def generate_ROC(ADVERSE_EVENT, plot_raw=False):
    f = open('scores/%s_drug_scores.csv' %(ADVERSE_EVENT),'rU')
    reader = csv.reader(f)
    reader.next()

    drug_list = []
    drugbank2name = dict()

    drug2score = defaultdict(dict)

    data_sets = []
    for drug_name,drugbank_id,training,mfpt,mfpt_target,bc,bc_target,sn,sn_target,isp,isp_target in reader:
        to_map = [training,mfpt,bc,sn,isp]
        data_sets.append(map(float,to_map))
        
        drug_list.append(drugbank_id)
        drugbank2name[drugbank_id] = drug_name
        
        drug2score['mfpt'][drugbank_id] = float(mfpt)
        drug2score['bc'][drugbank_id] = float(bc)
        drug2score['sn'][drugbank_id] = float(sn)
        drug2score['isp'][drugbank_id] = float(isp)
    f.close()

    examples = np.array([row[1:] for row in data_sets])
    labels = np.array([row[0] for row in data_sets])


    # Initialize random forest classifier
    clf_init = lambda: ensemble.RandomForestClassifier(n_estimators=100, max_features=1, max_depth=None, oob_score=True)

    # Random Forest
    clf = clf_init()
    clf.fit(examples, labels)

    rf_auroc = metrics.roc_auc_score(labels, clf.oob_decision_function_[:,1])
    fpr, tpr, thresholds = metrics.roc_curve(labels, clf.oob_decision_function_[:,1])
    rf_roc_curve = [fpr, tpr]
    rf_auroc = metrics.auc(fpr, tpr)

    print "Random Forest AUROC:", rf_auroc


    # Export OOB (out-of-bag) probabilities
    print "Exporting classifier probabilities to /probabilities"
    if not os.path.exists('probabilities'):
        os.makedirs('probabilities')
    outf = open('probabilities/%s_classifier_probs.csv' %ADVERSE_EVENT, 'w')
    writer = csv.writer(outf)
    writer.writerow(['drug_name', 'drugbank_id', 'ground_truth', 'SubNet'])
    for drugbank_id , val, prob in zip(drug_list,labels,clf.oob_decision_function_[:,1]):
        writer.writerow([drugbank2name[drugbank_id],drugbank_id,int(val),prob])
    outf.close()


    # Raw scores
    roc_curves = dict()
    aurocs = dict()

    all_metrics = ['mfpt', 'bc', 'sn', 'isp']

    scores = defaultdict(list)
    for drug in drug_list:
        for metric in all_metrics:
            scores[metric].append(drug2score[metric][drug])

    for metric in all_metrics:
        fpr, tpr, thresholds = metrics.roc_curve(labels, scores[metric])
        roc_curves[metric] = [fpr, tpr]
        aurocs[metric] = metrics.auc(fpr, tpr)
        #print metric,aurocs[metric]


    # Generate ROC plot
    font_large = {'family' : 'sans-serif',
                  'weight' : 'normal',
                  'size'   : 14}

    font_small = {'family': 'sans-serif',
                  'weight': 'normal',
                  'size': 11}

    colors = dict()
    colors['rf'] = '#ca0a37'
    colors['mfpt'] = 'k'
    colors['bc'] = 'blue'
    colors['sn'] = 'green'
    colors['isp'] = 'orange'

    fig = plt.figure(figsize=(5,5))

    plt.rc('font', **font_large)

    ax = plt.axes()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    ax.xaxis.set_ticks_position('bottom')
    ax.xaxis.set_label_position('bottom')
    ax.xaxis.grid(True)

    ax.yaxis.set_ticks_position('left')
    ax.yaxis.grid(True)

    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')

    p0 = plt.plot( rf_roc_curve[0], rf_roc_curve[1], lw=2, color=colors['rf'], label="Random Forest (%0.2f)" %(rf_auroc))

    if plot_raw == True:
        for metric in all_metrics:
           p0 = plt.plot( roc_curves[metric][0], roc_curves[metric][1], lw=2, color=colors[metric], label="%s (%0.2f)" %(metric.upper(),aurocs[metric]))

    plt.plot([0,1], [0,1], '--', linewidth=2, alpha=0.5, color='grey')
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    plt.legend(loc='lower right', frameon=False, prop={'size':11} )

    plt.title('MADSS Random Forest: %s' %ADVERSE_EVENT)

    if not os.path.exists('figs'):
        os.makedirs('figs')

    filename = '%s_MADSS_ROC' %ADVERSE_EVENT
    if plot_raw == True:
        filename += '+raw'
    plt.savefig('figs/%s.pdf' %filename)
