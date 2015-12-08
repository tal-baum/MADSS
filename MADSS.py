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
Main MADSS script. Run specifying phenotype of interest, for example:
python MADSS.py MI
to investigate acute myocardial infarction.

"""

import os
import sys
import csv
import re
import numpy as np
import networkx as nx
import cPickle as pickle
import matplotlib.pyplot as plt
from collections import defaultdict

# Locate necessary MADSS libraries
wd = os.getcwd()
sys.path.insert(0, wd+'/madss_libraries/')

import madss_credits
madss_credits.print_intro()

if len(sys.argv) < 2:
    sys.exit('Error: please specify an adverse event (e.g. "python MADSS.py MI"). Exiting.')

# Define phenotype to be evaluated
ADVERSE_EVENT = str(sys.argv[1]).upper() # "MI", "Gastro", "Liver", "Kidney", "LQTS"

# If only running one connectivity function
if len(sys.argv) > 2:
    metrics = [str(sys.argv[2])]
    print "Running %s only without scoring drugs" %metrics[0].upper()
else:
    metrics = ['mfpt', 'bc', 'sn', 'isp']
    print "Running all connectivity functions for %s" %ADVERSE_EVENT

# ------ Load interactome ------
import madss_interactome
ppi, node_list, H = madss_interactome.load_network()

# ------ Assign seeds ------
if ADVERSE_EVENT == 'MI':
    seeds = ['ENSP00000358301','ENSP00000206249','ENSP00000366307','ENSP00000304236','ENSP00000287936','ENSP00000290866','ENSP00000392858','ENSP00000361125','ENSP00000258743','ENSP00000229135']
elif ADVERSE_EVENT == 'Gastro':
    seeds = ['ENSP00000287641','ENSP00000356438','ENSP00000366506','ENSP00000334145','ENSP00000354612','ENSP00000308541','ENSP00000218099','ENSP00000039007']
elif ADVERSE_EVENT == 'Liver':
    seeds = ['ENSP00000258743','ENSP00000295897','ENSP00000263341','ENSP00000308541','ENSP00000264708','ENSP00000356438','ENSP00000206249','ENSP00000412237','ENSP00000353483','ENSP00000356771','ENSP00000334145','ENSP00000351671','ENSP00000356694','ENSP00000216117','ENSP00000253408']
elif ADVERSE_EVENT == 'Kidney':
    seeds = ['ENSP00000360997', 'ENSP00000295897', 'ENSP00000265132', 'ENSP00000340858', 'ENSP00000379204', 'ENSP00000328173', 'ENSP00000241052', 'ENSP00000356399', 'ENSP00000265689', 'ENSP00000384400', 'ENSP00000360541', 'ENSP00000366124', 'ENSP00000377192', 'ENSP00000378408', 'ENSP00000338082', 'ENSP00000222390', 'ENSP00000287936', 'ENSP00000348170', 'ENSP00000298556', 'ENSP00000280357', 'ENSP00000265023', 'ENSP00000277480', 'ENSP00000229794', 'ENSP00000352835', 'ENSP00000365663', 'ENSP00000355759', 'ENSP00000265970', 'ENSP00000364252', 'ENSP00000234347', 'ENSP00000164139', 'ENSP00000360519', 'ENSP00000272190', 'ENSP00000367102', 'ENSP00000264938', 'ENSP00000368727']
elif ADVERSE_EVENT == 'LQTS':
    seeds = ['ENSP00000262186','ENSP00000155840','ENSP00000337255','ENSP00000290310','ENSP00000266483','ENSP00000348573','ENSP00000217381','ENSP00000328968','ENSP00000243457','ENSP00000266376','ENSP00000341940','ENSP00000349588','ENSP00000322460']

print 'Number of seeds: ', len(seeds)

# ------ Prepare to calculate node connectivities ------
Sj_scores = defaultdict(dict)

if not os.path.exists('results'):
    os.makedirs('results')

def read_Sj_file(METRIC):
    results_file = open('results/%s_Sj_%s.txt' %(METRIC, ADVERSE_EVENT),'r')
    
    Sj_from_file = dict()
    neighborhood_count = 0
    for line in results_file:
        newline = re.split('\t|,',line)  # 0:Sj  1:protein_id
        Sj = float(newline[0])
        protein_id = newline[1].strip()
        Sj_from_file[protein_id] = Sj
        if Sj > 0:
            neighborhood_count += 1
    results_file.close()
    print METRIC.upper(), "neighborhood size:",neighborhood_count

    return Sj_from_file

# ------ Run connectivity functions ------
# Mean first passage time (MFPT)
if os.path.isfile('results/mfpt_Sj_%s.txt' %ADVERSE_EVENT) and 'mfpt' in metrics:
    print "\nLoading mean first passage time from results file"
    Sj_scores['mfpt'] = read_Sj_file('mfpt')
elif 'mfpt' in metrics:
    import madss_mfpt
    print "\nRunning mean first passage time"
    Sj_scores['mfpt'] = madss_mfpt.calc_mfpt_Sj(ppi, node_list, H, ADVERSE_EVENT, seeds)

# Betweenness centrality (BC)
if os.path.isfile('results/bc_Sj_%s.txt' %ADVERSE_EVENT) and 'bc' in metrics:
    print "\nLoading betweenness centrality from results file"
    Sj_scores['bc'] = read_Sj_file('bc')
elif 'bc' in metrics:
    import madss_bc
    print "\nRunning betweenness centrality"
    Sj_scores['bc'] = madss_bc.calc_bc_Sj(ppi, node_list, H, ADVERSE_EVENT, seeds)

# Shared neighbors (SN)
if os.path.isfile('results/sn_Sj_%s.txt' %ADVERSE_EVENT) and 'sn' in metrics:
    print "\nLoading shared neighbors from results file"
    Sj_scores['sn'] = read_Sj_file('sn')
elif 'sn' in metrics:
    import madss_sn
    print "\nRunning shared neighbors"
    Sj_scores['sn'] = madss_sn.calc_sn_Sj(ppi, node_list, H, ADVERSE_EVENT, seeds)

# Inverse shortest path (ISP)
if os.path.isfile('results/isp_Sj_%s.txt' %ADVERSE_EVENT) and 'isp' in metrics:
    print "\nLoading inverse shortest path from results file"
    Sj_scores['isp'] = read_Sj_file('isp')
elif 'isp' in metrics:
    import madss_isp
    print "\nRunning inverse shortest path"
    Sj_scores['isp'] = madss_isp.calc_isp_Sj(ppi, node_list, H, ADVERSE_EVENT, seeds)


# ------ Assign each drug to most highly connected target ------
if metrics != ['mfpt', 'bc', 'sn', 'isp']: # (i.e., if only seeking to calculate one function's scores)
    sys.exit("Only calcualated connectivity for %s. Exiting." %metrics[0].upper())

import madss_scoring
print '\n------------'
print "Scoring drugs"
# Open gold standard for adverse event (gt = ground truth)
gt_drugs, id2gt = madss_scoring.open_gold_standard(ADVERSE_EVENT)

# Collect drug targets from DrugBank
drug_list, drugbank_targets, drugbank2name, ensembl2gene = madss_scoring.get_drugbank_targets(ADVERSE_EVENT, gt_drugs, node_list)

# Score drugs
madss_scoring.score_drugs(ADVERSE_EVENT, drug_list, drugbank_targets, drugbank2name, ensembl2gene, id2gt, Sj_scores)


# ------ Train classifier and generate ROC plot ------
import madss_classifier
print '\n------------'
print "Training classifier"
madss_classifier.generate_ROC(ADVERSE_EVENT) #, plot_raw=True)

