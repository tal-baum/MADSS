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
For a given adverse event, this script loads a gold standard (list of drugs
that do/ do not cause the AE) as well as DrugBank targets for those drugs.
For each connectivity function, it then assigns each drug the score of the
most highly connected target as scored by the given connectivity function.
These results are stored as a matrix in /scores.

For any AE other than MI, a MySQL database of DrugBank will have to be
set up to obtain drug targets.

"""

import os
import sys
import csv
import numpy as np
import scipy as sp
import cPickle as pickle
from collections import defaultdict
import MySQLdb

# Collect AE-specific drugs
def open_gold_standard(ADVERSE_EVENT):
    print "\nLoading gold standard for %s" %ADVERSE_EVENT.upper()
    gt_drugs = []

    f = open('gold_standards/gold_standard_%s.csv' %ADVERSE_EVENT,'rU')
    reader = csv.reader(f)
    reader.next()

    id2gt = dict() # gt = ground truth
    for drug,drugbank_id,gt in reader:
        if drugbank_id != '':
            gt_drugs.append(drugbank_id)
            id2gt[drugbank_id] = int(gt)
    f.close()

    return gt_drugs, id2gt


# Gather DrugBank targets (e.g. targets, enzymes, transporters)
def get_drugbank_targets(ADVERSE_EVENT, gt_drugs, node_list):
    # Initialize dictionaries
    drug_list = []
    drugbank_targets = dict()
    drugbank2name = dict()
    ensembl2gene = dict()

    try:
        drugbank_targets = pickle.load(open("stored_vals/drugbank_targets_%s.p" %ADVERSE_EVENT, "rb"))
        for drugbank_id in gt_drugs:
            if drugbank_id in drugbank_targets:
                if not drugbank_id in drug_list:
                    drug_list.append(drugbank_id)

        drugbank2name = pickle.load(open("stored_vals/drugbank2name_%s.p" %ADVERSE_EVENT, "rb"))
        ensembl2gene = pickle.load(open("stored_vals/ensembl2gene_%s.p" %ADVERSE_EVENT, "rb"))
        print "DrugBank targets loaded from pickled file"

    except:
        print "Loading DrugBank targets from database"
        dbConnection = MySQLdb.connect(read_default_file='~/.my.cnf', db = 'compound_drugbank03')
        dbCursor = dbConnection.cursor()

        SQL = 'select distinct drugbank_id,drugname from drug_categories order by drugname'
        num_drugs = dbCursor.execute(SQL)
        drug_results = dbCursor.fetchall()

        for drugbank_id,drugname in drug_results:
            if drugbank_id in gt_drugs:
                if not drugbank_id in drug_list:
                    drug_list.append(drugbank_id)
                if not drugbank_id in drugbank_targets:
                    drugbank_targets[drugbank_id] = set()
                drugbank2name[drugbank_id] = drugname
        print "%d drugs with targets loaded from DrugBank" %len(drug_list)

        for i,drug in enumerate(drug_list):
            percent_complete = float(i+1) / len(drug_list) * 100
            print "%.2f" %percent_complete,
            
            # Targets
            SQL = '''SELECT DISTINCT drug_allTargets_human.drugbank_id,idmapping.id, partner_protein.Gene_Name
                     FROM compound_drugbank03.drug_allTargets_human,compound_drugbank03.partner_protein,protein_uniprot.idmapping
                     WHERE idmapping.id_type = 'Ensembl_PRO'
                     AND idmapping.uniprot = partner_protein.UniProt_ID
                     AND partner_protein.Partner_ID = drug_allTargets_human.Target_Partner
                     AND drug_allTargets_human.drugbank_id = '%s' ''' %drug
            dbCursor.execute(SQL)
            target_results = dbCursor.fetchall()
            
            for drugbank_id,target,name in target_results:
                if target in node_list:
                    drugbank_targets[drugbank_id].add(target)
                    ensembl2gene[target] = name
            
            # Enzymes
            SQL = '''SELECT DISTINCT drug_allEnzymes_human.drugbank_id,idmapping.id, partner_protein.Gene_Name
                     FROM compound_drugbank03.drug_allEnzymes_human,compound_drugbank03.partner_protein,protein_uniprot.idmapping
                     WHERE idmapping.id_type = 'Ensembl_PRO'
                     AND idmapping.uniprot = partner_protein.UniProt_ID
                     AND partner_protein.Partner_ID = drug_allEnzymes_human.Enzyme_Partner
                     AND drug_allEnzymes_human.drugbank_id = '%s' ''' %drug
            
            dbCursor.execute(SQL)
            enzyme_results = dbCursor.fetchall()
            
            for drugbank_id,enzyme,name in enzyme_results:
                if enzyme in node_list:
                    drugbank_targets[drugbank_id].add(enzyme)
                    ensembl2gene[enzyme] = name
            
            # Transporters
            SQL = '''SELECT DISTINCT drug_allTransporters_human.drugbank_id,idmapping.id, partner_protein.Gene_Name
                     FROM compound_drugbank03.drug_allTransporters_human,compound_drugbank03.partner_protein,protein_uniprot.idmapping
                     WHERE idmapping.id_type = 'Ensembl_PRO'
                     AND idmapping.uniprot = partner_protein.UniProt_ID
                     AND partner_protein.Partner_ID = drug_allTransporters_human.Transporter_Partner
                     AND drug_allTransporters_human.drugbank_id = '%s' ''' %drug
            
            dbCursor.execute(SQL)
            transporter_results = dbCursor.fetchall()
            
            for drugbank_id,transporter,name in transporter_results:
                if transporter in node_list:
                    drugbank_targets[drugbank_id].add(transporter)
                    ensembl2gene[transporter] = name

        dbCursor.close()
        dbConnection.commit()
        dbConnection.close()

        pickle.dump(drugbank_targets, open("stored_vals/drugbank_targets_%s.p" %ADVERSE_EVENT, "wb"))
        pickle.dump(drugbank2name, open("stored_vals/drugbank2name_%s.p" %ADVERSE_EVENT, "wb"))
        pickle.dump(ensembl2gene, open("stored_vals/ensembl2gene_%s.p" %ADVERSE_EVENT, "wb"))
        print '\r'

    return drug_list, drugbank_targets, drugbank2name, ensembl2gene


def find_best_target(Sj_list, drugbank_targets, drugbank2name, drugbank_id):
    best_score = float()
    best_target = ''

    set_init_best = True    
    
    if drugbank_targets[drugbank_id] == set([]):
        # print "\t",drugbank_id, drugbank2name[drugbank_id], "no targets found"
        return 'x','x'
    
    for target in drugbank_targets[drugbank_id]:        
        if set_init_best == True:
            best_score = float(Sj_list[target])
            best_target = target
            set_init_best = False

        if set_init_best == False:
            if float(Sj_list[target]) > best_score:
                best_score = float(Sj_list[target])
                best_target = target
    
    return best_score,best_target


def score_drugs(ADVERSE_EVENT, drug_list, drugbank_targets, drugbank2name, ensembl2gene, id2gt, Sj_scores):
    # Initialize score array
    all_scores = np.array([['drugbank_id']])
    for drugbank_id in drug_list:
        all_scores = np.vstack((all_scores,drugbank_id))

    # Add training labels
    drug_labels = np.array([['causes_AE']])
    for drugbank_id in all_scores[1:]:
        drugbank_id = drugbank_id[0]
        if drugbank_id in id2gt:
            drug_labels = np.vstack((drug_labels,id2gt[drugbank_id]))
        else:
            drug_labels = np.vstack((drug_labels,''))
    all_scores = np.hstack((all_scores,drug_labels))

    # Iterate over each connectivity function
    METRICS = ['mfpt', 'bc', 'sn', 'isp']
    for METRIC in METRICS:
        print 'Scoring', METRIC.upper()
        db_col = np.array([[METRIC.upper()]])
        target_col = np.array(['%s Target' %METRIC.upper()])
        
        Sj_list = Sj_scores[METRIC]
        
        best_scores = dict()
        best_targets = dict()

        for drugbank_id in drug_list:
            best_scores[drugbank_id], best_targets[drugbank_id] = find_best_target(Sj_list, drugbank_targets, drugbank2name, drugbank_id)
            
            db_col = np.vstack((db_col,str(best_scores[drugbank_id])))
            temp_target = str(best_targets[drugbank_id])
            target_col = np.vstack((target_col,ensembl2gene[temp_target] if temp_target in ensembl2gene else 'x'))
        
        all_scores = np.hstack((all_scores,db_col))
        all_scores = np.hstack((all_scores,target_col))

    # Save results
    if not os.path.exists('scores'):
        os.makedirs('scores')

    print "Writing results to /scores"
    outf = open('scores/%s_drug_scores.csv' %(ADVERSE_EVENT),'w')
    writer = csv.writer(outf)

    for drugbank_id,gt,mfpt,mtarget,bc,btarget,sn,starget,isp,itarget in all_scores:
        if mfpt != 'x':
            writer.writerow([drugbank2name[drugbank_id] if drugbank_id in drugbank2name else "drug_name",drugbank_id,gt,mfpt,mtarget,bc,btarget,sn,starget,isp,itarget])
          
    outf.close()