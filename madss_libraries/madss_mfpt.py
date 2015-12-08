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
Script to compute the mean first passage time as described by Kemeny and
Snell Finite Markov Chains (1960) pp 78-.

"""

import os
import sys
import csv
import numpy as np
import networkx as nx

def restart_line():
    sys.stdout.write('\r')
    sys.stdout.flush()

def calc_mfpt_Sj(ppi, node_list, H, ADVERSE_EVENT, seeds):
    #sorting all genes (to align with the adjacency matrix)
    numNodes = len(ppi)
    sorted_reachable = sorted(ppi)

    seeds = [gene for gene in seeds if gene in ppi]
    seedIndices = [sorted_reachable.index(gene) for gene in seeds]
    complementIndices = list(set(range(len(sorted_reachable))) - set(seedIndices))

    #calculate mfpt matrix
    A = np.zeros(shape=(numNodes,numNodes))
    for i,a in enumerate(sorted_reachable):
        for j,b in enumerate(sorted_reachable):
            if b in ppi[a]:
                A[i,j] = 1
                A[j,i] = 1

    if os.path.isfile('stored_vals/mfptMatrix.npy'):
        print "Loading mfpt matrix...",
        sys.stdout.flush()
        mfptMatrix = np.load('stored_vals/mfptMatrix.npy')
    else:
        # Note, for this to work (i.e. K be the correct shape) A must be a 2D array, not a matrix.
        print "Forming mfpt matrix...",
        K = np.diag(A.sum(1))
        B = np.linalg.inv((K-A + np.dot(K, np.ones(A.shape))))
        mfptMatrix = A.sum()*( np.dot( np.ones(A.shape),np.diag(np.diag(B)) ) - B )
        if not os.path.exists('stored_vals'):
            os.makedirs('stored_vals')
        np.save('stored_vals/mfptMatrix.npy', mfptMatrix)
    print 'done'

    Sjs = []
        
    for (geneIndex, geneName) in enumerate(sorted_reachable):
        sys.stdout.write('Calculating Sj for node %s (Index %d): ' %(geneName,geneIndex)),

        denominator = sum(mfptMatrix[:,geneIndex]) / numNodes

        inSet = sum(mfptMatrix[seedIndices,geneIndex])/len(seeds)

        complement = sum(mfptMatrix[complementIndices,geneIndex])/len(complementIndices)

        Sj = (complement - inSet) / denominator
        Sjs.append((Sj,sorted_reachable[geneIndex]))

        sys.stdout.write("%s" %Sj)
        sys.stdout.flush()
        restart_line()

    sorted_Sjs = sorted(Sjs)
    sorted_Sjs.reverse()

    print 'Sjs calculated.'
    # save Sjs for later use
    sj_file_handler = open('results/mfpt_Sj_%s.txt' %(ADVERSE_EVENT),'w')
    sj_writer = csv.writer(sj_file_handler)

    Sj_dict = dict()
    neighborhood = []
    for i, (sj, geneName) in enumerate(sorted_Sjs):
        sj_writer.writerow([sj, geneName])
        Sj_dict[geneName] = sj

        if sj > 0:
            if not geneName in neighborhood:
                neighborhood.append(geneName)

    #print "Number of Sj values:",len(sorted_Sjs)
    print "Neighborhood size:",len(neighborhood)

    sj_file_handler.close()

    return Sj_dict