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
Script to compute adapted shared neighbors connectivity function.

"""

import os
import sys
import numpy as np
import networkx as nx
import cPickle as pickle

def restart_line():
    sys.stdout.write('\r')
    sys.stdout.flush()

#----------------------------------------------------------------------------------------

def calc_Tc(H,i,j,neighbors):
    Tc = 0
    
    n_i = neighbors[i]
    n_j = neighbors[j]
    
    set_union = n_i | n_j         # union

    set_intersection = n_i & n_j  # intersection
    
    Tc = len(set_intersection) / float(len(set_union))
    
    return Tc

#----------------------------------------------------------------------------------------

def calc_Sj_stored(H,j,seeds_left,Tc_alls, neighbors, complement):
    Tc_seeds = 0
    #Tc_comp  = 0
    Tc_all   = Tc_alls[j]
    
    for i in seeds_left:
        Tc_temp = 0
        if i != j:
            Tc_temp = calc_Tc(H,i,j,neighbors)
            Tc_seeds += Tc_temp
    Tc_comp = Tc_all - Tc_seeds
    
    numSeeds = Tc_seeds / float(len(seeds_left))
    
    numComp = Tc_comp / float(len(complement))
    
    denominator = Tc_all / float(len(H))
    
    Sj = (numSeeds - numComp) / float(denominator)
    
    return Sj
#----------------------------------------------------------------------------------------


def calc_sn_Sj(ppi, node_list, H, ADVERSE_EVENT, seeds):
    seeds_left = []
    for node in node_list:
        for seed in seeds:
            if node == seed:
                if not node in seeds_left:
                    seeds_left.append(node)
    # print "Number of seeds left:",len(seeds_left)

    complement = []
    for node in node_list:
        if not node in seeds_left:
            if not node in complement:
                complement.append(node)

    # Try to load stored neighbors and Tcs
    try:
        neighbors = pickle.load(open("stored_vals/string700_neighbors.p", "rb") )
        Tc_alls = pickle.load(open("stored_vals/string700_Tc_alls.p" , "rb") )
        print "stored values loaded"

    except:
        #### Build neighbor dictionary
        print "building neighbor dictionary...",
        sys.stdout.flush()

        neighbors = dict()
        for node in H:
            neighbors[node] = set(nx.all_neighbors(H,node))
        pickle.dump(neighbors, open("stored_vals/string700_neighbors.p" , "wb"))
        print "done"
        sys.stdout.flush()

        #### Calculate Tc_all
        print "building Tc_all dictionary...",
        sys.stdout.flush()
        Tc_alls = dict()

        node_count = 0
        for j in H:
            Tc_all = 0

            for i in H:
                Tc_temp = 0
                if i != j:
                    Tc_temp = calc_Tc(H,i,j,neighbors)
                    Tc_all += Tc_temp
            
            Tc_alls[j] = Tc_all
            
            node_count += 1
            if node_count % 500 == 0:
                print node_count,
                sys.stdout.flush()
        pickle.dump(Tc_alls, open("stored_vals/string700_Tc_alls.p" , "wb"))
        print 'done',len(Tc_alls),'\n'
    ####

    Sj_dict = dict()

    #----- Calc Sj ---------
    for i, j in enumerate(H):
        sys.stdout.write('Calculating Sj for node: %s (Index %d): ' %(j,i)),

        Sj_dict[j] = calc_Sj_stored(H,j,seeds_left,Tc_alls, neighbors, complement)
        
        sys.stdout.write('%s' %Sj_dict[j])
        sys.stdout.flush()
        restart_line()

    Sj_file = open('results/sn_Sj_%s.txt' %ADVERSE_EVENT, 'w')
    neighborhood = []

    for j in H:
        Sj_file.write("%s\t%s\n" %(Sj_dict[j],j))
        if Sj_dict[j] > 0:
            if not j in neighborhood:
                neighborhood.append(j)

    #print "Number of Sj values:",len(Sj_dict)
    print "\nNeighborhood size:",len(neighborhood)

    Sj_file.close()
    return Sj_dict
