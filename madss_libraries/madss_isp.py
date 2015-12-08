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
Script to compute adapted inverse shortest path connectivity function.

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

def calc_Sj_stored(H,j,seeds_left,isp_alls,complement):
    isp_seeds = 0
    isp_all   = isp_alls[j]
            
    for i in seeds_left:
        isp_temp = 0
        if i != j:
            isp_temp = 1 / float(nx.shortest_path_length(H,i,j))
            isp_seeds += isp_temp
    isp_comp = isp_all - isp_seeds
    
    numSeeds = isp_seeds / float(len(seeds_left))
    
    numComp = isp_comp / float(len(complement))
    
    denominator = isp_all / float(len(H))
    
    Sj = (numSeeds - numComp) / float(denominator)
    
    #print j,Sj,isp_seeds,isp_comp,isp_all
    
    return Sj
#----------------------------------------------------------------------------------------

def calc_isp_Sj(ppi, node_list, H, ADVERSE_EVENT, seeds):
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

    # Try to load stored vals, otherwise create
    try:
        isp_alls = pickle.load(open("stored_vals/string700_isp_alls.p", "rb") )
        print "stored values loaded"
    except:
        # Calculate isp_all
        print "building ISP dictionary...",
        sys.stdout.flush()
        isp_alls = dict()

        node_count = 0
        for j in H:
            isp_all = 0

            for i in H:
                isp_temp = 0
                if i != j:
                    isp_temp = 1 / float(nx.shortest_path_length(H,i,j))
                    isp_all += isp_temp
            
            isp_alls[j] = isp_all
            
            node_count += 1
            if node_count % 500 == 0:
                print node_count,
                sys.stdout.flush()
        pickle.dump(isp_alls, open("stored_vals/string700_isp_alls.p", "wb"))
        print 'done',len(isp_alls),'\n'


    Sj_dict = dict()

    #----- Calc Sj ---------
    for i, j in enumerate(H):
        sys.stdout.write('Calculating Sj for node: %s (Index %d): ' %(j,i)),

        Sj_dict[j] = calc_Sj_stored(H,j,seeds_left,isp_alls,complement)

        sys.stdout.write('%s' %Sj_dict[j])
        sys.stdout.flush()
        restart_line()

    Sj_file = open('results/isp_Sj_%s.txt' %ADVERSE_EVENT, 'w')
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