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
Loads interactome (human protein-protein interaction network, PPIN) from
STRING v9.1 above a confidence threshold of 700 (out of 1000), either as
a pickled file or from a MySQL database. Processes network using NetworkX
to find the largest connected island and returns this graph (H).

"""

import os
import sys
import numpy as np
# import MySQLdb
import networkx as nx
import cPickle as pickle

def load_network():
    try:
        print "Loading STRING data from file"
        data = pickle.load(open("string700_data.p", "rb") )
    except:
        # Note: data must be in the form [('A','B'), ('B','C'), ('A','D'), ('B', 'D'), ('E', 'F')]
        print "Loading STRING data from database"
        
        # dbConnection = MySQLdb.connect(read_default_file='~/.my.cnf', db = 'pathway_string')
        # dbCursor = dbConnection.cursor()

        # SQL = 'SELECT distinct protein1,protein2,combined_score FROM protein_links_detailed_human WHERE combined_score >= 700'

        # list_length = dbCursor.execute(SQL)
        # print "Number of rows loaded:",list_length

        # data = []
        # for i in range(list_length):
        #     result = dbCursor.fetchone()
        #     data.append((result[0],result[1]))

        # pickle.dump(data, open("string700_data.p", "wb"))

        # dbCursor.close()
        # dbConnection.commit()
        # dbConnection.close()

    G = nx.Graph()
    G.add_edges_from(data)
    #print len(G)

    print "Using NetworkX to find the largest island"
    H = max(nx.connected_component_subgraphs(G), key=len)

    ppi = dict()
    node_list = []
    for a, b in H.edges(): 
        if not a in ppi:
            ppi[a] = set()
            node_list.append(a)
        if not b in ppi:
            ppi[b] = set()
            node_list.append(b)
        
        ppi[a].add(b)
        ppi[b].add(a)

    #sorting all genes (to align with the adjacency matrix)
    numNodes = len(ppi)
    print "Number of nodes:",len(node_list)

    return ppi, node_list, H