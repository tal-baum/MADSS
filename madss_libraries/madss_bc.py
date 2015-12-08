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
Script to compute adapted betweenness centrality connectivity function.

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
'''Adapted from NetworkX betweenness_centrality function
https://networkx.github.io/documentation/latest/_modules/networkx/algorithms/centrality/betweenness.html#betweenness_centrality

Copyright (C) 2004-2012, NetworkX Developers
Aric Hagberg <hagberg@lanl.gov>
Dan Schult <dschult@colgate.edu>
Pieter Swart <swart@lanl.gov>
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

  * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.

  * Redistributions in binary form must reproduce the above
    copyright notice, this list of conditions and the following
    disclaimer in the documentation and/or other materials provided
    with the distribution.

  * Neither the name of the NetworkX Developers nor the names of its
    contributors may be used to endorse or promote products derived
    from this software without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

def sub_betweenness_centrality(G, subset):
    print "Calculating betweenness for subset of length", len(subset)
    betweenness=dict.fromkeys(G,0.0) # b[v]=0 for v in G
    
    for i,node in enumerate(subset):
        sys.stdout.write('%d: s=%s' %(i,node))
        sys.stdout.flush()
        restart_line()
        # single source shortest paths
        S,P,sigma=_single_source_shortest_path_basic(G,node)
        # accumulation
        betweenness=_accumulate_basic(betweenness,S,P,sigma,node)

    return betweenness
    
def _single_source_shortest_path_basic(G,s):
    S=[]
    P={}
    for v in G:
        P[v]=[]
    sigma=dict.fromkeys(G,0.0)    # sigma[v]=0 for v in G
    D={}
    sigma[s]=1.0
    D[s]=0
    Q=[s]
    while Q:   # use BFS to find shortest paths
        v=Q.pop(0)
        S.append(v)
        Dv=D[v]
        sigmav=sigma[v]
        for w in G[v]:
            if w not in D:
                Q.append(w)
                D[w]=Dv+1
            if D[w]==Dv+1:   # this is a shortest path, count paths
                sigma[w] += sigmav
                P[w].append(v) # predecessors 
    return S,P,sigma

def _accumulate_basic(betweenness,S,P,sigma,s):
    delta=dict.fromkeys(S,0)
    while S:
        w=S.pop()
        coeff=(1.0+delta[w])/sigma[w]
        for v in P[w]:
            delta[v] += sigma[v]*coeff
        if w != s:
            betweenness[w]+=delta[w]
    return betweenness

#----------------------------------------------------------------------------------------

def calc_bc_Sj(ppi, node_list, H, ADVERSE_EVENT, seeds):
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

    # Calculate betweenness centralities
    seed_centralities = sub_betweenness_centrality(H, seeds_left)
    print "seed centralities calculated"

    try:
        centralities =  pickle.load(open("stored_vals/centralities_no_rescale.p", "rb") )
    except:
        centralities = sub_betweenness_centrality(H, node_list)
        if not os.path.exists('stored_vals'):
            os.makedirs('stored_vals')
        pickle.dump(centralities, open("stored_vals/centralities_no_rescale.p", "wb"))
    print "all centralities loaded"

    # Comp centralities are merely all_centralities - seed_centralities
    comp_centralities = dict()
    for node in H:
        comp_centralities[node] = centralities[node] - seed_centralities[node]
    print "comp centralities calculated"

    Sj_file = open('results/bc_Sj_%s.txt' %(ADVERSE_EVENT), 'w')

    Sj_dict = dict()
    neighborhood = []
    for node in H:
        if centralities[node] != 0:
            denominator = centralities[node]/( (len(node_list)-1) * (len(node_list)-2) )
            
            if node in seeds_left:
                node_in_seeds = 1
                node_in_comp = 0
            elif node in complement:
                node_in_seeds = 0
                node_in_comp = 1
            
            inSeed = seed_centralities[node]/( (len(seeds_left)-node_in_seeds) * (len(node_list)-2) )
            
            compVal = comp_centralities[node]/( (len(complement)-node_in_comp) * (len(node_list)-2) )
            
            Sj = (inSeed - compVal) / denominator
            
        else:
            inSeed = 0
            compVal = 0
            denominator = 0
            Sj = 0
            
        Sj_dict[node]= Sj
        Sj_file.write("%s\t%s\n" % (str(Sj),node))
        
        if Sj > 0:
            if not node in neighborhood:
                neighborhood.append(node)

    print "Neighborhood size:",len(neighborhood)

    Sj_file.close()

    return Sj_dict