#Script by Roberto Corso, 2020
#Made as part of the internship at Università degli studi di Palermo

import igraph
import pandas as pd
import numpy as np
import sys
import os
from scipy import stats

### PARAMETERS ###

if __name__ == '__main__':
    path = sys.argv[1] #path to edgelist of bipartite network
    tail = sys.argv[2] #which tail of the hypergeometric distribution should be used to validate the network
    if tail != 'over' and tail != 'under' and tail != 'both': #raises an error if an invalid argument is passed
        sys.exit('Error: invalid argument "tail"') 
    method = sys.argv[3] #which method should be used to validate the network, 'B' for Bonferroni or 'FDR' for False Discovery Rate
    if method != 'B' and method != 'FDR': #raises an error if an invalid argument is passed
        sys.exit('Error: invalid argument "method"')
    stat_threshold = float(sys.argv[4]) #threshold of the statistical significance chosen by the user
    if not 0 <= stat_threshold <= 1:
        sys.exit('Error: stat_threshold must be between 0 and 1') #raises an error if an invalid value is passed
    if len(sys.argv) == 5: #if no addition to the output name has been specified
        name_extension = '_validated' #default value of the argument
    else:
        name_extension = sys.argv[5] #extension of the default name of the output

    ### BIPARTITE NETWORK ###

    df = pd.read_csv(path, sep = ' ', header = None).dropna(axis = 'columns') #copies the edgelist into a dataframe and removes empty columns
    vs_shift = min(df[0]) 
    df[0] -= vs_shift #indices of the vertices of the first set must start from 0
    df[1] += -min(df[1]) + max(df[0]) +1 #indices of the vertices of the second set must start from the position after the last index of the first set
    n_vertices = len(set(df[0])) #number of vertices in the first set
    g = igraph.Graph(df.values.max()+1) #creates a graph object without edges
    g.add_edges([tuple(x) for x in df.values]) #adds the edges from the dataframe
    if sorted(set(df[0])) != list(range(min(df[0]), max(df[0])+1)): #checks if there are holes in the sequence of vertices of the first set
        g.delete_vertices(g.vs.select(_degree = 0)) #fixes potential holes in the vertex sequence by removing vertices added by Graph function
    g.vs[:n_vertices]['original_id'] = sorted(set(df[0]+vs_shift)) #stores the original id of each vertex of the first set

    ### PROJECTED NETWORK ###

    g.vs[0:n_vertices]['type'] = 0 #specifies the set of each vertex
    g.vs[n_vertices:]['type'] = 1
    g_proj = g.bipartite_projection(types = 'type', multiplicity = True, which = 0) #creates the projected network, 'multiplicity' stores edge weights, 'which' is the main set
    g_proj.vs['set_2_neighbors'] = g.degree(range(n_vertices)) #stores the degree of the verticces of the first set of the bipartite network
    g_proj.vs['original_id'] = g.vs[:n_vertices]['original_id'] #stores the original id of each vertex of the first set
    g_proj.delete_vertices(g_proj.vs.select(_degree = 0)) #removes all vertices with degree zero from the projected network

    ### STATISTIC VALIDATION ###

    n_vertices_set_2 = len(set(df[1])) #number of vertices in the second set
    n_edges = len(g_proj.es) #number of edges of the preojected network
    Nt = n_vertices*(n_vertices-1)/2 #number of hypotheses
    for i, e in enumerate(g_proj.es):
        g_proj.es[i]['hypergeom_param'] = tuple([g_proj.es(i)['weight'][0], g_proj.vs(e.source)['set_2_neighbors'][0], g_proj.vs(e.target)['set_2_neighbors'][0]]) #for each edge stores the edge weight and the degree of its vertices in the original bipartite network
    if tail == 'over' or tail == 'both': #survival function is used to evaluate over-expressions
        hypergeom_dict_over = {hypergeom_tuple:stats.hypergeom.sf(hypergeom_tuple[0]-1, n_vertices_set_2, hypergeom_tuple[1], hypergeom_tuple[2]) for hypergeom_tuple in set(g_proj.es['hypergeom_param'])} #creates a dictionary whose keys are the tuples of parameters of the hypergeometric function and whose values are the values of the survival function used to evaluate over-expressions
    if tail == 'under' or tail == 'both': #cumulative function is used to evaluate under-expressions
        hypergeom_dict_under = {hypergeom_tuple:stats.hypergeom.cdf(hypergeom_tuple[0], n_vertices_set_2, hypergeom_tuple[1], hypergeom_tuple[2]) for hypergeom_tuple in set(g_proj.es['hypergeom_param'])} #creates a dictionary whose keys are the tuples of parameters of the hypergeometric function and whose values are the values of the cumulative function used to evaluate under-expressions

    output = pd.DataFrame([g_proj.vs[e.source]['original_id'] for e in g_proj.es], columns = ['source']) #adds the column of source vertices to the edgelist
    output = pd.concat([output, pd.DataFrame([g_proj.vs[e.target]['original_id'] for e in g_proj.es], columns = ['target'])], axis = 1) #adds the column of target vertices to the edgelist
    output = pd.concat([output, pd.DataFrame(g_proj.es['weight'], columns = ['weight'])], axis = 1) #adds the column of edge weights to the edgelist
    if 'hypergeom_dict_over' in locals():
        output = pd.concat([output, pd.DataFrame([hypergeom_dict_over[g_proj.es[i]['hypergeom_param']] for i in range(n_edges)], columns = ['p-value_over'])], axis = 1) #adds the over-expression p-value of each edge to the edgelist
    if 'hypergeom_dict_under' in locals():
        output = pd.concat([output, pd.DataFrame([hypergeom_dict_under[g_proj.es[i]['hypergeom_param']] for i in range(n_edges)], columns = ['p-value_under'])], axis = 1) #adds the under-expression p-value of each edge to the edgelist
    
    if method == 'B': #Bonferroni
        B_value = stat_threshold/Nt #the corrected thershold to test edges
        if 'hypergeom_dict_over' in locals():
            output['test_over'] = 'fail' #adds the column to the output file with a default value
            output.loc[output['p-value_over'] < B_value, 'test_over'] = 'success' #updates the value for all the edges that passed the test
            if tail == 'both':
                output = output[['source', 'target', 'weight', 'p-value_over', 'test_over', 'p-value_under']] #reorders the columns for easier reading of the output
        if 'hypergeom_dict_under' in locals():
            output['test_under'] = 'fail' #adds the column to the output file with a default value
            output.loc[output['p-value_under'] < B_value, 'test_under'] = 'success' #updates the value for all the edges that passed the test
    
    elif method == 'FDR': #False Discovery Rate
        array_check = np.arange(1, n_edges+1)*stat_threshold/Nt #the threshold increases linearly with the number of edges
        if 'hypergeom_dict_over' in locals():
            p_values_sorted = np.sort(output['p-value_over']) #creates a sorted list of the over-expression p-values
            fdr_value = max(p_values_sorted[p_values_sorted <= array_check], default = 0) #finds the highest p-value that satisfies the condition
            output['test_over'] = 'fail' #adds the column to the output file with a default value
            output.loc[output['p-value_over'] <= fdr_value, 'test_over'] = 'success' #updates the value for all the edges that passed the test
            if tail == 'both':
                output = output[['source', 'target', 'weight', 'p-value_over', 'test_over', 'p-value_under']] #reorders the columns for easier reading of the output
        if 'hypergeom_dict_under' in locals():
            p_values_sorted = np.sort(output['p-value_under']) #creates a sorted list of the over-expression p-values
            fdr_value = max(p_values_sorted[p_values_sorted <= array_check], default = 0) #finds the highest p-value that satisfies the condition
            output['test_under'] = 'fail' #adds the column to the output file with a default value
            output.loc[output['p-value_under'] <= fdr_value, 'test_under'] = 'success' #updates the value for all the edges that passed the test
            
    output.to_csv(os.path.splitext(path)[0]+'_'name_extension+'.txt', sep = ' ', index = False) #writes the output edgelist into a file, info on splitext https://stackoverflow.com/questions/678236/how-to-get-the-filename-without-the-extension-from-a-path-in-python