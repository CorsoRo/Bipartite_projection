Script by Roberto Corso, 2020
Made as part of the internship at Università degli studi di Palermo

This is a script for the statistical validation of undirected unweighted bipartite networks.
A bipartite network is a network whose vertices can be categorized in 2 sets and in which each edge connects vertices of different sets.
The script creates a projection of the bipartite network by selecting all vertices of the first set and linking each pair if the two vertices have common neighbors in the second set.
The script then validates the new network against a null hypothesis by associating to each edge a p-value given by the probability of co-occurrences in the sets of neighbors of the two vertices obtained with the hypergeometric distribution.
Each p-value is then tested against a statistical threshold chosen by the user. The statistical threshold is corrected to account the multiple hypotheses testing. Two corrections are implemented in this script: the Bonferroni correction and the False Discovery Rate correction.
This analysis has been showcased in this article: https://arxiv.org/pdf/1008.1414.pdf.
It is also possible to evaluate both over-expressed and under-expressed co-occurencies. This method has been used in http://dx.doi.org/10.1080/14697688.2014.969889.

This script requires the following arguments:
path: the path to the edgelist of the bipartite network;
tail: can be "over" to evaluate over-expressed co-occurences, "under" to evaluate under-expressed co-occurencies, or "both";
method: the correction to account for multiple hypotheses testing, can be either "B" for Bonferroni correction or "FDR" for False Discovery Rate correction;
stat_threshold: the desired statistical significance for the validation, must be either int or float;
name_extension (optional): any string to add to the output name, dafault value is "_validated".

This script uses the libraries igraph, pandas, numpy, sys, os and scipy.

The input edgelist must contain the IDs of the two vertices of each link on columns separated by a single space without header; edge weights are not accounted for. The bipartite network is considered to be undirected
The output reports the edgelist of the projected network including edge weights, over-expression or under-expression p-values for each edge (or both) and the result of the test, either "fail" or "success". The IDs of the output edgelist are the same as those of the input edgelist.