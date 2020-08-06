# ENGSCI233: Project - Heuristic algorithm development
# Sooyong Kim - 573360750 - skim618

# PREPARATION:
# Completing the GETTING_STARTED exercies

from project_utils import *
auckland = read_network('network.graphml')
rest_homes = get_rest_homes('rest_homes.txt')

# We can use in-built networkx algorithms to compute shortest paths. For example,
# finding a shortest path between Everil Orr and Kumeu Village can be done as follows:

import networkx as nx
path = nx.shortest_path(auckland, 'Everil Orr', 'Kumeu Village', weight='weight')
distance = nx.shortest_path_length(auckland, 'Everil Orr', 'Kumeu Village', weight='weight')

# Here, we specify that we use the weight attribute to weight the edges for the shortest path algorithm.
# We can plot the path using the plot_path function from project_utils:

plot_path(auckland, path, save=False)

########################################

auckland['Auckland Airport']



