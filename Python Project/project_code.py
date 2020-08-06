# ENGSCI233: Project - Heuristic algorithm development
# Sooyong Kim - 573360750 - skim618
import networkx as nx
import numpy as np
from project_utils import *
from time import time

# ** below are implemented functions ** 
# -----------------------------------------------------------------------------------------
def shortest_dest(network, source, dest_path):
    """ Finds the node from the list that is closest to the starting node in the network

    Parameters
    ----------
    network : networkx.Graph
        The graph that contains the node and edge information
    source : string
        the name of starting node
    dest_path: list
        list of names of possible destination nodes

    Returns
    ----------
    shortest_node : string
        the name of node from dest_path that is closest to the source
    temp_length : float
        the shortest time (in hours) it takes to travel along the source to the shortest_node

    Notes
    ----------
    time between each pair is obtained through the shortest_path_length function from the networkx module

    """
    
    # initialise starting length and destination for comparison
    temp_length = nx.shortest_path_length(network, source, dest_path[0], weight='weight')
    shortest_node = dest_path[0]

    # iterate through each of the destination nodes
    for dest in dest_path:
        # calculate length of current node
        length = nx.shortest_path_length(network, source, dest, weight='weight')

        # if current length is lower than the inital, replace the inital with current and save the node
        if length < temp_length:        
            temp_length = length
            shortest_node = dest

    return shortest_node, temp_length

def nearest_neighbour(network, dest_path):
    """ Creates the shortest path from the source node to all the destination nodes

    Parameters
    ----------
    network : networkx.Graph
        The graph that contains the node and edge information
    dest_path: list
        list of names of possible destination nodes

    Returns
    ----------
    actual_path : list
        the actual path taken including the bus route from source to the destination nodes
    path : list
        the path taken node to node from the source (no bus route)
    totDist : float
        the total time taken (hours) for the trip

    Notes
    ----------
    The bus path between each pair is obtained through the shortest_path function from the networkx module.
    This function calls shortst_dest function to obtain the nearest node
    This function will assume the source and final destination is Auckland Airport
    """

    # initialise the path taken list, starting from the source
    source = 'Auckland Airport'
    path = [source]
    totDist = 0                                 # initialise the total time taken

    # keep iterating through until all of the destinations have been reached
    while len(dest_path) != 0:
        closest_node, distance = shortest_dest(network, source, dest_path)
        for i in range (len(dest_path)):
            if closest_node == dest_path[i]:
                pop = dest_path.pop(i)          # remove the closest node from the destination list
                break
        path.append(pop)                        # append the popped node to the path taken list
        totDist = totDist + distance            # cumulatively sum the distances
        source = pop                            # set the new source as the popped node

    # add the time taken to return to Auckland Airport
    totDist = totDist + nx.shortest_path_length(network, source, 'Auckland Airport', weight='weight')

    # append Auckland Airport at the end of the path taken list
    path.append('Auckland Airport')
    
    # obtain the bus route between nodes 
    actual_path = []
    # iterating through pairs of rest homes in the path 
    for rest1, rest2 in zip(path[:-1], path[1:]):
        bus_route = nx.shortest_path(network, rest1, rest2, weight='weight')
        actual_path.extend(bus_route[:-1]) 
    actual_path.append(path[-1])

    return actual_path, path, totDist

def partition_homes(network, rest_homes):
    """ Partitions the list of rest_homes into west and east categories

    Parameters
    ----------
    network : networkx.Graph
        The graph that contains the node and edge information
    rest_homes: list
        list of names of all possible destination nodes

    Returns
    ----------
    westLng : list
        list of nodes that are categoried in the west side
    eastLng : list
        list of nodes that are categoried in the east side

    Notes
    ----------
    The sum of the lengths of the two outputs should equal the length of the total list.
    The lists are categorised by a longitude limit
    """

    # initialise an empty array to store all the longitudes
    lng = []

    # iterate through all rest homes and store their longitudes
    for rest in rest_homes:
        lng.append(network.nodes[rest]['lng'])

    lng = np.sort(lng)                                  # sort the longitudes in ascending order
    midlng = lng[int(len(lng)/2)]                       # assign a limit 

    # initialise west and east lists
    westLng = []
    eastLng = []

    # if the longitude of the rest home is lower than that of the limit, append it to the west list.
    # else, append it to the east.
    for rest in rest_homes:
        if network.nodes[rest]['lng'] < midlng:
            westLng.append(rest)
        else:
            eastLng.append(rest)

    return westLng, eastLng

def partition_homes2(network, rest_homes):
    """ Partitions the list of rest_homes into 4 lists

    Parameters
    ----------
    network : networkx.Graph
        The graph that contains the node and edge information
    rest_homes: list
        list of names of all possible destination nodes

    Returns
    ----------
    destination_lists : list
        list of lists of homes that are categorised by their longitude and latitude
        [West(south), West(north), East(south), East(north)]

    Notes
    ----------
    The sum of the lengths of the four outputs should equal the length of the total list.
    The lists are categorised by calling the function partition_homes(), then splitting it according to their latitude
    """

    # obtain west and east homes by calling the funciton below
    west_homes, east_homes = partition_homes(network, rest_homes)

    # initialise the lists to store the latitudes
    latWest = []
    latEast = []

    # iterate through both list of homes, storing their latitudes
    for west in west_homes:
        latWest.append(network.nodes[west]['lat'])
    for east in east_homes:
        latEast.append(network.nodes[east]['lat'])

    # sort the latitudes in ascending order
    latWest = np.sort(latWest)
    latEast = np.sort(latEast)

    # assign a limit 
    midWest = latWest[int(2*len(latWest)/3)]
    midEast = latEast[int(len(latEast)/2)]

    # initialise the four lists based on latitudes
    West1 = []
    West2 = []
    East1 = []
    East2 = []

    # if the latitude of the rest home is lower than that of the limit, append it to the south list.
    # else, append it to the north.
    for west in west_homes:
        if network.nodes[west]['lat'] < midWest:
            West1.append(west)
        else:
            West2.append(west)
    for east in east_homes:
        if network.nodes[east]['lat'] < midEast:
            East1.append(east)
        else:
            East2.append(east)
    
    # store the lists in a final list
    destination_lists = [West1, West2, East1, East2]

    return destination_lists

# -----------------------------------------------------------------------------------------
# ** BELOW IS THE MAIN FUNCTION ** 

def main():
    """Reads in a network, and will return four ordered lists of network node names, 
        representing the paths taken by four couriers"""

    # read in the network, and the rest home data
    auckland = read_network('network.graphml')
    rest_homes = get_rest_homes('rest_homes.txt')

    # partition the rest_homes into four for the four couriers
    destinations = partition_homes2(auckland, rest_homes)

    # iterate through each of the lists and obtain their shortest path
    for i in range (len(destinations)):
        actualPath, justPath, distance = nearest_neighbour(auckland, destinations[i])
        # print('The time it takes for path_{:d} is {:.1f} hours'.format(i+1, distance))

        # save the path data to a text file
        fp = open('path_{:d}.txt'.format(i+1),'w')
        for node in justPath:
            fp.write(node + '\n')
        fp.close()

        # save the plot 
        plot_path(auckland, actualPath, save='path_{:d}.png'.format(i+1))

if __name__ == "__main__":
    t0 = time()  # Start time
    main()
    t1 = time() # End time
    # print('Time taken for this code to compute solutions is: {:.1f} seconds'.format(t1-t0))



