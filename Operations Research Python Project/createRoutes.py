# ENGSCI263: Group 19 Project - Warehouse vehicle routing project
# Project.py

# PURPOSE:
# This file contains the functions necassary for creating the vehicle routes

# import modules
import numpy as np
import pandas as pd
import itertools 
from pulp import *

# read in data
warehouseLocations = pd.read_csv('WarehouseLocations.csv')
warehouseDemand = pd.read_csv('demandestimate1.csv')
warehouseDurations = pd.read_csv('WarehouseDurations.csv', index_col = 0)

def create_routes(timeOfWeek, northClosed = False, simulation = False):
    ''' Creates set of feasible routes for a given time of the week. 
    
    Parameters:
    -------
    timeOfWeek : string
        Which time of the week to create routes for.
    northClosed : boolean
        Set whether the Northern distribution centre is closed or not.
    simulation : boolean
        Set whether to run the function in simulation mode

    Returns:
    --------
    feasibleRoutes : Pandas DataFrame
        Collection of feasible routes. 

    Notes:
    ------
    - 'Sunday' is not a valid parameter and 'Weekday' can be used to find routes assuming 
    each weekday has the same average demand value. 
    - This function saves the routes to a text file.
    '''

    # partition stores according to warehousePartition.png
    locationsNorth = []
    locationsSouth = []
    for i in range(len(warehouseLocations)):
        if warehouseLocations.iloc[i]['Long'] < warehouseLocations.iloc[i]['Lat']*1.1109890109890244 + 215.76870087912135:
            locationsNorth.append(warehouseLocations.iloc[i]['Store'])
        else:
            locationsSouth.append(warehouseLocations.iloc[i]['Store'])
    locationsNorth.remove('Distribution North')
    locationsSouth.remove('Distribution South')

    # remove Noel Leeming stores if day of the week is Saturday
    if timeOfWeek == 'Saturday':
        locationsNorth = [x for x in locationsNorth if x.startswith('The Warehouse')]
        locationsSouth = [x for x in locationsSouth if x.startswith('The Warehouse')]

    # determine all feasible routes for each partition
    feasibleRoutes = []
    totalDuration = []
    for i in range(2,7):
        # weekday routes cannot have more than 3 stores due to demand values
        if timeOfWeek != 'Saturday' and i > 4:
            break
        # find all possible combinations of stores that are size i
        all_combinations = [itertools.combinations(locationsNorth, i),itertools.combinations(locationsSouth, i)]        
        for half_combinations in all_combinations:
            # remove combinations which are excessively long (more than 5 hours)
            valid_combinations = [x for x in half_combinations if sum(warehouseDurations.loc[x[i]][x[i+1]] for i in range(len(x)-1)) < 18000]
            for combination in valid_combinations:
                # determine if combination does not exceed truck capacity
                pallets = sum(warehouseDemand.loc[warehouseDemand['Name'].isin([c for c in combination]), timeOfWeek])
                if pallets <= 20:
                    # if valid, find all possible permutations of the combination
                    all_permutations = itertools.permutations(combination)
                    duration = []
                    routes = []
                    for permutation in all_permutations:
                        permutation = list(permutation)
                        # add the relevant distribution centre start/ends to the permutation
                        if permutation[0] in locationsNorth and not northClosed:
                            permutation.insert(0,'Distribution North')
                            permutation.extend(['Distribution North'])
                        else:
                            permutation.insert(0,'Distribution South')
                            permutation.extend(['Distribution South'])
                        routes.append(permutation)
                        # determine the total duration of the permutation
                        duration.append(sum(warehouseDurations.loc[permutation[i]][permutation[i+1]] for i in range(len(permutation)-1)))
                    # add the most optimal permutation to the list of feasible routes in terms of duration
                    feasibleRoutes.append(routes[duration.index(min(duration))])
                    # save the duration for future reference (plus unloading times)
                    totalDuration.append(min(duration) + 600 * pallets)
                else:
                    continue

    # convert feasible routes list to pandas dataframe including a Route number and Cost column
    feasibleRoutes = pd.DataFrame(pd.Series(feasibleRoutes), columns = ['Stores'])
    feasibleRoutes.insert(0,'Route',range(1, 1 + len(feasibleRoutes)))
    # determine total cost of every route
    totalCosts = []
    for i in range(len(totalDuration)):
        if totalDuration[i] > 4*3600:
            totalCosts.append(4*175 + 250/3600 * (totalDuration[i] - 4*3600))
        else:
            totalCosts.append(totalDuration[i] * 175 / 3600)
    feasibleRoutes.insert(2,'Cost',totalCosts)
    
    # save feasible routes to text file
    if northClosed and not simulation:
        filename = 'northClosedFeasibleRoutes' + timeOfWeek + '.txt'
    else:
        filename = 'feasibleRoutes' + timeOfWeek + '.txt'
    with open(filename,'w') as outfile:
        feasibleRoutes.to_csv(outfile, index = False, sep = ' ')

    # print(feasibleRoutes)
    return feasibleRoutes
    
if __name__ == "__main__":
    # generate feasible routes (this takes approx 10 mins)
    # feasibleRoutesWeekday = create_routes('Weekday')
    # feasibleRoutesSaturday = create_routes('Saturday')
    feasibleRoutesWeekdayNorthClosed = create_routes('Weekday',True)
    feasibleRoutesSaturdayNorthClosed = create_routes('Saturday',True)