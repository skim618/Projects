# ENGSCI263: Group 19 Project - Warehouse vehicle routing project
# LinearProgram.py

# PURPOSE:
# This file contains the functions necassary for the Linear Program

# Import modules
import pulp 
import numpy as np
import pandas as pd 
from pulp import *  

# Read in data
feasibleRoutesWeekday = pd.read_csv('feasibleRoutesWeekday.txt', delim_whitespace = True)
feasibleRoutesSaturday = pd.read_csv('feasibleRoutesSaturday.txt', delim_whitespace = True)   
northClosedFeasibleRoutesWeekday = pd.read_csv('northClosedFeasibleRoutesWeekday.txt', delim_whitespace = True)
northClosedFeasibleRoutesSaturday = pd.read_csv('northClosedFeasibleRoutesSaturday.txt', delim_whitespace = True)

def inRoute(store, route):
    '''Determines whether a certain store is in a certain route.'''
    if (route.find(store) != -1):
        return True
    return False

def LinearProgram (feasibleRoutes, NorthernClosed = False, Saturday = False):
    ''' Solves a mixed integer program to find the least-cost routing schedule given a set
    of feasible routes.
    
    Parameters:
    -------
    feasibleRoutes : Pandas DataFrame
        Colletion of feasible routes.
    NorthernClosed: boolean
        Changes formulation to accomodate for closing of the Northern distribution centre.
    Saturday : boolean
        Changes formulation to accomodate for closing of Noel Leeming stores on Saturdays.

    Returns:
    --------
    prob.objective : float
        Objective function value.
    optimalRoutes : list
        List of routes that solve the integer program
    '''

    # read in and extract relevant data
    RouteNo = feasibleRoutes['Route']
    Routes = feasibleRoutes['Stores']
    Costs = feasibleRoutes['Cost']
    warehouseNames = list(pd.read_csv('WarehouseLocations.csv', delimiter=',', usecols=[2], skiprows=[1,2])['Store'])
    # remove Noel Leeming stores from the formulation if it is Saturday
    if Saturday:
        warehouseNames = [x for x in warehouseNames if x.startswith('The Warehouse')]

    # intialise and create the store matrix array
    StoreMatrix = np.zeros((len(warehouseNames),len(Routes)))
    for i in range(len(Routes)):
        for j in range(len(warehouseNames)):
            if (inRoute(warehouseNames[j], Routes[i])):
                StoreMatrix[j][i] = 1
    print(StoreMatrix)

    # The problem variable of the routes taken
    vars = LpVariable.dicts("Route", RouteNo, 0, None, LpInteger)

    # Create the objective function
    prob = LpProblem("Vehicle Routing Problem", LpMinimize)

    # The objective function is entered 
    prob += lpSum([vars[i]*Costs[i-1] for i in RouteNo])

    # Constraints are entered (Add more if missed)
    for i in range(len(warehouseNames)):
        prob += lpSum(vars[j] * StoreMatrix[i][j-1] for j in RouteNo) == 1, "Route {} Constraint".format(i)
    
    if NorthernClosed:
        prob += lpSum(vars[i] * 1) <= 30, "Truck limit"
    else:
        prob += lpSum(vars[i] * 1) <= 25, "Truck limit"

    # Write the problem data and solve
    prob.writeLP("VehicleRoutingProblem.lp")
    prob.solve()

    # Print the status
    print("Status:", LpStatus[prob.status])
    print("Delivery Costs = ", value(prob.objective))

    # Print the results
    optimalRoutes = []
    for v in prob.variables():
        if(v.varValue != 0):
            print(v.name, "=", v.varValue)
            optimalRoutes.append(v.name)

    # Return only the route number rather than the full variable name
    optimalRoutes = pd.Series([int(x.split('_')[1]) for x in optimalRoutes], name = 'Route')
    
    return value(prob.objective), optimalRoutes

if __name__ == "__main__":
    # solve linear program for different scenarios
    weekdayObjective, weekdayRoutes = LinearProgram(feasibleRoutesWeekday)
    saturdayObjective, saturdayRoutes = LinearProgram(feasibleRoutesSaturday,Saturday=True)
    northClosedWeekdayObjective, northClosedWeekdayRoutes = LinearProgram(northClosedFeasibleRoutesWeekday,NorthernClosed=True)
    northClosedSaturdayObjective, northClosedSaturdayRoutes = LinearProgram(northClosedFeasibleRoutesSaturday,NorthernClosed=True,Saturday=True)
   
    # convert solution into readable format
    weekdayRoutes = pd.merge(feasibleRoutesWeekday[['Route','Stores']], weekdayRoutes, how = 'inner', on = 'Route')
    saturdayRoutes = pd.merge(feasibleRoutesSaturday[['Route','Stores']], saturdayRoutes, how = 'inner', on = 'Route')
    northClosedWeekdayRoutes = pd.merge(northClosedFeasibleRoutesWeekday[['Route','Stores']], northClosedWeekdayRoutes, how = 'inner', on = 'Route')
    northClosedSaturdayRoutes = pd.merge(northClosedFeasibleRoutesSaturday[['Route','Stores']], northClosedSaturdayRoutes, how = 'inner', on = 'Route')

    # save solution as text file
    with open('weekdayRoutes.txt','w') as outfile:
        weekdayRoutes['Stores'].to_csv(outfile, index = False)

    with open('saturdayRoutes.txt','w') as outfile:
        saturdayRoutes['Stores'].to_csv(outfile, index = False)

    with open('northClosedWeekdayRoutes.txt','w') as outfile:
        northClosedWeekdayRoutes['Stores'].to_csv(outfile, index = False)

    with open('northClosedSaturdayRoutes.txt','w') as outfile:
        northClosedSaturdayRoutes['Stores'].to_csv(outfile, index = False)   

    # # display results
    print('\nCosts per day with NDC open:\nWeekdays ->\t'+str(weekdayObjective)+'\nNumber of routes ->\t',len(weekdayRoutes))
    print('Saturdays ->\t'+str(saturdayObjective)+'\nNumber of routes ->\t',len(saturdayRoutes))
    print('Total cost per month ->\t'+ str(4*(5*weekdayObjective+saturdayObjective)))
    print('\nCosts per day with NDC closed:\nWeekdays ->\t'+str(northClosedWeekdayObjective)+'\nNumber of routes ->\t',len(northClosedWeekdayRoutes))
    print('Saturdays ->\t'+str(northClosedSaturdayObjective)+'\nNumber of routes: \t', len(northClosedSaturdayRoutes))
    print('Total cost per month ->\t'+ str(4*(5*northClosedWeekdayObjective+northClosedSaturdayObjective)))
    print('\nOptimal Routes with NDC Open \n------------------------------')
    print('\nWeekdays: ')
    [print(str(i)+'\t' + weekdayRoutes['Stores'][i][1:(len(weekdayRoutes['Stores'][i])-23)]) for i in range(len(weekdayRoutes))]
    print('\nSaturdays:')
    [print(str(i)+'\t' + saturdayRoutes['Stores'][i][1:(len(saturdayRoutes['Stores'][i])-23)]) for i in range(len(saturdayRoutes))]
    print('\nOptimal Routes with NDC Closed \n------------------------------')
    print('\nWeekdays:')
    [print(str(i)+'\t' + northClosedWeekdayRoutes['Stores'][i][23:(len(northClosedWeekdayRoutes['Stores'][i])-23)]) for i in range(len(northClosedWeekdayRoutes))]
    print('\nSaturdays:')
    [print(str(i)+'\t' + northClosedSaturdayRoutes['Stores'][i][23:(len(northClosedSaturdayRoutes['Stores'][i])-23)]) for i in range(len(northClosedSaturdayRoutes))]
    print('\n')