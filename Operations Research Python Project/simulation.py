# ENGSCI263: Group 19 Project - Warehouse vehicle routing project
# simulation.py

# PURPOSE:
# Evaluate the quality of your schedule by creating a simulation 
# to estimate the actual cost of satisfying all pallet demand at 
# every store. Your simulation should take into variations in demand 
# and sensibly approximate the effect of traffic. Hence, give an 
# estimate of the cost of operating your proposed routing schedule, 
# with and without the Northern distribution centre.

# import modules
import matplotlib.pyplot as plt
import statistics as st
import statsmodels.stats.weightstats as sms
import numpy as np
import pandas as pd
from scipy import stats
from createRoutes import *
import random
import seaborn as sns
import warnings
warnings.filterwarnings("ignore", category = FutureWarning)

# these are the objective function values from solving the LP
weekdayObjective = 10672.590470525001
saturdayObjective = 3203.883531159722
northClosedWeekdayObjective = 11767.394668587502
northClosedSaturdayObjective = 3524.309717427778
# read in data
weekdayRoutes = pd.read_csv('weekdayRoutes.txt', converters = {'Stores' : eval})
saturdayRoutes = pd.read_csv('saturdayRoutes.txt', converters = {'Stores' : eval})
northClosedWeekdayRoutes = pd.read_csv('northClosedWeekdayRoutes.txt', converters = {'Stores' : eval})
northClosedSaturdayRoutes = pd.read_csv('northClosedSaturdayRoutes.txt', converters = {'Stores' : eval})

# create distribution of demands
weekdayDemands = np.zeros((40,100))
for i in range(100):
    for j in range(len(warehouseDemand['Name'])):
        store = warehouseDemand['Name'].iloc[j]
        selectedWeekdayDemands = warehouseDemand.loc[warehouseDemand['Name']==store,'Monday':'Friday']
        sampled = selectedWeekdayDemands.iloc[0,:].sample(n = 5, replace = True)
        sampled = list(sampled)
        weekdayDemands[j,i] = np.mean(sampled)

# calculate distribution of costs/month with northern distribution open
northOpenObjective = [4*(5*weekdayObjective+saturdayObjective)] * 100
for i in range(100):
    warehouseDemand['Weekday'] = weekdayDemands[:,i]
    # loop through weekday routes
    for route in weekdayRoutes['Stores']:
        routeFiltered = [x for x in list(route) if not x.startswith('Distribution')]
        pallets = sum(warehouseDemand.loc[warehouseDemand['Name'].isin([x for x in routeFiltered]), 'Weekday'])
        # take into account traffic, assume a delay anywhere between 0 and 30 minutes
        duration = sum(warehouseDurations.loc[route[i]][route[i+1]] for i in range(len(route)-1)) + random.randint(0,1800) + 600*pallets
        # in the case that there are more than 20 pallets we send another truck
        if pallets > 20:
            warehouseDemand = pd.read_csv('demandestimate1.csv')
            original = sum(warehouseDemand.loc[warehouseDemand['Name'].isin([x for x in routeFiltered]), 'Weekday'])
            # the truck should only carry the difference between the new and the original number of pallets required
            duration -= 600*pallets
            duration += 600*(pallets - original)
            if duration > 4*3600:
                cost = 4*175 + 250/3600 * (duration - 4*3600)
                northOpenObjective[i] += cost
            else:
                cost = duration * 175 / 3600
                northOpenObjective[i] += cost
        # if there are less than 20 pallets we don't send another truck but we alter the cost of the original route
        else:
            # add the duration of the new number of pallets
            if duration > 4*3600:
                cost = 4*175 + 250/3600 * (duration - 4*3600)
                northOpenObjective[i] += cost
            else:
                cost = duration * 175 / 3600
                northOpenObjective[i] += cost
            # subtract the duration of the original number of pallets
            warehouseDemand = pd.read_csv('demandestimate1.csv')
            pallets = sum(warehouseDemand.loc[warehouseDemand['Name'].isin([x for x in routeFiltered]), 'Weekday'])
            duration = sum(warehouseDurations.loc[route[i]][route[i+1]] for i in range(len(route)-1)) + 600*pallets
            if duration > 4*3600:
                    cost = 4*175 + 250/3600 * (duration - 4*3600)
                    northOpenObjective[i] -= cost
            else:
                cost = duration * 175 / 3600
                northOpenObjective[i] -= cost

    # loop through saturday routes
    for route in saturdayRoutes['Stores']:
        routeFiltered = [x for x in list(route) if not x.startswith('Distribution')]
        pallets = sum(warehouseDemand.loc[warehouseDemand['Name'].isin([x for x in routeFiltered]), 'Saturday'])
        duration = sum(warehouseDurations.loc[route[i]][route[i+1]] for i in range(len(route)-1)) + 600*pallets
        # subtract the original duration without traffic
        if duration > 4*3600:
            cost = 4*175 + 250/3600 * (duration - 4*3600)
            northOpenObjective[i] -= cost
        else:
            cost = duration * 175 / 3600
            northOpenObjective[i] -= cost
        # add the new duration with traffic
        duration += random.randint(0,1800)
        if duration > 4*3600:
            cost = 4*175 + 250/3600 * (duration - 4*3600)
            northOpenObjective[i] += cost
        else:
            cost = duration * 175 / 3600
            northOpenObjective[i] += cost


# plot distribution
sns.distplot(northOpenObjective)
plt.xlabel('Cost of Operation / Month ($)')
plt.show()
print('\nNorthern Centre Open \n---------------------')
print("Mean cost with northern open: {}".format(st.mean(northOpenObjective)))
print("Northern open 95% Percentile: {}".format(sms.DescrStatsW(northOpenObjective).tconfint_mean(alpha = 0.05)))

# calculate distribution of costs/month with northern distribution closed
northClosedObjective = [4*(5*northClosedWeekdayObjective+northClosedSaturdayObjective)] * 100
for i in range(100):
    warehouseDemand['Weekday'] = weekdayDemands[:,i]
    # loop through weekday routes
    for route in northClosedWeekdayRoutes['Stores']:
        routeFiltered = [x for x in list(route) if not x.startswith('Distribution')]
        pallets = sum(warehouseDemand.loc[warehouseDemand['Name'].isin([x for x in routeFiltered]), 'Weekday'])
        # take into account traffic, assume a delay anywhere between 0 and 30 minutes
        duration = sum(warehouseDurations.loc[route[i]][route[i+1]] for i in range(len(route)-1)) + random.randint(0,1800) + 600*pallets
        # in the case that there are more than 20 pallets we send another truck
        if pallets > 20:
            warehouseDemand = pd.read_csv('demandestimate1.csv')
            original = sum(warehouseDemand.loc[warehouseDemand['Name'].isin([x for x in routeFiltered]), 'Weekday'])
            # the truck should only carry the difference between the new and the original number of pallets required
            duration -= 600*pallets
            duration += 600*(pallets - original)
            if duration > 4*3600:
                cost = 4*175 + 250/3600 * (duration - 4*3600)
                northClosedObjective[i] += cost
            else:
                cost = duration * 175 / 3600
                northClosedObjective[i] += cost
        # if there are less than 20 pallets we don't send another truck but we alter the cost of the original route
        else:
            # add the duration of the new number of pallets
            if duration > 4*3600:
                cost = 4*175 + 250/3600 * (duration - 4*3600)
                northClosedObjective[i] += cost
            else:
                cost = duration * 175 / 3600
                northClosedObjective[i] += cost
            # subtract the duration of the original number of pallets
            warehouseDemand = pd.read_csv('demandestimate1.csv')
            pallets = sum(warehouseDemand.loc[warehouseDemand['Name'].isin([x for x in routeFiltered]), 'Weekday'])
            duration = sum(warehouseDurations.loc[route[i]][route[i+1]] for i in range(len(route)-1)) + 600*pallets
            if duration > 4*3600:
                    cost = 4*175 + 250/3600 * (duration - 4*3600)
                    northClosedObjective[i] -= cost
            else:
                cost = duration * 175 / 3600
                northClosedObjective[i] -= cost

    # loop through saturday routes
    for route in northClosedSaturdayRoutes['Stores']:
        routeFiltered = [x for x in list(route) if not x.startswith('Distribution')]
        pallets = sum(warehouseDemand.loc[warehouseDemand['Name'].isin([x for x in routeFiltered]), 'Saturday'])
        duration = sum(warehouseDurations.loc[route[i]][route[i+1]] for i in range(len(route)-1)) + 600*pallets
        # subtract the original duration without traffic
        if duration > 4*3600:
            cost = 4*175 + 250/3600 * (duration - 4*3600)
            northClosedObjective[i] -= cost
        else:
            cost = duration * 175 / 3600
            northClosedObjective[i] -= cost
         # add the new duration with traffic
        duration += random.randint(0,1800)
        if duration > 4*3600:
            cost = 4*175 + 250/3600 * (duration - 4*3600)
            northClosedObjective[i] += cost
        else:
            cost = duration * 175 / 3600
            northClosedObjective[i] += cost

# plot distribution
sns.distplot(northClosedObjective)
plt.xlabel('Cost of Operation / Month ($)')
plt.show()
print('\nNorthern Centre Closed \n---------------------')
print("Mean cost with northern closed: {}".format(st.mean(northClosedObjective)))
print("Northern closed 95% Percentile: {}".format(sms.DescrStatsW(northClosedObjective).tconfint_mean(alpha = 0.05)))
print('\n')
