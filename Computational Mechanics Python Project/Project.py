# ENGSCI263: Group 15 Project - Onehunga Aquifer
# Project.py

# PURPOSE:
# This file contains the functions necassary for creating the computer model evaluating
# the pressure and concentration of the Ongehunga aquifer

# PREPARATION:
# Review the CALIBRATION notes, Wairakei model notebook and functions discussed in the Google doc.
# Your group should be (1) deriving the equations, (2) implementing the equations as a model, (3) calibrating the model.

# SUBMISSION:
# Aim to finish functions below by 29th of August

# imports
import numpy as np
import random
from matplotlib import pyplot as plt
from matplotlib import cm
import time
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm as colmap
from copy import copy
from scipy.stats import multivariate_normal

def ODE_model(t, p, c, q, a, b, d, Mo, p1):
    ''' 
    Return the derivative dp/dt and dc/dt at time, t, for given parameters.

    Parameters:
    ----------
    t : float
        Independent variable.
    p : float
        Dependent variable (Pressure).
    c : float
        Dependent variable (Concentration).
    q : float
        Extraction rate.
    a : float
        Source/sink strength parameter.
    b : float
        Recharge strength parameter.
    d : float
        Surface leeching parameter.
    Mo : float
        Total system mass.
    p0 : float
        Pressure at low pressure boundary.
    p1 : float
        Pressure at how pressure boundary.
    Csrc : float
        Concentration at the surface?

    Returns:
    ----------
    dpdt : float
        Derivative of pressure with respect to independent variable.
    dcdt : float
        Derivative of concentration with respect to independent variable.

    Notes:
    ----------
    Concentration depends on direction flow
    '''

    # Concentration depends on direction flow
    if (p > 0.0):
        c_dash = c
    else:
        c_dash = 0
    
    # derivative equations for pressure and concentration respectively
    dpdt = -a*q - 2*b*p + b*p1
    dcdt = (((-b/a)*(p)*(c_dash-c)) + ((b/a)*(p-p1)*c) - d*(p-((p1)/2.))) / Mo

    return dpdt, dcdt

def improved_euler_step(f, tk, Pk, Ck, h, q, *pars):
    '''
    Compute a single Euler step.
	
    Parameters:
	----------
	f : callable
		Derivative function.
	tk : float
		Independent variable at beginning of step.
	Pk : float
		Solution for pressure at beginning of step.
   	Ck : float
		Solution for concentration at beginning of step.
	h : float
		Step size.
    q : float 
        Extraction rate
    *pars: array-like
        An array of parameters used in the ODE model

	Returns:
	----------
	xk1 : float
		Solution at end of the improved Euler step.
	'''	
    #Euler step
    dPdt,dCdt = f(tk,Pk,Ck,q,*pars)
    P1Euler=  Pk + h*dPdt
    C1Euler=  Ck + h*dCdt

    #Improved Euler step
    Pincrement, Cincrement = f(tk+h, P1Euler, C1Euler, q,*pars)
    P1 = Pk + 0.5*h*(Pincrement + dPdt)
    C1 = Ck + 0.5*h*(Cincrement + dCdt)
    return P1,C1

def model_ODE(f, t, p0, c0, qex, qex_times, pars):
    '''
    This function takes in the ODE model, a timeframe, as well as conditions and parameters, returning a vector
    of copper concentration values

    Parameters:
    ----------
    f: callable
        derivative function 
    t: array-like
        Array of time values to generate the appropriate copper conentrations, assumed to be of equal intervals
    p0: float
        The inital pressure inside the aquifer
    c0: float
        The initial concentration of copper inside the aquifer
    qex: float / varying
        The rate of extraction of water from the aquifer
    qex_times: array-like
        The list of times matching to the extraction rates provided
    *pars: array-like
        An array of parameters used in the ODE model

    Returns:
    ----------
    P: array-like
        An array of aquifer pressures corresponding to the input time frame
    C: array-like
        An array of copper concentrations corresponding to the input time frame
    '''
    nx = len(t)		# compute number of Euler steps to take
    dt = t[1]- t[0] #computing distance between time intervals
  
    C = np.zeros(nx)					# array to store solution
    C[0] = c0
    P = np.zeros(nx)
    P[0] = p0

    for k in range(0,nx-1):
        P[k+1],C[k+1] = improved_euler_step(f, t[k], P[k], C[k], dt, qex[k], *pars) 
    return P, C

def load_data():
    ''' 
    Loads the experimental data from the .csv file into respective lists.

    Returns:
    ----------
    TimeCu : list
        Data entry for the first column of copper data
    DataCu : list
        Data entry for the second column of copper data
    TimeP : list
        Data entry for the first column of pressure data
    DataP : list
        Data entry for the second column of pressure data
    TimeQ : list
        Data entry for the first column of extraction data
    DataQ : list
        Data entry for the second column of extraction data

    Notes:
    ----------
    This function is hard coded to load in the data from:
        ac_cu.csv, ac_p.csv, ac_q.csv
    '''
    # read in the copper data, and extract respective columns
    data = np.genfromtxt('ac_cu.csv', skip_header=1, delimiter=',')
    TimeCu = list(data[0:,0])
    DataCu = list(data[0:,1])
    
    # read in the pressure data, and extract respective columns
    data = np.genfromtxt('ac_p.csv', skip_header=1, delimiter=',')
    TimeP = list(data[0:,0])
    DataP = list(data[0:,1])

    # read in the extraction data, and extract respective columns  
    data = np.genfromtxt('ac_q.csv', skip_header=1, delimiter=',')
    TimeQ = list(data[0:,0])
    DataQ = list(data[0:,1])

    # Years, CU, Years, P, Years, Q
    return TimeCu,DataCu,TimeP,DataP,TimeQ,DataQ

def solve_LPM_pressure(t, a, b, d, Mo, p1):
    '''
    Returns the modelled pressure at a time t

    Parameters:
    ----------
    t : array-like
        Times to find the pressure at
    *pars : list
        Parameters of the function

    Returns:
    ----------
    p_vals : list
        The pressures corresponding to the inputted times
    '''

    print("Trying a combination")

    TimeCu, DataCu, TimeP, DataP, TimeQ, DataQ = load_data()
    p, c = model_ODE(ODE_model, TimeQ, DataP[0], DataCu[0], DataQ, TimeQ, [a, b, d, Mo, p1])
    
    p_vals = []

    for t_val in t:
        for i in range(1979, 2020, 1):
            if t_val < i:
                top_time = i
                bottom_time = (i - 1)
                break
    
        # Linearly interpolate to find p at t
        proportion = t_val - bottom_time
        difference = p[int(top_time - 1980)] - p[int(bottom_time - 1980)]
        p_val = p[int(bottom_time - 1980)] + proportion * difference
        p_vals.append(p_val)

    return p_vals

def solve_LPM_concentration(t, a, b, d, Mo, p1):
    '''
    Returns the model concentration at a specified time t

    Parameters:
    ----------
    t : array-like
        Times to find the pressure at
    *pars : list
        Parameters of the function

    Returns:
    ----------
    p_vals : list
        The pressures corresponding to the inputted times
    '''

    print("Trying a combination")

    TimeCu, DataCu, TimeP, DataP, TimeQ, DataQ = load_data()
    p, c = model_ODE(ODE_model, TimeQ, DataP[0], DataCu[0], DataQ, TimeQ, [a, b, d, Mo, p1])
    
    c_vals = []

    for t_val in t:
        for i in range(1979, 2020, 1):
            if t_val < i:
                top_time = i
                bottom_time = (i - 1)
                break
    
        # Linearly interpolate to find p at t
        proportion = t_val - bottom_time
        difference = c[int(top_time - 1980)] - c[int(bottom_time - 1980)]
        c_val = c[int(bottom_time - 1980)] + proportion * difference
        c_vals.append(c_val)

    return c_vals

def plot_ODE(futuretimes, futureQ, futureQtimes, return_year = False, show_plot = True, pars = False):
    ''' Plot the given experimental data against the built model
    '''
    # load the data
    TimeCu, DataCu, TimeP, DataP, TimeQ, DataQ = load_data()

    # Calibrate the data with pressure
    if pars is False:
        guess_pars = [2.5052275623378775e-05, 0.050834759990897566, 44.22, 31014.71811346603, 22.55619824383838]
        bounds = (0.0, np.inf)
        pars, covar1 = curve_fit(solve_LPM_pressure, TimeP, DataP, guess_pars, maxfev = 20000, bounds=bounds)
    
        guess_pars = pars
        bounds = ([pars[0] * 0.9999, pars[1] * 0.9999, 0.0, 0.0, pars[4] * 0.9999], [pars[0] * 1.0001, pars[1] * 1.0001, np.inf, np.inf, pars[4] * 1.0001])
        pars, covar2 = curve_fit(solve_LPM_concentration, TimeCu, DataCu, guess_pars, maxfev = 20000, bounds=bounds)

        # Calculate the stds
        var1 = np.diag(covar1)
        var2 = np.diag(covar2)

        var = [var1[0], var1[1], var2[2], var2[3], var2[4]]
        std = [np.sqrt(i) for i in var]

        print(pars)
        print(std)


    # Append the future times and Q values
    Times = TimeQ
    for i in futuretimes:
        Times.append(i)

    for k in futureQ:
        DataQ.append(k)

    p, c = model_ODE(ODE_model, Times, DataP[0], DataCu[0], DataQ, TimeQ, pars)
    if return_year:
        for time in range(len(c)):
            if c[time] > 2.0:
                break

        return time, std, pars
    

    if show_plot:
        # Plotting functions
        fig,ax = plt.subplots(1,2)
        ax[0].plot(TimeCu, DataCu, 'rx', label = 'Copper data')
        ax[1].plot(TimeP, DataP, 'kx', label = 'Pressure data')
        ax[1].plot(Times, p, 'b-', label='model')
        ax[0].plot(Times, c, 'b-', label = 'model')

        ax[1].legend()
        ax[0].legend()
        ax[0].title.set_text('Concentration of Copper vs Time')
        ax[1].title.set_text('Pressure in Aquifer vs Time')

        ax[0].set_ylabel('Copper concentration (mg/litre)')
        ax[0].set_xlabel('Year')
        ax[1].set_ylabel('Pressure (MPa)')
        ax[1].set_xlabel('Year')

        plt.show()

    else:
        return Times, p, c

def plot_whatIF(futureTimes):
    ''' Plot 4 possible outcomes of the problem
         - Constant extraction of water
         - Double extraction of water 
         - Half extraction of water
         - No extraction
    '''
    # Load the data
    TimeCu, DataCu, TimeP, DataP, TimeQ, DataQ = load_data()

    # Calibrate the data with pressure
    pars, covar = curve_fit(solve_LPM_pressure, TimeP, DataP, [0.3,0.2,0.3,48000,10000,30000,160])

    Times = TimeQ
    for i in futureTimes:
        Times.append(i)

    # Case 1 : Constant extraction
    temp1 = np.array([20.0 * 1000.0 / 86400.0] * 32)
    DataQ1 = np.hstack((DataQ,temp1))

    # Case 2 : Double extraction
    temp2 = [20.0 * 1000.0 / 86400.0*2] * 32
    DataQ2 = np.hstack((DataQ,temp2))

    # Case 3 : Half extraction
    temp3 = [(20.0 * 1000.0 / 86400.0)/2] * 32
    DataQ3 = np.hstack((DataQ,temp3))

    # Case 4 : No extraction
    temp4 = [0] * 32
    DataQ4 = np.hstack((DataQ,temp4))

    p1, c1 = model_ODE(ODE_model, Times, DataP[0], DataCu[0], DataQ1, TimeQ, pars)
    p2, c2 = model_ODE(ODE_model, Times, DataP[0], DataCu[0], DataQ2, TimeQ, pars)
    p3, c3 = model_ODE(ODE_model, Times, DataP[0], DataCu[0], DataQ3, TimeQ, pars)
    p4, c4 = model_ODE(ODE_model, Times, DataP[0], DataCu[0], DataQ4, TimeQ, pars)

    # Convert data to nice units from SI units before displaying:
    # Convert concentration from kg/kg to mg/litre
    c_data_units = []
    for j in DataCu:
        c_data_units.append(j * 1000.0)
    DataCu = c_data_units

    c_units = []
    for a in c1:
        c_units.append(a * 1000.0)
    c1 = c_units[37:70]
    bestfit = c_units[0:38]

    c_units = []
    for b in c2:
        c_units.append(b * 1000.0)
    c2 = c_units[37:70]

    c_units = []
    for c in c3:
        c_units.append(c * 1000.0)
    c3 = c_units[37:70]

    c_units = []
    for d in c4:
        c_units.append(d * 1000.0)
    c4 = c_units[37:70]

    # Plotting functions
    plt.plot(TimeCu, DataCu, 'ro', label = 'Copper data')
    plt.plot(TimeQ[0:38], bestfit, 'k-', label = 'Best-fit')
    plt.plot(TimeQ[37:70], c1, 'r-', label = 'Double extraction')
    plt.plot(TimeQ[37:70], c2, 'b-', label = 'Constant extraction')
    plt.plot(TimeQ[37:70], c3, 'g-', label = 'Half extraction')
    plt.plot(TimeQ[37:70], c4, color = "orange", label = 'No extraction')
    plt.plot(TimeQ,[2.0]*71 , '--', color = "gray" , label = 'Safe Limit')
    plt.title('Concentration of Copper vs Time')
    plt.ylabel('Copper concentration (mg/litre)')
    plt.xlabel('Year')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    # Commands for predictions
    if False:
        plot_ODE(np.arange(2019,2051,1),[20.0] * 32, np.arange(2019,2051,1))

    # Commands for plotting different cases
    if False:
        plot_whatIF(np.arange(2019,2051,1))

    # Commands for plotting residuals
    if False:
        time, covar, pars = plot_ODE([], [], [], return_year=True)
        time, p, c = plot_ODE([], [], [], show_plot=False)
        TimeCu, DataCu, TimeP, DataP, TimeQ, DataQ = load_data()

        cres = []
        pres = []
        for i in range(len(TimeCu)):
            cres.append((DataCu[i] - solve_LPM_concentration([TimeCu[i]], *pars)[0]))

        for i in range(len(TimeP)):
            pres.append(DataP[i] - solve_LPM_pressure([TimeP[i]], *pars)[0])

        fig,ax = plt.subplots(2,2)
        ax[1][0].plot(TimeCu, cres, 'kx', label = 'Copper Residuals')
        ax[1][1].plot(TimeP, pres, 'kx', label = 'Pressure Residuals')
        ax[1][0].plot(TimeCu, [0 for i in TimeCu], 'k--')
        ax[1][1].plot(TimeP, [0 for i in TimeP], 'k--')

        ax[1][1].legend()
        ax[1][0].legend()

        ax[1][0].set_ylabel('Copper concentration (mg/litre)')
        ax[1][0].set_xlabel('Year')
        ax[1][1].set_ylabel('Pressure (MPa)')
        ax[1][1].set_xlabel('Year')

        ax[0][0].plot(TimeCu, [i for i in DataCu], 'rx', label = 'Copper data')
        ax[0][1].plot(TimeP, DataP, 'kx', label = 'Pressure data')
        ax[0][1].plot(TimeQ, p, 'b-', label='model')
        ax[0][0].plot(TimeQ, c, 'b-', label = 'model')

        ax[0][1].legend()
        ax[0][0].legend()
        ax[0][0].title.set_text('Concentration of Copper vs Time')
        ax[0][1].title.set_text('Pressure in Aquifer vs Time')

        ax[0][0].set_ylabel('Copper concentration (mg/litre)')
        ax[0][0].set_xlabel('Year')
        ax[0][1].set_ylabel('Pressure (MPa)')
        ax[0][1].set_xlabel('Year')

        plt.show()
    
    # Commands for plotting unsafe year graph
    if False:
        times = []
        for i in range (15,41,1):
            time, covar, pars = plot_ODE(np.arange(2019,2101,1), [i] * 82, np.arange(2019,2101,1), return_year=True)
            times.append(time)

        fig,ax = plt.subplots(1,1)
        ax.plot(np.arange(15,41,1), [j + 1980 for j in times], 'rx', label = 'Year of unsafe concentration')
        ax.legend()
        ax.title.set_text('Extraction Rate vs Year of unsafe concentration')

        ax.set_ylabel('Year of unsafe concentration')
        ax.set_xlabel('Extraction Rate (millions of litres per day)')

        plt.show()

    # Commands for resampling
    if True:
        time, std, pars = plot_ODE([],[],[],return_year=True)
        
        TimeCu,DataCu,TimeP,DataP,TimeQ,DataQ = load_data()
        
        DataCu2 = []
        for i in DataCu:
            DataCu2.append(i * 1000.0)
        DataCu = DataCu2

        pars_un = [[]]
        for i in range (100):
            for j in range (5):
                numb = np.random.randn() * std[j] + pars[j]
                pars_un[i].append(numb)
            if i != 99:
                pars_un.append([])

        times, p, c = plot_ODE([], [], [], pars=pars, show_plot=False)

        fig,ax = plt.subplots(1,1)
        ax.plot(TimeCu, DataCu,'bo',label='data')
        ax.plot(times, c, 'r-', label='best-fit')

        for pari in pars_un:
            times, p, c = plot_ODE([], [], [], pars=pari, show_plot=False)
            ax.plot(times, c, 'k-', alpha=0.2, lw=0.5)

        ax.plot([], [], 'k-', lw=0.5, label='posterior samples')
        ax.legend()

        plt.show()


    
        
    