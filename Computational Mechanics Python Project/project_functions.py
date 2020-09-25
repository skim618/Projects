# ENGSCI263: Group 15 Project - Onehunga Aquifer
# project_functions.py

# PURPOSE:
# To develop a computer model for the Onehunga Aquifer

# SUBMISSION:
# Submission is due on the 21st September 12PM

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
from functools import partial

def ODE_model(t, p, c, q, a, b, d, Mo, p0p1, p0):
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
    p0p1 : float
        Sum of pressure at low and high boundary.
    p0 : float
        Pressure at low pressure boundary.

    Returns:
    ----------
    dpdt : float
        Derivative of pressure with respect to independent variable.
    dcdt : float
        Derivative of concentration with respect to independent variable.

    Notes:
    ----------
    Concentration depends on direction flow.
    '''

    # concentration depends on direction flow
    if (p > p0):
        c_dash = c
    else:
        c_dash = 0
    
    # derivative equations for pressure and concentration respectively
    dpdt = -a*q - 2*b*p + b*p0p1                      
    dcdt = (((-b/a)*(p - p0)*(c_dash-c)) + ((b/a)*(p-(p0p1-p0))*c) - d*(p-((p0p1)/2.))) / Mo

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
        Extraction rate.
    *pars: array-like
        An array of parameters used in the ODE model.

	Returns:
	----------
	P1 : float
		Solution for pressure at end of the improved Euler step.
    C1 : float
		Solution for concentration at end of the improved Euler step.

    Notes:
    ----------
    Length of *pars must be equal to the length of paramenters in the called derivative function.
    *pars is assumed to be in the order of [a, b, d, Mo, p0p1, p0]
	'''	
    # euler step
    dPdt,dCdt = f(tk,Pk,Ck,q,*pars)
    P1Euler=  Pk + h*dPdt
    C1Euler=  Ck + h*dCdt

    # improved Euler step
    Pincrement, Cincrement = f(tk+h, P1Euler, C1Euler, q,*pars)
    P1 = Pk + 0.5*h*(Pincrement + dPdt)
    C1 = Ck + 0.5*h*(Cincrement + dCdt)
    return P1,C1

def model_ODE(f, t, p0, c0, qex, qex_times, pars):
    '''
    This function numerically solves the ODEs in the derivative function 
    and returns a vector of solutions at given times. 

    Parameters:
    ----------
    f: callable
        Derivative function.
    t: array-like
        Array of time values to generate the appropriate copper conentrations, 
        assumed to be of equal intervals.
    p0: float
        The inital pressure inside the aquifer.
    c0: float
        The initial concentration of copper inside the aquifer.
    qex: float / varying
        The rate of extraction of water from the aquifer.
    qex_times: array-like
        The list of times matching to the extraction rates provided.
    *pars: array-like
        An array of parameters used in the ODE model.

    Returns:
    ----------
    P: array-like
        An array of aquifer pressures corresponding to the input time frame.
    C: array-like
        An array of copper concentrations corresponding to the input time frame.

    Notes:
    ----------
    Length of *pars must be equal to the length of paramenters in the called derivative function.
    *pars is assumed to be in the order of [a, b, d, Mo, p0p1, p0]
    '''
    nx = len(t)		                    # compute number of Euler steps to take
    dt = t[1]- t[0]                     # computing distance between time intervals

    # initialise arrays to store solutions
    C = np.zeros(nx)				
    C[0] = c0
    P = np.zeros(nx)
    P[0] = p0

    # use improved Euler to solve at each step
    for k in range(0,nx-1):
        P[k+1],C[k+1] = improved_euler_step(f, t[k], P[k], C[k], dt, qex[k], *pars) 
    return P, C

def load_data():
    ''' 
    Loads the experimental data from the .csv file into respective lists.

    Returns:
    ----------
    TimeCu : list
        Data entry for the first column of copper data.
    DataCu : list
        Data entry for the second column of copper data.
    TimeP : list
        Data entry for the first column of pressure data.
    DataP : list
        Data entry for the second column of pressure data.
    TimeQ : list
        Data entry for the first column of extraction data.
    DataQ : list
        Data entry for the second column of extraction data.

    Notes:
    ----------
    This function is hard coded to load in the data from:
        ac_cu.csv, ac_p.csv, ac_q.csv
    '''
    # read in the copper data, extract respective columns, and convert into a list
    data = np.genfromtxt('ac_cu.csv', skip_header=1, delimiter=',')
    TimeCu = list(data[0:,0])
    DataCu = list(data[0:,1])
    
    # read in the copper data, extract respective columns, and convert into a list
    data = np.genfromtxt('ac_p.csv', skip_header=1, delimiter=',')
    TimeP = list(data[0:,0])
    DataP = list(data[0:,1])

    # read in the copper data, extract respective columns, and convert into a list  
    data = np.genfromtxt('ac_q.csv', skip_header=1, delimiter=',')
    TimeQ = list(data[0:,0])
    DataQ = list(data[0:,1])

    return TimeCu,DataCu,TimeP,DataP,TimeQ,DataQ

def solve_LPM_pressure(t, a, b, d, Mo, p0p1, p0):
    '''
    Returns a list of modelled pressures at a given time.

    Parameters:
    ----------
    t : array-like
        Times at which the pressure is to be solved for.
    a : float
        Source/sink strength parameter.
    b : float
        Recharge strength parameter.
    d : float
        Surface leeching parameter.
    Mo : float
        Total system mass.
    p0p1 : float
        Sum of pressure at low and high boundary.
    p0 : float
        Pressure at low pressure boundary.

    Returns:
    ----------
    p_vals : list
        List of modelled pressures at the given times.
    '''

    # load in the experimental data
    TimeCu, DataCu, TimeP, DataP, TimeQ, DataQ = load_data()

    # solve for pressure and concentration with given parameters
    p, c = model_ODE(ODE_model, TimeQ, DataP[0], DataCu[0], DataQ, TimeQ, [a, b, d, Mo, p0p1, p0])
    
    # initialise a list to store the pressure values
    p_vals = []

    # nested for loop to determine top and bottom times for interpolation
    for t_val in t:
        for i in range(1979, 2020, 1):
            if t_val < i:
                top_time = i
                bottom_time = (i - 1)
                break
    
        # linearly interpolate to find the pressure at the given time, and append it to the list
        proportion = t_val - bottom_time
        difference = p[int(top_time - 1980)] - p[int(bottom_time - 1980)]
        p_val = p[int(bottom_time - 1980)] + proportion * difference
        p_vals.append(p_val)

    return p_vals

def solve_LPM_concentration(t, a, b, d, Mo, p0p1, p0):
    '''
    Returns a list of modelled concentrations at a given time.

    Parameters:
    ----------
    t : array-like
        Times at which the pressure is to be solved for.
    a : float
        Source/sink strength parameter.
    b : float
        Recharge strength parameter.
    d : float
        Surface leeching parameter.
    Mo : float
        Total system mass.
    p0p1 : float
        Sum of pressure at low and high boundary.
    p0 : float
        Pressure at low pressure boundary.

    Returns:
    ----------
    c_vals : list
        List of modelled concentrations at the given times.
    '''

    # load in the experimental data
    TimeCu, DataCu, TimeP, DataP, TimeQ, DataQ = load_data()

    # solve for pressure and concentration with given parameters
    p, c = model_ODE(ODE_model, TimeQ, DataP[0], DataCu[0], DataQ, TimeQ, [a, b, d, Mo, p0p1, p0])
    
    # initialise a list to store the concentration values
    c_vals = []

    # nested for loop to determine top and bottom times for interpolation
    for t_val in t:
        for i in range(1979, 2020, 1):
            if t_val < i:
                top_time = i
                bottom_time = (i - 1)
                break
    
        # linearly interpolate to find the pressure at the given time, and append it to the list
        proportion = t_val - bottom_time
        difference = c[int(top_time - 1980)] - c[int(bottom_time - 1980)]
        c_val = c[int(bottom_time - 1980)] + proportion * difference
        c_vals.append(c_val)

    return c_vals

def plot_ODE(futuretimes, futureQ, futureQtimes, return_year = False, show_plot = True, pars = False):
    ''' 
    Plots the given experimental data against the built model.

    Parameters:
    ----------
    futuretimes : list
        A list of future times to solve for.
    futureQ : list
        A list of future extraction rates.
    futureQtimes : list
        A list of times corresponding to future extraction rates.
    return_year : boolean
        If true, returns the time, covariance matrix, and calibrated parameters.
    show_plot : boolean
        If ture, displays the plot.
    pars : boolean
        If false, calibrate parameters.

    Returns:
    ----------
    exceed : int
        index at which the concentration exceeds the safe level (2.0mg/litre).
    Times : list
        List of times at which solutions are solved for.
    covar1 : array-like
        Covariance matrix of parameters.
    pars : list
        List of calibrated parameters.
    p: array-like
        An array of aquifer pressures corresponding to the input time frame.
    c: array-like
        An array of copper concentrations corresponding to the input time frame.

    Notes:
    ----------
    When return_year is true, this function will return time, covar1, and pars.
    Else, it will return Times, p and c. 
    pars is assumed to be in the order of [a, b, d, Mo, p0p1, p0]
    '''
    # load in the experimental data
    TimeCu, DataCu, TimeP, DataP, TimeQ, DataQ = load_data()
    covar1 = []

    if pars is False:
        # initialise uncertainty in pressure data, guess parameters, and parameter bounds
        sigma = [0.25] * len(DataP)
        guess_pars = [2.5052275623378775e-05, 0.050834759990897566, 3.22, 31014.71811346603, 5.5619824383838, 0]
        bounds = (0.0, np.inf)

        # automatically calibrate parameters to pressure data using curve_fit()
        pars, covar1 = curve_fit(solve_LPM_pressure, TimeP, DataP, guess_pars, maxfev = 20000, bounds=bounds, sigma=sigma, absolute_sigma=False)

        # create a new function called C3 which replicates solve_LPM_concentration but with 
        # only the concentration dependent parameters
        C2 = lambda a, b, p0p1, t, d, Mo, p0: solve_LPM_concentration(t,a,b,d,Mo,p0p1,p0)
        C3 = partial(C2, pars[0], pars[1], pars[4])

        # initialise uncertainty in concentration data
        sigma = [0.05] * len(DataCu)

        # automatically calibrate concentration dependent parameters to concentration data using curve_fit()
        pars2, covar2 = curve_fit(C3, TimeCu, DataCu, [0, 30000, 0], maxfev = 20000, bounds=bounds, sigma=sigma, absolute_sigma=True)

        # replace concentration dependent parameters in the initial parameters
        pars[2] = pars2[0]
        pars[3] = pars2[1]
        pars[5] = pars2[2]

    # append future times and extraction rates
    Times = TimeQ
    for i in futuretimes:
        Times.append(i)
    for k in futureQ:
        DataQ.append(k)

    # numerically solve ODEs using model_ODE with the calibrated parameters
    p, c = model_ODE(ODE_model, Times, DataP[0], DataCu[0], DataQ, TimeQ, pars)

    if return_year:
        # determine the index at which the concentration exceeds the safe level
        for exceed in range(len(c)):
            if c[exceed] > 2.0:
                break
        return exceed, covar1, pars
    
    if show_plot:
        # create two subplots, one for pressure and one for concentration
        fig, ax = plt.subplots(1,2)

        # plot concentration data against concentration model, with respective labels and titles
        ax[0].plot(TimeCu, DataCu, 'rx', label = 'Copper data')
        ax[0].plot(Times, c, 'b-', label = 'model')
        ax[0].title.set_text('Concentration of Copper vs Time')
        ax[0].set_ylabel('Copper concentration (mg/litre)')
        ax[0].set_xlabel('Year')
        ax[0].legend()

        # plot pressure data against pressure model, with respective labels and titles
        ax[1].plot(TimeP, DataP, 'kx', label = 'Pressure data')
        ax[1].plot(Times, p, 'b-', label='model')
        ax[1].title.set_text('Pressure in Aquifer vs Time')
        ax[1].set_ylabel('Pressure (MPa)')
        ax[1].set_xlabel('Year')
        ax[1].legend()

        # show plot
        plt.tight_layout()
        plt.show()

    else:
        return Times, p, c



    
        
    