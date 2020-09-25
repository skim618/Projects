# ENGSCI263: Group 15 Project - Onehunga Aquifer
# Project.py

# PURPOSE:
# This file contains all the functions necassary for creating the computer model evaluating
# the pressure and concentration of the Onehunga aquifer.

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
        plt.show()

    else:
        return Times, p, c

def main(experimental = False, predictions = False, whatIFs = False, residuals = False, unsafeYears = False, resampling = False):
    '''
    Main function of Onehunga Aquifer model. Returns plots of given commands.

    Parameters:
    ----------
    experimental : boolean
        If true, plots and displays experimental data only.
    predictions : boolean
        If true, plots and displays model predictions for a given future time list.
    whatIFs : boolean
        If true, plots and displays model predictions for different scenarios.
    residuals : boolean
        If true, plots and displays residuals of the model.
    unsafeYears : boolean
        If true, plots and displays a graph of unsafe years against extraction rates.
    resampling : boolean
        If true, plots and displays uncertainty of the model.

    Notes:
    ----------
    There are three scenarios if whatIFs is called.
        Case 1: No extraction
        Case 2: Continued extraction
        Case 3: Doubled extraction
    '''

    if experimental:
        # load in the experimental data
        TimeCu, DataCu, TimeP, DataP, TimeQ, DataQ = load_data()

        # create subplots with twin axes, and append to separate list
        fig, ax1 = plt.subplots(2, 1)
        ax2 = []
        ax2.append(ax1[0].twinx())
        ax2.append(ax1[1].twinx())

        # plot copper data with extraction rate
        ax1[0].plot(TimeCu, DataCu, 'b-', label = 'Copper data')
        ax2[0].plot(TimeQ, DataQ, 'r--', label = 'Extraction Rate')

        # plot pressure data with extraction rate
        ax1[1].plot(TimeP, DataP, 'g-', label = 'Pressure data')
        ax2[1].plot(TimeQ, DataQ, 'r--', label = 'Extraction Rate')

        # apply legends
        ax1[0].legend(loc=2)
        ax1[1].legend(loc=2)
        ax2[0].legend()
        ax2[1].legend()

        # apply axes labels and titles
        ax1[0].set_ylabel('Copper Concentration (mg/L)')
        ax1[0].set_xlabel('Time (yr)')
        ax2[0].set_ylabel('Extraction Rate (ML/day)')
        ax1[1].set_ylabel('Pressure (MPa)')
        ax1[1].set_xlabel('Time (yr)')
        ax2[1].set_ylabel('Extraction Rate (ML/day)')
        ax1[0].title.set_text("Copper Concentration vs Extraction Rate")
        ax1[1].title.set_text("Pressure vs Extraction Rate")

        # display plots
        plt.show()

    if predictions:
        # input future time into plot_ODE() to obtain future predictions
        plot_ODE(np.arange(2019,2051,1),[20.0] * 32, np.arange(2019,2051,1))

    if whatIFs:
        # load in the experimental data
        TimeCu, DataCu, TimeP, DataP, TimeQ, DataQ = load_data()

        # obtain prediction solutions of pressure and concentration for different cases
        times, p0, c0 = plot_ODE(np.arange(2019,2051,1),[0.0] * 32, np.arange(2019,2051,1), show_plot=False)        # case 1: No extraction
        times, p20, c20 = plot_ODE(np.arange(2019,2051,1),[20.0] * 32, np.arange(2019,2051,1), show_plot=False)     # case 2: Continued extraction
        times, p40, c40 = plot_ODE(np.arange(2019,2051,1),[40.0] * 32, np.arange(2019,2051,1), show_plot=False)     # case 3: Doubled extraction

        # plot predictions on concenctration model
        fig, ax = plt.subplots(1,2)
        ax[0].plot(TimeCu, DataCu, 'kx', label = 'Copper data')
        ax[0].plot(times, c0, 'g-', label = 'No Extraction')
        ax[0].plot(times, c20, 'y-', label = 'Continued Extraction')
        ax[0].plot(times, c40, 'r-', label = 'Doubled Extraction')
        ax[0].plot(times, [2.0] * len(times), 'k--')                    # copper safe level indicator

        # plot predictions on pressure model
        ax[1].plot(TimeP, DataP, 'kx', label = 'Pressure data')
        ax[1].plot(times, p0, 'g-', label='No Extraction')
        ax[1].plot(times, p20, 'y-', label='Continued Extraction')
        ax[1].plot(times, p40, 'r-', label='Doubled Extraction')

        # apply legends to both plots, and add respective titles
        ax[1].legend()
        ax[0].legend()
        ax[0].title.set_text('Concentration of Copper vs Time')
        ax[1].title.set_text('Pressure in Aquifer vs Time')

        # add axes labels
        ax[0].set_ylabel('Copper concentration (mg/litre)')
        ax[0].set_xlabel('Year')
        ax[1].set_ylabel('Pressure (MPa)')
        ax[1].set_xlabel('Year')

        plt.show()
        
    if residuals:
        # obtain covariance matrix, calibrated parameters and solutions for pressure and concentration
        time, covar, pars = plot_ODE([], [], [], return_year=True)
        time, p, c = plot_ODE([], [], [], show_plot=False)

        # load in the experimental data
        TimeCu, DataCu, TimeP, DataP, TimeQ, DataQ = load_data()

        # initialise arrays to store residuals
        cres = []
        pres = []

        # calculate concentration residuals by taking the difference, and append accordingly
        for i in range(len(TimeCu)):
            cres.append((DataCu[i] - solve_LPM_concentration([TimeCu[i]], *pars)[0]))

        # calculate pressure residuals by taking the difference, and append accordingly
        for i in range(len(TimeP)):
            pres.append(DataP[i] - solve_LPM_pressure([TimeP[i]], *pars)[0])

        # create a 2x2 subplots 
        fig,ax = plt.subplots(2,2)

        # plot residuals on the bottom row
        ax[1][0].plot(TimeCu, cres, 'kx', label = 'Copper Residuals')
        ax[1][1].plot(TimeP, pres, 'kx', label = 'Pressure Residuals')
        ax[1][0].plot(TimeCu, [0 for i in TimeCu], 'k--')
        ax[1][1].plot(TimeP, [0 for i in TimeP], 'k--')
        ax[1][1].legend()                                               # apply legends
        ax[1][0].legend()
        ax[1][0].set_ylabel('Copper concentration (mg/litre)')          # add axes labels
        ax[1][0].set_xlabel('Year')
        ax[1][1].set_ylabel('Pressure (MPa)')
        ax[1][1].set_xlabel('Year')

        # plot model with data on the first row
        ax[0][0].plot(TimeCu, [i for i in DataCu], 'rx', label = 'Copper data')
        ax[0][1].plot(TimeP, DataP, 'kx', label = 'Pressure data')
        ax[0][1].plot(TimeQ, p, 'b-', label='model')
        ax[0][0].plot(TimeQ, c, 'b-', label = 'model')
        ax[0][1].legend()                                               # apply legends
        ax[0][0].legend()
        ax[0][0].title.set_text('Concentration of Copper vs Time')      # add titles
        ax[0][1].title.set_text('Pressure in Aquifer vs Time')
        ax[0][0].set_ylabel('Copper concentration (mg/litre)')          # add axes labels
        ax[0][0].set_xlabel('Year')
        ax[0][1].set_ylabel('Pressure (MPa)')
        ax[0][1].set_xlabel('Year')

        # display plots
        plt.show()
    
    if unsafeYears:
        # initialise empty list to store unsafe years
        times = []
        error = []

        # load in the experimental data
        TimeCu, DataCu, TimeP, DataP, TimeQ, DataQ = load_data()

        # loop over varing extraction rates and append the unsafe year
        for i in range (15,41,1):
            exceed, covar, pars = plot_ODE(np.arange(2019,2101,1), [i] * 82, np.arange(2019,2101,1), return_year=True)

            # rerun samples to generate error bars
            ps = np.random.multivariate_normal(pars, covar, 10)

            # Create a list to store all the exceeding times
            exceeds = []

            for pi in ps:
                # create a new function called C3 which replicates solve_LPM_concentration but with 
                # only the concentration dependent parameters
                C2 = lambda a, b, p0p1, t, d, Mo, p0: solve_LPM_concentration(t,a,b,d,Mo,p0p1,p0)
                C3 = partial(C2, pi[0], pi[1], pi[4])
                
                # initialise uncertainty in concentration data
                sigma = [0.05] * len(DataCu)

                # automatically calibrate concentration dependent parameters to concentration data using curve_fit()
                pars2, covar2 = curve_fit(C3, TimeCu, DataCu, [45, 30000, 0], maxfev = 20000, bounds=(0, np.inf), sigma=sigma, absolute_sigma=True)

                # concentration dependent parameter samples from posterior
                pn = np.random.multivariate_normal([pars2[0], pars2[1], pars2[2]], covar2, 10)
                
                # combine the possibilities
                p_combined = np.array([pi,pi,pi,pi,pi,pi,pi,pi,pi,pi])
                for j in range(len(p_combined)):
                    p_combined[j, 2] = pn[j, 0]
                    p_combined[j, 3] = pn[j, 1]
                    p_combined[j, 5] = pn[j, 2]

                # append the unsafe year for each non-negative posterior sample in the combined list
                for p_comb in p_combined:
                    if any(p_comb < 0):
                        continue
                    exceed_uncert, covar, pars = plot_ODE(np.arange(2019,2101,1), [i] * 82, np.arange(2019,2101,1), return_year=True, pars=p_comb)
                    exceeds.append(exceed_uncert)

            # sort the years in ascending order
            exceeds.sort()

            # append the lower and upper bounds of sampled unsafe years to act as error bars
            error.append(-exceeds[4] + exceed)

            # and append the unsafe year to the time array
            times.append(exceed)

        # plot the unsafe years against extraction rates, with respective titles and lables
        fig,ax = plt.subplots(1,1)
        ax.errorbar(np.arange(15,41,1), [j + 1980 for j in times], yerr=error, fmt = '-o')      # apply error bars 
        ax.legend()                                                     # apply legend
        ax.title.set_text('Extraction Rate vs Year of unsafe concentration')
        ax.set_ylabel('Year of unsafe concentration')
        ax.set_xlabel('Extraction Rate (millions of litres per day)')

        # display plots
        plt.show()

    if resampling:
        # obtain covariance matrix, calibrated parameters and solutions for pressure and concentration
        exceed, covar, pars = plot_ODE([],[],[],return_year=True)  

        # load in the experimental data   
        TimeCu,DataCu,TimeP,DataP,TimeQ,DataQ = load_data()

        # obtain prediction solutions of pressure and concentration for different cases
        times, p0, c0 = plot_ODE(np.arange(2019,2051,1), [0.0] * 32, np.arange(2019,2051,1), pars=pars, show_plot=False)
        times, p20, c20 = plot_ODE(np.arange(2019,2051,1), [20.0] * 32, np.arange(2019,2051,1), pars=pars, show_plot=False)
        times, p40, c40 = plot_ODE(np.arange(2019,2051,1), [40.0] * 32, np.arange(2019,2051,1), pars=pars, show_plot=False)

        # append extrapolated times
        times_ex = copy(times)
        for i in (np.arange(2019,2051,1)):                        
            TimeQ.append(i)

        # append doubled extraction rates
        DataQ0 = copy(DataQ)
        DataQ20 = copy(DataQ)
        DataQ40 = copy(DataQ)
        for i in ([40.0] * 32):                        
            DataQ40.append(i)
            DataQ0.append(0.0)
            DataQ20.append(i / 2.0)

        # plot different cases of concentration and pressure models with their respective data
        fig,ax = plt.subplots(1,2)
        ax[0].plot(TimeP, DataP,'bo',label='data')
        ax[0].plot(times, p0, 'g-', label='best-fit (no extraction)')
        ax[0].plot(times, p20, 'y-', label='best-fit (continued extraction)')
        ax[0].plot(times, p40, 'r-', label='best-fit (doubled extraction)')
        ax[1].plot(TimeCu, DataCu,'bo',label='data')
        ax[1].plot(times, c0, 'g-', label='best-fit (no extraction)')
        ax[1].plot(times, c20, 'y-', label='best-fit (continued extraction)')
        ax[1].plot(times, c40, 'r-', label='best-fit (doubled extraction)')

        # parameter samples from posterior
        ps = np.random.multivariate_normal(pars, covar, 10)

        for pi in ps:
            # create a new function called C3 which replicates solve_LPM_concentration but with 
            # only the concentration dependent parameters
            C2 = lambda a, b, p0p1, t, d, Mo, p0: solve_LPM_concentration(t,a,b,d,Mo,p0p1,p0)
            C3 = partial(C2, pi[0], pi[1], pi[4])
            
            # initialise uncertainty in concentration data
            sigma = [0.05] * len(DataCu)

            # automatically calibrate concentration dependent parameters to concentration data using curve_fit()
            pars2, covar2 = curve_fit(C3, TimeCu, DataCu, [45, 30000, 0], maxfev = 20000, bounds=(0, np.inf), sigma=sigma, absolute_sigma=True)

            # concentration dependent parameter samples from posterior
            pn = np.random.multivariate_normal([pars2[0], pars2[1], pars2[2]], covar2, 10)
            
            # combine the possibilities
            p_combined = np.array([pi,pi,pi,pi,pi,pi,pi,pi,pi,pi])
            for j in range(len(p_combined)):
                p_combined[j, 2] = pn[j, 0]
                p_combined[j, 3] = pn[j, 1]
                p_combined[j, 5] = pn[j, 2]

            # plot the solutions for each non-negative posterior sample in the combined list
            for p_comb in p_combined:
                if any(p_comb < 0):
                    continue
                p0, c0 = model_ODE(ODE_model, times_ex, DataP[0], DataCu[0], DataQ0, TimeQ, p_comb)
                p20, c20 = model_ODE(ODE_model, times_ex, DataP[0], DataCu[0], DataQ20, TimeQ, p_comb)
                p40, c40 = model_ODE(ODE_model, times_ex, DataP[0], DataCu[0], DataQ40, TimeQ, p_comb)
                ax[0].plot(times_ex, p0, 'g-', alpha=0.2, lw=0.75)
                ax[1].plot(times_ex, c0, 'g-', alpha=0.2, lw=0.75)
                ax[0].plot(times_ex, p20, 'y-', alpha=0.2, lw=0.75)
                ax[1].plot(times_ex, c20, 'y-', alpha=0.2, lw=0.75)
                ax[0].plot(times_ex, p40, 'r-', alpha=0.2, lw=0.75)
                ax[1].plot(times_ex, c40, 'r-', alpha=0.2, lw=0.75)

        # add labels to posterior samples and apply legends
        ax[0].plot([], [], 'g-', lw=0.5, label='posterior samples (no extraction)')
        ax[1].plot([], [], 'g-', lw=0.5, label='posterior samples (no extraction)')
        ax[0].plot([], [], 'y-', lw=0.5, label='posterior samples (continued extraction)')
        ax[1].plot([], [], 'y-', lw=0.5, label='posterior samples (continued extraction)')
        ax[0].plot([], [], 'r-', lw=0.5, label='posterior samples (doubled extraction)')
        ax[1].plot([], [], 'r-', lw=0.5, label='posterior samples (doubled extraction)')
        ax[0].legend()
        ax[1].legend()

        # display plots
        plt.show()

if __name__ == "__main__":
    main(experimental = False, predictions = False, whatIFs = False, residuals = False, unsafeYears = True, resampling = False)


    
        
    