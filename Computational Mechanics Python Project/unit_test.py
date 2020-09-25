from project_functions import *
from matplotlib import pyplot as plt
import numpy as np

# This file contains the unit tests and benchmark functions for project_copy.py

def test_ODE_model():
    '''
    This function tests the functionally of the ODE_model() function in project_copy.py file.

    #####ODE_model(t, p, c, q, a, b, d, Mo, p0p1, p0)#####
    ODE_model() function should return correct values of dp/dt and dc/dt based on the formula given in
    https://canvas.auckland.ac.nz/courses/46529/pages/project-formulations?module_item_id=947913.
    '''

    # Case 1: When p > p0, C' = c
    dpdt, dcdt = ODE_model(1980, 5000, 400, 20, 0.3,0.2,0.3, 50 ,0.07,0.03) 
    assert abs(dpdt - (-2005.986)) < 1e-3
    assert abs(dcdt - (26636.45354)) < 1e-5
    
    # Case 2: When p < p0, C' = 0
    dpdt, dcdt = ODE_model(1980, 0.03, 400, 20, 0.3,0.2,0.3, 50 ,11000,5000)
    assert abs(dpdt - (2193.988))< 1e-3
    assert abs(dcdt - (-58633.34685)) <1e-5

    # Case 3: When p = p0, C' = 0
    dpdt, dcdt = ODE_model(1980, 5000, 400, 20, 0.3,0.2,0.3, 50 ,11000,5000) 
    assert dpdt == 194
    assert abs(dcdt - (-5330.333333)) < 1e-6

def test_improved_euler_step():
    '''
    This function tests the functionally of the improved_euler_step() function in project_copy.py file.

    #####improved_euler_step(f, tk, Pk, Ck, h, q, *pars)#####
    improved_euler_step() function should return correct values of p and c after a single imporved euler step.
    '''
    # When p > p0, C' = c
    p,c = improved_euler_step(ODE_model,1980,5000,400,0.5,20,0.3,0.2,0.3, 50, 0.07, 0.03)
    assert abs(p - (4097.3063)) < 1e-4
    assert abs(c - (189824.1169))< 1e-4

    # When p < p0, C' = 0
    p,c = improved_euler_step(ODE_model,1980,0.03,400,0.5,20,0.3,0.2,0.3, 50 ,11000, 5000)
    assert abs(p - (987.3246)) < 1e-4
    assert abs(c - (834544.395)) < 1e-3

    # When p = p0, C' = 0
    p,c = improved_euler_step(ODE_model,1980,5000,400,0.5,20,0.3,0.2,0.3, 50 ,11000, 5000)
    assert abs(p - (5087.3)) < 1e-1
    assert abs(c - (5886.17283)) < 1e-5

    # When p > p0 -> p < p0,  C' = c -> C' = 0
    p,c = improved_euler_step(ODE_model,1980,0.5,400,5,20,0.3,0.2,0.3, 50, 1, 0.45)
    assert abs(p - (0.51)) < 1e-2
    assert abs(c - (-397.55))< 1e-2

    # When p < p0 -> p > p0,  C' = 0 -> C' = c 
    p,c = improved_euler_step(ODE_model,1980,0.45,400,5,20,0.001,0.2,0.3, 50, 1.1,0.5)
    assert abs(p - (0.451)) < 1e-3
    assert abs(c - (200.1))< 1e-1

    # When p = p0 -> p > p0,  C' = 0 -> C' = c
    p,c = improved_euler_step(ODE_model,1980,0.45,400,5,20,0.001,0.2,0.3, 50, 1.1,0.45)
    assert abs(p - (0.451)) < 1e-3
    assert abs(c - (799.99851))< 1e-5

    # When step size = 0
    p,c = improved_euler_step(ODE_model,1980,5000,400,0,20,0.3,0.2,0.3, 50, 0.04, 0.03)
    assert p == 5000
    assert c == 400

def test_load_data():
    '''
    This function tests the functionally of the load_data() function in project_copy.py file.
    load_data() function should retreive correct values from csv files as arrays.
    '''

    TimeCu,DataCu,TimeP,DataP,TimeQ,DataQ = load_data()

    # Get the time period of copper concentrations in years
    assert(TimeCu[0] == 1980)
    assert(TimeCu[-1] == 2015)

    # Get the copper concentration data  during the time period above in mg/litre
    assert(DataCu[0] == 0.0000000e+00)
    assert(DataCu[-1] == 9.7636318e-01)

    # Get the time period of pressures in years
    assert(TimeP[0] == 1980)
    assert(TimeP[-1] == 2016)

    # Get the pressure data during the time period above in MPa
    assert(DataP[0] == 3.1290113e-02)
    assert(DataP[-1] == -2.9783908e-02)

    # Get the time period of extractions in years
    assert(TimeQ[0] == 1980)
    assert(TimeQ[-1] == 2018)

    # Get the extraction data during the time period above in 10^6 litre/day
    assert(DataQ[0] == 3.15e+01)
    assert(DataQ[-1] == 1.72e+01)

def benchmark_model_ODE():
    '''
    This function plots numerical solution and analytical solution of the pressure ODE
    to test if mode_ODE function in project_copy.py file has been implemented correctly.

    Parameters used
    ---------------
    P(0) = 0.03 MPa
    c0 = 0 mg/Litre
    q = 31.5 ML/day
    a = 0.3 
    b = 0.2 
    d = 0.3
    M0 = 50000
    P0 = 0.01 MPa
    P1 = 0.03 MPa
    '''
    # Load data
    TimeCu, DataCu, TimeP, DataP, TimeQ, DataQ = load_data()

    # Initialise subplots
    fig,ax = plt.subplots(1,2)

    # Plot numerical solution for pressure
    p, c = model_ODE(ODE_model, np.linspace(1980,2016,100), 0.03, 0, [31.5]*100, np.linspace(1980,2016,100), [0.3,0.2,0.3,50000,0.04,0.01])
    ax[0].plot(np.linspace(1980,2016,100), p, 'kx', label='Numerical Pressure')

    # Plot analytical solution for pressure
    analytical_p = 23.635 * np.exp(-0.4 * np.linspace(0.,36,100)) - 23.605
    ax[0].plot(np.linspace(1980,2016,100), analytical_p, 'g', label='Analytical Pressure')

    ax[0].set_title('Pressure ODE benchmark')
    ax[0].set_xlabel('Time (Years)')
    ax[0].set_ylabel('Pressure (MPa)')
    ax[0].legend()

    # Plot numerical solution for concentration
    p, c = model_ODE(ODE_model, np.linspace(1980,2016,100), 0.03, 10.0, [31.5]*100, np.linspace(1980,2016,100), [1.0,0.0,1.0,100000,1.0,1.0])
    ax[1].plot(np.linspace(1980,2016,100), c, 'kx', label='Numerical Concentration')

    # Plot analytical solution for concentration
    analytical_c = 0.0001575 * (np.linspace(0.0,36,100))**2 + 10.0
    ax[1].plot(np.linspace(1980,2016,100), analytical_c, 'g', label='Analytical Concentration')

    ax[1].set_title('Concentration ODE benchmark')
    ax[1].set_xlabel('Time (Years)')
    ax[1].set_ylabel('Concentration (mg/Litre)')
    ax[1].legend()
    
    save_figure = True
    if not save_figure:
        plt.show()
    else:
	    plt.savefig('benchmark.png',dpi=300, bbox_inches = "tight")

def test_correct_time_interval():
    '''
    This function plots solution with current time interval and solutions with different time intervals to test
    if an appropariate time interval has been chosen for model_ODE function in project_copy.py file.

    Parameters used
    ---------------
    P(0) = 0.03 MPa
    c0 = 0 mg/Litre
    q = 31.5 ML/day
    a = 0.3 
    b = 0.2 
    d = 0.3
    M0 = 50000
    P0 = 0.01 MPa
    P1 = 0.03 MPa
    '''
    #close any open figures
    plt.close()
    # Load data
    TimeCu, DataCu, TimeP, DataP, TimeQ, DataQ = load_data()


    # Plot numeric solutions with varying time intervals
    for i in range(10,100):
        Times = np.linspace(1980,2016,i)
        q = np.interp(Times,TimeQ,DataQ)
        p, c = model_ODE(ODE_model, Times, 0.01, 0, [0]*i*10, q, [0.3,0.2,0.3,50000,0.04,0.01])
        if i == 37:
            plt.plot(Times, p, 'kx',lw=0.5, alpha = 1, label = 'Chosen time interval for our model (dt=1)')
        else:
            plt.plot(Times, p, 'b-',lw=0.1, alpha = 1)

    # Plot numerical solution against analytical solution
    plt.title('Pressure ODE solution with different time intervals')
    plt.xlabel('Years')
    plt.ylabel('Pressure')
    plt.legend()
    
    save_figure = True
    if not save_figure:
	    plt.show()
    else:
	    plt.savefig('timestep.png',dpi=300, bbox_inches = "tight")

if __name__ == "__main__":
    # Unit test for ODE_model in project_copy.py file
    if True:
        test_ODE_model()

    # Unit test for improved_euler_step in project_copy.py file
    if True:
        test_improved_euler_step()

    # Unit test for load_data in project_copy.py file
    if True:
        test_load_data()

    # Benchmark for pressure ODE
    if True:
        benchmark_model_ODE()

    # Test for appropariate time interval for pressure ODE
    if True:
        test_correct_time_interval()
