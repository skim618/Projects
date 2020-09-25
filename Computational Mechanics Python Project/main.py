# ENGSCI263: Group 15 Project - Onehunga Aquifer
# main.py

# PURPOSE:
# To investigate how pressure and concentration behaves with varying volumes of
# water extraction in the Onehunga Aquifer. 

# SUBMISSION:
# Submission is due on the 21st September 12PM

# import modules and functions
import numpy as np
from project_functions import *

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
    This main function uses functions imported from project_functions.
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

        save_figure = True
        if not save_figure:
            plt.tight_layout()
            plt.show()
        else:
	        plt.savefig('experimental.png',dpi=300, bbox_inches = "tight")

    if predictions:
        # uncomment commands to display desired plots
        # case 1: No extraction
        plot_ODE(np.arange(2019,2051,1),[0.0] * 32, np.arange(2019,2051,1))         

        # case 2: Continued extraction
        plot_ODE(np.arange(2019,2051,1),[20.0] * 32, np.arange(2019,2051,1))        

        # case 3: Doubled extraction
        plot_ODE(np.arange(2019,2051,1),[40.0] * 32, np.arange(2019,2051,1))        

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

        save_figure = True
        if not save_figure:
            plt.tight_layout()
            plt.show()
        else:
	        plt.savefig('whatIF_scenario.png',dpi=300, bbox_inches = "tight")
        
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

        save_figure = True
        if not save_figure:
            plt.tight_layout()
            plt.show()
        else:
	        plt.savefig('residuals.png',dpi=300, bbox_inches = "tight")
    
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

        save_figure = True
        if not save_figure:
            plt.tight_layout()
            plt.show()
        else:
	        plt.savefig('unsafeYears.png',dpi=300, bbox_inches = "tight")

    if resampling:
        # obtain covariance matrix, calibrated parameters and solutions for pressure and concentration
        exceed, covar, pars = plot_ODE([],[],[],return_year=True)  

        # load in the experimental data   
        TimeCu,DataCu,TimeP,DataP,TimeQ,DataQ = load_data()

        # obtain prediction solutions of pressure and concentration for different cases
        times, p0best, c0best = plot_ODE(np.arange(2019,2051,1), [0.0] * 32, np.arange(2019,2051,1), pars=pars, show_plot=False)
        times, p20best, c20best = plot_ODE(np.arange(2019,2051,1), [20.0] * 32, np.arange(2019,2051,1), pars=pars, show_plot=False)
        times, p40best, c40best = plot_ODE(np.arange(2019,2051,1), [40.0] * 32, np.arange(2019,2051,1), pars=pars, show_plot=False)

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
        ax[0].plot(times, p0best, 'g-', label='best-fit (no extraction)')
        ax[0].plot(times, p20best, 'y-', label='best-fit (continued extraction)')
        ax[0].plot(times, p40best, 'r-', label='best-fit (doubled extraction)')
        ax[1].plot(TimeCu, DataCu,'bo',label='data')
        ax[1].plot(times, c0best, 'g-', label='best-fit (no extraction)')
        ax[1].plot(times, c20best, 'y-', label='best-fit (continued extraction)')
        ax[1].plot(times, c40best, 'r-', label='best-fit (doubled extraction)')

        # parameter samples from posterior
        ps = np.random.multivariate_normal(pars, covar, 10)

        # Initialise empty list that will contain the concentration at 2050 for each sample
        conc0 = []
        conc20 = []
        conc40 = []

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

                # Append the concentration value at 2050 to the appropriate conc list
                conc0.append(c0[-1])
                conc20.append(c20[-1])
                conc40.append(c40[-1])

        # Print a 90% confidence interval for the concentration at 2050 for each
        # of the three scenarios
        conc0.sort()
        conc20.sort()
        conc40.sort()
        print("No extraction:",c0best[-1],"+/-",(conc0[-5] - c0best[-1]))
        print("Continued extraction:",c20best[-1],"+/-",(conc20[-5] - c20best[-1]))
        print("Doubled extraction:",c40best[-1],"+/-",(conc40[-5] - c40best[-1]))


        # add labels to posterior samples and apply legends
        ax[0].plot([], [], 'g-', lw=0.5, label='posterior samples (no extraction)')
        ax[1].plot([], [], 'g-', lw=0.5, label='posterior samples (no extraction)')
        ax[0].plot([], [], 'y-', lw=0.5, label='posterior samples (continued extraction)')
        ax[1].plot([], [], 'y-', lw=0.5, label='posterior samples (continued extraction)')
        ax[0].plot([], [], 'r-', lw=0.5, label='posterior samples (doubled extraction)')
        ax[1].plot([], [], 'r-', lw=0.5, label='posterior samples (doubled extraction)')
        ax[0].legend()
        ax[1].legend()

        # add axes labels and titles
        ax[0].title.set_text('Scenario forecasts with uncertainty: Pressure')
        ax[1].title.set_text('Scenario forecasts with uncertainty: Concentration')
        ax[0].set_ylabel('Pressure (MPa)')
        ax[0].set_xlabel('Time (years)')
        ax[1].set_ylabel('Concentration (mg/litre)')
        ax[1].set_xlabel('Time (years)')

        save_figure = False
        if not save_figure:
	        plt.show()
        else:
	        plt.savefig('resampling.png',dpi=300, bbox_inches = "tight")

if __name__ == "__main__":
    main(experimental = False, predictions = False, whatIFs = False, residuals = False, unsafeYears = False, resampling = True)