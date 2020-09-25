Ongehunga Aquifer Model
Nick Wright, Biyuan Chen, Hyunbin Ko, Sooyong Kim

What is this project used for?
Auckland’s drinking water is sourced from large reservoirs in the Waitākere and Hūnua
Ranges, the Waikato River, and groundwater aquifers beneath the city. In particular, the
Onehunga Aquifer, an old fractured lava flow that contains fresh water in its pore space, has
been permitted to supply up to 20 million litres per day since 2000 (Auckland demand is
about 400 million litres per day).

Due to the growing drought in the Auckland area – an increase uptake from the Ongehunga
Aquifer is being considered. However, there are concerns that an increase in Water Uptake
will result in lower pressures within the aquifer, and therefore drawing surrounding waters
with higher copper concentration into the aquifer.

This project contains code designed to model the Onehunga Aquifer in terms of the overall
pressure and copper concentration.

What does the project contain?
The project contains three files: main.py, project_functions.py and unit_test.py. project_functions.py
contains the functions required to simulate the aquifer. The following is a list of functions contained
in the file and their purpose:
 ODE_model: This function takes a given time, concentration, pressure, and parameters. It
then returns the rate of change of both the concentration and the pressure, calculated
through the following ODEs shown in the project report. For further information as to how
these ODEs are derived, see the project report.

 improved_euler_step: This function takes all parameters and variables at a given time and
utilizing the ODE_model function and improved euler method, returns the predicted
concentration and pressure at a specified time step in the future.

 model_ODE: This function takes an array of values representing the extraction rate over
time and the respective times, along with initial values for copper concentration and
pressure, and all the parameters. It also required an arbitrary array of times for which to
solve the ODE. It returns the estimated concentration and pressure values at these times,
based on the improved euler method.

 Load_data: This function requires the files, ac_q.txt, ac_p.txt, and ac_cu.txt to be
downloaded and be in the working repository. This function is responsible for reading in
these sets of data into lists, in order to execute calibration. It creates a list of copper
concentration, pressure, and extraction rate from the year 1980 to 2018.

 solve_LPM_pressure(concentration): This function returns a set of pressure (concentration)
values only, depending on the inputted times. This is important when using curve_fit
calibration (see plot_ODE)

 plot_ODE: This function is designed to plot the data for concentration and pressure, along
with a calibrated model. This function has several Boolean variables that can be turned on
and off when called in the main function. Automatically, the function will simply show the
plots of these two sets of data. If return_year is True, then the function will instead return 3
outputs: the year at which the copper data exceeds 2mg/L, the covariance matrix of the
parameters and the parameters that were calibrated to the data. If show_plot is False, then
the values for times, pressures, and concentrations will be returned. The boolean pars will 
calibrate the parameters if set to false.

The main.py is responsible for running functions created in project_functions.py and acts as the main 
simulation of the model. It contains one main function: 

 main: This function controls the different plottings of the model according to various boolean inputs:
The range of actions which can be taken in the main function are controlled through turning the respective
boolean to true:

- experimental: Outputs plots outlining the relationship between pressure, concentration and extraction rate 
using the handout data.
- predictions: plots the model predictions using an input time list.
- whatIFs: plots the model predictions for the four outlined scenarios of the resource consent.
- residuals: plots the residual errors of the best fit copper and pressure models. 
- unsafeYears: plots the year where the copper concentration exceeds the unsafe limit of 2mg/L against varying 
extraction rates.
- resampling: plots the best fit concentration and pressure models along with its posterior samples against the
handout data.

All unit tests and benchmarking are implemented in unit_test.py. 

