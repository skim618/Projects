# Projects

## Group Python Project: Truck Scheduling Project (2020)

The Warehouse Group operates The Warehouse and Noel Leeming stores. Each store needs to receive goods daily to ensure their shelves are fully-stocked daily. They operate a fleet of 25 trucks in order to move these goods from their two distribution centre to their stores around Auckland.
On each day, each store receives pallets of goods from a distribution centre based on store sales. Therefore, the number of pallets shipped to each store differs each day. For this model, we will work in units of pallets, and we will not differentiate between different product categories.

Each truck can carry up to 20 pallets of goods, and operates on a trip schedule that will have each truck deliver goods to a selection of stores, and return to the warehouse. Once at the store, a pallet takes on average 10 minutes to unload. Current policy requires each scheduled trip take no more than four hours, on average, to complete; this includes both driving time and unloading time. Each truck costs $175 per hour to operate and can operate two (approximately) four-hour shifts per day. You may assume that the two shifts start at 8am or 2pm, and that each store only receives one delivery per day.
However, traffic conditions on Auckland roads are not always ideal, so the driving time required may well be longer or shorter depending on the time of day. This means some trucks may take more than four hours to complete their trip. In such cases, the extra time costs Foodstuffs $250 per hour.

On days where there are not sufficient trucks to satisfy all demand, either because of a shortage of truck time or an excess in store demand for pallets, additional trucks can be ‘wet-leased’ (vehicle rental that includes a driver) from Mainfreight for a cost of $1500 for every four hours of on-duty time, charged in four-hour blocks.
The Warehouse Group would like to determine a suitable truck logistics plan such that costs are minimised. 

### How do I get set up? ###

* createRoutes.py creates a set of text files that contain all potentially feasible routes. Run first. 
* linearProgram.py is the mixed integer program that solves and finds the least costing routes.
* routeVisualisation.ipynb is a Python Notebook that contains the visualisation of the least costing routes.
* simulation.py is the simulation of potential variations within the previously mentioned least-cost routes.
* warehousePartition.R visualises how the routes were partitioned.

## Group Python Project: Resource Consent Project for the Onehunga Aquifer (2020)

The Onehunga Aquifer has been a significant source of water in the Auckland region, providing up to 20 million litres per day since the year 2000. However, due to the growing drought in the Auckland region, Watercare has proposed to increase the water extraction from the aquifer to satisfy increased demands. There are potential environmental impacts concerned with this proposal:
- Pressure decline in the aquifer: As the rate of water extraction increases, the pressure in the aquifer decreases as water flows across the aquifer under a pressure gradient.
- Contamination of water: Contaminants, such as copper, dissolve into stormwater at the surface and then leaches down to pressure lows in ground water. The contamination correlates to the pressure decline as a result of increased water extraction. Water extracted can therefore be unsafe to drink if the level of copper exceeds a certain threshold.

Watercare have applied to the Auckland Regional Council (ARC) for resource consent to double their maximum take from its current limit of 20 million litres per day. They claim that the aquifer can support the increased take without copper concentrations exceeding safe levels.

This project undertakes a computer modelling study on groundwater quality of the Onehunga Aquifer addressing both pressure and copper concentration changes. A coupled lumped parameter model of concentration and copper is formulated and solved numerically using improved Euler's method. The parameters are then calibrated automatically using Python's built-in curve_fit function. The model is then used to predict into the future, for four different possible outcomes of the consent hearing. This includes approval of the expanded consent (a maximum of 40 million litres per day), no change, reduced consent (to a recommended safe level) or, more drastically, an indefinite moratorium on aquifer usage. A recommendation is then given in the final report. 

### Source files
The project contains three files: main.py, project_functions.py and unit_test.py. 
* project_functions.py contains the functions required to simulate the aquifer.
* The main.py is responsible for running functions created in project_functions.py and acts as the main simulation of the model.
* All unit tests and benchmarking are implemented in unit_test.py. 


## Python Project: Vaccine Distribution Project: (2020)

At long last, after many months of New Zealand’s isolation from the world, a
vaccine for the novel coronavirus has been developed! A first batch will soon be
arriving at Auckland Airport from the overseas facility that produced it, and this project
addresses the distribution of the vaccines to the most vulnerable population in
Auckland: the elderly living in residential care settings. A team of
four couriers available to do this distribution, and the job is to get the vaccine
delivered from the airport to all 143 rest homes in the city, as quickly as possible.
The project's aim is to figure out how to deliver the vaccines using the public transportation system.

This project uses a nearest-neighbour algorithm that returns four
ordered lists of network node names, representing the paths taken by the four
couriers. 

### Source files
The project contains three files: project_code.py, project_skim618.py and project_utils.py. 
* project_code.py contains the functions and algorithms required to partition the routes
* project_skim618.py contains the main function using the functions in project_code.py. 
* project_utils.py contains useful funcitosn for reading and writing out data in .txt format. 


## C Project (2019)

This project requires the writing of C code that will create a game based on the theme of 
a Warehouse. The aim of the game is to move the boxes in the warehouse to designated spots by
controlling a worker that will push the box upon contact. There are three key files:
- project.c: This is the source file. In this source file is all the funcitons written to run the game.
- project.h :This is the header file that contains the prototype declarations for the
functions you have to write. Both source files (project.c and test_project.c) 
include this header file.
- test_project.c: This is the source file that contains the main() function. This file is used to
launch the game.

## MATLAB Project (2019)
This project requires the writing of MATLAB code that will allow us to create two different kinds of images
from an arbitrary sequence of images that feature someone (or something) moving through a static
scene. These images could be from a series of photos or a selection of frames taken from a movie.

One kind of image that will be created is known as an “Action shot”, which is a single image that features
the moving person or object shown in different locations. These kinds of images are sometimes
used in extreme sport magazines, in order to show people moving over a jump or performing some
other kind of action. See https://en.wikipedia.org/wiki/ActionShot for a few examples. Action shots
are also popular in astronomy to show the motion of celestial objects such as the moon, as they
track across the sky, as shown in the eclipse scene below.

The other kind of image that will be created is one where the action has been removed, eliminating the
moving person (or object) from the image, leaving only the static scene visible. A static scene could
be useful if we want to remove a pesky pedestrian from an image.


