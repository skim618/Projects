# Projects

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

## Python Project (2020)

At long last, after many months of New Zealand’s isolation from the world, a
vaccine for the novel coronavirus has been developed! A first batch will soon be
arriving at Auckland Airport from the overseas facility that produced it, and this project
addresses the distribution of the vaccines to the most vulnerable population in
Auckland: the elderly living in residential care settings. A team of
four couriers available to do this distribution, and the job is to get the vaccine
delivered from the airport to all 143 rest homes in the city, as quickly as possible.
There is a catch, however: due to the severe economic impact of this period
of isolation, private travel for your team will not be possible. There’s just no
money in the budget for petrol, unfortunately. The project's aim is to figure 
out how to deliver the vaccines using the public transportation system.

This project uses a nearest-neighbour algorithm that returns four
ordered lists of network node names, representing the paths taken by the four
couriers. 

## Computational Mechanics: Python Resource Consent Project for the Onehunga Aquifer (2020)

The Onehunga Aquifer has been a significant source of water in the Auckland region, providing up to 20 million litres per day since the year 2000. However, due to the growing drought in the Auckland region, Watercare has proposed to increase the water extraction from the aquifer to satisfy increased demands. There are potential environmental impacts concerned with this proposal:
- Pressure decline in the aquifer: As the rate of water extraction increases, the pressure in the aquifer decreases as water flows across the aquifer under a pressure gradient.
- Contamination of water: Contaminants, such as copper, dissolve into stormwater at the surface and then leaches down to pressure lows in ground water. The contamination correlates to the pressure decline as a result of increased water extraction. Water extracted can therefore be unsafe to drink if the level of copper exceeds a certain threshold.

Watercare have applied to the Auckland Regional Council (ARC) for resource consent to double their maximum take from its current limit of 20 million litres per day. They claim that the aquifer can support the increased take without copper concentrations exceeding safe levels.

This project undertakes a computer modelling study on groundwater quality of the Onehunga Aquifer addressing both pressure and copper concentration changes. A coupled lumped parameter model of concentration and copper is formulated and solved numerically using improved Euler's method. The parameters are then calibrated automatically using Python's built-in curve_fit function. The model is then used to predict into the future, for four different possible outcomes of the consent hearing. This includes approval of the expanded consent (a maximum of 40 million litres per day), no change, reduced consent (to a recommended safe level) or, more drastically, an indefinite moratorium on aquifer usage. A recommendation is then given in the final report. 


