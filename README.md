# covid-age

This repository contains the C++ code used to run a stochastic SEIR model for COVID-19. Specifically, the model is fitted to Chicago's data for the first year of the pandemic. The model is modified from Northwestern University Malaria Modelling team's [early work](https://github.com/numalariamodeling/covid-chicago) during the COVID-19 pandemic. The code in this repository is used in the following publication:
- Toh KB, Runge M, Richardson RAK, Hladish TJ and Gerardin J. Design of effective outpatient sentinel surveillance for COVID-19 decision-making: a modeling study. _Currently under review in BMC Infectious Disease_. Available as a preprint in [medRxiv](https://www.medrxiv.org/content/10.1101/2022.10.21.22281330v1.full).

## General structure

`src` folder contains the "backbone" of the event driven algorithm and the SEIR model.
`exp` folder contains the individual folder to each "experiment". `chicago_yr1` subfolder contains code and parameters to reproduce first year of pandemic in Chicago. `1pop` subfolder is used in the main experiment conducted for the sentinel surveillance paper; `2pop` subfolder is used to model two partially mixing subpopulations (0 to 39 years old and > 40 years old). `2pop` allows us to investigate the behaviour of surveillance indicator when two subpopulations experience different transmission intensity.
See `covid-age-description.txt` for detailed description of the parameters, links to the model diagram and descriptions, and how they are specified in the code.

## Usage

Use the `Makefile` (`make all`) in `src` to compile the model, then use `Makefile` in the individual experiment folder to compile the code specific to the experiment via `make model`. Then you can do `./model 1234` to run the model with a random number seed of 1234.
