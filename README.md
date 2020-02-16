# MATLAB-choreography

## Installation

Download or clone this repository to the folder on your computer where you want to run your pipeline. Open MATLAB and navigate to the newly downloaded folder. 

Set up the data and figure directories by running the function:
```matlab
initialise_directories()
```
## Set up configuration

Duplicate the 'default_config.yaml' file in the config folder and modify the choreography parameters to suit your rig set up. Edit the 'configfile' variable in 'main_choreography.m' and 'main_plots.m' to match the desired config file name.

## Adding data

Copy output data from Multi-Worm Tracker (MWT) into /data/mwt_data.

## Run Choreography

Run the command:
```matlab
main_choreography
```
## Generate Plots

Run the command:
```matlab
main_plots
```


