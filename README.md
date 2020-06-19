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



## Additional packages (included)

- Colormaps
- append_pdfs
- yamlmatlab



## Changelog



### 2020-06-15

- Ridgelines now plot the full duration that objects are tracked for, rather than the first 30 seconds (time over which average metrics are calculated)
- Added 'curve' to the list of Choreography features output (config file)
- Added 'curve_smooth' parameter (moving mean, window = 21 frames)
- Added 'average curve' boxplot
- Changed the 'path' figures to be print all paths across a multi-page pdf, each page displaying 20 objects.
- Add 'additional_plots.m' - work in progress plots
  - Curve violins
  - Curve ridgeline - red highlights > 10 deg
  - Curve traces - grey lines: all animals, red line: mean trace, red shaded area: standard deviation
  - Area plots - fluctuation in average object area against time tracked.

