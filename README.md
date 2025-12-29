# Diel vertical migration parameter estimation using acoustics backscatter data

This repository contains a generic recipe to estimate the upper and lower boundaries of sound scattering layer as a function of underwater PAR (photosythetically active radiation) to represent dial vertical migration of organisms (mesozooplankton in this example) to be coded in respective biogeochemical models.

## Data source

The repository includes the jupyter notebook, ML_DVM.with.ECOSMOv2025.ipynb,that:

1) Processes ADCP backscatter density at ATWAIN mooring station (82N 31E). The raw data has been retrieved from Norwegian Marine Data Center (https://metadata.nmdc.no/metadata-api/landingpage/df62e3e2931b185b23542f1081183282). Density anomaly is calculated for each depth level.
2) Loads OSISAF sea ice concentration corresponding to ATWAIN location.
3) Ice concentration is used to modify the shortwave radiation retrieved from ERA5 reanalysis.
4) Weekly averages for each hour (e.g. hour = 00:00 - 01:00, hour = 01:00 - 02:00, etc) are calculated both for the ADCP data, and atmospeheric forcing data.
5) For documentation purposes, creates the temperature and salinity profiles for ATWAIN location obtained from Copernicus Arctic Reanalysis Model ( https://doi.org/10.48670/moi-00007 ). These profiles were used to nudge the 1D GOTM-ECOSMO model setup for the ATWAIN station to perform a biogeochemistry simulation to produce underwater PAR data.
6) Loads variables for the 1D GOTM-ECOSMO simulation and performs a similar weekly averaging for the hourly model output.
7) Normalises (value range 0 - 1) ADCP backscatter anomaly to each individual profile, locates the center of mass (value = 1.0) and assigns the locations for the upper and lower boundaries (value = 0.8) above and below the center of mass respectively. Handles missing data points.
8) The forcing data, retrieved boundary locations and their corresponding PAR values are stored as a dataset for use in decision tree algorithms.
9) The code divides the processed dataset to four cases: 1. Polar night, 2. Polar day, 3. Day cycle, 4. Night cycle. Each case performs a decision tree approach to estimate a PAR value with respect to environmental conditions.
10) The final output of the decision tree are the input values (e.g. number of daylight hours, surface light, integrated food etc.) and their corresponding PAR values as outputs. These input and output values then can be coded as if/else cases in respective biogeochemical models.   

The provided Jupyter notebook is intended to be a generic code and includes an example workflow for the ATWAIN station in the Barents Sea. Results are time and location specific and each case will differ. Increasing sample size for the ADCP data will improve the codes generic value. The presented case here is a recipe to perform other experiments.

## Baseline, validation procedure and metrics

The intended use of this repository is to estimate underwater PAR (photosythetically active radiation) corresponding to the upper and lower boundaries of high concentration band of vertically migrating organisms. As such, the output of the machine learning algorithm is used as input (e.g. environemental driver to activate movement and locate the position) for diel vertical migration (DVM) of a migrating state variable in biogeochemical models. 

In this regard, a biogeochemical model simulation without DVM is considered in this work as the baseline. The inputs that activate DVM in biogeochemical models were intentionally calculated using the least possible decision tree depth to increase their genericness across model regions and usage of different biogeochemical model. Thus, the validation of these input values alone may not be the appropriate approach, but the validation of the actual biogeochemical simulations in a 3D configuration, and their comparison to simulations without DVM would yield better estimation of the success of this approach. Initial comparisons of different biogeochemical models with and without DVM module has been documented in this report: https://doi.org/10.5281/zenodo.14926411 (D5.2). This document also showcases simple prescribed DVM parameters, though the simple vs ML DVM has not been compared specifically. The comprehensive validation procedure are unique to each biogeochemical model and its physical configuration. The references to the validations for the respective models in D5.2 will be included here in future as they get puclically available.   

## List of dependencies

The Jupyter notebook requires the following python libraries:
- matplotlib
- numpy
- datetime
- pandas
- netCDF4
- scipy
- pickle
- os
- math

The modelinput.py file is included in this repository.
