# Berkeley-COVID-testing

This repository is the companion to "Optimizing COVID-19 control with asymptomatic testing in a university environment", for both reproduction of simulations and figures shown in the manuscript and a brief guide for using the model with a custom parameterization.

## Contents
The repository is split into three main directories
- Final-Figs
  - Contains R scripts to generate each figure in the manuscript as well as the .Rdata files which are called upon
- all-model-runs
  - Contains self contained subdirectories for each simulation performed which can be directly copied out to a server or cluster, with all dependency files included in each
- model-sandbox
  - A simple directory with containing a script for running the model and the basic parameter files. This can be the jumping off point for those looking to adapt the model to their community

## Code Structure
The code builds up to a single function `replicate.epidemic()` which runs a number of realizations of the model and returns a list of 4 objects.
- `mean.dat` : Named list of means and CI bounds relevant epidemic information
- `mean.long` : Named list of epidemic trajectories for each model realization
- `pop.mat.chain` : List containing an item for each model realization with information on individual times of infection, viral dynamics, etc.
- `mean.R.mat` : Named list containing information about Reff distributions

`replicate.epidemic` takes in a set of parameters with no default parameters
Parameter Name | Description
`n.reps` | Number of model realizations
`virus.par` | Data frame with information about Reff distribution and generation time distribution (loaded from `"virus.par.9.4.Rdata"` in model-sandbox)
`input.par` | Data frame containing parameters for interventions (loaded from `"pop.par.base.Rdata"` in model-sandbox)
`titer.dat` | Data frame containing individual viral titer trajectories, number of individuals must match total population size (loaded from `"titer.dat.20K.Rdata"` in model-sandbox)
`times` | Number of days to simulate
`bay.area.prev` | Prevalence (as a fraction) in the surrounding community
`bay.area.R` | Reff of imported cases from the surrounding community
`burnin` | Duration of burnin
`serial_mean` | Mean of serial interval distribution for disease
`serial_sd` | Standard deviation of serial interval for disease
`length_timestep` | length of timestep in simulation (in days)
`superspreader` | Boolean, `TRUE` will implement superspreading dynamics
`LOD` | Viral titer needed for positive test
`input.pop` | List of sizes of each subpopulation
`n.init.exposed.vector` | List of number of individuals initially exposed in each subpopulation
`employ.id.vector` | List of employee ids marking the start of a new subpopulation (i.e. if the population is a single group this is simply `c(1)`. If the subpopulations are 5000 apiece this is `c(1,5001)`)

## model-sandbox basics
As stated above this is a basic directory containing the key pieces to set up customized model runs, including examples of how to change NPI intervention parameters (e.g. the parameters stored in `"pop.par.base.Rdata"`). The current script `model-sandbox.R` is parameterized to run for a single population of size 650, with group size limits of 15 under 90% compliance, one day lags for testing turnaround, contact tracing, and symptomatic isolation, with biweekly testing. These are declared in the head of the script, as global variables along with several other useful parameters, which are then loaded into the parameter data frame if required before calling `replicate.epidemic()`. We do not set up every parameter sent to `replicate.epidemic()` in this way for brevity, but other parameters could be changed in the same manner if they are of interest.

For the parameterization given, `model-sanbox.R` can run locally in modest time. However do note that for larger populations or higher numbers of realizations it will obviously take longer, benchmarking a run with a single iteration will allow you to roughly approximate how much longer run 10x,50x,100x, etc. epidemics may take.

`model-sandbox.R` also contains some rudimentary plotting functions as an example of how to work with the `replicate.epidemic()` output, creating and saving a 4 panel plot with the name given to `OUTPUT_PLOT_NAME` alongside the `.Rdata` file containing the output structure. This plot shows mean daily and cumulative incidence and isolations for the epidemic. ![Example Plot](/model-sandbox/example.png)
