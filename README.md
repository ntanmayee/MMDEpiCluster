# MMDEpiCluster
## Description
- This repository contains files required to run and replicate the MMDiff association pattern clustering experiment.
- Instructions are provided further down this document.
## Instructions
### Prerequisites 
- Check out this repository.
- Preferred operating system: Linux / OSX.
- Functioning anaconda https://www.anaconda.com package manager.
### Environment Setup
- Import the conda environment by running `conda env create -f Conda-env/environment.yml` from this project folder.
- Activate the conda environment by running `source activate experiment_env`.
### Install MMDiff
- Check out the MMDiff3 package from github *[NEEDS LINK]*.
- Start R.
- From R, load the package `devtools` by running `library("devtools")`.
- Install MMDiff by running `install("/Path/To/MMDiff/")`.
### Generate Configuration Files
- Navigate to `Experiment_NB/Exp_conf`
- Run all configuration generation scripts located in `Simulation_Script/Configs` from within `Experiment_NB/Exp_conf`
### Data Generation.
- Switch directories to `Experiment_NB`.
- From within `Experiment_NB`, start up Jupyter by running `jupyter notebook`.
- In Jupyter, open the notebook *Validation Experiment* and execute all paragraphs.
### Data Analysis
- In Jupyter, open the `Data Analysis` notebook
- In paragraph 3, set `sim_wd` variable to point to the simulation directory
- In paragraph 3, set `exp_wd` variable to point to the experiment directory to be analysed
- Execute all paragraphs



