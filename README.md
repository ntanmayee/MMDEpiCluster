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
### Run Setup & Generate Configuration Files
- Navigate to `MMDEpiCluster` directory
- Change permission of `setup.sh` to an executable file, i.e. `chmod 700 ./setup.sh`
- Run `./setup.sh`. This creates the working directory `Experiment_WD`, and generates configuration files for the experiment. These will be located in  `Experiment_WD/Exp_conf`
- In order to generate testing configuration files navigate to `Experiment_WD/Exp_conf_test`. From within this directory run `python3 ../../Simulation_Script/generate_config_test1.py` & `python3 ../../Simulation_Script/generate_config_test2.py` 
### Data Generation.
- Switch directories to `Experiment_WD`.
- From within `Experiment_WD`, start up Jupyter by running `jupyter notebook`.
- In Jupyter, open the notebook *Validation Experiment* and execute all paragraphs.
### Data Analysis
- In Jupyter, open the *Data Analysis* notebook
- In paragraph 3, set `sim_wd` variable to point to the simulation directory
- In paragraph 3, set `exp_wd` variable to point to the experiment directory to be analysed
- Execute all paragraphs



