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
- Install MMDiff by running `install("/Path/To/MMDiff/)`.
### Data Generation.
- Switch directories to `Experiment_NB`.
- From within `Experiment_NB`, start up Jupyter by running `jupyter notebook`.
- In Jupyter, open the notebook *Validation Experiment* and execute all paragraphs..



