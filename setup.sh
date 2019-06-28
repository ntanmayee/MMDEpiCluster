#!/bin/sh

#  setup.sh
#  
#
#  Created by David Helekal on 28/06/2019.
#  
mkdir Experiment_WD
cp -a Experiment_NB/. Experiment_WD/

cd Experiment_WD
mkdir Exp_conf
mkdir Exp_conf_test

cd Exp_conf

for CONFIG_GEN in ../../Simulation_Script/Configs/*.py
do
    python3 $CONFIG_GEN
done
