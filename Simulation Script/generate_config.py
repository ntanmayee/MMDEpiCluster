import json
import os
import random

random.seed(1234567890)

for i in range(0,5):
    
    data = {'config':{
        'name':'Experiment_HiDense_'+str(i),
        'random_seed':random.randint(1, 1234567890),
        'n_patterns':8,
        'n_histones':10,
        'n_modifications':10,
        'locations_per_pattern':300,
        'sampling_depth_foreground':60,
        'sampling_depth_background':20
    },
    'cluster':{
        'prob':0.4
    },
    'histone':{
        'width':300,
        'sigma':150
    },
    'foreground_prior_beta':{
        'alpha': 0.5,
        'beta': 0.5
    },
    'background_prior_beta':{
        'alpha': 1,
        'beta': 1
    }}
    jstr = json.dumps(data, indent=4)
    with open("configHiDense"+str(i)+".json", "w") as text_file:
        print(jstr, file=text_file)


