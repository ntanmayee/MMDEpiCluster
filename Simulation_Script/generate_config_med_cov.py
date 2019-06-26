import json
import os
import random

random.seed(1234567890)

for i in range(0,5):
    
    data = {'config':{
        'name':'Experiment_Med_Cov_'+str(i),
        'group':'med',
        'random_seed':random.randint(1, 1234567890),
        'n_patterns':4,
        'n_histones':10,
        'n_modifications':10,
        'locations_per_pattern':300,
        'sampling_depth_foreground':240,
        'sampling_depth_background':60
    },
    'cluster':{
        'prob':0.4,
        'intensity_alpha':0.8,
        'intensity_beta':0.4
    },
    'histone':{
        'width':300,
        'sigma':80
    },
    'foreground_prior_beta':{
        'alpha': 0.3,
        'beta': 0.2
    },
    'background_prior_beta':{
        'alpha': 1.0,
        'beta': 1.0
    }}
    jstr = json.dumps(data, indent=4)
    with open(data['config']['name']+"conf"+".json", "w") as text_file:
        print(jstr, file=text_file)


