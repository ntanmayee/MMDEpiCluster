import json
import os
import random

random.seed(1234567890)

for i in range(0,6):
    
    data = {'config':{
        'name':'Experiment_Med_Cov_'+str(i),
        'group':'dense_med',
        'random_seed':random.randint(1, 1234567890),
        'n_patterns':8,
        'n_histones':20,
        'n_modifications':10,
        'locations_per_pattern':300,
        'sampling_depth_foreground':160,
        'sampling_depth_background':40
    },
    'cluster':{
        'prob':0.4,
        'intensity_alpha':0.4,
        'intensity_beta':0.2
    },
    'histone':{
        'width':200,
        'sigma':75
    },
    'foreground_prior_beta':{
        'alpha': 0.7,
        'beta': 0.9
    },
    'background_prior_beta':{
        'alpha': 1.0,
        'beta': 1.0
    }}
    jstr = json.dumps(data, indent=4)
    with open(data['config']['name']+"_"+data['config']['group']+"_"+"conf"+".json", "w") as text_file:
        print(jstr, file=text_file)


