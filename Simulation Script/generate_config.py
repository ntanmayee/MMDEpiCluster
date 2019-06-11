import json
import os

data = {'config':{
    'name':'Experiment One',
    'random_seed':123456,
    'n_patterns':10,
    'n_histones':10,
    'n_modifications':10,
    'locations_per_pattern':10,
    'sampling_depth_foreground':60,
    'sampling_depth_background':40
},
'cluster':{
    'prob':0.5
},
'histone':{
    'width':300,
    'sigma':150
},
'foreground_prior_beta':{
    'alpha': 0.3,
    'beta': 0.3
},
'background_prior_beta':{
    'alpha': 1,
    'beta': 1
},
'mixing_prior_beta':{
    'alpha': 4,
    'beta': 0.6
}}

jstr = json.dumps(data, indent=4)
with open("config.json", "w") as text_file:
    print(jstr, file=text_file)