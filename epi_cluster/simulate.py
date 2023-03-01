import numpy as np
import tensorflow as tf
import tensorflow_probability as tfp
import argparse 
from pathlib import Path
from os.path import join
import json
import pandas as pd
import logging


class Params(object):
    pass


def draw_sample(tensor, index, mask):
    tensor_f = tf.reshape(tensor, [-1])
    index_f = tf.reshape(index, [-1])
    mask_f = tf.reshape(mask, [-1])
    bool_mask_f = tf.greater(mask_f, 0)
    
    index_masked = tf.boolean_mask(index_f, bool_mask_f)
    tensor_masked = tf.boolean_mask(tensor_f, bool_mask_f)

    return (tf.stack([tf.cast(tensor_masked,dtype=tf.int32), index_masked]))


class Simulation(object):
    def __init__(self, params):
        # possibly learn these parameters later with STAN
        self.name = params.name
        self.output_dir = params.output_dir
        self.n_patterns = params.n_patterns
        self.n_histones = params.n_histones
        self.n_modifications = params.n_modifications
        self.locations_per_pattern = params.locations_per_pattern
        self.sampling_depth_foregroud = params.sampling_depth_foreground
        self.sampling_depth_backgroud = params.sampling_depth_background
        self.cluster_prob = params.cluster_prob
        self.cluster_intensity_alpha = params.cluster_intensity_alpha
        self.cluster_intensity_beta = params.cluster_intensity_beta
        self.histone_width = params.histone_width
        self.histone_sigma = params.histone_sigma
        self.foreground_prior_beta = params.foreground_prior_beta
        self.foreground_prior_alpha = params.foreground_prior_alpha
        self.background_prior_beta = params.background_prior_beta
        self.background_prior_alpha = params.background_prior_alpha
        
        self.patterns = self.generate_patterns()
        self.session = tf.InteractiveSession()
        self.save_config()

    def run_simulation(self):
        logging.info('Started run_simulation()')
        # create cluster for each pattern, sample foreground reads
        for i in range(self.n_patterns):
            for j in range(len(self.patterns[i])):
                confounding_vector = tf.constant(self.patterns[i][j], dtype=tf.float32)
                confounding_vector *= tfp.distributions.Beta(self.cluster_intensity_alpha, self.cluster_intensity_beta).sample()

                cluster, _, confounder_distribution = self.create_cluster(confounding_vector)
                with tf.name_scope('summaries'):
                    tf.summary.tensor_summary("Association "+str(i)+" Cluster "+str(j), confounder_distribution.probs)
                if j > 0:
                    clusters = tf.concat([cluster, clusters],axis=1)
                else:
                    clusters = cluster
            
            # sample background reads
            draws_background = self.create_nucleosome_draws(type='background')
            driver_background, _ = self.create_random_drivers()
            draws_background = self.flatten_reshape(draws_background)
            driver_background = self.flatten_reshape(driver_background)
            index_background = self.create_mod_index()
            background_samples = draw_sample(draws_background, index_background, driver_background)

            # run and write to file
            for j in range(self.locations_per_pattern):
                logging.info('Simulating reads...')
                _ = tf.summary.merge_all()
                out = self.session.run(
                    {'background': tf.transpose(background_samples, [1,0]), 
                     'foreground': tf.transpose(clusters, [1,0])})
            
                background_df = pd.DataFrame({'position': out['background'][:,0], 'modification': out['background'][:,1]})
                foreground_df = pd.DataFrame({'position': out['foreground'][:,0], 'modification': out['foreground'][:,1]})

                df = pd.concat([background_df, foreground_df])
                self.write_reads_as_bed(df, f'peak_{j}')

    def write_reads_as_bed(self, df, name):
        logging.info('Writing reads to %s' % name)
        df.to_csv(join(self.output_dir, 'data', name), sep='\t', index=False)

    def create_random_drivers(self):
        logging.info('Creating random drivers...')
        nucleosome_probabilities = tfp.distributions.Beta(self.background_prior_alpha, self.background_prior_beta).sample([self.n_modifications, self.n_histones])
        nucleosome_distribution = tfp.distributions.Bernoulli(probs=nucleosome_probabilities)
        samples = nucleosome_distribution.sample(self.sampling_depth_backgroud)
        return (samples, nucleosome_distribution) ###[N x self.n_modifications x self.n_histones]

    def generate_patterns(self):
        logging.info('Generating enrichment patterns...')
        perm = np.random.permutation(self.n_modifications)
        size_sum = 0
        clusters = []
        while size_sum < self.n_modifications:
            clust_sz = np.random.binomial(self.n_modifications-size_sum-1, self.cluster_prob, 1) + 1
            clust_vec = np.pad(np.ones(clust_sz), (size_sum, self.n_modifications-size_sum-clust_sz), 'constant', constant_values=(0,0))
            clust_vec = clust_vec[perm]
            
            clusters.append(clust_vec)
            size_sum = size_sum+clust_sz[0]
        return(clusters)
    
    def save_config(self, indent=3):
        logging.info('Saving config/patterns to file...')
        Path(join(self.output_dir, 'data')).mkdir(parents=True, exist_ok=True)

        patterns_to_save = [list(p) for p in self.patterns]
        with open(join(self.output_dir, 'patterns.json'), 'w') as fp:
            json.dump(patterns_to_save, fp, indent=indent)
        with open(join(self.output_dir, 'config.txt'), 'w') as fp:
            import datetime
            fp.write(f'Time of start: {datetime.datetime.now()}\n\n')
            fp.write(json.dumps(vars(self), indent=indent))

    def create_cluster(self, confounding_vector):
        logging.info('Creating a cluster of modifications...')
        foreground_draws = self.create_nucleosome_draws(type='foreground')
        driver, confounder_distribution, nucleosome_distribution = self.create_confounding_drivers(confounding_vector)

        foreground_draws = self.flatten_reshape(foreground_draws)
        driver = self.flatten_reshape(driver)

        index = self.create_mod_index()
        samples = draw_sample(foreground_draws, index, driver)
        
        return(samples, nucleosome_distribution, confounder_distribution)
    
    def create_mod_index(self):
        logging.info('In create_mod_index()...')
        return tf.broadcast_to(
            tf.range(self.n_modifications), 
            [self.sampling_depth_foreground * self.n_histones, self.n_modifications])  ####[N*self.n_histones, self.n_modifications]
    
    def create_confounding_drivers(self, confounding_vector):
        logging.info('Creating confounding drivers...')
        nucleosome_probabilities = tfp.distributions.Beta(self.foreground_prior_alpha, self.foreground_prior_beta).sample([self.n_histones]) ###[self.n_histones]
        nucleosome_probabilities = tf.expand_dims(nucleosome_probabilities, 0) ###[1 x self.n_histones]
        nucleosome_probabilities = tf.tile(nucleosome_probabilities, multiples=[self.n_modifications, 1]) ###[self.n_modifications x self.n_histones]
        
        confounder = tf.expand_dims(confounding_vector, 1) ###[self.n_modifications x 1]
        confounder = tf.tile(confounder, multiples=[1, self.n_histones]) ###[self.n_modifications x self.n_histones]
        
        nucleosome_distribution = tfp.distributions.Bernoulli(probs=nucleosome_probabilities) ###[self.n_modifications x histnum]
        confounder_distribution = tfp.distributions.Bernoulli(probs=confounder) ###[self.n_modifications x histnum]

        nucleosome_draws = nucleosome_distribution.sample(self.sampling_depth_foreground) ###[self.sampling_depth_foreground x self.n_modifications x histnum]
        confounder_draws = confounder_distribution.sample(self.sampling_depth_foreground) ###[self.sampling_depth_foreground x self.n_modifications x histnum]
        
        draw = nucleosome_draws * confounder_draws
            
        return (draw, confounder_distribution, nucleosome_distribution)
    
    def flatten_reshape(self, tensor):
        tensor_t = tf.transpose(tensor, [0, 2, 1])
        tensor_r = tf.reshape(tensor_t, [-1, self.n_modifications])
        return(tensor_r)
    
    def create_nucleosome_draws(self, type):
        logging.info('Creating nucleosome draws for %s' % type)
        if type == 'foreground':
            N = self.sampling_depth_foregroud
        elif type == 'background':
            N = self.sampling_depth_backgroud
        else:
            raise ValueError(f'Unknown read type: {type}. Should be either foreground/background')
        
        offset = (self.histone_width / 2) - (2 * self.histone_sigma)
        means = tf.reshape(
            tf.tile(
                tf.linspace(
                    self.histone_width / 2, (self.n_histones - 1) * self.histone_width + self.histone_width / 2, self.n_histones), 
                [self.n_modifications]), 
            [self.n_modifications, self.n_histones])
        truncated_normal_draws = tf.truncated_normal([N, self.n_modifications, self.n_histones], mean=0.0, stddev=self.histone_sigma, dtype=tf.float32)
        location_draws = tf.add(truncated_normal_draws, means)
        location_draws = location_draws - offset
        return (location_draws) ###[N x self.n_modifications x self.n_histones]


def main(args):
    logging.info('Args to main(): %s' % str(args))
    args = json.load(open(args.params))
    params = Params()
    for key in args:
        setattr(params, key, args[key])
    
    simulation_object = Simulation(params)
    simulation_object.run_simulation()


if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    parser = argparse.ArgumentParser(prog = 'epi_cluster', description = 'Simulate Chip-Seq for correlated histone modifications')
    parser.add_argument('params', required=True)
    args = parser.parse_args()
    main(args)
