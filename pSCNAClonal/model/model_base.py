# -*- coding: utf-8 -*-
'''
Abstract classes which all models for classifying paired data should sub-class. 

Created on 2011-03-31

@author: Andrew Roth

JointSNVMix-0.6.2
joint_snv_mix.classification.base.EMModelTrainer
joint_snv_mix.classification.base.PriorParser

================================================================================

Modified on 2013-07-29

@author: Yi Li

pyloh.model.model_base

================================================================================

Modified on 2014-04-15

@author: Yi Li
'''
import sys
import pickle as pkl

from ConfigParser import ConfigParser
from abc import abstractmethod

import numpy as np

from pSCNAClonal.model.utils import *
from pSCNAClonal.preprocess.data.pools import stripePool

class ProbabilisticModel(object):
    def __init__(self, max_copynumber, subclone_num, baseline_thred):
        self.max_copynumber = max_copynumber
        self.subclone_num = subclone_num
        self.baseline_thred = baseline_thred
        self.priors_parser = PriorParser()
        self._init_components()

    @abstractmethod
    def read_priors(self, priors_filename):
        raise NotImplementedError

    #restore stripePool from pkl file
    def read_stripePool(self, pkl_path):
        infile = open(pkl_path,'rb')
        # self.data = pkl.load(infile)
        self.stripePool = pkl.load(infile) # 实际上是stripePool
        infile.close()


    @abstractmethod
    def preprocess(self):
        raise NotImplementedError
        
    def run(self, max_iters, stop_value):
        trainer = self.model_trainer_class(self.priors, self.stripePool, self.max_copynumber,
                                           self.subclone_num, max_iters, stop_value)
        
        trainer.train()
        trainer.predict()
        
        self.trainer = trainer
        
    def write_results(self, output_filename_base):
        results_file_name = output_filename_base + '.pSCNAClonal.output.pkl'
        outfile = open(results_file_name, 'wb')
        pkl.dump(self.trainer, outfile, protocol=2)
        
        outfile.close()

    @abstractmethod
    def _init_components(self):
        raise NotImplemented


#JointSNVMix
class ModelTrainer(object):
    def __init__(self, priors, stripePool, max_copynumber, subclone_num, max_iters, stop_value):
        self.priors = priors
        
        self.stripePool = stripePool
        
        self.max_copynumber = max_copynumber
        
        self.subclone_num = subclone_num
        
        self.max_iters = max_iters
        
        self.stop_value = stop_value
        
        self.iters = 0
        
        self.ll = 0
            
        self._init_components()

    @abstractmethod
    def train(self):
        raise NotImplemented

    @abstractmethod
    def _print_running_info(self, new_log_likelihood, old_log_likelihood, ll_change):
        raise NotImplemented  

    @abstractmethod
    def _init_components(self):
        raise NotImplemented

class ConfigParameters(object):
    def __init__(self, max_copynumber, subclone_num):
        self.max_copynumber = max_copynumber
        self.subclone_num = subclone_num
        
        self._init_components()

    @abstractmethod
    def _init_components(self):
        raise NotImplemented

class ModelParameters(object):
    def __init__(self, priors, stripePool, config_parameters):
        self.priors = priors
        self.stripePool = stripePool
        self.config_parameters = config_parameters
        
        self._init_parameters()

    @abstractmethod
    def _init_parameters(self):
        raise NotImplemented
    
    
class LatentVariables(object): 
    def __init__(self, stripePool, config_parameters):
        self.stripePool = stripePool
        self.config_parameters = config_parameters

    @abstractmethod
    def _init_components(self):
        raise NotImplemented


class ModelLikelihood(object):
    def __init__(self, priors, stripePool, config_parameters):
        self.priors = priors
        self.stripePool = stripePool
        self.config_parameters = config_parameters

    @abstractmethod
    def get_log_likelihood(self, parameters, priors):
        raise NotImplemented


#JointSNVMix
class PriorParser(object):
    def __init__(self):                
        self.priors = {}
        
    def read_priors(self, priors_filename, max_copynumber):
        self.parser = ConfigParser() # configuration file ?
        self.parser.read(priors_filename)

        copynumber = get_copynumber(max_copynumber)
        copynumber_num = get_copynumber_num(max_copynumber)
        
        self.priors['omega'] = np.zeros(copynumber_num) # generate a copynumrt_num * 1 arrays filled with 0

        for i, cn in enumerate(copynumber):
            self.priors['omega'][i] = self.parser.getfloat('omega', str(cn))


