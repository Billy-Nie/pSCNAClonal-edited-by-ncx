#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: run.py
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2017-04-10 16:50:28
#       History: Author : Andrew Roth
# =============================================================================
'''
import argparse
from pSCNAClonal.preprocess.run_preprocess import process as run_process
from pSCNAClonal.model.run_model import run_model
from pSCNAClonal.postprocess.run_postprocess import run_postprocess

parser = argparse.ArgumentParser(description="Run pSCNAClonal to ",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter) # TODO : describe the perpose of this software
parser.add_argument('--version', action='version', version='PSSP-0.0.0')
subparsers = parser.add_subparsers()

###################
#   preprocess    #
###################


parserPreprocess = subparsers.add_parser("preprocess",help = "out put preprocess format")

parserPreprocess.add_argument("--nBamName",help="BAM file for normal sample")

parserPreprocess.add_argument("--tBamName",help="BAM file for tumor sample")

parserPreprocess.add_argument("--bedName",help="BED file for segments sample")

parserPreprocess.add_argument("--refFaName",help="FASTA file for reference genome")

parserPreprocess.add_argument("--pathPrefix",default="",help="Base name of the preprocessed input file to be created")

parserPreprocess.add_argument("--subcloneNumber",default=2,type=int,help="Set the subclone number")

parserPreprocess.add_argument("--coverage",type=int,default=30,help="Set the coverage number")

parserPreprocess.add_argument("--maxCopyNumber",type=int,default=6,help="Set the maximum copy number")

parserPreprocess.add_argument("--baselineThredLOH",default=0.3,type=float,help="baseline Thred of LOH")

parserPreprocess.add_argument("--baselineThredAPM",default=0.01,type=float,help="baseline Thred of APM")

parserPreprocess.add_argument("--minDepth",default=20,type=int,help="Minimum reads depth required for both normal and tumor samples. Default is 20")

parserPreprocess.add_argument("--minBqual",default=10,type=int,help="Minimum base quality required. Default is 10")

parserPreprocess.add_argument("--minMqual",default=10,type=int,help="Minimum mapping quality required.Default is 10")

parserPreprocess.add_argument("--processNum",default=1,type=int,help="Number of processes to launch for preprocessing. Default is 1.")

parserPreprocess.add_argument("--bedCorrectedPath",default="",help="The name of corrected BICseq result file")

parserPreprocess.add_argument("--pklPath",default="",help="Load the pkl path")

parserPreprocess.add_argument("--gcCorrectionMethod",default="auto",help="The gc correction method, one of auto and visual")

parserPreprocess.add_argument("--pklFlag",default=False,type=bool,help="The pkl flag")

parserPreprocess.set_defaults(func = run_process)



########################
#   model process      #
########################
parserModel = subparsers.add_parser('model',help="Output model parameters format")

parserModel.add_argument('--pklPath',help="Path to the stripePool.pkl generated by preprocess")
parserModel.add_argument('--max_copynumber',default=6,type = int,help = "Maximum copy number of each stripe allows to take. Default is 6.")
parserModel.add_argument('--subclone_num',default=0,type=int,help="Number of subclones within the tumor sample. If not provided, \
                                go through [1, 5] and select the best model. Default is None.")
parserModel.add_argument('--baseline_thred',default=0.16,type=float,help="The threshold of LOH SNP sites fraction within each segment to \
                            define the segment as baseline. Default is 0.16.")
parserModel.add_argument('--max_iters',default=30,type=int,help="Maximum number of iterations for training. Default is 30.")
parserModel.add_argument('--stop_value',default=1e-6,type=float,help="Stop value of EM iterations")
parserModel.add_argument('--output_filename_base',default="model_result",help="the name base of the output file")
#parser_run_model.add_argument('--priors', default=None, type=str,
#                             help='''File of the prior distribution. If not provided,
#                                use uniform prior. Default is None.''')

parserModel.set_defaults(func=run_model)

###########################
#    post process         #
###########################
parser_postprocess = subparsers.add_parser('postprocess',help='Extract various result files from the outputfile')

parser_postprocess.add_argument('--output_file_base', default='postprocess_output', help='Base name of the output file to be created')

parser_postprocess.set_defaults(func=run_postprocess)

############################
#          run             #
############################

args = parser.parse_args()
args.func(args)


