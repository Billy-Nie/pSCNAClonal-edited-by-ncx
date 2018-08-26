# -*- coding: utf-8 -*-
'''
Created on 2018-08-01

@author


'''

import os
import sys
import pickle as pkl

import numpy as np
import scipy as sp

PLT_AVAIL = True
try:
    from matplotlib import pyplot as plt
except:
    PLT_AVAIL = False

from pSCNAClonal import constants
from pSCNAClonal.preprocess.data.pools import stripePool
from pSCNAClonal.model import model_base,joint_model



def run_postprocess(args):

    file_name = args.output_file_base + ".pSCNAClonal.output.pkl"
    infile = open(file_name,'rb')

    trainer = pkl.load(infile)
    stripePool = trainer.stripePool
    extract_postprocess_table(stripePool.stripes, stripePool.segPool.segments, True)
    infile.close()


def extract_postprocess_table(stripes, segments, simulation_flag):
    if not os.path.exists("postprocess"):
        os.makedirs("postprocess")
    outputFile = open("postprocess/postprocess_table.txt", "wr")
    if simulation_flag:
        outputFile.write(
            "chrom\tstart\tend\tnReadNum\ttReadNum\tcopy_number_predict\tphi_predict\tstripe id\tcopy_number_real\tphi_real\n")
        for i in range(len(segments)):
            segment = segments[i]
            outputFile.write("{0}\t{1}\t{2}\t{3}\t{4}\t".format(segment.chromName, segment.start, segment.end, segment.nReadNum, segment.tReadNum))
            baseline = True
            for j in range(len(stripes)):
                if i in stripes[j].segsIdxL:
                    outputFile.write("{0}\t{1}\t{2}\t".format(stripes[j].copyNumber, stripes[j].phi,stripes[j].sid))
                    baseline = False
            if baseline:
                outputFile.write("{0}\t{1}\tbaseline\t".format(2, 1))

            if segment.start >= 210682863 and segment.end <= 210782863:
                outputFile.write("3\t0.9\n")
            elif segment.start >= 152728665 and segment.end <= 152828665:
                outputFile.write("0\t0.8\n")
            elif segment.start >= 152829765 and segment.end <= 153029765:
                outputFile.write("2\t1\n")
    else:
        outputFile.write("chrom\tstart\tend\tnReadNum\ttReadNum\tcopy_number_predict\tphi_predict\n")
        for i in range(len(segments)):
            segment = segments[i]
            outputFile.write("{0}\t{1}\t{2}\t{3}\t{4}\t".format(segment.chromName, segment.start, segment.end, segment.nReadNum, segment.tReadNum))
            for j in range(len(stripes)):
                if i in stripes[j].segsIdxL:
                    outputFile.write("{0}\t{1}\t".format(stripes[j].copyNumber, stripes[j].phi))

    outputFile.close()


"""
def extract_paired_counts(stripePool, output_file_base):
    counts_file_name = output_file_base + ".pSCNAClonal.counts"
    outfile = open(counts_file_name,'w')
    stripes = stripePool.stripes

    print "Extracting paired counts file..."
    sys.stdout.flush()

    #TODO : 以下的内容需要和初师兄确认
    outfile.write('\t'.join(['#stripe_name','normal_A','normal_B','tumor_A','tumor_B','chrom','pos']) + '\n')

    for j in range(0,stripePool.stripeNum):
        for i in range(0,stripes[j].pairedCounts.shape[0]):
            outfile.write(stripes[j].name+'\t'+'\t'.join(map(str,stripes[j].pairedCounts[i])) + '\n')

    outfile.close()

def extract_stripes(stripePool, output_file_base):
    stripe_file_name = output_file_base + 'pSCNAClonal.stripes'
    outfile = open(stripe_file_name,'w')
    stripes = stripePool.stripes

    print "Extracting stripes file"
    sys.stdout.flush()

    #TODO:以下内容需要和初师兄确认,tag到底是干啥😂
    outfile.write('\t'.join(['#stripe_name','stripe_id','normal_read_num','tumor_read_num','tag','copy_number','geno_type','phi','subclone_cluster']))

    for j in range(0,stripePool.stripeNum):
        outfile.write('\t'.join(map(str,[stripes[j].name,stripes[j].sid,stripes[j].nReadNum,stripes[j].tReadNum,
                                         stripes[j].tag,stripes[j].copyNumber,stripes[j].genotype,stripes[j].phi,
                                         stripes[j].subclone_cluster])) + '\n')

    outfile.close()


def extract_BAFheatmap(stripePool,output_file_base):
    BAF_counts_min = constants.BAF_COUNTS_MIN
    BAF_counts_max = constants.BAF_COUNTS_MAX

    outheatmap_dir_name = output_file_base + 'pSCNAClonal.heatmap'
    if not os.path.exists(outheatmap_dir_name):
        os.mkdir(outheatmap_dir_name)

    seg_num = len(stripePool.segPool.segments)

    #TODO : stripe里面没有BAFCounts属性，从它里面的segPools里面拿了
    for j in range(0,seg_num):
        BAF_counts_j = stripePool.segPool.segments[j].BAFCounts
        seg_name_j = stripePool.segPool.segments[j].name
        BAF_counts_sub = BAF_counts_j[BAF_counts_min:BAF_counts_max, BAF_counts_min:BAF_counts_max]
        color_max_j = BAF_counts_sub.max()

        print "plotting segment {0}...".format(seg_name_j)
        sys.stdout.flush()

        plt.figure(figsize=(8, 8), dpi=150)
        plt.xlim((0, 100))
        plt.ylim((0, 100))
        plt.xticks(sp.linspace(0, 100, 11), sp.linspace(0, 1, 11))
        plt.yticks(sp.linspace(0, 100, 11), sp.linspace(0, 1, 11))
        plt.xlabel('Tumor sample B allele frequency')
        plt.ylabel('Normal sample B allele frequency')
        plt.imshow(BAF_counts_j, vmin=0, vmax=max(1, color_max_j))
        cbar = plt.colorbar(ticks=[0, color_max_j], orientation='vertical', shrink=0.78)
        cbar.ax.set_yticklabels(['0', '>= ' + str(int(color_max_j))])
        plt.savefig('./' + outheatmap_dir_name + '/' + seg_name_j, bbox_inches='tight')



def extract_stripe_plot(stripePool, output_file_base):
    subclone_plot_file_name = output_file_base + '.pSCNAClonal.stripePlot.png'

    subclone_prev_lst = []
    copynumber_lst = []
    stripe_num = stripePool.stripeNum

    print "Extracting stripe plot file"
    sys.stdout.flush()

    #TODO:这里这个baseline以及subclone_prev要怎么改？
    for j in range(0,stripe_num):
        if stripePool.segPool.segments[j].baselineLabel == True or stripePool.stripes[j].genotype == 'PM':
            continue
        subclone_prev_lst.append(stripePool.stripes[j].phi)
        copynumber_lst.append(stripePool.stripes[j].copyNumber)

    X = len(subclone_prev_lst)

    plt.figure(figsize=(8,8), dpi=150)
    plt.plot(range(1,X+1),subclone_prev_lst,'o')
    plt.xlim(0,X+1)
    plt.ylim(0,1)
    plt.xlabel('Copy number')
    plt.ylabel('Subclonal cellular prevalence')
    plt.xticks(sp.linspace(1,X,X),copynumber_lst)
    plt.yticks(sp.linspace(0,1,11),['0%', '10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%', '100%'])
    plt.savefig(subclone_plot_file_name,bbox_inches='tight')
"""






def extract_summary(model_parameters, config_parameters, ll, output_file_base):
    summary_file_name = output_file_base + '.pSCNAClonal.summary'
    outfile = open(summary_file_name, 'w')

    phi = model_parameters.parameters['phi']
    subclone_cluster = '\t'.join(map(str, range(1, config_parameters.subclone_num + 1)))
    subclone_prev = '\t'.join(map("{0:.3f}".format, phi.tolist()))

    print "Extracting summary file..."
    sys.stdout.flush()

    outfile.write("Model : joint\n")
    outfile.write("Maximum copy number : %s\n" % (config_parameters.max_copynumber))
    outfile.write("Subclone number : %s\n" % (config_parameters.subclone_num))
    outfile.write("Subclone cluster :             %s\n" % subclone_cluster)
    outfile.write("Subclonal cellular prevalence : %s\n" % subclone_prev)
    outfile.write("Optimum log-likelihood : %s\n" % (ll))

    outfile.close()


