#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
# =============================================================================
#      FileName: converter.py
#          Desc:
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2018-01-31 14:07:22
#       History:
# =============================================================================
'''
import os
import pickle as pkl
import sys
from multiprocessing import Pool

import numpy as np
import pysam

import pSCNAClonal.constants as constants
import matplotlib.pyplot as plt
from pSCNAClonal.preprocess.data.pools.segmentPool import SegmentPool
from pSCNAClonal.preprocess.data.pools.stripePool import StripePool
from pSCNAClonal.preprocess.iofun import (PairedCountsIterator,
                                            PairedPileupIterator)
from pSCNAClonal.preprocess.mcmc import MCMCLM
from pSCNAClonal.preprocess.plotGC import GCStripePlot, GCStripePoolPlot
from pSCNAClonal.preprocess.utils import (get_BAF_counts, AnswerIndex,
                                            dump_seg_to_txt,
                                            normal_heterozygous_filter)



np.set_printoptions(threshold=np.inf)




class BamConverter:

    def __init__(self,
                 nBamName,
                 tBamName,
                 bedName,
                 refFaName,
                 pathPrefix="",
                 subcloneNumber=2,
                 coverage=30,
                 maxCopyNumber=6,
                 baselineThredLOH=0.3,
                 baselineThredAPM=0.01,
                 minDepth=20,
                 minBqual=10,
                 minMqual=10,
                 processNum=1,
                 bedCorrectedPath="",
                 pklPath=""):
        self._nBamName = nBamName
        self._tBamName = tBamName
        self._bedName = bedName
        self._refFaName = refFaName

        self.__pathPrefix = pathPrefix

        self.__subcloneNumber = subcloneNumber
        self.__coverage = coverage
        self.__maxCopyNumber = maxCopyNumber
        self.__baselineThredLOH = baselineThredLOH
        self.__baselineThredAPM = baselineThredAPM

        self.__minDepth = minDepth
        self.__minBqual = minBqual
        self.__minMqual = minMqual

        self.__processNum = processNum
        self.__bedCorrectedPath=bedCorrectedPath
        self.__pklPath = pklPath

        self._segPool = None

    def convert(self, readFromBed=True, method="auto", pklFlag=False):
        if not pklFlag:
            self._load_segs(readFromBed)
            # self._correct_bias(method)
            self._load_allele_counts()
            self._dump(self._segPool, ".temp.segPool")

            pklFile = open(self.__pklPath, 'wb')
            pkl.dump(self._segPool, pklFile, protocol=2)
            pklFile.close()
        else:
            pklFile = open(self.__pklPath, 'rb')
            self._segPool = pkl.load(pklFile )
            pklFile .close()


        blSegsL = self._get_baseline()

        stripePool = self._generate_stripe()


        stripePool.stripeNum = len(stripePool.stripes)

        self._dump(stripePool, "stripePool.pkl")
        self._dump_txt(stripePool, "stripePool.txt")

        self._dump(self._segPool, "segPool.pkl")
        self._dump_txt(self._segPool, "segPool.txt")

        if self.__answerFilePath != "":
            self.__dump_seg_to_txt()
        # self.__generate_ROC_table(None,stripePool,self._segPool.segments)

    def __dump_seg_to_txt(self):
        """
        output table for R, draw figures
        """
        dump_seg_to_txt(self._segPool, self.__pathPrefix)


    def __generate_ROC_table(self,blSegs,stripePool,segPool):

        #根据每一个segment的开始和结束的位置来确定变异的类型
        #具体变异类型由answers/rd30stage0.seg.txt中第二列的括号里的东西给定
        def decide_veriation_id(segment):
            if segment.start >= 210682863 and segment.end <= 210782863:
                return 2
            elif segment.start >= 152728665 and segment.end <= 152828665:
                return 1
            elif segment.start >= 209873187 and segment.end <= 209973187:
                return 3
            elif segment.start >= 152829765 and segment.end <= 153029765:
                return "base_line"

        #统计每一个stripe里面的segment的变异类型，选择变异类型最多的那一个作为这个stripe的变异类型
        def decide_stripe_variation_id(stripe,segPool):
            variation_1 = 0
            variation_2 = 0
            variation_3 = 0
            variation_baseline = 0
            for segsidx in stripe.segsIdxL:
                seg_variation_id = decide_veriation_id(segPool[segsidx])
                if seg_variation_id == 2:
                    variation_2 += 1
                elif seg_variation_id == 1:
                    variation_1 += 1
                elif seg_variation_id == 3:
                    variation_3 += 1
                elif seg_variation_id == "base_line":
                    variation_baseline += 1
            max_number = max(variation_3,variation_baseline,variation_1,variation_2)
            if max_number == variation_2:
                return 2
            elif max_number == variation_1:
                return 1
            elif max_number == variation_3:
                return 3
            elif max_number == variation_baseline:
                return "base_line"

        #计算某一个stripe的TP,FP,FN,TN(True Positive,False Positive)
        def calculate_TP_FP_FN_TN(stripeSegs,segPool,variation_id):
            TP = 0
            FP = 0
            FN = 0
            TN = 0
            for seg in segPool:
                true_variation_id = decide_veriation_id(seg)
                if seg in stripeSegs:
                    if variation_id == true_variation_id:
                        TP += 1
                    else:
                        FP += 1
                else:
                    if variation_id != true_variation_id:
                        TN += 1
                    else:
                        FN += 1
            return TP,FP,FN,TN

        def get_stripe_segs(stripe,segPool):
            stripe_seg = []
            for segID in stripe.segsIdxL:
                stripe_seg.append(segPool[segID])

            return stripe_seg
        def write_plot_data_to_file(sid_l,PPV_l,type_l):
            if not os.path.exists("plot_data/PPV_data"):
                os.makedirs("plot_data/PPV_data")
            outfile = open("plot_data/PPV_data/ppv_data.txt",'wr')
            outfile2 = open("plot_data/PPV_data/ppv_type.txt", "wr")
            outfile2.write("#stripe id\ttype\n")
            outfile.write("#stripe id\tPPV\n")
            for i in range(0,len(sid_l)):
                outfile.write("{0}\t{1}\n".format(sid_l[i],PPV_l[i]))
                outfile2.write("{0}\t{1}\n".format(sid_l[i],type_l[i]))
            outfile.close()
            outfile2.close()


        type_l = []
        stripe_id_l = []
        PPV_l = [] # positive predictive value

        for stripe in stripePool.stripes:
            stripe_variation_id = decide_stripe_variation_id(stripe,segPool)
            TP,FP,FN,TN = calculate_TP_FP_FN_TN(get_stripe_segs(stripe,segPool),segPool,stripe_variation_id)
            print stripe.sid
            print stripe_variation_id
            PPV = float(TP)/(float(TP) + float(FP))
            print "PPV: {:.2%} ".format(PPV)
            PPV_l.append(PPV)
            stripe_id_l.append(stripe.sid)
            type_l.append(stripe_variation_id)

        write_plot_data_to_file(stripe_id_l,PPV_l,type_l)




    def _dump_txt(self, Pool, outFilePath):
        """
        out put Pool into plain txt
        """
        Pool.output_txt(self.__pathPrefix + "/" +  outFilePath)

    def _load_allele_counts(self):
        self._get_counts(self._tBamName, self._segPool)

    def _dump(self, data, dpFileName):
        fileName = self.__pathPrefix + dpFileName
        outFile = open(fileName, 'wb')
        pkl.dump(data, outFile, protocol=2)
        outFile.close()

    def _generate_stripe(self):
        """
        generate stripe from segs
        """
        # 此处应该对每一个tag进行单独的聚类
        # 也有可能无法确定有多少个条带
        #
        # 另一种方案应该是对所有的条带进行聚类，然后从中挑选出目标条带，分解为新
        # 条带。

        yDown = constants.YDOWN
        yUp = constants.YUP
        # 此处应该近似为最大拷贝数与亚克隆数的乘积，作为手工输入也可以
        #stripeNum = constants.STRIPENUML[self._segPool[-1].idx]
        #noiseStripeNum = constants.NOISESTRIPENUML[self._segPool[-1].idx]
        stripeNum = constants.STRIPENUM
        noiseStripeNum = constants.NOISESTRIPENUM

        tempSP = StripePool(self._segPool, self._segPool.baseline, yDown, yUp,
                            stripeNum, noiseStripeNum)
        tempSP.get(byTag = False, plot= True)

        # 如上一步中添加了baseline选项
        # 那么这里的排序需要先排除baseline选型
        gspp = GCStripePoolPlot(tempSP)
        gspp.plot()

        tempSP.stripes.sort(key = lambda item: int(item.tag))

        for idx, sp in enumerate(tempSP.stripes):
            tempSP.stripes[idx].id = idx

        return tempSP


    def _load_segs(self, readFromBed=True):
        """
        load segments for each tumor sample
        """

        # for tBamName, bedName, coverage, subcloneNumber, i in zip(self._tBamNameL,
            # self._bedNameL, self.__coverageL, self.__subcloneNumberL, range(len(self._tBamNameL))):
            # print >> sys.stdout, 'Loading segments from bam file:\n{0}\n'.format(tBamName)
            # print >> sys.stdout, 'and bed file with gc:\n{0}\n'.format(bedName)
        self._segPool = SegmentPool(self.__maxCopyNumber,self.__coverage)
        if not readFromBed:
            print "****************************"
            print "Load read count from BAM file"
            print "***************************"
            nBam = pysam.Samfile(self._nBamName, 'rb')
            tBam = pysam.Samfile(self._tBamName, 'rb')
            self._segPool.load_seg_bam(nBam, tBam, self._bedName)
            nBam.close()
            tBam.close()
        else:
            self._segPool.load_seg_bed(self._bedName)

    def _correct_bias(self, method="auto"):
        """
        correct bias of each tumor sample
        """
        if "auto" == method:
            self._MCMC_GC_C(self._segPool, self.__subcloneNumber)
        elif "visual" == method:
            self._V_GC_C(self._segPool, len(self._segPool.segments))

    def _get_baseline(self):
        """
        get the baseline segments
        calculate baseline of each SegmentPool
        return: the baseline segments list of each SegmentPool
        """

        #tempBL = segPool.get_baseline(self.__maxCopyNumber,
        #                             self.__subcloneNumber,
        #                              self.__baselineThredLOH,
        #                              self.__baselineThredAPM,
        #                              isPreprocess=True)

        tempBL = self._segPool.get_baseline(self.__maxCopyNumber,
                                            self.__subcloneNumber,
                                            self.__baselineThredLOH,
                                            self.__baselineThredAPM,
                                            isPreprocess=True)

        return tempBL

    def _MCMC_GC_C(self, segPool, subcloneNumber):
        """
        The interception is irrelevant for correction, set as median
        MCMCLM only returns the m and c, then correct the segPool here
        """

        mcmclm = MCMCLM(segPool, 0, subcloneNumber, self.__maxCopyNumber)
        m, c = mcmclm.run()
        print "MCMC slope = {}".format(m)

        x = np.array(map(lambda seg: seg.gc, segPool.segments))
        y = np.array(map(lambda seg: np.log(seg.tReadNum + 1) -
                         np.log(seg.nReadNum + 1), segPool.segments))

        yCorrected = self._correct(x, y, m, c)

        for i in range(len(yCorrected)):
            segPool.segments[i].tReadNum = np.exp(
                yCorrected[i] +
                np.log(segPool.segments[i].nReadNum + 1)
            ) - 1
            segPool.segments[i].log_ratio = np.log(
                (yCorrected[i] + 1.0) /
                (segPool.segments[i].nReadNum + 1.0)
            )

        print "gc corrected, with slope = {0}, intercept = {1}".\
            format(m, c)

    def _correct(self, x, y, slope, intercept):
        K = np.percentile(y, 50)
        A = slope * x + intercept
        return y - A + K

    def visualize(self, segPool):
        gsp = GCStripePlot(segPool.segments, len(segPool.segments))
        print "total number: {}".format(segPool.segNum)
        gsp.plot()
        x, y, m, c = gsp.output()
        print "x, y, m, c"
        print x, y, m, c

    def _V_GC_C(self, segPool, sampleNumber=10000):
        gsp = GCStripePlot(segPool.segments, sampleNumber)
        print >> sys.stdout, "total number: {}".format(len(segPool.segments))
        gsp.plot()
        print >> sys.stdout, "x, y, m, c"
        print >> sys.stdout, gsp.output()

        x = np.array(map(lambda seg: seg.gc, segPool.segments))
        y = np.array(map(lambda seg: np.log(seg.tReadNum + 1) -
                         np.log(seg.nReadNum + 1), segPool.segments))
        yCorrected = self._correct(x, y, m, c)

        for i in range(len(yCorrected)):
            segPool.segments[i].tReadNum = np.exp( yCorrected[i] +
                np.log(segPool.segments[i].nReadNum + 1)) - 1
            segPool.segments[i].log_ratio = np.log(
                (yCorrected[i] + 1.0) /
                (segPool.segments[i].nReadNum + 1.0)
            )

        print "gc corrected, with slope = {0}, intercept = {1}".\
            format(slope, intercept)

    def _get_counts(self, tBamName, segPool):
        """
        get allele counts of target bam file
        save the counts into segPool
        """

        segNum = len(segPool.segments)
        processNum = self.__processNum
        print "processNum = {}".format(processNum)

        if processNum > segNum:
            processNum = segNum

        pool = Pool(processes=processNum)

        argsL = []

        for j in range(0, segNum):
            segName = segPool.segments[j].name
            chromName = segPool.segments[j].chromName
            chromIdx = segPool.segments[j].chromIdx
            start = segPool.segments[j].start
            end = segPool.segments[j].end

            argsT = (
                segName,
                chromName,
                chromIdx,
                start,
                end,
                self._nBamName,
                tBamName,
                self._refFaName,
                self.__minDepth,
                self.__minBqual,
                self.__minMqual)

            argsL.append(argsT)

        countsTL = pool.map(process_by_segment, argsL)

        for j in range(0, segNum):
            pairedCountsJ, BAFCountsJ = countsTL[j]
            segPool.segments[j].pairedCounts = pairedCountsJ
            # segPool.segments[j].BAFCounts = BAFCountsJ

# ===============================================================================
#  Function
# ===============================================================================


def process_by_segment(argsT):
    segName, chromName, chromIdx, start, end, nBamName,\
        tBamName, refFaName, minDepth, minBqual,\
        minMqual = argsT

    print 'Preprocessing segment {0}...'.format(segName)
    sys.stdout.flush()

    nBam = pysam.Samfile(nBamName, 'rb')
    tBam = pysam.Samfile(tBamName, 'rb')
    refFasta = pysam.Fastafile(refFaName)

    normalPileupIter = nBam.pileup(chromName, start, end)
    tumorPileupIter = tBam.pileup(chromName, start, end)

    pairedPileupIter = PairedPileupIterator(
        normalPileupIter, tumorPileupIter, start, end)
    pairedCountsIter = PairedCountsIterator(
        pairedPileupIter,
        refFasta,
        chromName,
        chromIdx,
        minDepth,
        minBqual,
        minMqual)

    pairedCountsJ, BAFCountsJ = iterator_to_counts(pairedCountsIter)
    countsTuple_j = (pairedCountsJ, BAFCountsJ)

    nBam.close()
    tBam.close()
    refFasta.close()

    return countsTuple_j


def iterator_to_counts(pairedCountsIter):
    buff = 100000

    pairedCountsJ = np.array([[], [], [], [], [], []], dtype=int).transpose()
    BAFCountsJ = np.zeros((100, 100))
    buffCounts = []
    i = 0

    for counts in pairedCountsIter:
        buffCounts.append(counts)
        i = i + 1

        if i < buff:
            continue

        buffCounts = np.array(buffCounts)

        if buffCounts.shape[0] != 0:
            BAFCountsBuff = get_BAF_counts(buffCounts)
            BAFCountsJ += BAFCountsBuff

        buffCountsFiltered = normal_heterozygous_filter(buffCounts)

        if buffCountsFiltered.shape[0] != 0:
            pairedCountsJ = np.vstack((pairedCountsJ, buffCountsFiltered))

        buffCounts = []
        i = 0

    buffCounts = np.array(buffCounts)

    if buffCounts.shape[0] != 0:
        BAFCountsBuff = get_BAF_counts(buffCounts)
        BAFCountsJ += BAFCountsBuff

    buffCountsFiltered = normal_heterozygous_filter(buffCounts)

    if buffCountsFiltered.shape[0] != 0:
        pairedCountsJ = np.vstack((pairedCountsJ, buffCountsFiltered))

    return (pairedCountsJ, BAFCountsJ)
