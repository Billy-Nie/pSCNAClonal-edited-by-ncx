#!/bin/sh
#=============================================================================
#     FileName: run_origin.sh
#         Desc: pipeline
#       Author: Chu Yanshuo
#        Email: chu@yanshuo.name
#     HomePage: http://yanshuo.name
#      Version: 0.0.1
#   LastChange: 2016-12-05 21:20:06
#      History:
#=============================================================================

# software path
#pSCNAClonal_PATH=/data/yschu/projects/p-SCNAClonal/p-scnaclonal/
pSCNAClonal_PATH=/home/cxnie/p-scnaclonal-pscnaclonal/


#result_folder=/data/yschu/projects/p-SCNAClonal/pipeline/
result_folder=/home/cxnie/p-scnaclonal-pscnaclonal/pipeline/
normal_bam=/data/yschu/projects/phy-SCNAClonal/subsimtree/rd30/normal_30.bam
tumor_bam=/data/yschu/projects/phy-SCNAClonal/subsimtree/rd30/tumor_30.stage0.bam
segments_bed=/data/yschu/projects/phy-SCNAClonal/subsimtree/rd30/sh_result/rd30stage0.seg.splited.gc.txt #bedName

reference_genome=/data/yschu/data/fasta/reference/chr1.fa
pathPreFix=$result_folder

BICseq_bed_corrected=$segments_bed.gccorrected
pkl_path=$result_folder.temp.pkl
echo "pkl path"
echo $pkl_path

#default 20
min_depth=1
#default
min_base_qual=10
#default
min_map_qual=10
#default = 1
process_num=$1
echo "process_num"
echo $process_num

#default = 6
max_copynumber=8
# parser_run_model.add_argument('--subclone_num', default=0, type=int, help='''Number of subclones within the tumor sample. If not provided, go through [1, 5] and select the best model. Default is None.''')
subclone_num=3
#parser_run_model.add_argument('--baseline_thred', default=0.09, type=float, help='''The threshold of LOH SNP sites fraction within each segment to define the segment as baseline. Default is 0.09.''')
baseline_thred_LOH=0.16
#parser_run_model.add_argument('--priors', default=None, type=str, #                             help='''File of the prior distribution. If not provided, #                                use uniform prior. Default is None.''')
#priors=None
#parser_run_model.add_argument('--max_iters', default=30, type=int, help='''Maximum number of iterations for training. Default is 30.''')
baseline_thred_APM=0.6

max_iters=30
#parser_run_model.add_argument('--stop_value', default=1e-6, type=float, help='''Stop value of the EM algorithm for training. If the change of log-likelihood is lower than this value, stop training. Default is 1e-6.''')
gcCorrectionMethod="auto"
coverage=30

##test operation##



####
###run###
python2.7 $pSCNAClonal_PATH/run.py preprocess \
     $normal_bam \
     $tumor_bam \
     $segments_bed \
     $reference_genome \
     --pathPrefix $pathPreFix \
     --subcloneNum $subclone_num \
     --coverage $coverage \
     --maxCopyNumber $max_copynumber \
     --baselineThredLOH $baseline_thred_LOH \
     --baselineThredAPM $baseline_thred_APM \
     --minDepth $min_depth  \
     --minBqual $min_base_qual \
     --minMqual $min_map_qual \
     --processNum $process_num \
     --bedCorrectedPath $BICseq_bed_corrected \
     --pklPath $pkl_path\
    --gcCorrectionMethod $gcCorrectionMethod \
