'''
# =============================================================================
#      FileName: data.py
#        Author: Chu Yanshuo
#         Email: chu@yanshuo.name
#      HomePage: http://yanshuo.name
#       Version: 0.0.1
#    LastChange: 2017-04-24 15:30:17
#       History: YI lI
# =============================================================================
'''

class Segment:

    def __init__(self):
        self.name = ""
        self.chromIdx = -1
        self.chromName = ""
        self.start = -1
        self.end = -1
        self.nReadNum = -1
        self.tReadNum = -1
        self.gc = -1

        self.LOHFrac = -1
        self.LOHStatus = 'NONE'
        self.APMFrac = -1
        self.APMStatus = 'NONE'

        self.pairedCounts = None
        # self.BAFCounts = None

        self.baselineLabel = 'FALSE'
        self.tag = '0'
        self.stripeIdx = -1
        self.stripeID = ''
        self.alleleType = 'NONE'

        self.copyNumber = -1

    def toName(self):
        return "chrom\t\
                start\t\
                end\t\
                gc\t\
                nReadNum\t\
                tReadNum\t\
                baselineLabel\t\
                tag\t\
                stripeID\t\
                copyNumber\
                "

    def toString(self):
        return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\
                    \t{8}\t{9}".format(
            str(self.chromName), str(self.start),
            str(self.end), str(self.gc),
            str(self.nReadNum), str(self.tReadNum), str(self.baselineLabel),
            str(self.tag), str(self.stripeID),
            str(self.copyNumber))