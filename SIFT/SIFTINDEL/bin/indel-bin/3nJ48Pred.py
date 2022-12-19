# This program use the if-else rules from decision tree to make
# predictions for 3n indels.
#
# Usage: python 3nJ48Pred.py featurefile predfile
# Jing Hu @ Franklin & Marshall College
# 2/4/2013
#
##################################################################


import sys

def predJ48(fracPFDomainsAffected, isRepeat, isDisorder, dnaConservScoreLeft):
    if fracPFDomainsAffected <= 0 :
        if isRepeat == 1 :
            if isDisorder == 1 :
                if dnaConservScoreLeft <= 0.301 :
                    predEffect = "neutral"
                    confidence = 0.696
                    ClassificationPath = "Classification_path_1"
                if dnaConservScoreLeft > 0.301 :
                    predEffect = "damaging"
                    confidence = 0.667
                    ClassificationPath = "Classification_path_2"
            if isDisorder == 0 :
                predEffect = "damaging"
                confidence = 0.826
                ClassificationPath = "Classification_path_3"
        if isRepeat == 0 :
            if isDisorder == 1 :
                predEffect = "neutral"
                confidence = 0.918
                ClassificationPath = "Classification_path_4"
            if isDisorder ==0 :
                if dnaConservScoreLeft <= 1.405 :
                    predEffect = "neutral"
                    confidence = 0.720
                    ClassificationPath = "Classification_path_5"
                if dnaConservScoreLeft > 1.405 :
                    predEffect = "damaging"
                    confidence = 0.765
                    ClassificationPath = "Classification_path_6"
    if fracPFDomainsAffected > 0 :
        if isDisorder == 1 :
            if isRepeat == 1 :
                predEffect = "damaging"
                confidence = 0.943
                ClassificationPath = "Classification_path_7"
            if isRepeat == 0 :
                if dnaConservScoreLeft <= 2.058 :
                    predEffect = "neutral"
                    confidence = 0.636
                    ClassificationPath = "Classification_path_8"
                if dnaConservScoreLeft > 2.058 :
                    predEffect = "damaging"
                    confidence = 0.846
                    ClassificationPath = "Classification_path_9"
        if isDisorder == 0 :
            predEffect = "damaging"
            confidence = 0.894
            ClassificationPath = "Classification_path_10"

    return (predEffect, confidence, ClassificationPath)


def main():
    fin = open(sys.argv[1], 'r')
    fout = open(sys.argv[2], 'w')

    count = {}

    line = fin.readline()
    for line in fin:
        line = line.split('\n')[0]
        info = line.split('\t')
        indel = info[0]
        ensgID = info[1]
        isRepeat = int(info[2])
        if info[3]=='?' or info[4]=='?' or info[5]=='?':
            continue
        dnaConservScoreLeft = float(info[3])
        fracPFDomainsAffected = float(info[4])
        avgDisorder = float(info[5])
        if avgDisorder > 0.5:
            isDisorder = 1
        else:
            isDisorder = 0
        pred = predJ48(fracPFDomainsAffected, isRepeat, isDisorder, dnaConservScoreLeft)
        
        fout.write(indel+'\t'+ ensgID+'\t'+ pred[0]+'\t'+str(pred[1])+'\t3N_'+pred[2]+'\n')
        

    fout.close()
    fin.close()

main()
sys.exit(0)
    
