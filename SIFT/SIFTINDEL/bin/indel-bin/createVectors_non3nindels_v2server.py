'''
Copyright 2011 Jing Hu. All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, 
	this list of
      conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, 
      this list of conditions and the following disclaimer in the documentation 
      and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY FRANKLIN & MARSHALL COLLEGE ''AS IS'' 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL COPYRIGHT HOLDER OR CONTRIBUTORS 
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are 
those of the authors and should not be interpreted as representing
official policies, either expressed or implied, of Franklin & Marshall 
College.
'''
###############################################################################
#
# Updated : 5/10/2011
#
# This program extract feature vectors for each indel at gene lfevel.
# 1) AvgNumExons: average number of exons of the gene
#
# 2) num_of_tx_unaffected = absolute number of transcripts unaffected.
# Expect if this is high, then more likely to be neutral.
#
# 3) fract_of_tx_unaffected =  ( num_of_tx_unaffected / total number of transcripts).
# Ranges from 0 to 1.  Expect if this is high, more likely to be neutral.
#
# 4) avg_relative_indel_location, max_relative_indel_location, min_relative_indel_location :
# From the 3rd line of output for each transcript, calculate location of indel as a percentage
# (position of indel/total number of residues). You will have a value for each transcript.
# To get a single value for a gene, calculate the average of all values over all affected
# transcripts.  Please also calculate the max and min over all affected transcripts. I
# expect this value to capture indels at the end of transcripts (i.e. indels at the end
# will have values close to 1).  Ranges from 0 to 1. Expect if this is close to 1, then more
# likely to be neutral (indel is at end of protein).  Probably only one of the three values
## calculated (avg, max,min) will end up being used by the classifier.
#
# 5) Max_relative_indel_location  (*******chosen)
#
# 6) min_relative_indel_location
#
# 7) centered_avg_relative_indel_location:  Modification of avg_relative_indel_location
# to capture indels at the beginning and ends of transcripts.  Define
# centered_avg_relative_indel_location =  abs (0.5- abs_relative_indel_location).
# What does this do? I hope this will capture indels at the _beginning_ and end of transcripts.
# If an indel is at the end of the transcript , this new value will be close to  0.5
# (e.g. avg_relative_indel_location ~ 1, and centered_avg_relative_indel_location = 0.5 = abs (0.5 - 1)).
# If an indel is at the beginning of the transcript, centered_avg_relative_indel_location will also
# be close to 0.5 (abs (0.5-0)) . Ranges from 0 to 0.5, where if value is high,
# I expect the indel to be neutral.
#
# 8) avg_overlap_of_protein_1, Is overlap very different from position of indel? I don't think
# it is. If it is, then for each transcript, calculate fraction of overlap (overlap/total protein length),
# and then calculate the average over all transcripts.
#
# 9) avg_overlap_of_protein_2
#
# 10) avg_overlap_of_protein_mostDomains
#
# 11) fraction_with_mRNA_Decay	 : (# of transcripts that have NMD = yes / total number of transcripts).
# If this there is no NMD and/or the number of unaffected transcripts is high, then fraction_with_NMD
# will be low. This ranges from 0 to 1. If the value is low, I expect the indel
# to be neutral (less transcripts degraded).
#
# 12) fract_AllDomain_affected: ranges from 0 to 1. If the value is low, expect the indel to be neutral.
#
# 13) fract_PfamDomain_affected
#
# 14) fract_SuperfamDomain_affected
#
# 15) fract_SignalPeptideDomain_affected
#
# 16) Fract_conserved_affected : ranges from 0 to 1. The lower the fraction, the fewer conserved
# regions have been affected. If the value is low, expect the indel to be neutral.
#  (*******chosen)
#
# 17)  min_disExonBoundary: minimum distance of indel to the exon boundary  (*******chosen)
#  
# 18) num_Paralog
#
# 19) ka/ks: (added in this version)
#     Total entries: 701650, 642254 of which has dN entires, 9298 of which has dS = 0
#
# 20) fract_conservedAminoAcid_affected: percentage of conserved amino acid positions lost due to indel
#   (*******chosen)
# Usage: python createVectors_2011_4_3.py inputFile chrNum outputFile
#
# Copyright: Jing Hu & Pauline Ng
# Franklin & Marshall College, Singapore Genome Institute
# Date: April 3, 2011
#
##############################################################################

import sys
import os
from time import *

##########################################################################
#
# if indel happens at the beginging of the gene (i.e., in the first 75 DNA bases),
# , or within < 5th percentile of length,  then look for the next in-frame ATG,
# that there is only one protein 1.
# otherwise, start translating from the very begining.
#
###########################################################################
BEGINING_OFFSET = 75  


###############################################################
#
# This function extract all features about each gene
#
###############################################################
def extractFeatures(affectedEnstList, chrNum, enstInfoDict, \
                    DNAConvervationScores,fileName, codinginfoPath, buildVersion):
    fout = open(fileName, 'a')
    fout.write(affectedEnstList[0][0]+'\t'+affectedEnstList[0][1]+'\t')
    
##    numExons = []
##    for t in affectedEnstList:
##        #numExons.append(t[-2])
##        numExons.append(t[15])
##    avgNumExons = float(sum(numExons))/len(numExons)
##    fout.write(str(avgNumExons)+'\t')   #1) AvgNumExons: average number of exons of the gene
##        
##    
##    fout.write(str(affectedEnstList[0][11])+'\t') #2) num_of_tx_unaffected 
##    
##    fout.write(str(float(affectedEnstList[0][11])/affectedEnstList[0][12])+'\t') #3) fract_of_tx_unaffected
##
    #4,5,6) avg_relative_indel_location, max_relative_indel_location, min_relative_indel_location
    tempList = []
    for t in affectedEnstList:
        tempList.append(float(t[3])/t[4])
    avg_relative_indel_location = sum(tempList)/len(tempList)
    max_relative_indel_location = max(tempList)
    min_relative_indel_location = min(tempList)
##    fout.write(str(avg_relative_indel_location)+'\t'+str(max_relative_indel_location)+'\t'+\
##               str(min_relative_indel_location)+'\t')
    fout.write(str(max_relative_indel_location)+'\t')
    

##    #feature 7: centered_avg_relative_indel_location
##    centered_avg_relative_indel_location = abs(0.5 - avg_relative_indel_location)
##    fout.write(str(centered_avg_relative_indel_location)+'\t')
##    
##    #8) avg_overlap_of_protein_1
##    overlaps = []
##    for t in affectedEnstList:
##        #overlaps.append(t[5])
##        overlaps.append(float(t[5])/t[4])
##    fout.write(str(float(sum(overlaps))/len(overlaps))+'\t')
##    
##    #9) avg_overlap_of_protein_2
##    overlaps = []
##    for t in affectedEnstList:
##        #overlaps.append(t[6])
##        overlaps.append(float(t[6])/t[4])
##    fout.write(str(float(sum(overlaps))/len(overlaps))+'\t')
##    
##    #10) avg_overlap_of_protein_mostDomains: avg overlap of protein
##    #    (1 or 2) with most number of protein domains with the original protein
##    overlaps = []
##    for t in affectedEnstList:
##        if t[8]>t[7]:
##            #overlaps.append(t[6])
##            overlaps.append(float(t[6])/t[4])
##        else:
##            #overlaps.append(t[5])
##            overlaps.append(float(t[5])/t[4])
##    fout.write(str(float(sum(overlaps))/len(overlaps))+'\t')
##
##    #11) fraction_with_mRNA_Decay
##    num_mRNA_decay = 0
##    for t in affectedEnstList:
##        if t[13] != 'NONE':
##            num_mRNA_decay += 1
##    fract_num_mRNA_decay = float(num_mRNA_decay)/affectedEnstList[0][12]
##    fout.write(str(fract_num_mRNA_decay)+'\t')

    #12) fract_AllDomain_affected
    #13) fract_PfamDomain_affected
    #14) fract_SuperfamDomain_affected
    #15) fract_SignalPeptideDomain_affected
    indel = affectedEnstList[0][0]
    idENSG = affectedEnstList[0][1]

    #find the affected regions, also all coding exons of all transcripts of this gene
    affectedRegions, codingExonsList, strandList, idENSPList = \
                     findAffectedRegions(indel, idENSG, chrNum, enstInfoDict, affectedEnstList)

    #faeture 16: frac_conserved_affected (fraction of conserved regions affected)
    fractConservedAffected = calFractDNAConservedAffected(affectedRegions, codingExonsList, DNAConvervationScores[chrNum])
    fout.write(str(fractConservedAffected))
    fout.write('\t')

    #feature 17: minimum distance of indel to the exon boundary
    tempList = []
    for t in affectedEnstList:
        #tempList.append(t[-1])
        tempList.append(t[16])
    fout.write(str(min(tempList))+'\t')

##    #feature 18: the number of paralogs
##    id_ENSG = affectedEnstList[0][1]
##    if id_ENSG not in paralogDict:
##        fout.write('0\t')
##    else:
##        fout.write(str(len(paralogDict[id_ENSG]))+'\t')
##    
##    
##
##    #feature 19: average ka/ks
##    id_ENSG = affectedEnstList[0][1]
##    if id_ENSG not in kaks:
##        fout.write('-1\t')
##    else:
##        fout.write(str(kaks[id_ENSG])+'\t')

    #codingExonsList, strandList
    #checkFrames(codingExonsList, strandList)
##    a = assignFrames(codingExonsList, strandList)
##    if not a:
##        print('ALARM!!!!!!!!')
##    else:
##        print('a='+str(len(a)))

        
    #feature 20: fract_conservedAminoAcid_affected: percentage of conserved amino acid positions lost due to indel
    proteins = []
    for t in affectedEnstList:
        proteins.append((t[17],t[2],t[18])) #idENSP, protein sequence, protien 1 sequence
    
        #print((t[17],t[2],t[18]))

   # print('\n\n\n\n\n***************')
    #print(calcFractProteinConservedAffected(proteins))
   
    proteinLostFractions = calcFractProteinConservedAffected(proteins, codinginfoPath, buildVersion)
##    for item in proteinLostFractions:
##        fout.write(str(item)+'\t')
    fout.write(str(proteinLostFractions[0]))
        
    fout.write('\n')
    
    
    

    fout.close()

def findConservationFile(idENSP, codinginfoPath, buildVersion):
    '''
    Given a protein ID, finds if there is a conservation score file.
    If yes, returns a tuple of (proteinSeq, conservation score list, >medianFlag,>25%Flag, >75%Flag)
    '''
    path = codinginfoPath + '/indelFiles_V66/Percentiles_V66/'
    if buildVersion == 'hs38':
        path = codinginfoPath + '/indelFiles_hg38_v83/Percentiles/'
    fileName = idENSP+'.fasta.query.globalX.info.percentiles'

    findFileFlag = False
    if os.path.exists(path+'inconsistent_proteins_JH_032311/'+ fileName):
        fullFileName = path+'inconsistent_proteins_JH_032311/'+ fileName
        findFileFlag = True
    else:
# Oct 31, 2014 Pauline bug fix
        for i in range(10):
            if os.path.exists(path+str(i)+'/'+ fileName):
                fullFileName = path+str(i)+'/'+ fileName
                findFileFlag = True
                break

    if findFileFlag:    #find the file
        aseq = ''
        indexes = []
        cScores = []
        _50Flags = ''
        _25Flags = ''
        _75Flags = ''
        fin = open(fullFileName, 'r')
        for line in fin:
            line = line.split()
            indexes.append(int(line[0][:-1]))
            aseq += line[0][-1]
            cScores.append(float(line[1]))
            _50Flags += line[2]
            _25Flags += line[3]
            _75Flags += line[4]
        fin.close()
        pTuple = (aseq, indexes, cScores, _50Flags, _25Flags, _75Flags)
            
    if  findFileFlag:
        return pTuple
    else:
        return None
                
            

def overlapRange(originalProteinSeq, protein1):
    '''
    Find the overlapping sequence range between two proteins
    '''
    if len(protein1)==0:
        return None
    mPos = []
    for i in range(len(originalProteinSeq)):
        if originalProteinSeq[i]==protein1[0]:
            mPos.append(i)

    maxOverlap = 0
    overlapRegion= None
    for s in mPos:
        tempSeq = originalProteinSeq[s:]
        i = 0
        while i<min(len(tempSeq), len(protein1)):
            if tempSeq[i]!=protein1[i]:
                break
            i+=1
        if i>maxOverlap:
            maxOverlap = i
            overlapRegion = (s+1, s+i)

    #print(overlapRegion)
    return overlapRegion


def percentageLost(pTuple, overlapresidues):
    '''
    This function calculates what percentage of conserved residues
    are lost for the protein due to indel
    '''
    _75totalConservedResidues = 0
    _50totalConservedResidues = 0
    _25totalConservedResidues = 0
    _75lostConservedResidues = 0
    _50lostConservedResidues = 0
    _25lostConservedResidues = 0
    for i in range(len(pTuple[1])):
        if pTuple[5][i] == '1':
            _75totalConservedResidues+= 1
            if pTuple[1][i]<overlapresidues[0] or pTuple[1][i]>overlapresidues[1]:
                _75lostConservedResidues+=1
        if pTuple[4][i] == '1':
            _25totalConservedResidues+= 1
            if pTuple[1][i]<overlapresidues[0] or pTuple[1][i]>overlapresidues[1]:
                _25lostConservedResidues+=1
        if pTuple[3][i] == '1':
            _50totalConservedResidues+= 1
            if pTuple[1][i]<overlapresidues[0] or pTuple[1][i]>overlapresidues[1]:
                _50lostConservedResidues+=1

    if _75totalConservedResidues == 0:
        _75lostFraction = 0
    else:
        _75lostFraction = float(_75lostConservedResidues)/_75totalConservedResidues

    if _50totalConservedResidues == 0:
        _50lostFraction = 0
    else:
        _50lostFraction = float(_50lostConservedResidues)/_50totalConservedResidues

    if _25totalConservedResidues == 0:
        _25lostFraction = 0
    else:
        _25lostFraction = float(_25lostConservedResidues)/_25totalConservedResidues
                
    return (_25lostFraction, _50lostFraction, _75lostFraction)

def calcFractProteinConservedAffected(proteins, codinginfoPath, buildVersion):
    '''
    This function calculates the fraction of lost conserved amino acid residues
    '''
    lostFractions = [[], [], []]
    for p in proteins:
        idENSP = p[0]
        if p[1][-1]=='*':
            proteinSeq = p[1][:-1]
        else:
            proteinSeq = p[1]
        if p[2][-1]=='*':
            p1Seq = p[2][:-1]
        else:
            p1Seq = p[1]

        pTuple = findConservationFile(idENSP, codinginfoPath, buildVersion) #aseq, indexes, cScores, _50Flags, _25Flags, _75Flags
        
        if pTuple!=None:
            if pTuple[0] == proteinSeq:
                overlapResidues = overlapRange(proteinSeq, p1Seq)
                lostPercentage = percentageLost(pTuple, overlapResidues)
                #lostFractions.append(lostPercentage)
                for i in range(len(lostFractions)):
                    lostFractions[i].append(lostPercentage[i])
            else:
                for i in range(len(lostFractions)):
                    lostFractions[i].append(0)
        else:
            for i in range(len(lostFractions)):
                lostFractions[i].append(0)

    

    return (max(lostFractions[0]), max(lostFractions[1]), max(lostFractions[2]))

        
    


###############################################################
#
# This function calculates the fract_domain_affected
#
###############################################################
def calFractDNAConservedAffected(affectedRegions, codingExonsList, conservations):
    numTotalConservedPositions = 0
    numAffectedConservedPositions = 0
    for region in affectedRegions:
        for ind in range(region[0], region[1]+1):
            if ind in conservations and conservations[ind]>=1:
                numAffectedConservedPositions+=1

    allIndexList = []
    unionExons = getSuperTranscript(codingExonsList)
    for exon in unionExons:
        for ind in range(exon[0], exon[1]+1):
            if ind in conservations and conservations[ind]>=1:
                numTotalConservedPositions+=1
                
##    print 'numTotalConservedPositions'
##    print numTotalConservedPositions
##    print 'numAffectedConservedPositions'
##    print numAffectedConservedPositions
##    print

    if numTotalConservedPositions==0:
        return 0
    else:
        return float(numAffectedConservedPositions)/numTotalConservedPositions
                
            
 


###############################################################
#
# This function see if a position falls in a functional domains
#
###############################################################
def inFuncDomains(domainList, proteinInd):
    if len(domainList) == 0: #if no functional domains for this protein
        return (False, False, False, False)

    inDomainFlag = False   #it is in a protein functional domain
    inPFDomainFlag = False #it is in a Pfamily domain
    inSMDOmainFlag = False #it is in a Super family domain
    inSPDomainFlag = False #it is in a Signal Peptide family domain
    for funcDomain in domainList:
        ind1 = funcDomain[0]
        ind2 = funcDomain[1]
        funcDomainType = funcDomain[2]
        if proteinInd >= ind1 and proteinInd <=ind2:
            inDomainFlag = True     #this position is in functional domain
            if funcDomainType[:2] == 'PF':
                inPFDomainFlag = True  #this position is in functional domain but also PFam domains
            if funcDomainType[:2] == 'SM':
                inSMDOmainFlag = True
            if funcDomainType == 'Sigp':
                inSPDomainFlag = True
    return (inDomainFlag, inPFDomainFlag, inSMDOmainFlag, inSPDomainFlag)
    

###############################################################
#
# This function calculates the fract_domain_affected
#
###############################################################
def calFractDomainAffected(affectedRegions, codingExonsList, ProteinDomains, strandList, idENSPList):
    #dictionary marking if a position is in functional site (any functional site or PFam site)
    funcDomainMark = {}
    for exonList in codingExonsList:
        print 'coding Exon'
        print exonList
        for exon in exonList:
            for i in range(exon[0], exon[1]+1):
                #the first item markes if a position in any func site
                #the second item markes if a position in only Pfam site
                if i not in funcDomainMark:
                    funcDomainMark[i] = [False, False, False, False]

    for ind in range(len(codingExonsList)):
        id_ENSP = idENSPList[ind]
        if id_ENSP not in ProteinDomains:
            continue
        exonList = codingExonsList[ind]
        cdnaInd = 0
        for exon in exonList:
            if strandList[ind]=='1':
                for i in range(exon[0], exon[1]+1):
                    cdnaInd += 1
                    proteinInd = (cdnaInd - 1) // 3 + 1
                    domainFlag = inFuncDomains(ProteinDomains[id_ENSP], proteinInd)
                    funcDomainMark[i] = domainFlag
                    
            else:
                for i in range(exon[1], exon[0]-1, -1):
                    cdnaInd += 1
                    proteinInd = (cdnaInd - 1) // 3 + 1
                    domainFlag = inFuncDomains(ProteinDomains[id_ENSP], proteinInd)
                    funcDomainMark[i] = domainFlag
##        print id_ENSP
##        print exonList
##        print ProteinDomains[id_ENSP]

    numAllDomains = 0           #total number of positions in functional domains
    numAllPFDomains = 0         #total number of positions in Pfam domains
    numAllSFDomains = 0         #total number of positions in Superfamily domains
    numAllSPDomains = 0         #total number of positions in Signal Peptides
    
    numAffectedDomains = 0      #number of affected positions in functional domains
    numAffectedPFDomains = 0    #number of affected positions in Pfam domains
    numAffectedSFDomains = 0    #number of affected positions in Superfamily domains
    numAffectedSPDomains = 0    #number of affected positions in signal peptides
    

    

    for ind in funcDomainMark:
        if funcDomainMark[ind][0]:
            numAllDomains+=1
        if funcDomainMark[ind][1]:
            numAllPFDomains+=1
        if funcDomainMark[ind][2]:
            numAllSFDomains+=1
        if funcDomainMark[ind][3]:
            numAllSPDomains+=1

    

    for region in affectedRegions:
        for ind in range(region[0], region[1]+1):
            if funcDomainMark[ind][0]:
                numAffectedDomains += 1
            if funcDomainMark[ind][1]:
                numAffectedPFDomains += 1
            if funcDomainMark[ind][2]:
                numAffectedSFDomains += 1
            if funcDomainMark[ind][3]:
                numAffectedSPDomains += 1

##    print 'number of residues: '+str(len(funcDomainMark.keys()))
##    print 'numAllDomains: '+str(numAllDomains)
##    print 'numAffectedDomains: '+str(numAffectedDomains)
##    print
##    print 'numAllPFDomains: '+str(numAllPFDomains)
##    print 'numAffectedPFDomains: '+str(numAffectedPFDomains)
##    print
##    print 'numAllSFDomains: '+str(numAllSFDomains)
##    print 'numAffectedSFDomains: '+str(numAffectedSFDomains)
##    print
##    print 'numAllSPDomains: '+str(numAllSPDomains)
##    print 'numAffectedSPDomains: '+str(numAffectedSPDomains)
##    
##    print
##    print

    
    affectedDomainPercentage = [None, None, None, None]
    if numAllDomains != 0:
        affectedDomainPercentage[0] = float(numAffectedDomains)/numAllDomains
    if numAllPFDomains!=0:
        affectedDomainPercentage[1] = float(numAffectedPFDomains)/numAllPFDomains
    if numAllSFDomains!=0:
        affectedDomainPercentage[2] = float(numAffectedSFDomains)/numAllSFDomains
    if numAllSPDomains!=0:
        affectedDomainPercentage[3] = float(numAffectedSPDomains)/numAllSPDomains
        
    return affectedDomainPercentage

    


###########################################################################
#
# This function finds the information about coding exons:
# such as how many coding exons, which coding exon indel happens,
# and the exact starting and ending position of each coding exon (not
# including UTRs.
#
###########################################################################
def findCodingExonInfo(startPos, endPos, UTR5_Starts, UTR5_Ends, UTR3_Starts, UTR3_Ends,\
                       exons, ranks):
    exonStarts = exons[0].split(';')
    exonEnds = exons[1].split(';')
    ranks = ranks.split(';')
    #get exons in order (including UTRs) 
    exons = [None]*(len(ranks)+1)
    for i in range(len(ranks)):
        ind = int(ranks[i])
        exons[ind] = [int(exonStarts[i]), int(exonEnds[i])]
    
    numUTR5 = 0
    numUTR3 = 0

   
    #if UTR5_Starts !='':
    UTR5_Starts = UTR5_Starts.split(';')
    UTR5_Ends = UTR5_Ends.split(';')

    #if UTR3_Starts !='':
    UTR3_Starts = UTR3_Starts.split(';')
    UTR3_Ends = UTR3_Ends.split(';')

    if UTR5_Starts[0] != '':
        for i in range(len(UTR5_Starts)):
            ind1 = int(UTR5_Starts[i])
            ind2 = int(UTR5_Ends[i])
            numUTR5+=ind2-ind1+1

    if UTR3_Starts[0] != '':
        for i in range(len(UTR3_Starts)):
            ind1 = int(UTR3_Starts[i])
            ind2 = int(UTR3_Ends[i])
            numUTR3+=ind2-ind1+1

    
    #find coding exons and their starting and ending positions
    codingExons = []
    for exon in exons[1:]:
        trimedExon = exon
        flag1 = True
        flag2 = True
        if UTR5_Starts[0] != '':
            #splice UTR5'
            for i in range(len(UTR5_Starts)):
                ind1 = int(UTR5_Starts[i])
                ind2 = int(UTR5_Ends[i])
                #the whole exon is UTR
                if ind1==trimedExon[0] and ind2==trimedExon[1]:
                    flag1 = False
                    break
                else:
                    #if strand=='1': #if positve chain and UTR5 covers part of exon
                    if ind1==trimedExon[0] and ind2<trimedExon[1]:
                        #codingExons.append([ind2+1, exon[1]])
                        trimedExon = [ind2+1, trimedExon[1]]
                        #flag1 = False
                        break
                    #else:           #if negative chain and UTR5 covers part of exon
                    elif ind1>trimedExon[0] and ind2==trimedExon[1]:
                        #codingExons.append([exon[0], ind1-1])
                        trimedExon = [trimedExon[0], ind1-1]
                        #flag1 = False
                        break
        
        if UTR3_Starts[0] != '':
            #splice UTR3'
            for i in range(len(UTR3_Starts)):
                ind1 = int(UTR3_Starts[i])
                ind2 = int(UTR3_Ends[i])
                #the whole exon is UTR
                if ind1==trimedExon[0] and ind2==trimedExon[1]:
                    flag2 = False
                    break
                else:
                    #if strand=='1': #if positve chain and UTR3 covers part of exon
                    if ind1>trimedExon[0] and ind2==trimedExon[1]:
                        #codingExons.append([exon[0], ind1-1])
                        trimedExon = [trimedExon[0], ind1-1]
                        #flag2 = False
                        break
                    #else:           #if negative chain and UTR3 covers part of exon
                    elif ind2<trimedExon[1] and ind1==trimedExon[0]:
                        #codingExons.append([ind2+1, exon[1]])
                        trimedExon = [ind2+1, trimedExon[1]]
                        #flag2 = False
                        break
       
        if flag1 and flag2:       #it is not trimmed, or it is only partial trimed
            codingExons.append(trimedExon)

    #get the coding exon # where indel sits
    whichExon = 0
    if startPos!=endPos:
        startPos+=1
    for i in range(len(codingExons)):
        exon = codingExons[i]
        if startPos>=exon[0] and startPos<=exon[1] and endPos>=exon[0] \
           and endPos<=exon[1]:
            whichExon = i + 1
            break

    
    return codingExons,  whichExon

def subExon(exon1, exon2):
    '''
    Substract exon1 by exon2
    '''
    if exon2[0]>exon1[0] and exon2[1]<exon1[1]:
        return [exon1[0], exon2[0]-1], [exon2[1]+1, exon1[1]]
    elif exon2[0]==exon1[0] and exon2[1]<exon1[1]:
        return [exon2[1]+1, exon1[1]], None
    elif exon2[0]>exon1[0] and exon2[1]==exon1[1]:
        return [exon1[0], exon2[0]-1], None
    elif exon2[0]==exon1[0] and exon2[1]==exon1[1]:
        return None, None
    elif exon2[0]<exon1[0]:
        if exon2[1]>=exon1[0] and exon2[1]<exon1[1]:
            return [exon2[1]+1, exon1[1]], None
        else:
            return None, None
    elif exon2[1]>exon1[1]:
        if exon2[0]>exon1[0] and exon2[0]<=exon1[1]:
            return [exon1[0], exon2[0]-1], None
        else:
            return None, None
        
###############################################################
#
# This function substracts transcript2 from transcript 1, which
# means to find regions which are on transcript 1, but not on
# transcript 2.
#
###############################################################
def substractTranscript(transcript1, transcript2):
    if len(transcript2)==0:
        return transcript1

    blist = transcript1[:]

    for exon2 in transcript2:
        for exon1 in blist:
            #print('exon1 '+str(exon1))
            #print('exon2 '+str(exon2))
            if isExonOverlap(exon1, exon2):
                flag = True
                e1, e2 = subExon(exon1, exon2)
                if e1!=None:
                    blist.append(e1)
                if e2!=None:
                    blist.append(e2)
                #print (e1, e2)
                blist.remove(exon1)
            
            
    blist.sort()      

    return blist


##
##def substractTranscript(transcript1, transcript2):
##    if len(transcript2)==0:
##        return transcript1
##    substractedTranscript = []
##    for exon in transcript1:
##        for ind in range(exon[0], exon[1]+1):
##            if inExon(ind, [transcript1]) and not inExon(ind, [transcript2]):
##                if not inExon(ind-1, [transcript1]) or inExon(ind-1, [transcript2]):
##                    ind1 = ind
##                elif not inExon(ind+1, [transcript1]) or inExon(ind+1, [transcript2]):
##                    ind2 = ind
##                    substractedTranscript.append((ind1,ind2))
##
##    return substractedTranscript

def inExon(ind, transcriptList):
    for exons in transcriptList:
        for exon in exons:
            if ind>=exon[0] and ind<=exon[1]:
                return True
    return False
            
           


def isSubRange(exon1, exon2):
    '''
    Return True if exon1 is part of exon2
    '''
    if exon1[0]>=exon2[0] and exon1[1]<=exon2[1]:
        return True
    else:
        return False

    
def isExonOverlap(exon1, exon2):
    '''
    Return true if exon1 and exon2 has overlap
    '''
    if exon2[0]>=exon1[0] and exon2[0]<=exon1[1] or \
       exon2[1]>=exon1[0] and exon2[1]<=exon1[1]:
        return True
    else:
        return False

    
def mergeTwoExons(exon1, exon2):
    return [min(exon1[0], exon2[0]), max(exon1[1], exon2[1])]
    
###############################################################
#
# This function gets a super transcript by combining the
# exon regions of many transcripts
#
###############################################################
def getSuperTranscript(transcriptList):
    if len(transcriptList)==0:
        return []

    alist = []
    for transcript in transcriptList:
        for exon in transcript:
            if exon not in alist:
                alist.append(exon)

    alist.sort()

    blist = alist[:]
    
    Flag = True
    while Flag:
        Flag = False
        clist = []
        for i in range(len(blist)-1):
            exon1 = blist[i]
            exon2 = blist[i+1]
            if isExonOverlap(exon1, exon2):
                Flag = True
                newExon = mergeTwoExons(exon1, exon2)
                for j in range(len(blist)):
                    if j!=i and j!=i+1:
                        clist.append(blist[j])
                clist.append(newExon)
                clist.sort()
                
                
                blist = clist[:]
                break

        

    return blist
        

    
    
   
            
    
#################################################################
###
### This function gets a super transcript by combining the
### exon regions of many transcripts
###
#################################################################
##def getSuperTranscript(transcriptList):
##    if len(transcriptList)==0:
##        return []
##    flag = True
##    for t in transcriptList:
##        if len(t)>0:
##            flag = False
##    if flag:        #transcriptList is just a list of empty transcripts
##        return []
##    boundaryList = []
##    for transcript in transcriptList:
##        if len(transcript)==0:
##            continue
##        curminInd = min([x[0] for x in transcript])
##        curmaxInd = max([x[1] for x in transcript])
##        boundaryList.append(curminInd)
##        boundaryList.append(curmaxInd)
##
##    minInd = min(boundaryList)
##    maxInd = max(boundaryList)
##
##    superList = []
##    for ind in range(minInd, maxInd+1):
##        if inExon(ind, transcriptList) and not inExon(ind-1, transcriptList):
##            ind1 = ind
##        elif inExon(ind, transcriptList) and not inExon(ind+1, transcriptList):
##            ind2 = ind
##            superList.append((ind1, ind2))
##
##    return superList
##    
###############################################################
#
# This function tells if a position in in exon region of
# any transcripts of a gene
#
###############################################################
def inExon(ind, transcriptList):
    for exons in transcriptList:
        for exon in exons:
            if ind>=exon[0] and ind<=exon[1]:
                return True
    return False
    
###############################################################
#
# This function find all transcripts belong to the gene which
# contains the indel. Then merge all transcripts to create a
# single universal transcript that contains all coding exons.
# Then it removes the unaffected regions and keep the affected
# regions. It then calcualtes the frac_domain_affected and
# fract_conserved_affected
#
###############################################################
def findAffectedRegions(indel, idENSG, chrNum, enstInfoDict, affectedEnstList):
    universalTranscript = [None]
    unaffectedRegions = [None]
    affectedRegions = [None]
    
    indel = indel.split(',')
    if chrNum!='X' and chrNum!='Y':
        chrNum = int(chrNum)
    #chrNum = int(indel[0])
    startPos = int(indel[1])
    endPos = int(indel[2])
    allel = indel[3]

    
    transcriptExonDict = dict() #dictionary with the key as ID_ENST and value as the exon list

    enstInfoList = enstInfoDict[chrNum]
   
    allTranscripts = []  #all transcripts of a gene
    nonAffectedTranscripts = [] #non-affected transcripts of a gene
    strandList = [] #a list shows strands of each transcript
    idENSPList = [] #a list of idENSP of each transcript
    for enst in enstInfoList:
        info = enst.split()
        title = info[0]
        title = title.split('|')

        id_ENSG = title[0]      #gene ID
        if id_ENSG != idENSG:
            continue
        
        id_ENST = title[1]      #transcript ID
        chroNum = title[2]      #chorosome number
        id_ENSP = title[3]      #protein ID
        CDs = title[4:6]        #CDs
        UTR5_Starts = title[6]  #UTR5' starts  
        UTR5_Ends = title[7]    #UTR5' ends
        UTR3_Starts = title[8]  #UTR3' starts
        UTR3_Ends = title[9]    #UTR3' ends
        exons = title[10:12]    #exons
        ranks = title[13]       #ranks of exons
        strand = title[12]      #strand

        codingExons, whichExon = findCodingExonInfo(startPos, endPos, \
                                                    UTR5_Starts, UTR5_Ends, UTR3_Starts, UTR3_Ends,exons, ranks)
        
        transcriptExonDict[id_ENST] = codingExons
        
        allTranscripts.append(codingExons)
        if whichExon == 0:
            nonAffectedTranscripts.append(codingExons)
        strandList.append(strand)
        idENSPList.append(id_ENSP)
            

##    print 'allTranscripts'
##    print allTranscripts
    #step 1: Creates the universal transcripts by taking the union of all transcripts. 
    universalTranscript = getSuperTranscript(allTranscripts)  #universal transcript, which is a union of all transcripts
    #step 2: Create a union of unaffected transcripts. 
    unionUnaffected = getSuperTranscript(nonAffectedTranscripts)  #super transcript which is not affected

    #universalTranscript = [(1,100), (200,300), (400,500), (600,700)]
    #unaffectedRegions = [(40,80), (200,250), (450,500), (600,620), (680,700)]
        
##    print 'UniversalTranscript'
##    print universalTranscript
##    print
##    print 'union of unaffected transcripts'
##    print unionUnaffected
##    print

    #step 3: If there are unaffected transcripts but without mRNA decay, then we have two conditions
    # Condition 1: If the indel is after first 50 bases, then functional regions are those bases before the indel
    # Condition 2: otherwise the functional regions are those bases starting from the next
    #              in-frame starting codon.

    
    noMRNAdecayTranscriptList = []    #a list of affected transcripts without mRNA decay
    for t in affectedEnstList:
        if t[13]=='NONE':
            #print '---->'+str(t[14])
            exons = transcriptExonDict[t[14]]
            #each item has the following information:
            # 1) exons; 2) indel position in the original protein; 3)original protein size;
            # 4)overlap of protein 1 with original protein
            noMRNAdecayTranscriptList.append([exons, t[3], t[4], t[5]])  
##    print 'noMRNAdecayTranscriptList'
##    print noMRNAdecayTranscriptList
##    print

    unionFunctional = getFunctional(noMRNAdecayTranscriptList, strand, startPos, endPos)
##    print 'unionFunctional'
##    print unionFunctional
##    print

    
    #step 4: take the union of unionFunctional and unionUnaffected
    if len(unionUnaffected)!=0 and len(unionFunctional)!=0:
        F = getSuperTranscript([unionUnaffected,unionFunctional])
    elif len(unionUnaffected)!=0 and len(unionFunctional)==0:
        F = getSuperTranscript([unionUnaffected])
    elif len(unionUnaffected)==0 and len(unionFunctional)!=0:
        F = getSuperTranscript([unionFunctional])
    else:
        F = []
##    print 'All functional: '
##    print F
##    print

    #step 5: substract universalTranscript by F, we got all the affected regions
    affectedRegions = substractTranscript(universalTranscript, F)
##    print 'affectedRegions'
##    print affectedRegions
##    print
##        
    return affectedRegions, allTranscripts, strandList, idENSPList


###############################################################
#
# This function gets functional regions of the affected
# transcripts without mRNA decay
#
###############################################################
def getFunctional(noMRNAdecayTranscriptList, strand, startPos, endPos):
    functionalList = []  #a list of functional regions of each transcript
    for t in noMRNAdecayTranscriptList:
        funcRegions = []
        #within the first BEGINING_OFFSET positions, all positions
        #after the next in-frame starting codon are functional
        if t[1]*3 - 3 <= BEGINING_OFFSET:
            originalProteinSize = t[2]
            overlap = t[3]
            offset = (originalProteinSize - overlap)*3  #the missed regions at 5'
            exons = t[0]
            for i in range(len(exons)):
                offset = offset - (exons[i][1]-exons[i][0]+1)
                if offset <0:
                    break
            offset += (exons[i][1]-exons[i][0]+1)
            if strand == '1':
                funcRegions.append((exons[i][0]+offset, exons[i][1]))
            else:
                funcRegions.append((exons[i][0], exons[i][1]-offset))
            for exon in exons[i+1:]:
                funcRegions.append(exon)
        else:
            exons = t[0]
            for exon in exons:
                if strand=='1':
                    if exon[1]<startPos:
                        funcRegions.append(exon)
                    elif exon[0]<startPos and startPos<=exon[1]:
                        funcRegions.append((exon[0], startPos-1))
                else:
                    if exon[0]>endPos:
                        funcRegions.append(exon)
                    elif exon[0]<=endPos and endPos<exon[1]:
                        funcRegions.append((endPos+1, exon[1]))
        functionalList.append(funcRegions)
        #print 'funcRegions'
        #print funcRegions

        #get the union of all functional regions of all
##    print 'functionalList'
##    print functionalList
    unionFunctional = getSuperTranscript(functionalList)
    return unionFunctional

    
   


###############################################################
#
# This function gets the intersection of two region lists
#
###############################################################
def getIntersection(t1, t2):
    alist = []
    for exon in t1:
        for ind in range(exon[0], exon[1]+1):
            if inExon(ind, [t2]):
                alist.append(ind)
    blist = []
    for ind in alist:
        if ind in alist and ind-1 not in alist:
            ind1 = ind
        elif ind in alist and ind+1 not in alist:
            ind2 = ind
            blist.append((ind1, ind2))

    return blist
    
    
###############################################################
#
# This function reads the transcripts from corresponding
# chromsome file
#
###############################################################
def readExons(chrNum, codinginfoPath, buildVersion):
    #chromsomeList = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']
    chromsomeList = [chrNum]
    enstInfoDict = {}
    for chrom in chromsomeList:
        enstInfoList = []
        if buildVersion == 'hs37':
            fin = open(codinginfoPath + '/indelFiles_V66/chromosomes_V66_update/chr'+str(chrom)+'_valid_2013_1_17.txt', 'r')
        else:
            fin = open(codinginfoPath + '/indelFiles_hg38_v83/chromosomes_2016/chr'+str(chrom)+'_valid_2016_8_16.txt', 'r')
        lines = fin.read()
        lines = lines.split('>')[1:]
        for info in lines:
            enstInfoList.append(info)
        fin.close()
        enstInfoDict[chrom] = enstInfoList
    return enstInfoDict


#######################################################
# This function read all the conservation scores for
# each coding position of all chromsomes.
#
#######################################################
def readConservationScores(chrNum, codinginfoPath, buildVersion):
    #chromsomeList = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']
    #chromsomeList = [1]
    chromsomeList = [chrNum]
    consvInfoDict = {}
    for chrom in chromsomeList:
        conservations = {}
        if buildVersion == 'hs37':
            fin = open(codinginfoPath + '/indelFiles_V66/conservationScores/chr'+str(chrom)+'_cds_conservations.txt', 'r')
        else:
            fin = open(codinginfoPath + '/indelFiles_hg38_v83/conservationScores/chr'+str(chrom)+'_cds_conservations.txt', 'r')
        for line in fin:
            line = line.split()
            if line[1]=='fixedStep':
                continue
            conservations[int(line[0])] = float(line[1])
        fin.close()
        consvInfoDict[chrom] = conservations
    return consvInfoDict

################################################################################
### This function read the files in domains directory, and extrat domain
### positions for each ENSP
###
#################################################################################
##def getDomains():
##    if not os.path.exists('../domains/proteinDomains_6_10.txt'):
##        print 'ERROR: domains/proteinDomains_6_10.txt not found!'
##        return False
##    fin = open('../domains/proteinDomains_6_10.txt','r')
##    domains = {}
##    for line in fin:
##        line = line.split()
##        id_ENSP = line[0]
##        domains[id_ENSP] = []
##        if len(line)<3:
##            continue
##        #id_ENSP = line[0]
##        #domains[id_ENSP] = []
##        for d in line[2:]:
##            d = d.split('-')
##            domains[id_ENSP].append((int(d[0]), int(d[1]), d[2]))
##    fin.close()
##    return domains
##
##def readParalogs(threshold):
##    '''
##    This function reads the paralog information of for each gene
##    '''
##    if not os.path.exists('data/paralogs.txt'):
##        print 'ERROR: data/paralogous.txt not found!'
##        return False
##    paralogDict = {}
##    fin = open('data/paralogs.txt', 'r')
##    line = fin.readline()
##    while True:
##        line = fin.readline()
##        if not line:
##            break
##        info = line.split()
##        if len(info)<4:
##            continue
##        idGene = info[0]
##        paralogID = info[2]
##        percenIdentity = float(info[3])
##        if percenIdentity>=threshold:
##            if idGene not in paralogDict:
##                paralogDict[idGene] = [paralogID]
##            else:
##                if paralogID not in paralogDict[idGene]:
##                    paralogDict[idGene].append(paralogID)
##
##    return paralogDict
##
##def readKaks():
##    '''
##    This function reads the Ka/ks (dN/dS) information
##    '''
##    if not os.path.exists('data/paralogs.txt'):
##        print 'ERROR: data/paralogous.txt not found!'
##        return False
## 
##    kaks = dict()
##    fin = open('data/kaks.txt', 'r')
##    for line in fin:
##        line = line.split()
##        if len(line)==4:
##            id_ENSG = line[0]
##            id_ENST = line[1]
##            dN = float(line[2])
##            dS = float(line[3])
##            if dS==0:
##                continue
##            if id_ENSG not in kaks:
##                kaks[id_ENSG] = {}
##            else:
##                if id_ENST not in kaks[id_ENSG]:
##                    kaks[id_ENSG][id_ENST] = [dN/dS]
##                else:
##                    kaks[id_ENSG][id_ENST].append(dN/dS)
##    fin.close()
##
##    
##
##    geneKaks = {}
##    for id_ENSG in kaks:
##        avgKaks = []
##        for key in kaks[id_ENSG]:
##            avgKaks.append(min(kaks[id_ENSG][key]))
##        geneKaks[id_ENSG] = sum(avgKaks)/len(avgKaks)
##
## 
##    return geneKaks
    
 

def main():

    start = time()
##    print'\n--------------------------------------\n'
##    print 'Indel Information Extraction Program!'
##    print 'Author: Jing Hu & Pauline Ng\n'
##    print '--------------------------------------\n'

    fin = open(sys.argv[1], 'r')
    chrNum = sys.argv[2]
    #Gina added
    codinginfoPath = sys.argv[4]
    #Herty added for build options
    buildVersion = sys.argv[5]
    if chrNum!='X' and chrNum!='Y':
        chrNum = int(chrNum)
##    print 'Reading all transcripts...'
    enstInfoDict = readExons(chrNum, codinginfoPath, buildVersion) #read all transcripts
##    print '.....................Done!\n'

##    #protein domains
##    print '-->Reading protein domains...'
##    ProteinDomains = getDomains()
##    print '........................Done!\n'

##    for x in ProteinDomains:
##        print ProteinDomains[x]
##
##    return
    
    #conservation scores
    #print '-->Reading conservation scores of all exons...'
    DNAConvervationScores = readConservationScores(chrNum, codinginfoPath, buildVersion)
    #print '.........................................Done!\n'

    end1 = time()
    

##    #paralogs
##    print '-->Reading paralogs...'
##    paralogDict = readParalogs(50)
##    print '.........................................Done!\n'
##
##    #ka/ks
##    print '-->Reading ka/ks information...'
##    kaks = readKaks()
##    print '.........................................Done!\n'
        
    outputFileName = sys.argv[3]
##    fout = open(outputFileName, 'w')
##    fout.write('Indel\tid_ENSG\tmax_relative_indel_location\tfract_conservedDNABase_affected\tmin_disExonBoundary\t_25_fract_conservedAminoAcid_affected\n')
##    fout.close()
    
    affectedEnstList = []
    while True:
        line = fin.readline()
        if not line:
            break
        if len(line.split()) == 0: #finished reading the transcripts of one gene
            #do processing
            if len(affectedEnstList)==0:
                continue
##            print '\n\n---->Processing '+ affectedEnstList[0][1]
##            affectedRegions = extractFeatures(affectedEnstList, chrNum, enstInfoDict, \
##                                              ProteinDomains, DNAConvervationScores,\
##                                              paralogDict, kaks,  outputFileName)
            affectedRegions = extractFeatures(affectedEnstList, chrNum, enstInfoDict, \
                                              DNAConvervationScores,outputFileName, codinginfoPath, buildVersion)
            
            affectedEnstList = []
            
        elif line[0]=='>': #new transcript
            line = line.split()
            indel = line[0][3:]                #[0]indel  
            idENSG = line[1]                #[1] gene ID
            idENST = line[2]                #[14] transcript ID
            idENSP = line[3]                #[17] protein ID
            numExons = int(line[5])         #[15] number of exons
            for i in range(4):
                line = fin.readline()
            line = fin.readline()   
            origProtein = line.split()[2]   #[2]original protein
            #print(origProtein)
            line = fin.readline()   
            line = line.split() 
            indelPos = int(line[0])         #[3]indel position  
            origProSize = int(line[1])      #[4]original protein size
            line = fin.readline()
            protein1 = line.split()[2]      #[18] protein1 sequence due to indel
            line = fin.readline()   
            overlapPro1 = int(line.split()[1]) #[5]overlap of protein 1 with original protein
            line = fin.readline()
            line = fin.readline()
            line = line.split()
            if len(line)==2:
                overlapPro2 = int(line[1]) #[6]overlap of protein 2 with original protein
            else:
                overlapPro2 = 0
            line = fin.readline()
            line = line.split()
            numFuncDomain1 = int(line[0])           #[7]number of functional domains in protein 1
            numFuncDomain2 = int(line[1])           #[8]number of functional domains in protein 2
            numAllFuncDomains = int(line[2])        #[9]number of all functional domains in the original protein
##            for i in range(3):
##                line = fin.readline()
            line = fin.readline()
            line = fin.readline().split()
            numAffectedTrans = int(line[1])         #[10]number of affected transcripts
            numUnaffectedTrans = int(line[4])       #[11]number of unaffected transcripts
            numAllTrans = int(line[7])              #[12]number of all transcripts of the gene
            #flagNMD = fin.readline().split()[1]     #[13]if this indel is at NMD of this transcript
            nmdValue = fin.readline()[12:].split('\n')[0] #[13]the notation of mRNAdecay (none, NMD, and nonstop mRNA decay)
            #print idENST

            
            line = fin.readline()
            line = fin.readline()
            disExonBoundary = int(line.split()[-1])  #16  distance to the boundary
            line = fin.readline()
            
            
            affectedEnstList.append((indel, idENSG, origProtein, indelPos, origProSize, overlapPro1, overlapPro2,\
                                    numFuncDomain1, numFuncDomain2,  numAllFuncDomains, numAffectedTrans,\
                                     numUnaffectedTrans, numAllTrans, nmdValue, idENST, numExons, \
                                     disExonBoundary, idENSP, protein1))
            affectedEnstList.append((indel,                #0\
                                     idENSG,               #1\
                                     origProtein,          #2\
                                     indelPos,             #3\
                                     origProSize,          #4\
                                     overlapPro1,          #5\
                                     overlapPro2,          #6\
                                    numFuncDomain1,        #7\
                                     numFuncDomain2,       #8\
                                     numAllFuncDomains,    #9\
                                     numAffectedTrans,     #10\
                                     numUnaffectedTrans,   #11\
                                     numAllTrans,          #12\
                                     nmdValue,             #13\
                                     idENST,               #14\
                                     numExons,             #15\
                                     disExonBoundary,      #16\
                                     idENSP,               #17\
                                     protein1))             #18\
##            affectedEnstList.append((indel,                #0\
##                                     idENSG,               #1\
##                                     origProtein,          #2\
##                                     indelPos,             #3\
##                                     origProSize,          #4\
##                                     overlapPro1,          #5\
##                                     None,          #6\
##                                    None,        #7\
##                                     None,       #8\
##                                     None,    #9\
##                                     None,     #10\
##                                     None,   #11\
##                                     None,          #12\
##                                     nmdValue,             #13\
##                                     idENST,               #14\
##                                     None,             #15\
##                                     disExonBoundary,      #16\
##                                     idENSP,               #17\
##                                     protein1))             #18\
        
            
                                
            
            

    fin.close()
    end2 = time()
##    print('Time spent in reading files:'+str(end1-start))
##    print('Time spent in total: '+str(end2-start))

main()
                  
    
    
