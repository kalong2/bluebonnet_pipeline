'''
Copyright 2013 Jing Hu. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of
      conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list
      of conditions and the following disclaimer in the documentation and/or other materials
      provided with the distribution.

THIS SOFTWARE IS PROVIDED BY Franklin & Marshall College ''AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those of the
authors and should not be interpreted as representing official policies, either expressed
or implied, of Franklin & Marshall College.
'''
###############################################################################
#
# Updated : 28/3/2013
#
# This program extract feature vectors for each indel at gene level.
# 
#
# Copyright: Jing Hu & Pauline Ng
# Franklin & Marshall College, Singapore Genome Institute
# Date: 28/3/2013
#
##############################################################################

import sys
import os
from time import *
#import time
##########################################################################
#
# if indel happens at the beginging of the gene (i.e., in the first 75 DNA bases),
# , or within < 5th percentile of length,  then look for the next in-frame ATG,
# that there is only one protein 1.
# otherwise, start translating from the very begining.
#
###########################################################################
BEGINING_OFFSET = 75  

def avgFeature(aseq, adict):
    '''
    It calculates the average feature of aseq
    '''
    if len(aseq)==0:
        return '?'
    totalValue = 0
    num = 0
    for ch in aseq:
        if ch in adict:
            totalValue+=adict[ch]
            num+=1

    if num == 0:
        return 0
    return float(totalValue)/num

def hasHydro(astr):
    if len(astr)==0:
        return '?'
    for aa in 'AQE':
        if aa in astr:
            return '1'
    return '0'

def hasStrBreak(astr):
    if len(astr)==0:
        return '?'
    for aa in 'PGDS':
        if aa in astr:
            return '1'
    return '0'


###############################################################
#
# This function extract all features about each gene
#
###############################################################
def extractFeatures(affectedEnstList, chrNum, enstInfoDict, ProteinDomains,\
                    DNAConvervationScores, disorder, fileName):
    fout = open(fileName, 'a')
    fout.write(affectedEnstList[0][0]+'\t'+affectedEnstList[0][1]+'\t')
    
    tempList = []
    for t in affectedEnstList:
        tempList.append(float(t[3])/t[4])

    
    #Feature 1: DNA conservation scores on left side of indel position
    indel = affectedEnstList[0][0]
    loc1 = int(indel.split(',')[1])
    loc2 = int(indel.split(',')[2])
    dnaConservations = DNAConvervationScores[chrNum]
    fout.write(str(dnaConservations.get(loc1, '?'))+'\t')
 #              str(dnaConservations.get(loc2+1, '?'))+'\t')
    
    idENSG = affectedEnstList[0][1]

    
    proteins = []
    for t in affectedEnstList:
        proteins.append((t[17],t[2],t[18])) #idENSP, protein sequence, protien 1 sequence

 
    affectedDomainPercentage, avgDScorePerGene = \
                         calcProteinEffects(proteins, ProteinDomains,  disorder)

    #Feature 2: Precentage of PFDomain affected
    fout.write(str(affectedDomainPercentage[1])+'\t')


    #Feature 3: average disorder scores in the affected region
    fout.write(str(avgDScorePerGene))

    fout.write('\n')
    
    fout.close()

def findConservationFile(idENSP,codinginfoPath):
    '''
    Given a protein ID, finds if there is a conservation score file.
    If yes, returns a tuple of (proteinSeq, conservation score list, >medianFlag,>25%Flag, >75%Flag)
    '''
    #gina
    path = codinginfoPath + '/indelFiles_V66/Percentiles_V66/'
    fileName = idENSP+'.fasta.query.globalX.info.percentiles'

    findFileFlag = False
    if os.path.exists(path+'inconsistent_proteins_JH_032311/'+ fileName):
        fullFileName = path+'inconsistent_proteins_JH_032311/'+ fileName
        findFileFlag = True
    else:
        for i in range(9):
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


    if overlapRegion==None:
        return None
    
    if  overlapRegion[1]==len(originalProteinSeq):
        #print('\n\nprotein overlapping regions')
        #print(overlapRegion)
        
        return [overlapRegion]

    for i in range(-1, -1*len(originalProteinSeq), -1):
        if abs(i)>min(len(originalProteinSeq), len(protein1)) or originalProteinSeq[i]!=protein1[i] \
           or len(originalProteinSeq)+i+1==overlapRegion[1] or len(protein1)+i+1==overlapRegion[1]:
            break

    overlapRegion2 = []

    if i<-1:
        overlapRegion2 = ((len(originalProteinSeq)+i+2, len(originalProteinSeq)))
                   

    #print('\n\nprotein overlapping regions')
    #print(overlapRegion, overlapRegion2)
    
    if len(overlapRegion2)!=0:
        return [overlapRegion, overlapRegion2]
    else:
        return [overlapRegion]


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

def findConsScoreOnBothSide(pTuple, overlapResidues):
    '''
    The function returns the conservation score of surrounding amino acids,
    one amino acid on either side of the indel.
    '''
    conservationScores = pTuple[2]
    indexes = pTuple[1]
    leftIndex = overlapResidues[0][1]
    if len(overlapResidues)==2:
        rightIndex = overlapResidues[1][0]
    else:
        rightIndex = leftIndex
        

    if leftIndex in indexes:
        leftScore = conservationScores[indexes.index(leftIndex)]
    else:
        leftScore = '?'


    if rightIndex in indexes:
        rightScore = conservationScores[indexes.index(rightIndex)]
    else:
        rightScore = '?'
        
    return (leftScore, rightScore)
        

def findProteinFlankingSeq(overlapResidues, proteinSeq, p1Seq):
    #print(overlapResidues)
    if len(overlapResidues)==2:
        if overlapResidues[1][0]-overlapResidues[0][1]>1:  #deletion
            flankingSeq = proteinSeq[max(0, overlapResidues[0][1]-5):overlapResidues[0][1]]+'-'
            flankingSeq += proteinSeq[overlapResidues[0][1]:overlapResidues[1][0]-1]+'-'
            flankingSeq += proteinSeq[overlapResidues[1][0]-1:min(overlapResidues[1][0]+4, len(proteinSeq))]
        else:   #insertion
            #print('I am here')
            flankingSeq = proteinSeq[max(0, overlapResidues[0][1]-5):overlapResidues[0][1]]+'-'
            index1 = overlapResidues[0][1]
            #print(overlapResidues[1][0]-1,overlapResidues[1][1])
            index2 = overlapResidues[1][1]-overlapResidues[1][0]+1
            #print(index1, -index2)
            flankingSeq += p1Seq[index1:-1*index2]+'-'
            flankingSeq += proteinSeq[overlapResidues[1][0]-1:min(overlapResidues[1][0]+4, len(proteinSeq))]

    else:
        if len(proteinSeq)>len(p1Seq):  #deletion
            flankingSeq = proteinSeq[max(0, overlapResidues[0][1]-5):overlapResidues[0][1]]+'-'
            flankingSeq+= proteinSeq[overlapResidues[0][1]:]+'-'
        else:                           #insertion
            flankingSeq = proteinSeq[max(0, overlapResidues[0][1]-5):overlapResidues[0][1]]+'-'
            flankingSeq+= p1Seq[overlapResidues[0][1]:]+'-'

    #print(flankingSeq)
    return flankingSeq
            
    

def calcProteinEffects(proteins, ProteinDomains,  disorder):
    '''
    This function calculates effects on proteins, such information includes
    conservation scores on both sides of the indel,
    fraction of functional domains affected.
    '''
  #  conservationScore2Sides = [[], [],[]]

    numAllDomains = 0           #total number of positions in functional domains
    numAllPFDomains = 0         #total number of positions in Pfam domains
    numAllSFDomains = 0         #total number of positions in Superfamily domains
    numAllSPDomains = 0         #total number of positions in Signal Peptides
    
    numAffectedDomains = 0      #number of affected positions in functional domains
    numAffectedPFDomains = 0    #number of affected positions in Pfam domains
    numAffectedSFDomains = 0    #number of affected positions in Superfamily domains
    numAffectedSPDomains = 0    #number of affected positions in signal peptides

    avgDisorderScoreList = []
    numEList = []
    numHList = []
    numCList = []
    
    
    
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

 #       pTuple = findConservationFile(idENSP) #aseq, indexes, cScores, _50Flags, _25Flags, _75Flags

        overlapResidues = overlapRange(proteinSeq, p1Seq)
        #print('overlapResidues',overlapResidues)
        if overlapResidues==None:
            continue
        

        if idENSP in ProteinDomains:
            domains = ProteinDomains[idENSP]

            for domain in domains:
                domainType = domain[2]
                
                if domainType[:2] == 'PF':
                    numAllPFDomains += 1
                if domainType[:2] == 'SM':
                    numAllSFDomains+=1
                if domainType == 'Sigp':
                    numAllSPDomains+= 1
                numAllDomains += 1


            affectedResidues = []
            if len(overlapResidues)==2:
                if overlapResidues[0][1]==overlapResidues[1][0]:
                    affectedResidues = [overlapResidues[0][1]]
                else:
                    if overlapResidues[0][1]== overlapResidues[1][0] -1:
                        affectedResidues = [overlapResidues[0][1] + 0.5]
                    else:
                        for ind in range(overlapResidues[0][1]+1, overlapResidues[1][0]):
                            affectedResidues.append(ind)
            else:
                affectedResidues = [overlapResidues[0][1]+1]

            for domain in domains:
                ind1 = domain[0]
                ind2 = domain[1]
                domainType = domain[2]
                affectThisDomain = False
                for ind in affectedResidues:
                     if ind>=ind1 and ind<=ind2:
                         affectThisDomain = True
                         break
                if not affectThisDomain:
                    continue
                if domainType[:2] == 'PF':
                    numAffectedPFDomains += 1
                if domainType[:2] == 'SM':
                    numAffectedSFDomains+=1
                if domainType == 'Sigp':
                    numAffectedSPDomains+= 1
                numAffectedDomains += 1

        #get disorder Information
        indelResidues = []
        if len(overlapResidues)==2:
            if overlapResidues[0][1]==overlapResidues[1][0] or overlapResidues[0][1]== overlapResidues[1][0] -1:
                indelResidues = [overlapResidues[0][1]]
            else:
                for ind in range(overlapResidues[0][1]+1, overlapResidues[1][0]):
                    indelResidues.append(ind)
        else:
            indelResidues = [overlapResidues[0][1]+1]

        #print('ENSP', idENSP)
        if idENSP in disorder:
            dscores = disorder[idENSP]
            indelDScores = []
            #print(indelResidues)
            alist = [indelResidues[0]-1]+indelResidues+[indelResidues[-1]+1]
            for ind in alist:
                #print(ind)
                if ind-1 < len(dscores):
                    indelDScores.append(dscores[ind-1])
            if len(indelDScores)>0:
                avgDisorderScoreList.append(sum(indelDScores)/len(indelDScores))
           
                


  


    avgDScorePerGene = '?'
    if len(avgDisorderScoreList) != 0:
        avgDScorePerGene = sum(avgDisorderScoreList) / len(avgDisorderScoreList)

  

         
    affectedDomainPercentage = [0,0,0,0]
        

 
    
    if numAllDomains != 0:
        affectedDomainPercentage[0] = float(numAffectedDomains)/numAllDomains
    if numAllPFDomains!=0:
        affectedDomainPercentage[1] = float(numAffectedPFDomains)/numAllPFDomains
    if numAllSFDomains!=0:
        affectedDomainPercentage[2] = float(numAffectedSFDomains)/numAllSFDomains
    if numAllSPDomains!=0:
        affectedDomainPercentage[3] = float(numAffectedSPDomains)/numAllSPDomains




    return    affectedDomainPercentage, avgDScorePerGene
    


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
##    print indel
##    print idENSG
##    print chrNum
##    print enstInfoDict
##    print affectedEnstList

    
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

 
    noMRNAdecayTranscriptList = []    #a list of affected transcripts without mRNA decay
    for t in affectedEnstList:
        if t[13]=='NONE':
            #print '---->'+str(t[14])
            exons = transcriptExonDict[t[14]]
            #each item has the following information:
            # 1) exons; 2) indel position in the original protein; 3)original protein size;
            # 4)overlap of protein 1 with original protein
            noMRNAdecayTranscriptList.append([exons, t[3], t[4], t[5]])  

    unionFunctional = getFunctional(noMRNAdecayTranscriptList, strand, startPos, endPos)
    


    #step 4: take the union of unionFunctional and unionUnaffected
    if len(unionUnaffected)!=0 and len(unionFunctional)!=0:
        F = getSuperTranscript([unionUnaffected,unionFunctional])
    elif len(unionUnaffected)!=0 and len(unionFunctional)==0:
        F = getSuperTranscript([unionUnaffected])
    elif len(unionUnaffected)==0 and len(unionFunctional)!=0:
        F = getSuperTranscript([unionFunctional])
    else:
        F = []

    #step 5: substract universalTranscript by F, we got all the affected regions
    affectedRegions = substractTranscript(universalTranscript, F)
          
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
def readExons(chrNum,codinginfoPath):
    #chromsomeList = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']
    #Added by Herty
    buildVersion = sys.argv[5]
    chromsomeList = [chrNum]
    enstInfoDict = {}
    for chrom in chromsomeList:
        enstInfoList = []
        if buildVersion == 'hs37':
            fin = open(codinginfoPath + '/indelFiles_V66/chromosomes_V66_update/chr'+str(chrom)+'_valid_2013_1_17.txt', 'r')
        else:
            fin = open(  codinginfoPath + '/indelFiles_hg38_v83/conservationScores/chr' + str(chrom) + '_cds_conservations.txt', 'r')
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
def readConservationScores(chrNum,codinginfoPath):
    #chromsomeList = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']
    #chromsomeList = [1]
    #Added by Herty
    buildVersion = sys.argv[5]
    chromsomeList = [chrNum]
    consvInfoDict = {}
    for chrom in chromsomeList:
        conservations = {}
        if buildVersion == 'hs37':
            fin = open(codinginfoPath + '/indelFiles_V66/conservationScores/chr'+str(chrom)+'_cds_conservations.txt', 'r')
        else:
            fin = open(codinginfoPath + '/indelFiles_hg38_v83/conservationScores/chr' + str(chrom) + '_cds_conservations.txt', 'r')
        for line in fin:
            line = line.split()
            if line[1]=='fixedStep':
                continue
            conservations[int(line[0])] = float(line[1])
        fin.close()
        consvInfoDict[chrom] = conservations
    return consvInfoDict

##############################################################################
# This function read the files in domains directory, and extrat domain
# positions for each ENSP
#
###############################################################################
def getDomains(codinginfoPath):
    print "codinginfoPath: ", codinginfoPath
    #Added by Herty
    buildVersion = sys.argv[5]
    if buildVersion == 'hs37':
        if not os.path.exists(codinginfoPath + '/indelFiles_V66/domains/proteinDomains_2013_2_12.txt'):
            print ('ERROR: /indelFiles_V66/domains/proteinDomains.txt not found!')
            return False
    else:
        if not os.path.exists(codinginfoPath + '/indelFiles_hg38_v83/domains/proteinDomains_2016_5_8.txt'):
            print ('ERROR: /indelFiles_hg38_v83/domains/proteinDomains_2016_5_8.txt not found!')
            return False

    if buildVersion == 'hs37':
        fin = open(codinginfoPath+ '/indelFiles_V66/domains/proteinDomains_2013_2_12.txt','r')
    else:
        fin = open(codinginfoPath + '/indelFiles_hg38_v83/domains/proteinDomains_2016_5_8.txt', 'r')

    domains = {}
    for line in fin:
        line = line.split()
        id_ENSP = line[0]
        domains[id_ENSP] = []
        if len(line)<3:
            continue
        #id_ENSP = line[0]
        #domains[id_ENSP] = []
        for d in line[2:]:
            d = d.split('-')
            domains[id_ENSP].append((int(d[0]), int(d[1]), d[2]))
    fin.close()
    return domains


##############################################################################
# This function read the disorder information
# of proteins
#
###############################################################################
def readDisorderAA(chrNum,codinginfoPath):
    disorder = {}
    #Added by Herty
    buildVersion = sys.argv[5]
    if buildVersion == 'hs37':
        fin = open(codinginfoPath + '/indelFiles_V66/enspDisorder/ensp_disorder_'+str(chrNum), 'r')
    else:
        fin = open(codinginfoPath + '/indelFiles_hg38_v83/enspDisorder/ensp_chr' + str(chrNum) + '_disorder.txt', 'r')

    while True:
        line = fin.readline()
        if not line:
            break
        if line[0]=='>':
            enspID = line.split()[0][1:].split(':')[0]
            aseq = fin.readline().split()[0]
            alist = fin.readline().split(',')[:-1]
            disScores = []
            for item in alist:
                disScores.append(float(item))
            disorder[enspID] = disScores
    fin.close()
    
    return disorder



def main():
    start = time()

    fin = open(sys.argv[1], 'r')
    chrNum = sys.argv[2]
    #gina
    codinginfoPath = sys.argv[4]

    if chrNum!='X' and chrNum!='Y':
        chrNum = int(chrNum)
##    print 'Reading all transcripts...'
    enstInfoDict = readExons(chrNum,codinginfoPath) #read all transcripts
##    print '.....................Done!\n'

#    print '-->Reading protein domains...'
    ProteinDomains = getDomains(codinginfoPath)
#    print '........................Done!\n'
    
    #conservation scores
 #   print '-->Reading conservation scores of all exons...'
    DNAConvervationScores = readConservationScores(chrNum,codinginfoPath)
 #   print '.........................................Done!\n'

    #disorder scores and secondary structure
 #   print '-->Reading disorder and secondary structure info'
    disorder = readDisorderAA(chrNum,codinginfoPath)
    
    end1 = time()
        
    outputFileName = sys.argv[3]
    
    affectedEnstList = []
    while True:
        line = fin.readline()
        if not line:
            break
        if len(line.split()) == 0: #finished reading the transcripts of one gene
            #do processing
            if len(affectedEnstList)==0:
                continue
 #           print '\n\n---->Processing '+ affectedEnstList[0][1]

            affectedRegions = extractFeatures(affectedEnstList, chrNum, enstInfoDict, ProteinDomains, \
                                              DNAConvervationScores, disorder, outputFileName)
            
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
            
            
##            affectedEnstList.append((indel, idENSG, origProtein, indelPos, origProSize, overlapPro1, overlapPro2,\
##                                    numFuncDomain1, numFuncDomain2,  numAllFuncDomains, numAffectedTrans,\
##                                     numUnaffectedTrans, numAllTrans, nmdValue, idENST, numExons, \
##                                     disExonBoundary, idENSP, protein1))
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
                                
            
            

    fin.close()
    end2 = time()
 

main()
                  
    
    
