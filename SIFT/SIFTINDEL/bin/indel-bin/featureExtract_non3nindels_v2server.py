'''
Copyright 2011 Jing Hu. All rights reserved.

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

#############################################################################
# This program extract features for each indel
#
# Features to extract includes"
#   0) Indel input (chr, location, strand, bases) - directly from the input
#   1) ENSG id -< What gene it's in
#   2) ENST id <- what transcript it's in
#   3) ENSP id <- what protein it's in
#   4) Exon # <- what exon number it's in, relative to transcript (including UTR)
#   5) total # of coding exons in transcript (this will let you see
#       if it is the last or 2nd to last exon). <- Including UTR
#   6) coding exon #  <- Not including exons with UTRs
#   7) total # of coding exons <- Not including exons with UTRs
#   8) chromosome coordinates of exon indel lies in (chr: start-end)\
#   9) Strand
#  10) old bases and new bases at the indel potisiont
#  11) old coding sequence
#  12) new coding sequence after indel
#  13) old codon and new codon
#  14) original protein sequence
#  15) What amino acid position indel starts
#  16) total # of amino acids in corresponding normal protein
#  17) newly translated protein sequence (sequence 1). If there is indel, we consider two possible protiens: protein 1 and protein 2
#  (a new protein starts from the new in-frame starting codon after the stopping codon of protein 1)
#  18) overlap of protein 1 with original protein
#  19) newly translated protein 2 (if no such protein, then empty)
#  20) overlap of protein 2 with original protein.
#  21) Flanking sequence
#  22) For the gene, how many of them are affected, how many are not, and what is total number of transcripts per gene.
#  23) Whether it is mRNA decay.
#  24) Influenced exon coordinates
#  25) Wheater this indel can be found at ensembl database. If yes, what is the Ensembl variation ID.
#  
# 
# Usage: python FeatureEtract_12_14.py indelFile chrNum outputfile
# Example: python FeatureEtract_12_14.py test.txt 1 out.txt
#
# Copyright: Jing Hu & Pauline Ng
# Franklin & Marshall College, Singapore Genome Institute
# Date: May 9, 2011
# Version: 1.0
# Log:
#   1) Removed the calculation of overlapped functional regions
#
##############################################################################

import sys
import os
from math import *
from time import *

###########################################################################
#
#Codon Table mapping
#
###########################################################################
codon = {}
codon['TTT'] = 'F'; codon['TTC'] = 'F'; codon['TTA'] = 'L'; codon['TTG'] = 'L';
codon['TCT'] = 'S'; codon['TCC'] = 'S'; codon['TCA'] = 'S'; codon['TCG'] = 'S';
codon['TAT'] = 'Y'; codon['TAC'] = 'Y'; codon['TAA'] = '*'; codon['TAG'] = '*';
codon['TGT'] = 'C'; codon['TGC'] = 'C'; codon['TGA'] = '*'; codon['TGG'] = 'W';

codon['CTT'] = 'L'; codon['CTC'] = 'L'; codon['CTA'] = 'L'; codon['CTG'] = 'L';
codon['CCT'] = 'P'; codon['CCC'] = 'P'; codon['CCA'] = 'P'; codon['CCG'] = 'P';
codon['CAT'] = 'H'; codon['CAC'] = 'H'; codon['CAA'] = 'Q'; codon['CAG'] = 'Q';
codon['CGT'] = 'R'; codon['CGC'] = 'R'; codon['CGA'] = 'R'; codon['CGG'] = 'R';

codon['ATT'] = 'I'; codon['ATC'] = 'I'; codon['ATA'] = 'I'; codon['ATG'] = 'M';
codon['ACT'] = 'T'; codon['ACC'] = 'T'; codon['ACA'] = 'T'; codon['ACG'] = 'T';
codon['AAT'] = 'N'; codon['AAC'] = 'N'; codon['AAA'] = 'K'; codon['AAG'] = 'K';
codon['AGT'] = 'S'; codon['AGC'] = 'S'; codon['AGA'] = 'R'; codon['AGG'] = 'R';

codon['GTT'] = 'V'; codon['GTC'] = 'V'; codon['GTA'] = 'V'; codon['GTG'] = 'V';
codon['GCT'] = 'A'; codon['GCC'] = 'A'; codon['GCA'] = 'A'; codon['GCG'] = 'A';
codon['GAT'] = 'D'; codon['GAC'] = 'D'; codon['GAA'] = 'E'; codon['GAG'] = 'E';
codon['GGT'] = 'G'; codon['GGC'] = 'G'; codon['GGA'] = 'G'; codon['GGG'] = 'G';


###########################################################################
#
#Find the candidate transcript and to see which exon it is located
#if a positive integer is returned, means found the corresponding exon
#in the transcript;
#if a 0 is returned, means found the transcript, but not on the exon region;
#if -1 is returend, means not the right transcript.
#
###########################################################################
def findTranscript(exons, ranks, startPos, endPos):
    #we are using space-based coordinates.
    #0A1C2G3T4, therefore (1,3) represents CG, (1,1) represents position 1
    if startPos!=endPos:
        startPos = startPos + 1
    ranks = ranks.split(';')
    exonStarts = exons[0].split(';')
    exonEnds = exons[1].split(';')
    for i in range(len(exonStarts)):
        if startPos>=int(exonStarts[i]) and startPos<=int(exonEnds[i]) \
           or endPos>=int(exonStarts[i]) and endPos<=int(exonEnds[i]): #modifed May 14, 2013. Changed 3rd "and" to "or"
            return ranks[i], (exonStarts[i], exonEnds[i])  #found the exon
    minIndex = int(exonStarts[0])
    maxIndex = int(exonStarts[0])
    for i in range(len(exonStarts)):
        if min((int(exonStarts[i]), int(exonEnds[i])))<minIndex:
            minIndex = min((int(exonStarts[i]), int(exonEnds[i])))
        if max((int(exonStarts[i]), int(exonEnds[i])))<minIndex:
            maxIndex = max((int(exonStarts[i]), int(exonEnds[i])))

##    if id_ENST == 'ENST00000368138':
##        print id_ENST
##        print exonStarts
##        print exonEnds
##        print startPos
##        print endPos
##        print minIndex
##        print maxIndex
    
    if (startPos>=minIndex and startPos<=maxIndex) \
       or (endPos>=minIndex and endPos<=maxIndex):
        return 0, None   #indel not in the exon region, but in the transcript
    else:
        return -1, None  #indel not in the transcript



        
###########################################################################
#
# This function finds the information about coding exons:
# such as how many coding exons, which coding exon indel happens,
# and the exact starting and ending position of each coding exon (not
# including UTRs.
#
###########################################################################
def findCodingExonInfo(startPos, endPos, UTR5_Starts, UTR5_Ends, UTR3_Starts, UTR3_Ends,\
                       exons, ranks, strand):
    exonStarts = exons[0].split(';')
    exonEnds = exons[1].split(';')
    ranks = ranks.split(';')
    #get exons in order (including UTRs) 
    exons = [None]*(len(ranks)+1)
    for i in range(len(ranks)):
        ind = int(ranks[i])
        exons[ind] = [int(exonStarts[i]), int(exonEnds[i])]
##    print 'original exons'
##    print exons
    #print UTR5_Starts
    #print len(UTR5_Starts)

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

##    print 'coding exons'
##    print codingExons
    #get the coding exon # where indel sits
    whichExon = 0
    #the distance to exon boundary. If the indel is in the middle exon, then the distance is the smaller of the distances
    #to both boundaries. If the indel is in the first or last exon, then the distance is the distance to the non-UTR boundary
    disExonBoundary = None   
    if startPos!=endPos:
        startPos+=1
    for i in range(len(codingExons)):
        exon = codingExons[i]
	#modified at May 14, 2013 by Jing HU. changed 2nd "and" to "or"
        if startPos>=exon[0] and startPos<=exon[1] or endPos>=exon[0] \
           and endPos<=exon[1]:
            whichExon = i + 1
            #calculate the distance to boundary
            if i==0:
                disExonBoundary = exon[1]-endPos
            elif i==len(codingExons)-1:
                disExonBoundary = startPos-exon[0]
            else:
                disExonBoundary = min([exon[1]-endPos, startPos-exon[0]])
            break

##    print 'number of UTR5 ' +str(numUTR5)
##    print 'number of UTR3 '+str(numUTR3)
    #returns coding Exons, number of coding exons, which coding exon it is located,
    #exons in list format, and number of UTR5s and UTR3s
    return codingExons, len(codingExons), whichExon, disExonBoundary, exons, (numUTR5, numUTR3)

##########################################################################
#
# if indel happens at the beginging of the gene (i.e., in the first 75 DNA bases),
# , or within < 5th percentile of length,  then look for the next in-frame ATG,
# that there is only one protein 1.
# otherwise, start translating from the very begining.
#
###########################################################################
def translateStartPos(ind1,ind2, base, seqLength):
    if ind1<75 or float(ind1)/seqLength <=0.05:
        if ind1==ind2:  #insertion
            if ind2%3==0:
                return ind2+len(base)
            else:
                return ind2+len(base) + (3 - ind2 % 3)
        elif ind2>ind1 and base == '/':  #deletion
            if ind2 % 3 ==0:
                return ind1
            else:
                return ind1 + (3 - ind2%3)                                  
        else:                           #block substition, including true indels and SNPs (e.g., 21,10,12,ATGT)
            if ind2 % 3 == 0:
                return ind2 + (len(base) - (ind2-ind1))
            else:
                return (ind2//3+1)*3 + (len(base) - (ind2-ind1))
    else:
        return 0

    
            
    
###########################################################################
#
# This function translate a new CD sequence into proteins
#
###########################################################################
def translation(ind1, ind2, base, oldCDSeq, newCDSeq, UTR3string,functionDomains,  nmrOffset):
    extendedSequence = newCDSeq + UTR3string
##    print 'extended sequence'
##    print extendedSequence
    proteinSeq = ''
    proteinSeq2 = ''
    NMD_Flag = ''
    #start translating protein 1
    i = translateStartPos(ind1, ind2, base, len(oldCDSeq)) #get the starting position in the coding sequence. 
##    print '>>>new starting position->'+str(i)
    if i==0:
        aa = codon[extendedSequence[i:i+3]]
    else:
        #find start codon ATG or Met/M
        aa = codon[extendedSequence[i:i+3]]
        while aa!='M':
            i = i + 3
            if i+3>len(newCDSeq):
##                print 'Protien 1: No starting codon'
                return None, None
            
            aa = codon[extendedSequence[i:i+3]]

    newStartCodonPos = i  #if indel happens at first 75 bases, then we skip for next in-frame
                          #starting codon. This is the position of the new starting codon
    
    proteinSeq += aa
    i = i + 3
    #translate until stop codon and pass indel position
    boundaryFlag = False #test to see if all coding sequence has been consumed
    if i+3<=len(newCDSeq):
        aa = codon[extendedSequence[i:i+3]]
        while aa!='*':
            proteinSeq += aa
            i = i + 3
            if i+3>=len(newCDSeq)-3:
                boundaryFlag = True  #pass the CDs boundary
            if i+3>len(extendedSequence):
##                print 'Protien 1: No stopping codon'
                break    
            aa = codon[extendedSequence[i:i+3]]
        if aa=='*':
            proteinSeq += aa
            stopCodingPosition = i+3  #position of stopping codon on the coding sequence for protein 1
            if nmrOffset==None:
                NMD_Flag = 'NONE'
            elif nmrOffset>0:
                nonNMD_Range = [len(newCDSeq)-nmrOffset, len(newCDSeq)]
                if stopCodingPosition>=nonNMD_Range[0]:
                    NMD_Flag = 'NONE'
                else:
                    NMD_Flag = 'NMD'
        else:
            NMD_Flag = 'nonstop mRNA decay'
            
                  
    
    #if there is still CDs left, continue to see if protein 2 exsits
    if not boundaryFlag:   
        i = i + 3
        decrease = ind2 - ind1
        if base!='/':
            increase = len(base)
        else:
            increase = 0
        change = increase - decrease
        change = change % 3
        i = i + change #update i so that it can find next starting codon in-frame

        #try to find start codon ATG or Met/M
        hasProtein2 = True
        if i+3<=len(newCDSeq):
            aa = codon[newCDSeq[i:i+3]]
            while aa!='M':
                if i+3>len(newCDSeq)-3:
                    hasProtein2 = False
                    break
                i = i + 3
                aa = codon[extendedSequence[i:i+3]]
            if hasProtein2:
                proteinSeq2 += aa
                i = i + 3
                #translate until stop codon (It should find the last stop codon since it is in frame)

                if i+3<=len(newCDSeq):
                    aa = codon[newCDSeq[i:i+3]]
                    while aa!='*':
                        proteinSeq2 += aa
                        i = i + 3
                        if i+3 > len(newCDSeq):
                            #print 'Protein 2: No stopping codon. Something wrong!'
                            proteinSeq2=''
                            break
                        aa = codon[newCDSeq[i:i+3]]
                    if aa=='*':
                        proteinSeq2 += aa

    
        
    
    
    #translate the original exon into original protein    
    i = 0
    originalProteinSeq = ''
    aa = codon[oldCDSeq[i:i+3]]
    while aa!='M':
        i = i + 3
        if i+3>len(oldCDSeq):
##            print 'Original protein: No starting codon'
            return None, None
        aa = codon[oldCDSeq[i:i+3]]
       
    i = i + 3
    originalProteinSeq += aa
    #translate until stop codon and pass indel position
    #boundaryFlag = False #test to see if all coding sequence has been consumed
    while i+3<=len(oldCDSeq):
        aa = codon[oldCDSeq[i:i+3]]
        originalProteinSeq += aa
        i = i+3
        
##    if i+3<=len(oldCDSeq):
##        aa = codon[oldCDSeq[i:i+3]]
##        while aa!='*':
##            originalProteinSeq += aa
##            i = i + 3
##            if i+3>len(oldCDSeq):
##                print 'Original protein: No stopping codon'
##                #originalProteinSeq = ''?
##                break
##        
##            aa = codon[oldCDSeq[i:i+3]]
##        if aa=='*':
##            originalProteinSeq += aa

    indelPositionOnProtein = ind1 // 3 + 1    #indel position of the original protein
    #originalProteinSize = len(oldCDSeq) // 3  #position of indel on the original protein
    originalProteinSize = len(originalProteinSeq)  #size of the original protein


    #get the overlap residues  review comments: need to consider the new in-frame starting codon
    overlaps = overlapResidues(originalProteinSeq, proteinSeq, proteinSeq2, oldCDSeq,newCDSeq)
##    print overlaps
##    print '++++++'

    if overlaps[0]==None:
        overlap1 = 0
    else:
        overlap1 = 0
        for item in overlaps[0]:
            if item==None:
                continue
            overlap1 += item[1]-item[0]+1

    if overlaps[1]==None:
        overlap2 = 0
    else:
        overlap2 = 0
        for item in overlaps[1]:
            if item==None:
                continue
            overlap2 += item[1]-item[0]+1

##    print overlap1
##    print overlap2
    
    DomainsInfo = getNumDomains(functionDomains, overlaps)
    
    return [originalProteinSeq, proteinSeq, proteinSeq2, indelPositionOnProtein,\
            originalProteinSize, overlap1, overlap2,  NMD_Flag, DomainsInfo, len(functionDomains)], newStartCodonPos

###########################################################################
#This function get the number of function domains in each protein.
#
###########################################################################
def getNumDomains(domains, overlaps):
    
    #print domains
    adict={}
    for domain in domains:
        adict[str(domain[0])+'-'+str(domain[1])+'-'+domain[2]] = False
        
    if overlaps[0]==None or overlaps[0]==[None]:
        numDomain1 = 0
    else:
        numDomain1 = 0
        for domain in domains:
            for region in overlaps[0]:
                if domain[0]>=region[0] and domain[1]<=region[1]:
                    adict[str(domain[0])+'-'+str(domain[1])+'-'+domain[2]] = True
                    numDomain1 += 1
##            if domain[0]>=overlaps[0][0] and domain[1]<=overlaps[0][1]:
##                adict[str(domain[0])+'-'+str(domain[1])+'-'+domain[2]] = True
##                numDomain1 += 1

##    lostDomains1 = []  #lost domains after indel, only count domains kept by protein 1
##    for key in adict:
##        if adict[key]==False:
##            key = key.split('-')
##            lostDomains1.append((key[0]+'-'+key[1], key[2]))
##            
    if overlaps[1]==None or overlaps[1]==[None]:
        numDomain2 = 0
    else:
        numDomain2 = 0
        for domain in domains:
            for region in overlaps[1]:
                if domain[0]>=region[0] and domain[1]<=region[1]:
                    adict[str(domain[0])+'-'+str(domain[1])+'-'+domain[2]] = True
                    numDomain2 += 1
####                    
####            if domain[0]>=overlaps[1][0] and domain[1]<=overlaps[1][1]:
####                adict[str(domain[0])+'-'+str(domain[1])+'-'+domain[2]] = True
####                numDomain2 += 1
##
##    #create a list of lost domains
##    lostDomains2 = []  #lost domains after indels, count domains kept by both protein 1 and 2
##    for key in adict:
##        if adict[key]==False:
##            key = key.split('-')
##            lostDomains2.append((key[0]+'-'+key[1], key[2]))
##            

    return (numDomain1, numDomain2)
    



###########################################################################
# This function find the overlap between new proteins with original protein.
#
#
#!##########################################################################
def overlapResidues(originalProteinSeq, protein1, protein2, oldCDSeq,newCDSeq ):
    #non-3n indels or subsitutions
    overlapRegion1=None
    overlapRegion2=None
    if (len(newCDSeq)-len(oldCDSeq)) % 3 !=0:    
        #find the overlap of protein 1
        if len(protein1)==0:
            overlapRegion1=None
        else:
            mPos = []
            for i in range(len(originalProteinSeq)):
                #if originalProteinSeq[i]=='M':
                if originalProteinSeq[i]==protein1[0]:
                    mPos.append(i)
            maxOverLap = 0
            for s in mPos:
                tempSeq = originalProteinSeq[s:]
                i = 0
                while i<min(len(tempSeq), len(protein1)):
                    if tempSeq[i]!=protein1[i]:
                        break
                    i+=1
                if i>maxOverLap:
                    maxOverLap = i
                    overlapRegion1 = (s+1, s+i)
        #find the overlap of protein 2
        if len(protein2)==0:
            overlapRegion2 = None
        else:
            i = -1
            while abs(i)<=min(len(originalProteinSeq), len(protein2)):
                if originalProteinSeq[i]!=protein2[i]:
                    break
                i = i -1
            overlapRegion2 = (len(originalProteinSeq)+i+2, len(originalProteinSeq))

##        print 'overlapping'
##        print [(overlapRegion1), (overlapRegion2)]

        return [[overlapRegion1], [overlapRegion2]]
    else:  #3n subsitutions or indels
    #in this case, there is only protein 1, and protein 2 will be empty
        if len(protein1)==0:
            overlapRegion1=None
        else:
            mPos = []
            for i in range(len(originalProteinSeq)):
                #if originalProteinSeq[i]=='M':
                if originalProteinSeq[i]==protein1[0]:
                    mPos.append(i)
            maxOverLap = 0
            for s in mPos:
                #overlapping regions before mutation
                tempSeq = originalProteinSeq[s:]
                i = 0
                while i<min(len(tempSeq), len(protein1)):
                    if tempSeq[i]!=protein1[i]:
                        break
                    i+=1
                maxOverLap1 = i
                overlapRegionA = (s+1, s+i)
            
                #overlapping regions after mutation
                halfSeq = tempSeq[i:]   #sequence left for original sequence
                halfProtein = protein1[1:] #remaining of protein1

                halfSeq = list(halfSeq)
                halfSeq.reverse()
                halfProtein = list(halfProtein)
                halfProtein.reverse()
                i = 0
                while i<min(len(halfSeq), len(halfProtein)):
                    if halfSeq[i]!=halfProtein[i]:
                        break
                    i+=1
                maxOverLap2 = i
                #print '\n\nYYYYY'
                #print maxOverLap2
                overlapRegionB =(s + len(tempSeq)-i+1, s + len(tempSeq) )
                if (maxOverLap1+maxOverLap2) > maxOverLap:
                    maxOverLap = maxOverLap1+maxOverLap2
                    overlapRegion1 = (overlapRegionA,overlapRegionB)
                #print maxOverLap
                #print overlapRegions1
                

##            print 'overlapping'
##            print [overlapRegion1, None]
            return [overlapRegion1, None]



###########################################################################
#
# This function finds out how many coding bases are not in non-NMD regions:
# which is the last exon or <50 bases in the second to last exon
#
###########################################################################
def nmdRange(numUTR3,allExons):
    #allExons looks like [None, [15986364, 15988217]], so if length is 2, then only 1 exon
    if len(allExons)==2:
        return None
    offSet = 0
    #if there is no UTR3, then the last exon and the last
    #50 nucleotides in the second to last exon are all coding sequence
    if numUTR3 == 0: 
        offSet += max(allExons[-1])-min(allExons[-1])+1
        offSet += min(50, max(allExons[-2]) - min(allExons[-2]) + 1)
    else:
        #UTR3s in the last exon (case 1): only 1 UTR
        if numUTR3<=max(allExons[-1])-min(allExons[-1]) + 1:
            offSet += max(allExons[-1])-min(allExons[-1]) + 1 - numUTR3
            offSet += min(50, max(allExons[-2])-min(allExons[-2]) + 1)
        #UTR3s cover the whole last exon (>= 2 UTRs): case 2 and 3
        else:
            for i in range(-1, -1*len(allExons),-1):
                exon = allExons[i]
                numUTR3 = numUTR3 - (max(exon) - min(exon) + 1)
                if numUTR3 < 0:
                    curInd = i
                    break
            if numUTR3>0:
                return None
            numUTRatThisExon = (max(allExons[curInd]) - min(allExons[curInd]) + 1) + numUTR3 #num of UTRs which does not cover full exon
            if numUTRatThisExon<50:
                offSet = min(50, max(allExons[curInd]) - min(allExons[curInd]) + 1) - numUTRatThisExon
            else:
                offset = None
                
                
            
##        elif numUTR3 < max(allExons[-1])-min(allExons[-1]) + 1 + 50:
##            offSet += max(allExons[-1])-min(allExons[-1]) + 1 + 50  -numUTR3
##        #if UTR3 cover the last exon and >=50 bases in the second to last exon
##        else:
##            offSet += -100
    return offSet


    
###########################################################################
#
# This function finds if an indel has record in ensembl database.
# If it has, it returns a tuple of (True, rsID);
# otherwise, it returns a tuple of (False, None)
#
###########################################################################
def isInEnsembl(newCDSeq, chrNum,startPos, endPos, bases,ensemblDict):
    indel = str(chrNum)+','+str(startPos)+','+str(endPos)+','+bases
    for key in ensemblDict:
        if newCDSeq in ensemblDict[key] or indel in ensemblDict[key]:
            return (True, key)
    return (False, None)


###########################################################################
#
# This function finds the information about translated protein information:
#
###########################################################################
def findProteinInfo(chrNum, CDs, startPos, endPos, bases, codingExons, allExons, \
                    seq, strand, numUTRs, functionDomains, ensemblDict):
##    print 'exon seq'
##    print seq
    numUTR5 = numUTRs[0]
    numUTR3 = numUTRs[1]


    #updated on March 2013 by Jing HU
    CDs_Starts = CDs[0]
    CDs_Ends = CDs[1]

    numCDs = 0
    CDs_StartsIndex = []
    if CDs_Starts!='' and CDs_Ends!='':
        for ind in CDs_Starts.split(';'):
            CDs_StartsIndex.append(int(ind))
        CDs_EndsIndex = []
        for ind in CDs_Ends.split(';'):
            CDs_EndsIndex.append(int(ind))
        
        numCDs = max(CDs_EndsIndex)-min(CDs_StartsIndex) + 1


    if numUTR5!=0 and numUTR3!=0:
        if numCDs!=0:
            codingSeq = seq[numUTR5:numUTR5+numCDs]
        else:
            codingSeq = seq[numUTR5:-1*numUTR3] #the original CD sequence
    elif numUTR5!=0 and numUTR3==0:
        if numCDs!=0:
            codingSeq = seq[numUTR5:numUTR5+numCDs]
        else:
            codingSeq = seq[numUTR5:] #the original CD sequence
    elif numUTR5==0 and numUTR3!=0:
        if numCDs!=0:
            codingSeq = seq[:numCDs]
        else:
            codingSeq = seq[:-1*numUTR3] #the original CD sequence
    else:
        codingSeq = seq

            
    
##    if numUTR5!=0 and numUTR3!=0:
##        codingSeq = seq[numUTR5:-1*numUTR3] #the original CD sequence
##    elif numUTR5!=0 and numUTR3==0:
##        codingSeq = seq[numUTR5:] #the original CD sequence
##    elif numUTR5==0 and numUTR3!=0:
##        codingSeq = seq[:-1*numUTR3] #the original CD sequence
##    else:
##        codingSeq = seq

    
        
    numAA = len(codingSeq) // 3
    # if it is positive chain
    newCDSeq = ''       #new CD sequence after indel/substitution
    cdsPos = 0
    originalAASeq = ''
    newBases = bases
    if strand == '1':   #transcript is on positive chain
        for exon in codingExons:
            if startPos >=exon[0] and startPos<=exon[1]:
                cdsPos += startPos - exon[0] + 1
                break
            else:
                cdsPos += exon[1] - exon[0] + 1
        if startPos < endPos and bases=='/':       #deletion
            ind1 = cdsPos
            ind2 = (endPos-startPos) + ind1
            
            newCDSeq = codingSeq[0:ind1]+codingSeq[ind2:]
        elif startPos<endPos and bases!='/':        #indels
            ind1 = cdsPos
            ind2 = (endPos-startPos) + ind1
            newCDSeq = codingSeq[0:ind1] + bases + codingSeq[ind2:]
        else:                                       #insertion
            ind1 = cdsPos
            ind2 = cdsPos
            newCDSeq = codingSeq[0:ind1] + bases + codingSeq[ind2:]

        oldBases = codingSeq[ind1:ind2]
        
            

    else:   #transcript is on negative chain
        for exon in codingExons:
            if endPos>=exon[0] and endPos<=exon[1]:
                #print str(exon[0])+' --- '+str(exon[1])
                #cdsPos += exon[1] - endPos + 1
                cdsPos += exon[1] - endPos
                break
            else:
                cdsPos += exon[1] - exon[0] + 1
        if startPos < endPos and bases=='/':       #deletion
            ind1 = cdsPos
            ind2 = (endPos-startPos) + ind1
            newCDSeq = codingSeq[0:ind1]+codingSeq[ind2:]
        elif startPos<endPos and bases!='/':        #indels
            ind1 = cdsPos
            ind2 = (endPos-startPos) + ind1
            newBases = getComplements(bases)
            newCDSeq = codingSeq[0:ind1] + newBases + codingSeq[ind2:]
        else:                                       #insertion
            ind1 = cdsPos
            ind2 = cdsPos
            #print ind1
            #print ind2
            newBases = getComplements(bases)
            newCDSeq = codingSeq[0:ind1] + newBases + codingSeq[ind2:]

        oldBases = codingSeq[ind1:ind2]
        
##    print 'ind1'
##    print ind1
##    
##    print 'coding sequence'
##    print codingSeq
##    print 'new Coding sequence'
##    print newCDSeq
##    
    flankingSeq = getFlankSeq(codingSeq, ind1, ind2, bases)
##    print 'Flanking sequence:'
##    print flankingSeq
    
    
    #to see if the indel is in ensembl database by using newCDSeq information
    #inEnsemblFlag = isInEnsembl(newCDSeq, ensemblDict)
    inEnsemblInfo = isInEnsembl(newCDSeq, chrNum,startPos, endPos, bases, ensemblDict)
##    if inEnsemblInfo[0]:
##        print 'The indel is in Ensembl variation database'
##    else:
##        print 'The indel is NOT in Ensembl variation database'

    nmrOffset = nmdRange(numUTR3, allExons) #offset of non-NMR
    if numUTR3==0:
        proteinInfo, newStartCodonPos = translation(ind1, ind2, bases, codingSeq.upper(), newCDSeq.upper(), '',  functionDomains, nmrOffset)
        #return proteinInfo, flankingSeq, conservationSTA, inEnsemblInfo
    else:
        proteinInfo, newStartCodonPos = translation(ind1, ind2, bases, codingSeq.upper(), newCDSeq.upper(), seq[-1*numUTR3:].upper(), functionDomains, nmrOffset)

    #avgConservationScore = calAverageConservation(startPos, endPos, strand, codingExons, conservationList)
##    conservationSTA = calConservationSTA(startPos, endPos, strand, codingExons, conservationList, newStartCodonPos, ind1, ind2)
##    print 'Avg conservation score after indel: '+str(conservationSTA[0])
##    print 'Avg highly conservation score after indel: '+str(conservationSTA[1])
##    print 'Lost highly conserved positions: '+str(conservationSTA[2])+'; All conserved positions: '+str(conservationSTA[3])+'\n'
##    print 'Percentage of lost high conserved positions: '+str(conservationSTA[4])+'\n'
        
    #return proteinInfo, flankingSeq, conservationSTA, inEnsemblInfo
    #return proteinInfo, flankingSeq, inEnsemblInfo
    return proteinInfo, flankingSeq, inEnsemblInfo, codingSeq, newCDSeq, oldBases, newBases

##############################################################################
# This function get the flanking sequences in indel position
#
###############################################################################  
def getFlankSeq(codingSeq, ind1, ind2, bases):
    flankingSeq = ''
    wSize = 5
    if ind1<wSize:
        flankingSeq += codingSeq[:ind1]
    else:
        flankingSeq += codingSeq[ind1-wSize:ind1]
##    if ind1<ind2:   #deletion and block subsitution
##        toDelete = codingSeq[ind1:ind2].lower()
##        flankingSeq += toDelete
##    else:           #insertion
##        flankingSeq += bases.lower()

    if ind1<ind2:   #deletion and block subsitution
        toDelete = '['+codingSeq[ind1:ind2]+'/'
        if bases=='/':
            toDelete+='-]'
        else:
            toDelete+=bases+']'
        flankingSeq += toDelete
    else:           #insertion
        flankingSeq += '[-/'+bases+']'
        
    if len(codingSeq)-ind2<wSize:
        flankingSeq += codingSeq[ind2:]
    else:
        flankingSeq += codingSeq[ind2:ind2+wSize]

    return flankingSeq
        
    
            
##############################################################################
# This function returns the complement of a subsequence on the positive chain
# For example: 'ATATCG' on positive chain will be 'CGATAT' on negative chain
#
###############################################################################
def getComplements(bases):
    complement = {}
    complement['A'] = 'T'
    complement['T'] = 'A'
    complement['C'] = 'G'
    complement['G'] = 'C'
    newBases = ''
    for i in range(-1, -len(bases)-1, -1):
        newBases += complement[bases[i].upper()]
    return newBases
            
#################################################
# This function counts the number of transcripts
# per gene.
#
#################################################
def numTranscriptsPerGene(enstInfoList, id_ENSG):
    numTranscripts = 0

    for enst in enstInfoList:
        info = enst.split()
        title = info[0]        
        title = title.split('|')
        if title[0]==id_ENSG:
                numTranscripts+=1
    return numTranscripts

##############################################
# This function finds the following info:
#   ENSG ID, ENST ID, ENSP ID, exon number
#
##############################################
def extractIndelInfo(enstInfoDict, chrNum, startPos, endPos, bases,  outputFileName, \
                      ProteinDomains, ensembleVariations, errorFile1,errorFile2):
    indelAtExonFlag = False
    numTranscripts = 0  #number of transcripts which contain the indel (intron or exon regions)
    firstEncounterFlag = True
    numTranscriptsAffected = 0 #number of transcripts which have exons affected
    #chromsomeList = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y]
##    chromsomeList = [21]
##    enstInfoList = []   # a list to keep all transcript information
##    #read all chromsome files to get information
##    for chrom in chromsomeList:
##        fin = open('chr'+str(chrom)+'new.txt', 'r')
##        lines = fin.read()
##        lines = lines.split('>')[1:]
##        for info in lines:
##            enstInfoList.append(info)
##        fin.close()

    #transcript information
    enstInfoList = enstInfoDict[chrNum]
    #DNA conservation information of a chromesome
    #conservationList = DNAConvervationScores[chrNum]
    #ensembl variation information of a chromesome
    ensemblDict = ensembleVariations[chrNum]
    
  
    #all features extracted for each indel at each transcript which has the indel at exon region
    IndelFeatureList = []
    # a collection of all transcripts (in exons) which contain the indel
    # (does not have to be on exons). For example, 
    transcriptExonList = [] 

    for enst in enstInfoList:
        info = enst.split()
        title = info[0]
        seq = ''
        for s in info[1:]:
            seq = seq + s
        
        title = title.split('|')
        
        
        
        id_ENSG = title[0]      #gene ID
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


        #candidate transcript (we must find if the indel happens at exon too)
        #if foundCandidateTranscript(leftInd, rightInd, startPos, endPos) and strand==title[12]:
        #exonNum tells which exon rank (e.g., 2) it is on.
        #theExon tells which exon (e.g., (11230019, 11230079)) it is on
        exonNum, theExon = findTranscript(exons, ranks, startPos, endPos)
        
        #if exonNum>=0:
        #    numTranscripts+=1
            
        #print exonNum
        #print id_ENST
        if exonNum>=0: # and strand==title[12]:
            #numTranscripts+=1

            if firstEncounterFlag:
                numTranscripts = numTranscriptsPerGene(enstInfoList, id_ENSG)
                firstEncounterFlag = False
            
##            print '\n-----Processing transcript : '+id_ENST
            
            #get coding Exons, number of coding exons            
            codingExons, numCodingExons, whichCodingExon,disExonBoundary, allExons, (numUTR5, numUTR3) = findCodingExonInfo(startPos, endPos, \
                    UTR5_Starts, UTR5_Ends, UTR3_Starts, UTR3_Ends,exons, ranks, strand)

            #append to the list to create a list of coding exons of all transcripts which contains the indel
            #indels do not have to be on coding exons
            transcriptExonList.append((id_ENST, codingExons, whichCodingExon))
            
            if whichCodingExon==0:
                continue
            
            transcriptVector = []
            indelAtExonFlag = True #found the valid exon which has the indel
            numTranscriptsAffected+=1

            #if not sequence data skip
            if seq=='Sequenceunavailable' or id_ENSP=='':
##                print 'Transcript: ' + title[1]+' : sequence unavailable'
                ferror2 = open(errorFile2, 'a')
                ferror2.write(str(chrNum)+','+str(startPos)+','+ str(endPos)+','+bases+'\n')
                ferror2.close()
                continue
            
##            print title[1]+' '+str(exonNum)
            #fout = open(outputFileName, 'a')
            #feature 0
            #transcriptVector.append(str(chrNum)+','+str(startPos)+','+ str(endPos)+','+strand+','+bases)
            transcriptVector.append(str(chrNum)+','+str(startPos)+','+ str(endPos)+','+bases)
            #feature 1-3
            transcriptVector.append(id_ENSG)
            transcriptVector.append(id_ENST)
            transcriptVector.append(id_ENSP)
            #feature 4: which exon
            transcriptVector.append(exonNum)
            #feature 5: total # of exons including UTRs
            transcriptVector.append(len(exons[0].split(';')))            
            #Feature 6, 7:  which coding exon # indel sits            
            transcriptVector.append(whichCodingExon)
            transcriptVector.append(numCodingExons)            
            #Feature 8: chromsome coordinates (chr:start-end)
            transcriptVector.append('chr:'+str(chrNum)+':'+str(startPos)+'-'+str(endPos))

            

            #New features: protein and indel features       
##            if transcriptVector[3] not in ProteinDomains:
##                theDomains = []
##            else:
##                theDomains = ProteinDomains[transcriptVector[3]]
            theDomains = []
##            proteinInfo, flankingSeq =  findProteinInfo(startPos, endPos, bases, codingExons,\
##                                                        allExons, seq, strand,(numUTR5, numUTR3), ProteinDomains[id_ENSP])
##            proteinInfo, flankingSeq, inEnsemblInfo =  findProteinInfo(chrNum, startPos, endPos, bases, codingExons,\
##                                                        allExons, seq, strand,(numUTR5, numUTR3), ensemblDict)

            proteinInfo, flankingSeq, inEnsemblInfo,  codingSeq, newCDSeq, oldBases, newBases =  findProteinInfo(chrNum, CDs, startPos, \
                                                                                        endPos, bases, codingExons,\
                                                        allExons, seq, strand,(numUTR5, numUTR3), theDomains, ensemblDict)

                
                        
            if proteinInfo!=None:
##                print 'new protein 1'
##                print proteinInfo[1]
##                print 'new protein 2'
##                print proteinInfo[2]
                #fout.write(proteinInfo[0]+'\n')
                #feature 9-18
                for item in proteinInfo:
                    transcriptVector.append(item)
                
                transcriptVector.append(flankingSeq) #19
                #transcriptVector.append(conservationSTA)
                #transcriptVector.append(inEnsemblFlag)
                transcriptVector.append(inEnsemblInfo) #20
            else:
                #print 'XXXX'
                ferror2 = open(errorFile2, 'a')
                ferror2.write(str(chrNum)+','+str(startPos)+','+str(endPos)+','+ bases+'\n')
                ferror2.close()
                for i in range(12):
                #for i in range(10):
                    transcriptVector.append(None)
            #transcriptVector.append(theExon)
            transcriptVector.append(codingExons[whichCodingExon-1])  #21
            transcriptVector.append(strand)  #22
            transcriptVector.append((oldBases, newBases))  #23: old and new Bases which are to be mutated
            transcriptVector.append(codingSeq) #24: original coding sequence
            transcriptVector.append(newCDSeq)  #25: new coding sequence after mutation
            oldCodon, newCodon = findMutatedCodon(codingSeq, newCDSeq)
            transcriptVector.append((oldCodon, newCodon))  #26: old codon and new codon
            transcriptVector.append(disExonBoundary) #27: distance to the boundray
                
            #print findProteinInfo(startPos, endPos, codingExons, seq)
            #fout.close()
            IndelFeatureList.append(transcriptVector)  
            

            

            

    
    featureDict = sortFeaturesByGene(IndelFeatureList, enstInfoList)    #a dictionary of feature vector list

    if  indelAtExonFlag:
       # ferror1 = open(errorFile1,'a')
       # ferror1.write(str(chrNum)+','+str(startPos)+','+str(endPos)+','+ bases+'\n')
    #else:
        
        for key in featureDict:
            value = featureDict[key]
            IndelFeatureList = value[:-2]
            numAffectedTranscripts = value[-2]
            numberAllTranscripts = value[-1]
            
            fout = open(outputFileName, 'a')
            for transcriptVector in IndelFeatureList:
                if  transcriptVector[9]==None:
                    continue
                fout.write('>>>')
                for item in transcriptVector[0:9]:
                    fout.write(str(item)+'\t')
                fout.write('\n')
                fout.write('Strand: '+str(transcriptVector[22])+'\n')
                fout.write('Old Bases: '+transcriptVector[23][0]+'\t')
                fout.write('New Bases: '+transcriptVector[23][1]+'\n')
                fout.write('Old Coding Sequence: '.ljust(25)+transcriptVector[24]+'\n')
                fout.write('New Coding Sequence: '.ljust(25)+transcriptVector[25]+'\n')
                if transcriptVector[26][0]!='':
                    fout.write('Old Codon: '+transcriptVector[26][0]+'\t')
                    fout.write('New Codon: '+transcriptVector[26][1]+'\n')

                if transcriptVector[9][-1]!='*': #No stopping codon or stopping codon in the middle
                    if '*' in transcriptVector[9]:
                        fout.write('Ensembl gene annotation error-unable to process: Stopping codon in the middle!\n')
                    else:
                        fout.write('Ensembl gene annotation error-unable to process: No stopping codon!\n')
                    continue
                elif transcriptVector[9]==None:
                    fout.write('Ensembl gene annotation error-unable to process: No starting codon!\n')
                    continue
                fout.write('Original Protein: '.ljust(20)+transcriptVector[9]+'\n')  #original protein sequence
                fout.write(str(transcriptVector[12])+'\t'+str(transcriptVector[13])+'\n')  #indel position, and original protein length
                fout.write('Protein 1: '.ljust(20)+transcriptVector[10]+'\t') #protein sequence 1
                if transcriptVector[10][-1]!='*':
                    fout.write('Protein 1 no stopping codon!')
                fout.write('\n')
                fout.write('Overlap: '+str(transcriptVector[14])+'\n') #overlap of protein1 with original protein
                fout.write('Protein 2: '.ljust(20)+transcriptVector[11]+'\n') #protein sequence 2
                fout.write('Overlap: '+str(transcriptVector[15])+'\n') #overlap of protein2 with original protein
                #ovelapping domains on protein 1 and protein 2 and total number of original domains
                fout.write(str(transcriptVector[17][0])+'\t'+str(transcriptVector[17][1])+'\t'+str(transcriptVector[18])+'\n')
                #lost domains due to indel (only protein 1 is counted)
##                fout.write('Lost domains (only protein 1 is counted): ')
##                for item in transcriptVector[16][2]:
##                    fout.write(item[0]+'\t'+item[1]+'\t')
##                fout.write('\n')
##                #lost domains due to indel (protein 1 and 2 are counted)
##                fout.write('Lost domains (Both protein 1 and 2 are counted): ')
##                for item in transcriptVector[16][3]:
##                    fout.write(item[0]+'\t'+item[1]+'\t')
##                fout.write('\n')
                fout.write('Flaking sequence: '+transcriptVector[19]+'\n')  #flanking sequence
                fout.write('Affected: '+str(numAffectedTranscripts)+\
                           '\tNot Affected: '+str(numberAllTranscripts-numAffectedTranscripts)+\
                           '\t'+'Total TranScripts: '+str(numberAllTranscripts)+'\n')
                #indication of NMD (nonsense-mediated decay)
                if transcriptVector[16]!='':
                    fout.write('mRNA decay: '+ transcriptVector[16]+'\n')
                else:
                    fout.write('mRNA decay: NA\n')
                #fout.write('Influenced exon coordinates: '+str(transcriptVector[21][0])+\
                #           '-'+str(transcriptVector[21][1])+'\n')
                fout.write('Influenced exon coordinates: '+str(transcriptVector[21][0])+\
                           '-'+str(transcriptVector[21][1])+'\n')
                fout.write('Distance of the indel to the boundray of the exon: '+str(transcriptVector[27])+'\n')
            
                #print '$$$$$$$$$$$$$$$$$$$$$$$$$$'
                indelStartPos = int(transcriptVector[0].split(',')[1])
                indelEndPos = int(transcriptVector[0].split(',')[2])
                fout.write('Found in Ensembl variation database: '+str(transcriptVector[20][0])+'\t')
                fout.write('Ensembl Variation ID: '+str(transcriptVector[20][1])+'\n')

            fout.write('\n')
    

            fout.close()


def findMutatedCodon(oldSeq, newSeq):
    '''
    This function compares the old coding sequence with new coding sequence.
    If there are the same size, then find the original codon with new codon.
    Otherwise, just return None.
    '''
    oldCodon = ''
    newCodon = ''
    if len(oldSeq) != len(newSeq):
        return oldCodon, newCodon

    i = 0
    while i<=len(oldSeq)-3:
        if oldSeq[i:i+3]!=newSeq[i:i+3]:
            oldCodon+= oldSeq[i:i+3]
            newCodon+= newSeq[i:i+3]
        i+=3

    return oldCodon, newCodon
#############################################################################
#
# The function arrange transcripts by their gene names.
# It returns a dictionary.
# Key is the geneID, and the value is the list of features of each transcript,
# plus the number of affected transcripts and all transcripts per gene.
# Example: featureDict[ENSG00000243989] = [featurelist1, featurelist2,..,featurelistN, numAffected, totalTranscripts]
#
#############################################################################
def sortFeaturesByGene(IndelFeatureList, enstInfoList):
    featureDict = dict()
    for item in IndelFeatureList:
        id_ENSG = item[1]
        if id_ENSG not in featureDict:
            featureDict[id_ENSG] = [item]
        else:
            featureDict[id_ENSG].append(item)

    for id_ENSG in featureDict:
        featureDict[id_ENSG].append(len(featureDict[id_ENSG]))    #number of affected transcripts
        numTranscripts = numTranscriptsPerGene(enstInfoList, id_ENSG)
        featureDict[id_ENSG].append(numTranscripts)             #number of all transcripts

    return featureDict
              


#######################################################
# This function read all the transcript information
# of all chromsomes.
#
#######################################################
def readExons(chrNum, codinginfoPath, buildVersion):
    #chromsomeList = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']
    chromsomeList = [chrNum]
    enstInfoDict = {}
    for chrom in chromsomeList:
        enstInfoList = []
        #fin = open('/mnt1/indelFiles_V66/chromosomes_V66_update/chr'+str(chrom)+'_valid_2013_1_17.txt', 'r')
        if buildVersion == "hs37":
            fin = open(codinginfoPath + '/indelFiles_V66/chromosomes_V66_update/chr'+str(chrom)+'_valid_2013_1_17.txt', 'r')
        else:
            fin = open(codinginfoPath + '/indelFiles_hg38_v83/chromosomes_2016/chr' + str(chrom) + '_valid_2016_8_16.txt', 'r')
        lines = fin.read()
        lines = lines.split('>')[1:]
        for info in lines:
            enstInfoList.append(info)
        fin.close()
        enstInfoDict[chrom] = enstInfoList
    return enstInfoDict



#######################################################
# This function read the information about Ensembl
# non-3n (related to frameshift) mutations
#
#######################################################
def readEnsemblVariations(chrNum):
    #chromsomeList = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']
    chromsomeList = [chrNum]
    ensemblDict = {}
    for chrom in chromsomeList:
        ensemblChromDict = {}
##    	fin = open('/mnt1/indelFiles_V66/ensemblVariations/ensemblChr'+str(chrom)+'_info.txt','r')
##    	while True:
##    	    line = fin.readline()
##    	    if not line:
##    		break
##    	    if line[0]=='>':
##    		name = line.split()[0][3:]
##    		line = fin.readline()
##    		indel = line.split()[0]
##    		if name not in ensemblChromDict:
##    		    ensemblChromDict[name] = [indel]		
##                line = fin.readline()
##    		cdSeq = line.split()[0]
##    		ensemblChromDict[name].append(cdSeq)
##    	fin.close()
    	ensemblDict[chrom] = ensemblChromDict
    return ensemblDict
    	

    
def main():
    before = time()
##    print'\n--------------------------------------\n'
##    print 'Indel Information Extraction Program!'
##    print 'Author: Jing Hu & Pauline Ng\n'
##    print '--------------------------------------\n'
    codinginfoPath = sys.argv[4]
    #Herty added for build options
    buildVersion = sys.argv[5]
    chrNum = sys.argv[2]
    if chrNum!='X' and chrNum!='Y':
        chrNum = int(chrNum)
##    print 'Reading all transcripts...'
#    enstInfoDict = readExons(chrNum) #read all transcripts
    enstInfoDict = readExons(chrNum, codinginfoPath, buildVersion) #read all transcripts
##    print '.....................Done!\n'
    
    #indelInfo = sys.argv[1]

    #read the indel file
    #indel: e.g. 21,16365351,16365351,1,A
##    print '-->Reading indel list....'
    indelList = []
    inputFileName = sys.argv[1]
    fin = open(inputFileName,'r')
    for line in fin:
        line = line.split()
        if len(line)==0:
            continue
        indelList.append(line[0])
    fin.close()
##    print '...................Done!\n'

    errorFile1 = inputFileName + '.nonExonIndels'
    errorFile2 = inputFileName + '.missingIncompleteTranscripts'
    

    

    #protein domains
##    print '-->Reading protein domains...'
##    ProteinDomains = getDomains()
##    print '........................Done!\n'
    ProteinDomains = []
##    
##    #conservation scores
##    print '-->Reading conservation scores of all exons...'
##    DNAConvervationScores = readConservationScores(chrNum)
##    print '.........................................Done!\n'
##    
    #ensembl variations
##    print '-->Reading ensembl variations....'
    ensembleVariations = readEnsemblVariations(chrNum)
##    print '............................Done!\n'

    
    outFileName = sys.argv[3]

    for indelInfo in indelList:
        indelInfo = indelInfo.split(',')
        if len(indelInfo)!=4:
            print 'Incorrect Indel format'
            return
        chrNum = indelInfo[0]
        if chrNum!='X' and chrNum!='Y':
            chrNum = int(chrNum)
        indelStart = int(indelInfo[1])
        indelEnd = int(indelInfo[2])
        #strand = indelInfo[3]
        bases = indelInfo[3]
    
    
        #try to identify the indel type:
        #del: deletion; ins: insertion; snp: single nucletide polymorphism
        if bases == '/':
            indelType = 'del'
        elif len(bases.split('/'))==2:
            indelType = 'snp'
        else:
            indelType = 'ins'

        #extract indel information
        extractIndelInfo(enstInfoDict, chrNum, indelStart, indelEnd, bases, \
                         outFileName, ProteinDomains, ensembleVariations, errorFile1,errorFile2 )
##        print '\n\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n'

    after = time()

##    print('Time spent...'+str(after-before))


main()

    
    
