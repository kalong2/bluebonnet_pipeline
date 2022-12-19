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
#
# This program fix the indel so that all indels are at the left
# E.g.,  ATCATCatcCGT should be treated as atcATCATCCGT
#   To run: python inputfile flushedfile mappingfile
#
# Copyright: Jing Hu @ Franklin & Marshall College
# Date: Feb 4, 2013
#
###############################################################

import sys
import os

def findSmallestCandidate(aseq):
    for i in range(1, len(aseq)+1):
        seg = aseq[:i]
        if aseq == seg*(len(aseq)//i):
            return seg

def detectRepeat(aseq):
    '''
    Detects if there is repeats, if so, what is the repeat.
    '''
    #print(aseq)
    aseq = aseq.upper()
    parts = aseq.split('-')
    midSeq = parts[1]
    leftFlank = parts[0]
    rightFlank = parts[2]

    candidate = findSmallestCandidate(midSeq)

    repeatNum = len(midSeq)//len(candidate)

    rightRepeat = False
    leftRepeat = False
    #check if any repeats on the right
    for i in range(0, len(rightFlank), len(candidate)):
        left = i
        right = min(i+len(candidate), len(rightFlank))
        if rightFlank[left:right] == candidate:
            repeatNum += 1
            rightRepeat = True
        else:
            break
        

    for i in range(len(leftFlank),0, -1*len(candidate)):
        left = max(i-len(candidate), 0)
        right = i 
        if leftFlank[left:right]==candidate:
            repeatNum += 1
            leftRepeat = True
        else:
            break

    if rightRepeat or leftRepeat:
        return (candidate, repeatNum, leftRepeat, rightRepeat)
    else:
        return None

def main():

            
    
    fin = open(sys.argv[1], 'r')
    fout = open(sys.argv[2], 'w')  #flushed to left
    fout2 = open(sys.argv[3], 'w')  #remember the mapping between old and new coordinates

    for line in fin:
        #print(line)
        info = line.split('\n')[0].split('\t')
        indel = info[0]
        flankingSeq = info[1]
        upperFlankingSeq = flankingSeq.upper()
        
        
        repeatResult = detectRepeat(upperFlankingSeq)
      
        indelInfo = indel.split(',')
        chrNum = indelInfo[0]
        start = int(indelInfo[1])
        end = int(indelInfo[2])
        orientation = indelInfo[3]
        allele = indelInfo[4]

        flankingSeqInfo = upperFlankingSeq.split('-')
        leftFlankingSeq = flankingSeqInfo[0]
        midSeq = flankingSeqInfo[1]

        tempSeq = leftFlankingSeq+midSeq

        if repeatResult==None:
            fout.write(line)
            continue
            

        candidate = repeatResult[0]
        offset = 0
        while tempSeq[-1*len(midSeq)+offset:len(tempSeq)+offset]==midSeq:
            offset = offset - len(candidate)

        newStart = start + offset + len(candidate)
        newEnd = end + offset + len(candidate)
        newIndel = chrNum+','+str(newStart)+','+str(newEnd)+','+orientation+','+allele

        fout.write(newIndel+'\t'+flankingSeq+'\n')
        if indel!=newIndel:
            fout2.write(indel+'\t'+newIndel+'\t'+flankingSeq+'\n')
        
        

    fout.close()
    fout2.close()
    fin.close()


main()
