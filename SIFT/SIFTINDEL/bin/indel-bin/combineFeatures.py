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
    


def splitInParts(aseq):
    allele = ''
    for ch in aseq:
        if ord(ch)>=97 and ord(ch)<=122:
            allele += ch

    if len(allele)==0:
        return '', '', ''
    parts = aseq.split(allele)
    if len(parts)==0:
        return '', allele.upper(), ''
    else:
        return parts[0].upper(), allele.upper(), parts[1].upper()


def getIndex(alist, indel):
    for i in range(len(alist)):
        if alist[i][0]==indel:
            return i
    return None
            
        
def main():   
    indelList = []
    fileName = sys.argv[1]
    print "in mapping.py inputfilename", fileName
    fin = open(fileName, 'r')
    for line in fin:
##        if 'Warning: NCBI reference miscall!':
##            line = line.replace('Warning: NCBI reference miscall!', '')
        info = line.split('\n')[0].split('\t')
        indelList.append(info)
    fin.close()

    indelFeatures = {}
    fin = open(sys.argv[2], 'r')
    for line in fin:
        info = line.split()
        #gina
        if len(info)==0:
           continue
        indel = info[0]
        ensgID = info[1]
        indelFeatures[indel+'_'+ensgID] = info[2:]
        
    

    fin.close()

   
    fout = open(sys.argv[3], 'w')
    fout.write('Indel\tensgID\tisRepeat\tdnaConservScoreLeft\tfracPFDomainsAffected\tavgDisScore\n')
    print "in mapping.py, output filename", fout 
    adict = {}
    for indel in indelFeatures:
        adict[indel] = False

    for item in indelFeatures:
        indel = item.split('_')[0]
        ensgID = item.split('_')[1]
        
        index = getIndex(indelList, indel)
        if index == None:
            continue

        adict[item] = True

        repeatResult = detectRepeat(indelList[index][1])

        fout.write(indel+'\t'+ensgID+'\t')
        

        #output repeat features
        if repeatResult==None:
            fout.write('0\t')
        else:
            fout.write('1\t')

        for x in indelFeatures[item]:
            fout.write(str(x)+'\t')

        fout.write('\n')
        


main()
