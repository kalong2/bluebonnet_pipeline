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

#get the flanking sequence of each indel

import sys
import os
import random

def getComplement(aseq):
    newSeq = ''
    for i in range(len(aseq)-1, -1, -1):
        ch = aseq[i]
        if ch=='A':
            newSeq+='T'
        elif ch=='T':
            newSeq+='A'
        elif ch=='C':
            newSeq+='G'
        elif ch=='G':
            newSeq+='C'
        else:
            newSeq+=ch
    return newSeq
            

def main():
    fin = open(sys.argv[1], 'r')
    alist = []
    for line in fin:
        if line not in alist:
            alist.append(line)
    fin.close()
            
    fout = open(sys.argv[2], 'w')
    for line in alist:
        line = line.split('\n')[0].split('\t')
        flanking = line[1]
        info = line[0]
        indel = info.split('#')[0]
        if indel[-1]==',':
            indel = indel[:-1]

        info = indel.split(',')
        chrNum = info[0]
        startInd = info[1]
        endInd = info[2]
        orientation = info[3]
        allele = info[4]
        if '/' in allele:
            allele = '/'
        if orientation=='-1' and '/' not in allele:
            allele =  getComplement(allele)
        
        indel = chrNum+','+startInd+','+endInd+',1,'+allele
        fout.write(indel+'\t'+flanking+'\n')
  
    fout.close()
    


main()
