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
#  map 3 features back to original coordinates
#  Usage: python mapback.py flankSeqFile mapFile vectorFile mapbackVectorFile
#
#  2/4/2013
#
#######################################

import sys
import os

def main():
    mapping = {}
    fin = open(sys.argv[1], 'r')
    for line in fin:
        indel = line.split()[0]
        indel = indel.upper()
        mapping[indel] = indel
    fin.close()

    fin = open(sys.argv[2], 'r')
    for line in fin:
        line = line.split()
        oldIndel = line[0].upper()
        newIndel = line[-2].upper()

        if oldIndel in mapping:
            mapping[oldIndel] = newIndel        

    fin.close()

    adict = {}
    fin = open(sys.argv[3], 'r')

    for line in fin:
        indel = line.split()[0]
        info = indel.split(',')
        newIndel = info[0]+','+info[1]+','+info[2]+',1,'+info[3]
        adict[newIndel] = line      
    fin.close()


    count=0
    count2 = 0
    fout = open(sys.argv[4], 'w')

    for indel in mapping:
        mapIndel = mapping[indel]
 
        if mapIndel not in adict:
            #print(indel, mapIndel)
            count+=1
        else:
            count2+=1
            info = adict[mapIndel].split()[1:]
            fout.write(indel+'\t')
            for x in info:
                fout.write(x+'\t')
            fout.write('\n')
            
    fout.close()

    

    
    




main()
