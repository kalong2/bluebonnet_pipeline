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
#  Pipeline for predicting 3n indels
#
#  Usage:
#
#  Jing Hu @ Franklin & Marshall College
#  Date: 2/4/2013
#
##################################################

import os
import sys
#gina

binPath = os.path.dirname (os.path.abspath (__file__))
print binPath, " is the binPath" 
#Added by Jing Hu at July 15, 2013
def read3nSTOPindels(file):
    '''
    This function reads the file which contains early stop 3n indels.
    These indels should be predicted using non3n predictor.
    '''
    #return set()
    if not os.path.exists(file):
        print "read3nSTOPindels:", file," should exist\n";
        return set()
    aset = set()
    fin = open(file, 'r')
    for line in fin:
        aset.add(line.split('\t')[0])
    fin.close()
    return aset


def outputEmptyInfo(file):
    fout = open(file, 'w')
    fout.close()

def main():
    inputFile = sys.argv[1]
    print "*************", inputFile
    path = sys.argv[2]  #big folder for temporary files
    #Gina added for "." issue
#    temp1 = inputFile.split('.')[-1]
#    temp = inputFile.split('.'+temp1)[0]
#    temp = inputFile.split('.')[0]+'.'+inputFile.split('.')[1]

    # Fixing Gina's wrong fix to path issue
    (filepath, filename) = os.path.split(inputFile)
    (shortname, extension) = os.path.splitext(filename)
    temp = os.path.join(filepath, shortname)
    # We put into the sub folder

    print "TEMPORARY FILE root: ",temp

    #added by Jing Hu at July 15, 2013
    _3nStopFile = temp+'.3nSTOP'
    _3nStopIndels = read3nSTOPindels(_3nStopFile)
    
    file0 = temp+'.3n.flanking'
    print "this is file0:",file0
    file1 = temp+'.flankingSeq'
    file2 = temp+'.flanking.flushed'
    file3 = temp+'.mapping'
    file4 = temp+'.topipeline'
    file5 = temp+'.topipeline.vectors'
    print "this is file5:", file5
    file6 = temp+'.mapback.vectors'
    file7 = temp+'.3nfeatures'
    file8 = temp+'.3n.predictions'
    
    binPath = sys.argv[3] #    binPath = "/var/www/sift-bin/indel-bin"


    codinginfoPath = sys.argv[4]
    # Herty added for build options
    buildVersion = sys.argv[5]

    alist = []
    fin = open(inputFile, 'r')
    for line in fin:
        item = line.split('\n')[0].split('\t')
        flanking = item[1]
        info = item[0]
        indel = info.split('#')[0]
        if indel[-1]==',':
            indel = indel[:-1]

        #added by JIng Hu at July 15, 2013
        if indel in _3nStopIndels:
            continue

        info = indel.split(',')
        startInd = int(info[1])
        endInd = int(info[2])
        allele = info[4]
        if '/' in allele:
            allele = '/'
        if endInd!=startInd and (endInd-startInd)%3 == 0:
            alist.append(line)
        elif endInd==startInd and len(allele)%3==0:
            alist.append(line)
    fin.close()

    if len(alist)==0:
        outputEmptyInfo(file8)
        return

    #Step 1: output all 3nindels to a file
    fout = open(file0, 'w')
    for line in alist:
        fout.write(line)
    fout.close()

    #Step 2: preprocessing
    command = 'python '+binPath+'/flanking.py '+file0+' '+ file1
    print"STEP 2: PREPROCESSING: ",command,"\n";
    os.system(command)

    if not os.path.exists(file1):
        outputEmptyInfo(file8)
        return

    #Step 3: flushed to left
    command = 'python '+binPath+'/flushIndelLeft.py '+file1+' '+file2+' '+file3
    os.system(command+'\n')

    if not os.path.exists(file2):
        outputEmptyInfo(file8)
        return

    #Step 4: prepare for pipeline of extracting 3n indels
    command = 'python '+binPath+'/toPipeline.py '+ file2+' '+file4
    os.system(command)

    if not os.path.exists(file4):
        outputEmptyInfo(file8)
        return
    #Step 5: get 3 features by running the pipeline
    command = 'python '+binPath+'/prepare3Features3nServer.py '+file4+' '+ path +' '+ codinginfoPath +' '+ buildVersion
    os.system(command)
 
    if not os.path.exists(file5):
        outputEmptyInfo(file8)
        return

    #Step 6: map back to original coordiantes if we flushed left some of indels
    command = 'python '+binPath+'/mapback.py '+file1+ ' '+file3+' '+file5+' '+file6
    os.system(command)

    if not os.path.exists(file6):
        outputEmptyInfo(file8)
        return
    
    #Step 7: combine into 4 features
    command = 'python '+binPath+'/combineFeatures.py '+file1+' '+file6+ ' '+file7
    os.system(command)


    if not os.path.exists(file7):
        outputEmptyInfo(file8)
        return

    #Step 8: make predictions using 3n J48 rules
    command = 'python '+binPath+'/3nJ48Pred.py '+file7+' '+file8
    os.system(command)  
    
    
    
main()
