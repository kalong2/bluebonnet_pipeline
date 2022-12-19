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
# This is the pipeline for processing  3n indel files and generate 3 features:
#           1) dnaConservationLeft;
#           2) PFDomainAffected;
#           3) avgDisScore.
#       
# To run:
#       $python3.x prepare3Features3n.py inputfile path
# Programmer: Jing Hu @ Franklin & Marshall College
# Copyright:  Jing Hu & Pauline Ng
# Date:       January 31, 2013
#
##############################################################################

import sys
import os
import time
import random


def main():
    inputFile = sys.argv[1]
    path = sys.argv[2]  #the big folder for temporary files
    tempFolder = str(int(str(time.time()).split('.')[0])+ random.randint(1,100))
    path = path + '/'+tempFolder
    os.system('mkdir '+path)
    print "tempFolder ", path 
    codinginfoPath = sys.argv[3]
    #Herty added for build options
    buildVersion = sys.argv[4]
    binPath = os.path.dirname (os.path.abspath (__file__))
    print binPath
#    binPath = "/var/www/sift-bin/indel-bin"
    absoluteFileName = inputFile.split('/')[-1]
    print "absoluteFileName,",absoluteFileName

    #split the indel file into multiple indels file based
    #on their chromsome
    #fin = open(inputFile, 'r')
    try:
        fin = open(inputFile, 'r')
#        fin = open("test.txt", 'rw')
#    except IOError as (errno, strerror):
#        print "Raised an I/O error {0}: {1}\n".format(errno, strerror)
    except:
        print "unexpected error"
        
    for line in fin:
        if len(line.split())>0:
            indel = line.split()[0].upper()
            chrNum = indel.split(',')[0]
            fileName = path+'/'+absoluteFileName+'.'+str(chrNum)
            print "Absolute: ====================== absoluteFilename", fileName
            fout = open(fileName, 'a')
            fout.write(indel+'\n')
            fout.close()

    #for each indel file, convert into vectors
    hasVectorFile = False
    print "hasVectorFile = False"
    for f in os.listdir(path):
        chrNum = f.split('.')[-1]
        indelFileName = path+'/'+f
        featureFileName = path+'/'+f+'.features'
        vectorFileName = path+'/'+f+'.vectors'

        print "vectorFileName ", vectorFileName 

        print "featureFileName:", featureFileName        
        os.system('python '+binPath+'/'+'FeatureExtract_3nindels_v1server.py '+indelFileName+ ' ' \
                  +str(chrNum)+' '+featureFileName +' '+ codinginfoPath +' '+ buildVersion)
        if not os.path.exists(featureFileName):
	    continue
        
        os.system('python '+binPath+'/'+'createVectors_3nindels_v1server.py '+featureFileName+' ' \
                  +str(chrNum)+' '+vectorFileName +' '+ codinginfoPath +' '+ buildVersion)
        print "createVectors_3nindels_v1server.py", featureFileName, str(chrNum), vectorFileName, codinginfoPath, buildVersion
        hasVectorFile = True
        print "hasVectorFile = True"

 #   allVectorFile = path+'/'+absoluteFileName+'.vectors'
    allVectorFile = inputFile+'.vectors'
    print "input file ======================",inputFile,'.vectors'
    #in .vector file, columns are : indel, ensg, dnaConservationLeft
    #,PFdomainAffected, avgDisScore.
    if hasVectorFile:
        os.system('cat '+path+'/*.vectors >'+ allVectorFile)
        print "all vectors ======================",inputFile,'.vectors'


    
            

main()


