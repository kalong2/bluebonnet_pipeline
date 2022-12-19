###############################################################################
#
# This is the pipeline for processing indel files and generate prediction
# results.
# To run:
#       $python indelPredPipeline.py inputfile
# Programmer: Jing Hu @ Franklin & Marshall College
# Copyright: Jing Hu & Pauline Ng
# Date: May 23, 2011
#
##############################################################################

import sys
import os
import time
import stat
binPath = os.path.dirname (os.path.abspath (__file__))


def predJ48(max_relative_indel_location, fract_conservedDNABase_affected, \
                             min_disExonBoundary,_25_fract_conservedAminoAcid_affected):
    '''
    Prediction using J48 decision tree parsed from weka package
    '''
    if fract_conservedDNABase_affected <= 0.042735 :
        if fract_conservedDNABase_affected <= 0.011527 :
            predEffect = "neutral"
            confidence = 0.961
            ClassificationPath = "FS_Classification_path_1"
        if fract_conservedDNABase_affected > 0.011527 :
            if max_relative_indel_location <= 0.854749 :
                predEffect = "neutral"
                confidence = 0.915
                ClassificationPath = "FS_Classification_path_2"
            if max_relative_indel_location > 0.854749 :
                if _25_fract_conservedAminoAcid_affected <= 0.009494 :
                    predEffect = "neutral"
                    confidence = 0.814
                    ClassificationPath = "FS_Classification_path_3"
                if _25_fract_conservedAminoAcid_affected > 0.009494 :
                    predEffect = "disease"
                    confidence = 0.783
                    ClassificationPath = "FS_Classification_path_4"
    if fract_conservedDNABase_affected > 0.042735 :
        if _25_fract_conservedAminoAcid_affected <= 0 :
            predEffect = "neutral"
            confidence = 0.787
            ClassificationPath = "FS_Classification_path_5"
        if _25_fract_conservedAminoAcid_affected > 0 :
            if min_disExonBoundary <= 6 :
                if min_disExonBoundary <= 1 :
                    predEffect = "neutral"
                    confidence = 0.773
                    ClassificationPath = "FS_Classification_path_6"
                if min_disExonBoundary > 1 :
                    predEffect = "disease"
                    confidence = 0.529
                    ClassificationPath = "FS_Classification_path_7"
            if min_disExonBoundary > 6 :
                if max_relative_indel_location <= 0.086525 :
                    if max_relative_indel_location <= 0.004348 :
                        predEffect = "neutral"
                        confidence = 1.00
                        ClassificationPath = "FS_Classification_path_8"
                    if max_relative_indel_location > 0.004348 :
                        predEffect = "disease"
                        confidence = 0.579
                        ClassificationPath = "FS_Classification_path_9"
                if max_relative_indel_location > 0.086525 :
                    if fract_conservedDNABase_affected <= 0.062176 :
                        if min_disExonBoundary <= 108 :
                            predEffect = "disease"
                            confidence = 0.731
                            ClassificationPath = "FS_Classification_path_10"
                        if min_disExonBoundary > 108 :
                            predEffect = "neutral"
                            confidence = 0.818
                            ClassificationPath = "FS_Classification_path_11"
                    if fract_conservedDNABase_affected > 0.062176 :
                        predEffect = "disease"
                        confidence = 0.858
                        ClassificationPath = "FS_Classification_path_12"
    return (predEffect, confidence, ClassificationPath)

def getComplement(allele):
    newAllele = ''
    for ch in allele:
        if ch=='A':
            newAllele += 'T'
        elif ch=='T':
            newAllele += 'A'
        elif ch=='C':
            newAllele += 'G'
        elif ch=='G':
            newAllele += 'C'
    #fix the bug: Jing Hu 2/3/2013
    alist = list(newAllele)
    alist.reverse()
    newAllele = ''
    for x in alist:
        newAllele+=x
    return newAllele

    
def getNewIndelFormat(chrNum, startIdx, endIdx, direction, allele):
    '''
    The function returns the indel format in positive chain
    '''
    newIndelFormat = ''
    if direction=='1':
        newIndelFormat = chrNum+','+startIdx+','+endIdx+','+allele
    else:
        if startIdx!=endIdx:
            newIndelFormat = chrNum+','+startIdx+','+endIdx+','+allele
        else:
            newAllele = getComplement(allele)
            newIndelFormat = chrNum+','+startIdx+','+endIdx+','+newAllele
            
    #print(newIndelFormat)

    return newIndelFormat


#added on Feb 4, 2013
def is3nIndel(startIdx, endIdx, allele):
    '''
    Check to see if an indel is 3n or not
    '''
    startIdx = int(startIdx)
    endIdx = int(endIdx)
    if startIdx!=endIdx and (endIdx-startIdx) % 3 == 0:
        return True
    elif startIdx==endIdx and len(allele)>0 and len(allele) % 3 == 0:
        return True

    return False

def read3nSTOPindels(file):
    '''
    This function reads the file which contains early stop 3n indels.
    These indels should be predicted using non3n predictor.
    '''
    #return set()
    if not os.path.exists(file):
        return set()
    aset = set()
    fin = open(file, 'r')
    for line in fin:
        aset.add(line.split('\t')[0])
    fin.close()
    return aset
    
def main():
    start = time.time()
    inputFile  = sys.argv[1]
    binPath = sys.argv[3] #"/var/www/sift-bin/indel-bin"
    path = sys.argv[2]  #get temporary folder from system argument
#    print "PATH: ",path
    codinginfoPath = sys.argv[4]
    tempPath = sys.argv[5]
    #Herty added for build options
    buildVersion = sys.argv[6]
#    print "This is tempPath: ",tempPath
    #read the 3n indels which cause early stop
    _3nStopFile = ''
#    print "inputFile: ",inputFile
    (early_stop_filepath, filename) = os.path.split(inputFile)
    (shortname, extension) = os.path.splitext(filename)
    _3nStopFile = os.path.join(early_stop_filepath, shortname)
    _3nStopFile = _3nStopFile+'.3nSTOP'
#    print "_3nStopFile: ",_3nStopFile
    _3nStopIndels = read3nSTOPindels(_3nStopFile)
#    print "_3nStopIndels:",_3nStopIndels
#    segs = inputFile.split('.')
#    for item in segs[:-1]:
#        _3nStopFile += item+'.'
#    _3nStopFile = _3nStopFile+'3nSTOP'
#    _3nStopIndels = read3nSTOPindels(_3nStopFile)
    

    outputFile = inputFile+'.predictions'
    print "nonoutputFile",outputFile
    errorFile1 = inputFile+'.nonExonIndels'
    errorFile2 = inputFile+'.missingIncompleteTranscripts'

    absoluteFileName = inputFile.split('/')[-1]
#    print"absoluteFileName:",absoluteFileName

    try:
        fin = open(inputFile, 'r')
#        fin = open("test.txt", 'rw')
#    except IOError as (errno, strerror):
#        print "Raised an I/O error {0}: {1}\n".format(errno, strerror)
    except:
        print "unexpected error"

    for line in fin:
#	print "LINE: ",line,'\n'
#        print "_3nStopIndels",_3nStopIndels
    	#line = line.split('\n')[0]
       	if len(line.split(','))>3 and len(line.split())>0:
            #indel = line.split()[0].upper()
	    line = line.split('\n')[0]
	    line = line.split('\r')[0]
            fields = line.split(',')
            startIdx = ''
            endIdx = ''
            allele=''
  	    if len(fields)>=5:
            	chrNum = fields[0]
            	startIdx = fields[1]
            	endIdx = fields[2]
            	direction = fields[3]
            	allele = fields[4]
#                print "chrNum:",chrNum,",startIdx:",startIdx,",endIdx:",endIdx,",direction:",direction,",allele:",allele
            	#newIndelFormat = chrNum+','+startIdx+','+endIdx+','+allele
            	newIndelFormat = getNewIndelFormat(chrNum, startIdx, endIdx, direction, allele)
            else:
                #print "There you are\n"
            	chrNum = fields[0]
            	startIdx = fields[1]
		endIdx = fields[2]
		#direction = fields[3]
            	allele = fields[3]
            	newIndelFormat = chrNum+','+startIdx+','+endIdx+','+allele
            #print ('***processing '+newIndelFormat+'\n')
            if '/' in line:
                line = line.split(',')
                line = line[0]+','+line[1]+','+line[2]+','+line[3]+',-/'
            #print(line+'_______')                
            if not line in _3nStopIndels and is3nIndel(startIdx, endIdx, allele):
#                print "not line in _3nStopIndels",line
                continue
            fileName = tempPath+'/'+absoluteFileName+'.'+str(chrNum)
#            print "This is fileName: ",fileName
            fout = open(fileName, 'a')
	    fout.write(newIndelFormat+'\n')
	    #print("fileName: "+fileName)
	    #print("newIndelFormat: "+newIndelFormat)
            fout.close()
    

    #for each indel file, convert into vectors
    hasVectorFile = False
    for f in os.listdir(tempPath):
        chrNum = f.split('.')[-1]
        indelFileName = tempPath+'/'+f
        featureFileName = tempPath+'/'+f+'.features'
        vectorFileName = tempPath+'/'+f+'.vectors'
        os.system('python '+binPath+'/'+'featureExtract_non3nindels_v2server.py '+indelFileName+ ' ' \
                      +str(chrNum)+' '+featureFileName + ' ' + codinginfoPath + ' ' +  buildVersion)
	if not os.path.exists(featureFileName):
	    continue
        os.system('python '+binPath+'/'+'createVectors_non3nindels_v2server.py '+featureFileName+' ' \
                  +str(chrNum)+' '+vectorFileName + ' ' + codinginfoPath + ' ' +  buildVersion)
        hasVectorFile = True

    hasNonExonIndels = False
    hasMissingIncompleteTranscripts = False
    for f in os.listdir(tempPath):
        if f.split('.')[-1]=='nonExonIndels':
            hasNonExonIndels = True
        elif f.split('.')[-1]=='missingIncompleteTranscripts':
            hasMissingIncompleteTranscripts = True
    if hasNonExonIndels:
        os.system('cat '+tempPath+'/*.nonExonIndels >'+errorFile1)
    if hasMissingIncompleteTranscripts:
        os.system('cat '+tempPath+'/*.missingIncompleteTranscripts >'+ errorFile2)
            
    if not hasVectorFile:
        #if there is no vector (or valid indels or indels in exons of valid transcripts,
        #then 
        fout = open(outputFile, 'w')
        fout.close()
#       os.system('rm '+tempPath+'/*')
#    	os.system('rmdir '+tempPath)
    	return
    	

    allVectorFile = tempPath+'/'+absoluteFileName+'.vectors'
    #merge all vector files
    
    os.system('cat '+tempPath+'/*.vectors >'+ allVectorFile)

    indelsList = []

    fin = open(allVectorFile, 'r')
    for line in fin:
        indelVector = line.split()
        indelsList.append(indelVector)
    fin.close()

    fout = open(outputFile, 'w')
    for indel in indelsList:
        #print(indel)
        pred = predJ48(float(indel[2]), float(indel[3]),float(indel[4]),float(indel[5]))
        #fout.write(indel[2]+' '+indel[3]+' '+indel[4]+' '+indel[5])
        fout.write(str(indel[0])+' '+str(indel[1])+' '+pred[0]+' '+str(pred[1])+' '+pred[2]+'\n')
    
    fout.close()
    
 
    end = time.time()
       
            

main()


