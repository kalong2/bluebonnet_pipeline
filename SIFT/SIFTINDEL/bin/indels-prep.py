#!/usr/bin/env python
import sys

if len(sys.argv)==1 : 
  print 'Usage: ',sys.argv[0],'<input file>'
  sys.exit(1)

filename=sys.argv[1]

try:
  errFileName=filename+'.err'
  oFileErr=open(errFileName,'w')
except:
  print 'Cannot open %s for write!' % errFileName
  sys.exit(1)

try:
  iFile=open(filename)
except:
  print >>oFileErr,'Cannot open %s!' % filename
  sys.exit(1)

try:
  okFileName=filename+'.ok'
  oFileOk=open(okFileName,'w')
except:
  print >>oFileErr,'Cannot open %s for write!' % okFileName
  sys.exit(1)

chrTable=[str(l) for l in range(1,23)]+['X','Y']
coorTable=['A','T','G','C']
allelesTable=['/','-']+coorTable

def checkChr(chrNum):
  if chrNum not in chrTable:
    return True
  return False

def checkPos(pos):
  try:
     posInt=int(pos)
     return False
  except:
     return True

def checkCoor(coor):
  coors=coor.split('/')
  if len(coors)!=2:
    return True
  if (coors[0] not in coorTable) or (coors[1] not in coorTable):
    return True
  return False


def checkAlleles(alleles):
  if alleles=='': return True

  for c in alleles:
    if c not in allelesTable: return True

  slash=alleles.split('/')
  if len(slash)>2: return True
  return False

filename='input'
lineCount=0
error=0
for line in iFile:
  lineCount+=1
  line=line.strip()
  orgLine=line
  sharp=line.find('#')
  if sharp<0: sharpC=''
  else:
    sharpC=line[sharp+1:]
    line=line[:sharp].strip()

#  print line
  line=line.replace(' ','')  #strip all spaces
  line=line.replace('\t','') #strip all tabs
  
  cols=line.split(',')

  if len(cols)<5:
    if line=='': continue
    if line=='ONTENT-TYPE:APPLICATION/OCTET-STREAM':continue
    if line=='ONTENT-TYPE: TEXT/PLAIN':continue
    print >>oFileErr,'input:line %d:[%s] is not in the correct format.' % (lineCount,orgLine)
    error=1
    continue
  
  chrNum=cols[0]
  chrNum=chrNum.upper().replace('CHR','')
  if checkChr(chrNum):
    error=1
    print >>oFileErr,'%s:line %d:[%s] chrNum[%s] is not in the correct format.' % (filename,lineCount,orgLine,chrNum)
    continue

  pos1=cols[1]
  if checkPos(pos1):
    error=1
    print >>oFileErr,'%s:line %d:[%s] pos1[%s] is not in the correct format.' % (filename,lineCount,orgLine,pos1)
    continue

  pos2=cols[2]
  if checkPos(pos2):
    error=1
    print >>oFileErr,'%s:line %d:[%s] pos1[%s] is not in the correct format.' % (filename,lineCount,orgLine,pos2)
    continue

  #intPos1=int(pos1)
  #intPos2=int(pos2)
  #if intPos1>intPos2:
  #  error=1
  #  print >>oFileErr,'%s:line %d:[%s] pos1[%s] should not greater than pos2[%s]!' % (filename,lineCount,line,pos1,pos2)
  #  continue
  
  direction=cols[3]
  if direction!='1' and direction!='-1':
    error=1
    print >>oFileErr,'%s:line %d:[%s] direction[%s] is not in the correct format.' % (filename,lineCount,orgLine,direction)
    continue
  
  alleles=cols[4]
  if checkAlleles(alleles):
    error=1
    print >>oFileErr,'%s:line %d:[%s] alleles[%s] is not in the correct format.' % (filename,lineCount,orgLine,direction)
    continue

  comments=','.join(cols[5:])
  if sharpC!='': comments+='#'+sharpC
  if comments!='':
    print >>oFileOk,'%s,%s,%s,%s,%s,%s'%(chrNum,pos1,pos2,direction,alleles,comments)
  else:
    print >>oFileOk,'%s,%s,%s,%s,%s'%(chrNum,pos1,pos2,direction,alleles)

sys.exit(error)

