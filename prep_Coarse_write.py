###########################################################################
#imports
###########################################################################
import os
#from pathlib import Path
import numpy as np
###########################################################################
#module global variables
###########################################################################

PORO_NAME = 'PORO'

###########################################################################
#functions:
###########################################################################

###########################################################################

def mkDirIfExist(path):
  if not os.path.exists(path):
    os.mkdir(path)

###########################################################################

def mkSimpleInputFile(path, name, keyword, dataVec):
  fullNm = name + '_C' + '.ECLIN'
  fullFile = os.path.join(path, fullNm)
  form = '%f'
  delim = '\n'
  finish = '/'
  
  fileSave = open(fullFile, 'w')
  fileSave.write('%s\n' % keyword)
  np.savetxt(fileSave, dataVec,fmt = form, delimiter = delim)
  fileSave.write('%s' % finish)
  fileSave.close()

###########################################################################

def mkTransFiles(tranList, path, nameBase):
  addenums = {0:'X', 1:'Y', 2:'Z'}

  for coord in range(len(addenums)):
    name = nameBase + addenums[coord] 
    mkSimpleInputFile(path, name, name, tranList[coord])

###########################################################################

def wellType(index):
  if index:
    type = 'PROD'
  else:
    type = 'INJ'
  return type

def createWelspecs(wellLocs, fileStream):
  KEY = 'WELSPECS'

  fileStream.write('%s\n' % KEY)

  for type in range(2):
    title = wellType(type)
    if type == 0:
      FMT = '%s %s %d %d 1* WATER 2* SHUT YES/\n'
    else:
      FMT = '%s %s %d %d 1* OIL 2* SHUT YES/\n'
    for wellNum in range(len(wellLocs[type])):
      iDim = wellLocs[type][wellNum][0]
      jDim = wellLocs[type][wellNum][1]
      spec = title[0] + str(wellNum + 1)

      fileStream.write(FMT % (spec, title, iDim, jDim))

  fileStream.write('%s\n\n' % '/')

def createCompdat(wellLocs, Wis, rw, fileStream, loc_z):
  KEY = 'COMPDAT'
  FMT = '%s %d %d %d %d OPEN 1* %f 0.2 1* 0/\n'

  fileStream.write('%s\n' % KEY)

  for type in range(2):
    title = wellType(type)
    
    for wellNum in range(len(wellLocs[type])):
      iDim = wellLocs[type][wellNum][0]
      jDim = wellLocs[type][wellNum][1]
      spec = title[0] + str(wellNum + 1)
      
      point = [Wis[type][x][wellNum] for x in range(len(Wis[type]))]
      for perf in range(len(point)):
        fileStream.write(FMT % (spec, iDim, jDim, loc_z[type][wellNum][perf] + 1, loc_z[type][wellNum][perf] + 1,
                         point[perf]))

      if type == 0 or wellNum != (len(wellLocs[type]) - 1):
        fileStream.write('\n')

  fileStream.write('%s\n\n' % '/') 

def createWellstr(wellLocs, fileStream):
  KEY = 'WELLSTRE'
  FMT = '%s 1 0/\n'

  fileStream.write('%s\n' % KEY)
  for wellNum in range(len(wellLocs[0])):
    spec = 'I' + str(wellNum + 1)

    fileStream.write(FMT % (spec))

  fileStream.write('%s\n\n' % '/') 

def wellDescriptions(path, wellLocs, Wis, rw, count, loc_z):
  name = 'WELLDESC' + '_' + str(count) + '.txt'
  #fullFile = Path(path + '/' + name)
  fullFile = os.path.join(path, name)
  
  fileToSave = open(fullFile, 'w')

  createWelspecs(wellLocs, fileToSave)
  createCompdat(wellLocs, Wis, rw, fileToSave, loc_z)
  createWellstr(wellLocs, fileToSave)

  fileToSave.close()

###########################################################################

def wellConstraints(wellDat, path, count):
  fileName = 'WELLCONST'+ '_' + str(count) + '.txt'
  KEYS = ['WCONINJE', 'WCONPROD']
  FMTS = ['%s WATER OPEN BHP 2* %f /\n', '%s OPEN BHP 5* %f/\n']
  #fullFile = Path(path + '/' + fileName)  
  fullFile = os.path.join(path, fileName)

  fileToSave = open(fullFile, 'w')
  
  for type in range(2):
    title = wellType(type)
    head = KEYS[type]
    datUse = wellDat[type]


    fileToSave.write('%s\n' % head)
    for wellNum in range(len(datUse[0])):
      front = title + str(wellNum + 1)
      fileToSave.write(FMTS[type] % (front, datUse[0][0][0]))

    fileToSave.write('%s\n\n' % '/')

  fileToSave.close()

###########################################################################

def addDimensions(cGrids, fileStream):
  KEY = 'DIMENS'
  FMT = '%d %d %d\n'

  fileStream.write('%s\n' % KEY)
  fileStream.write(FMT % tuple(cGrids))
  fileStream.write('%s\n' % '/')

def addVolumes(cDims, fileStream):
  KEY = 'VOLUME'
  FMT = '%f /\n'

  volume = cDims[0]*cDims[1]*cDims[2]
  fileStream.write('%s\n' % KEY)
  fileStream.write(FMT % volume)
  fileStream.write('%s\n' % '/')

def addDepths(fact, cGrids, cDims, top, fileStream):
  KEY = 'DEPTH'
  FMT = '%d*%6.2f '

  gridInLayer = cGrids[0]*cGrids[1]
  dz = cDims[-1] / fact[-1]

  fileStream.write('%s\n' % KEY)
  for lyr in range(cGrids[-1]):
    fileStream.write(FMT % (gridInLayer, top + cDims[-1] * lyr))
  fileStream.write('\n%s\n' % '/')
  

def makeDimensionsInput(top, cGrids, cDims,fact, path, count):
  fileName = 'DIMENSIONS'+ '_' + str(count) + '.txt'
  #fullFile = Path(path + '/' + fileName) 
  fullFile = os.path.join(path, fileName)

  fileToSave = open(fullFile,'w')

  addDimensions(cGrids, fileToSave)
  addVolumes(cDims, fileToSave)
  addDepths(fact, cGrids, cDims, top, fileToSave)

  fileToSave.close()