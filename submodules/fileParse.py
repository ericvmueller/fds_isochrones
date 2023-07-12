"""
/***************************************************************************
 fileParse module

 This module contains a set of functions used to parse relevant data from
 various file types associated with FDS. The functionality is largely
 borrowed from the pyfdstools package:
 https://github.com/johodges/pyfdstools, but has been adapted for simplicity
 and to reduce requirements for external dependencies.

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

import numpy as np
import struct
from collections import defaultdict

class SLCT:
    def __init__(self,file):
        self.f=open(file,'rb')
        self.quantity=None
        self.sName=None
        self.uts=None
        self.iX=None
        self.eX=None
        self.iY=None
        self.eY=None
        self.iZ=None
        self.eZ=None
        self.times=None
        self.currentTime=None
        self.data=None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.f.close()

    # Return relevant SLCT header data
    def readHeader(self):
        self.f.seek(0)
        data = self.f.read(142)
        header = data[:110]
        size = struct.unpack('>iiiiii', data[115:139])
        tmp = header.split(b'\x1e')
        self.quantity = tmp[1].decode('utf-8').replace('\x00','').strip(' ')
        self.shortName = tmp[3].decode('utf-8').replace('\x00','').strip(' ')
        self.units = tmp[5].decode('utf-8').replace('\x00','').strip(' ')

        self.iX, self.eX, self.iY, self.eY, self.iZ, self.eZ = size

    # Return array of data times from SLCT file
    def readTimes(self):
        self.f.seek(0)
        self.readHeader()
        (NX, NY, NZ) = (self.eX-self.iX, self.eY-self.iY, self.eZ-self.iZ)
        data = self.f.read()
        if len(data) % 4 == 0:
            fullFile = np.frombuffer(data, dtype=np.float32)
        else:
            remainder = -1*int(len(data) % 4)
            fullFile = np.frombuffer(data[:remainder], dtype=np.float32)
        self.times = fullFile[2::(NX+1)*(NY+1)*(NZ+1)+5]
        self.f.seek(0)

    # Return a single time and data array from a slice file
    def readRecord(self):

        # if at the start of the file
        if self.f.tell()<142:
            self.readHeader()

        (NX, NY, NZ) = (self.eX-self.iX, self.eY-self.iY, self.eZ-self.iZ)
        shape = (NX+1, NY+1, NZ+1)
        _ = np.frombuffer(self.f.read(8), dtype=np.float32)
        self.currentTime = np.frombuffer(self.f.read(4), dtype=np.float32)
        _ = np.frombuffer(self.f.read(8), dtype=np.float32)
        try:
            self.data = np.frombuffer(self.f.read((NX+1)*(NY+1)*(NZ+1)*4),
                                 dtype=np.float32)
            self.data = np.reshape(self.data, shape, order='F')
        except:
            self.data = None

def parseOUT(file):
    lat,lon=-1e6,-1e6
    with open(file) as f:
        # while True:
        for i in range(257):
            line=f.readline()
            # found geographic info
            if line[1:11]=='Geographic':
                line=f.readline()
                line=f.readline()
                lat=float(line[20:])
                line=f.readline()
                lon=float(line[20:])
                break
            # no geographic info available, stop reading file
            if line[1:14]=='Miscellaneous':
                break

    return lat,lon

# parse SMV file, in this case to find grid and SLCT information
def parseSMV(file):
    with open(file) as f:
        linesSMV = f.readlines()

        grids=[]
        SLCTfiles=defaultdict(bool)
        for il in range(len(linesSMV)):
            if ('GRID' in linesSMV[il]):
                gridTRNX, gridTRNY, gridTRNZ = parseGRID(linesSMV, il)
                grids.append([gridTRNX.copy(),
                              gridTRNY.copy(),
                              gridTRNZ.copy()])

            if ('SLCT' in linesSMV[il]):
                file = '%s.sf'%(linesSMV[il+1][1:].split('.sf')[0])
                SLCTfiles[file] = defaultdict(bool)
                SLCTfiles[file]['QUANTITY'] = linesSMV[il+2].strip()
                SLCTfiles[file]['SHORTNAME'] = linesSMV[il+3].strip()
                SLCTfiles[file]['UNITS'] = linesSMV[il+4].strip()
                SLCTfiles[file]['LINETEXT'] = linesSMV[il]
                SLCTfiles[file]['MESH']=int(SLCTfiles[file]['LINETEXT'][5:10])
                SLCTfiles[file]['AGL']=float(SLCTfiles[file]['LINETEXT'][11:20])

    return SLCTfiles,grids

def parseGRID(lines, i):
    gridPts = [int(x) for x in lines[i+1].replace('\n','').split()]
    (gridTRNX, gridTRNY, gridTRNZ) = ([], [], [])
    (xind0, xind1) = (i+8, i+9+gridPts[0])
    (yind0, yind1) = (i+12+gridPts[0], i+13+gridPts[0]+gridPts[1])
    zind0 = i+16+gridPts[0]+gridPts[1]
    zind1 = i+17+gridPts[0]+gridPts[1]+gridPts[2]
    xlines = lines[xind0:xind1]
    ylines = lines[yind0:yind1]
    zlines = lines[zind0:zind1]
    for x in xlines:
        gridTRNX.append([float(y) for y in x.replace('\n','').split()])
    for x in ylines:
        gridTRNY.append([float(y) for y in x.replace('\n','').split()])
    for x in zlines:
        gridTRNZ.append([float(y) for y in x.replace('\n','').split()])
    gridTRNX = np.array(gridTRNX)
    gridTRNY = np.array(gridTRNY)
    gridTRNZ = np.array(gridTRNZ)
    return gridTRNX, gridTRNY, gridTRNZ


