"""
/***************************************************************************
 fileParse module

 This module contains a set of functions used to parse relevant data from
 various file types associated with FDS. Some parsing functionality is
 borrowed from pyfdstools: https://github.com/johodges/pyfdstools,
 but has been adapted for simplicity and to reduce requirements for
 external dependencies.
"""

import numpy as np
import struct
from collections import defaultdict

# class to control parsing of cartesian boundary files
class BNDF:
    def __init__(self,file):
        self.f=open(file,'rb')
        self.quantity=None
        self.shortName=None
        self.units=None
        self.nPatches=None
        self.times=None
        self.currentTime=None
        self.patchData={}

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.f.close()

    # Return relevant cartesian BNDF header data
    def readHeader(self):
        self.f.seek(0)
        buff = self.f.read(126)
        header = buff[:110]
        patchInfo = buff[110:]
        tmp = header.split(b'\x1e')
        
        self.quantity = tmp[1].decode('utf-8').replace('\x00','').strip(' ')
        self.shortName = tmp[3].decode('utf-8').replace('\x00','').strip(' ')
        self.units = tmp[5].decode('utf-8').replace('\x00','').strip(' ')

        buff = np.frombuffer(patchInfo, dtype=np.int32, count=4)
        self.nPatches = buff[2]

        for n in range(self.nPatches):
            self.patchData[n]={}
            buff = np.frombuffer(self.f.read(44), dtype=np.int32, offset=4)
            (dx,dy,dz) = (buff[1]-buff[0],buff[3]-buff[2],buff[5]-buff[4])
            if abs(buff[6]) == 1:
                dy = dy+1
                dz = dz+1
                data_count = dy*dz
            elif abs(buff[6]) == 2:
                dx = dx+1
                dz = dz+1
                data_count = dx*dz
            elif abs(buff[6]) == 3:
                dx = dx+1
                dy = dy+1
                data_count = dx*dy
            # nPts = np.max([dx,1])*np.max([dy,1])*np.max([dz,1])+2
            self.patchData[n]['dims']=[dx,dy,dz,
                    buff[0],buff[1],buff[2],
                    buff[3],buff[4],buff[5]]
            self.patchData[n]['count']=data_count
            self.patchData[n]['IOR']=buff[6]
            self.patchData[n]['OBST_ID']=buff[7]
            self.patchData[n]['OBST_NM']=buff[8]
        self.total_value_count=sum([value['count'] for value in self.patchData.values() if 'count' in value])


    # Return array of data times from cartesian BNDF file
    def getTimes(self):
        self.f.seek(0)
        self.readHeader()
        buff = self.f.read()
        if len(buff) % 4 != 0:
            buff = np.frombuffer(buff[:int(len(buff)/4)*4], dtype=np.float32)
        else:
            buff = np.frombuffer(buff, dtype=np.float32)
        padding = 2*self.nPatches+3
        self.times = buff[1::self.total_value_count+padding]
        self.f.seek(0)

    # Return a single time and data array from a cartesian BNDF file
    def getRecord(self):

        # if header is not read
        if self.nPatches is None:
            self.readHeader()
        # if at the start of the file
        if self.f.tell()<(126+44*self.nPatches):
            self.readHeader()

        _ = _ = self.f.read(4)
        self.currentTime = np.frombuffer(self.f.read(4), dtype=np.float32)
        _ = _ = self.f.read(4)
        for n in range(self.nPatches):
             _ = self.f.read(4)
             self.patchData[n]['values'] = np.frombuffer(self.f.read(self.patchData[n]['count']*4), dtype=np.float32)
             _ = self.f.read(4)
       
# class to control parsing of unstructured (GEOM) boundary files
class BNDE:
    def __init__(self,file):
        self.f=open(file,'rb')
        self.cc=True
        self.nVerts=None
        self.nFaces=None
        self.verts=None
        self.faceVerts=None
        self.faceCenters=None
        self.values=None
        self.times=None
        self.currentTime=None

        # always build geometry when creating class instance
        self.getGeometry()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.f.close()

    # get geometry for relevant unstructured BNDE file (vertex locations and face center locations)
    def getGeometry(self):
        gcfFile=open(self.f.name[:self.f.name.rfind('_')]+'.gcf','rb')
        buff = np.frombuffer(gcfFile.read(76), dtype=np.int32)
        self.nVerts = buff[15]
        self.nFaces = buff[16]
        if (self.nVerts>0) and (self.nFaces>0):
            _ = gcfFile.read(4)
            buff = np.frombuffer(gcfFile.read(self.nVerts*4*3),dtype=np.float32)
            _ = gcfFile.read(4)
            self.verts = np.reshape(buff,(self.nVerts,3))
            _ = gcfFile.read(4)
            buff = np.frombuffer(gcfFile.read(self.nFaces*4*3),dtype=np.int32)
            self.faceVerts = np.reshape(buff,(self.nFaces,3))
            _ = gcfFile.read(4)

            self.faceCenters = np.empty(self.faceVerts.shape)
            for iaxis in range(3):
                self.faceCenters[:,iaxis]=self.verts[self.faceVerts-1,iaxis].sum(axis=1)/3

    # get list of time stamps for unstructured BNDE data
    def getTimes(self):
        self.f.seek(0)
        buff = self.f.read()
        nvv,nfv = np.frombuffer(buff[:56], dtype=np.int32)[12:14]
        if len(buff) % 4 != 0:
            buff = np.frombuffer(buff[:int(len(buff)/4)*4], dtype=np.float32)
        else:
            buff = np.frombuffer(buff, dtype=np.float32)
        # assumes there are only EITHER vertex values or face values
        if nvv>0:
            padding = nvv+11
        else:
            padding = nfv+11
        self.times = buff[7::padding]
        self.f.seek(0)

    # Return a single time and data array from an unstructured BNDE file
    def getRecord(self):
        # if at start of file, skip over header data
        if self.f.tell()<(24):
            self.f.seek(0)
            self.f.read(24)

        _ = self.f.read(4)
        self.currentTime = np.frombuffer(self.f.read(4), dtype=np.float32)
        _ = self.f.read(8)
        nvv,nfv = np.frombuffer(self.f.read(16), dtype=np.int32)[2:4]
        _ = self.f.read(4)
        if nvv>0:
            self.cc = False
            _ = self.f.read(4)
            self.values = np.frombuffer(self.f.read(nvv*4), dtype=np.float32)
            _ = self.f.read(4)
        else:
            self.cc = True
            _ = self.f.read(4)
            self.values = np.frombuffer(self.f.read(nfv*4), dtype=np.float32)
            _ = self.f.read(4)

# function to parse CHID.out for ORIGIN_LAT and ORIGIN_LON
def parseOUT(file):
    # values if not defined
    lat,lon=-1e6,-1e6
    with open(file) as f:
        while True:
            line=f.readline()
            # found geographic info
            if line[1:11]=='Geographic':
                line=f.readline()
                line=f.readline()
                lat=float(line[20:])
                line=f.readline()
                lon=float(line[20:])
                break
            # no geographic info found, stop reading file now
            if line[1:14]=='Miscellaneous':
                break

    return lat,lon

# parse SMV file, in this case to find grid and SLCT information
def parseSMV(file):
    with open(file) as f:
        linesSMV = f.readlines()

        grids=[]
        BNDFfiles=defaultdict(bool)
        BNDEfiles=defaultdict(bool)
        for il in range(len(linesSMV)):
            if ('GRID' in linesSMV[il]):
                gridTRNX, gridTRNY, gridTRNZ = parseGRID(linesSMV, il)
                grids.append([gridTRNX.copy(),
                              gridTRNY.copy(),
                              gridTRNZ.copy()])

            if (('BNDF' in linesSMV[il]) or ('BNDC' in linesSMV[il])):
                file = '%s.bf'%(linesSMV[il+1][1:].split('.bf')[0])
                BNDFfiles[file] = defaultdict(bool)
                BNDFfiles[file]['QUANTITY'] = linesSMV[il+2].strip()
                BNDFfiles[file]['SHORTNAME'] = linesSMV[il+3].strip()
                BNDFfiles[file]['UNITS'] = linesSMV[il+4].strip()
                BNDFfiles[file]['LINETEXT'] = linesSMV[il]
                BNDFfiles[file]['MESH']=int(BNDFfiles[file]['LINETEXT'][5:10])
                if (BNDFfiles[file]['LINETEXT'][3]=='C'):
                    BNDFfiles[file]['CELL_CENTERD']=True
                else:
                    BNDFfiles[file]['CELL_CENTERD']=False

            if ('BNDE' in linesSMV[il]):
                file = '%s.be'%(linesSMV[il+1][1:].split('.be')[0])
                BNDEfiles[file] = defaultdict(bool)
                BNDEfiles[file]['QUANTITY'] = linesSMV[il+3].strip()
                BNDEfiles[file]['SHORTNAME'] = linesSMV[il+4].strip()
                BNDEfiles[file]['UNITS'] = linesSMV[il+5].strip()
                BNDEfiles[file]['LINETEXT'] = linesSMV[il]
                BNDEfiles[file]['MESH']=int(BNDEfiles[file]['LINETEXT'][5:10])
                # currently the only option
                BNDEfiles[file]['CELL_CENTERD']=True


    return BNDFfiles,BNDEfiles,grids

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
    # check for z stretching
    if (int(zlines[0][3:5])==-1):
        delta=(zind1-zind0-2)
        zlines = lines[zind0+delta:zind1+delta]
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


