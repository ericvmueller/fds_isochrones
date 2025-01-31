"""
/***************************************************************************
 processData module

 This module contains a set of functions used to process relevant data from
 FDS outputs and create GIS features. Some parsing functionality is
 borrowed from pyfdstools: https://github.com/johodges/pyfdstools,
 but has been adapted for simplicity and to reduce requirements for
 external dependencies.
"""

import processing
import numpy as np
from scipy import interpolate,signal
from .fileParse import BNDE,BNDF,parseSMV,parseOUT
from collections import defaultdict
from qgis.PyQt.QtCore import QVariant
from qgis.gui import QgsMapCanvas
from qgis.core import (
    QgsProcessing,
    QgsProcessingException,
    QgsProject,
    QgsFeature,
    QgsField,
    QgsGeometry,
    QgsPointXY,
    QgsVectorLayer)


def strip_ghost_cells(dat,d1,d2):
    dat=np.reshape(dat,(d1,d2), order='F')
    dat=dat[:-1,:-1]
    d1-=1
    d2-=1
    return dat.flatten(order='F'),d1,d2

def bndf2contour(
    feedback,
    CHID,
    fds_path,
    QUANTITY,
    threshold,
    t_step,
    dx_out,
    crs,
    xy_offset,
    dateTime,
    samplePoints):
    """

    Parameters
    ----------
    feedback : QgsProcessingFeedback
        class to provide feedback
    CHID : str
        FDS job name
    fds_path : str
        folder containing run data
    QUANTITY : str
        SCLT quantity
    threshold : float
        threshold value to denote fire arrival
    t_step : float
        time increment between isochrones
    t_step : float
        resolution to resample output
    crs: str
        code for coordinate reference system (e.g. EPSG:5070 for NAD83/Conus Albers)
    offset : QgisPoint
        [x_o, y_o] adjustment of domain to CRS units
    dateTime : QDateTime
        Date and time to be added for temporal animation
    samplePoints: bool
        save FDS sample points as temporary layer

    Returns
    -------
    contourOutput : QgsVectorLayer
        vector layer with isochrone contours as features
    maxTime : float
        maximum time of fire spread (used for symbology)

    """

    # parse SMV file for BNDF/E files and grid information
    [BNDFfiles,BNDEfiles,grids]=parseSMV(fds_path+'/'+CHID+'.smv')

    x=[]
    y=[]
    z=[]
    dx=[]
    dy=[]
    t_a=[]

    bndf_found=False

    for BNDFfile in BNDFfiles:
        if (BNDFfiles[BNDFfile]['QUANTITY']=='FIRE ARRIVAL TIME'):
            bndf_found=True
            feedback.pushInfo('Reading from '+BNDFfile+'...')
            with BNDF(fds_path+'/'+BNDFfile) as bndf:
                bndf.getTimes()
                NT = len(bndf.times)
                for i in range(NT):
                    bndf.getRecord()
                    
                for n in range(bndf.nPatches):

                    nm=bndf.patchData[n]['OBST_NM']
                    cc=BNDFfiles[BNDFfile]['CELL_CENTERD']
                    
                    xp=grids[nm-1][0][bndf.patchData[n]['dims'][3]:
                                      bndf.patchData[n]['dims'][4]+1][:,1]
                    yp=grids[nm-1][1][bndf.patchData[n]['dims'][5]:
                                      bndf.patchData[n]['dims'][6]+1][:,1]
                    zp=grids[nm-1][2][bndf.patchData[n]['dims'][7]:
                                      bndf.patchData[n]['dims'][8]+1][:,1]
                    lxp=len(xp)
                    lyp=len(yp)
                    lzp=len(zp)
            
                    tmp=bndf.patchData[n]['values']            
                                
                    # determine the patch orientation and build the point list
                    if (bndf.patchData[n]['IOR']==3):
                    # only worry about IOR=3 for now                                  
                        if cc:
                            tmp,lxp,lyp=strip_ghost_cells(tmp,lxp,lyp)
                            xp=xp[:-1]
                            yp=yp[:-1]

                        tmp_dx=(xp[1:]-xp[:-1])
                        dxp=0*xp
                        dxp[1:-1]=0.5*(tmp_dx[1:]+tmp_dx[:-1])
                        dxp[0]=0.5*tmp_dx[0]
                        dxp[-1]=0.5*tmp_dx[-1]

                        tmp_dx=(yp[1:]-yp[:-1])
                        dyp=0*yp
                        dyp[1:-1]=0.5*(tmp_dx[1:]+tmp_dx[:-1])
                        dyp[0]=0.5*tmp_dx[0]
                        dyp[-1]=0.5*tmp_dx[-1]

                        xp=np.tile(xp,lyp)
                        dxp=np.tile(dxp,lyp)
                        yp=np.repeat(yp,lxp)
                        dyp=np.repeat(dyp,lxp)
                        zp=np.repeat(zp,lxp*lyp)
                    
                        x.append(xp)
                        y.append(yp)
                        z.append(zp)
                        dx.append(dxp)
                        dy.append(dyp)
                        t_a.append(tmp)
                
    for BNDEfile in BNDEfiles:
        if (BNDEfiles[BNDEfile]['QUANTITY']=='FIRE ARRIVAL TIME'):
            bndf_found=True
            feedback.pushInfo('Reading from '+BNDEfile+'...')
            with BNDE(fds_path+'/'+BNDEfile) as bnde:
                bnde.getTimes()
                NT = len(bnde.times)
                for i in range(NT):
                    bnde.getRecord()
                    
                x.append(bnde.faceCenters[:,0])
                y.append(bnde.faceCenters[:,1])
                z.append(bnde.faceCenters[:,2])
                t_a.append(bnde.values)

    if not bndf_found:
        raise QgsProcessingException('ERROR: No relevant boundary files found')

    # strip no arrival values
    x=np.concatenate(x)
    y=np.concatenate(y)
    dx=np.concatenate(dx)
    dy=np.concatenate(dy)
    t_a=np.concatenate(t_a)
    x=x[t_a<9E5]
    y=y[t_a<9E5]
    dx=dx[t_a<9E5]
    dy=dy[t_a<9E5]
    t_a=t_a[t_a<9E5]
    maxTime = max(t_a)

    #create temporary layer containing sample points
    uri='Point?crs='+crs.authid()+'&field=id:integer&field=time:double&index=yes'
    pointLayer=QgsVectorLayer(uri, 'fds_sample_points', 'memory')
    pointLayer.startEditing()

    # do resampling
    if not dx_out==0:
        feedback.pushInfo('Resampling to output resolution...')
        x_i = np.arange(min(x),max(x),dx_out)
        y_i = np.arange(min(y),max(y),dx_out)
        [xxi,yyi]=np.meshgrid(x_i,y_i)
        t_ai=interpolate.griddata((x,y),t_a,(xxi,yyi))
        # box filter
        box = 1/9*np.ones((3,3))
        t_ai=signal.convolve2d(t_ai, box, boundary='symm', mode='same')
        x=xxi[~np.isnan(t_ai)]
        y=yyi[~np.isnan(t_ai)]
        t_ai=t_ai[~np.isnan(t_ai)]

        pointLayer = _addLayerPoints(feedback,x,y,t_ai,pointLayer,xy_offset)

    else:
        pointLayer = _addLayerPoints(feedback,x,y,t_a,pointLayer,xy_offset)


    pointLayer.commitChanges()
    # optional, add layer of fds sample points to map
    QgsProject.instance().addMapLayer(pointLayer,addToLegend=samplePoints)
    # calculate contours from point data
    contourOutput=_createContourLayer(pointLayer.source(),t_step)

    if not samplePoints:
        QgsProject.instance().removeMapLayer(pointLayer)

    contourOutput.startEditing()
    # contour remove unused attribute from contour plugin
    contourOutput.deleteAttribute(2)

    # add burned area estimation
    contourOutput = _addBurnedArea(feedback,contourOutput,dx,dy,t_a)

    # add datetime field for animation purposes
    if not dateTime.isNull():
        contourOutput = _addDateTime(feedback,contourOutput,dateTime)

    contourOutput.commitChanges()
    contourOutput.rollBack()

    return contourOutput,maxTime


# add to vector file of points from a 2D numpy array
def _addLayerPoints(feedback,x,y,data,pointLayer,xy_offset):

    # point used to populate layer
    point=QgsFeature()
    feedback.pushInfo('Creating virtual point layer...')
    total=len(data)
    count=0
    # loop through data array and get location from 'grid'
    for (xi,yi,ti) in zip(x,y,data):
        # Stop the algorithm if cancel button has been clicked
        if feedback.isCanceled():
            break
        count=count+1
        # shift FDS location relative to global origin
        point.setGeometry(QgsGeometry.fromPointXY(
            QgsPointXY(xi+xy_offset.x(),yi+xy_offset.y())))
        point.setAttributes([count,float(np.round(ti,3))])
        pointLayer.addFeatures([point])

        # Update the progress bar
        feedback.setProgress(int(100*count/total))

    return pointLayer

# get a vector file of points from a 2D numpy array
def _createContourLayer(inputSource,t_step):

    return  processing.run("contourplugin:generatecontours",
            {'InputLayer':inputSource,
            'InputField':'"time"',
            'DuplicatePointTolerance':0.001,
            'ContourType':0,
            'ExtendOption':None,
            'ContourMethod':3,
            'NContour':100000,
            'MinContourValue':None,
            'MaxContourValue':None,
            'ContourInterval':t_step,
            'ContourLevels':'',
            'LabelDecimalPlaces':-1,
            'LabelTrimZeros':False,
            'LabelUnits':'',
            'OutputLayer':'TEMPORARY_OUTPUT'})['OutputLayer']

# add datetime attribute to contour layer
def _addDateTime(feedback,contourLayer,dateTime):

    contourLayer.addAttribute(QgsField('datetime',QVariant.DateTime))
    for feature in contourLayer.getFeatures():
        # Stop the algorithm if cancel button has been clicked
        if feedback.isCanceled():
            break
        # add FDS time to igniton datetime
        feature['datetime'] = dateTime.addMSecs(round(1000*feature['time']))
        contourLayer.updateFeature(feature)

    return contourLayer


# add datetime attribute to contour layer
def _addBurnedArea(feedback,contourLayer,dx,dy,t_a):

    contourLayer.addAttribute(QgsField('burned_area',QVariant.Double))
    for feature in contourLayer.getFeatures():
        # Stop the algorithm if cancel button has been clicked
        if feedback.isCanceled():
            break
        # add FDS time to igniton datetime
        feature['burned_area'] = float(sum(dx[t_a<=feature['time']]*dy[t_a<=feature['time']]))
        contourLayer.updateFeature(feature)

    return contourLayer