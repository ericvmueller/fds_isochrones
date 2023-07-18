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
from .fileParse import SLCT,parseSMV,parseOUT
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


def slct2contour(
    feedback,
    CHID,
    fds_path,
    QUANTITY,
    threshold,
    t_step,
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

    # parse SMV file for SLCT and grid information
    [SLCTfiles,grids]=parseSMV(fds_path+'/'+CHID+'.smv')

    # scan SCLT files and extract relevant data
    append=False
    maxval=0
    slct_found=False
    for SLCTfile in SLCTfiles:
        with SLCT(fds_path+'/'+SLCTfile) as slct:
            slct.readHeader()

            # export if correct quantity
            if (slct.quantity == QUANTITY):
                slct_found=True
                feedback.pushInfo('Reading from '+SLCTfile+'...')
                # get data shape
                (NX, NY, NZ) = (slct.eX-slct.iX, slct.eY-slct.iY, slct.eZ-slct.iZ)
                shape = (NX+1, NY+1, NZ+1)
                # allocate array to store arrival time
                arrival_data = -1.0*np.ones(shape)
                # get output times
                slct.readTimes()
                NT = len(slct.times)
                time = np.zeros(NT)
                # for each time identify new data points which experience fire arrival
                for i in range(0, NT):
                    slct.readRecord()
                    tmp = np.reshape(slct.data, shape, order='F')
                    # use threshold to identify arrival
                    if not QUANTITY=='TIME OF ARRIVAL':
                        arrival_data[(tmp >= threshold) & (arrival_data<0)] = slct.currentTime

                # Time of arrival can be taken directly
                if QUANTITY=='TIME OF ARRIVAL':
                    arrival_data[:] = tmp

                # all SCLT should be 2D (x and y)
                arrival_data=np.squeeze(arrival_data)
                # anywhere that the fire does not arrive gets max val (to make nice contours)
                maxTime=slct.currentTime.item()
                arrival_data[arrival_data<0]=maxTime

                #create temporary layer containing sample points
                if not append:
                    uri='Point?crs='+crs.authid()+'&field=id:integer&field=time:double&index=yes'
                    pointLayer=QgsVectorLayer(uri, 'fds_sample_points', 'memory')
                    pointLayer.startEditing()
                    append=True

                pointLayer = _addLayerPoints(
                    feedback,arrival_data,grids[SLCTfiles[SLCTfile]['MESH']-1],pointLayer,xy_offset)

    if not slct_found:
        raise QgsProcessingException('ERROR: No relevant slice files found')

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
    # add datetime field for animation purposes
    if not dateTime.isNull():
        contourOutput = _addDateTime(feedback,contourOutput,dateTime)

    contourOutput.commitChanges()
    contourOutput.rollBack()

    return contourOutput,maxTime


# add to vector file of points from a 2D numpy array
def _addLayerPoints(feedback,data,grid,pointLayer,xy_offset):

    # point used to populate layer
    point=QgsFeature()
    feedback.pushInfo('Extracting SCLT data points...')
    total=len(grid[0])*len(grid[1])
    count=0
    # loop through data array and get location from 'grid'
    for ir in grid[0]:
        for jr in grid[1]:
            # Stop the algorithm if cancel button has been clicked
            if feedback.isCanceled():
                break
            count=count+1
            # shift FDS location relative to global origin
            point.setGeometry(QgsGeometry.fromPointXY(
                QgsPointXY(ir[1]+xy_offset.x(),jr[1]+xy_offset.y())))
            i,j=int(ir[0]),int(jr[0])
            time=data[i,j].item()
            point.setAttributes([count,time])
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