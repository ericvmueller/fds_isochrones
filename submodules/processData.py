import processing
import numpy as np
import pyfdstools as pfds
from collections import defaultdict
from qgis.gui import QgsMapCanvas
from qgis.core import (
    QgsProcessing,
    QgsProject,
    QgsFeature,
    QgsGeometry,
    QgsPoint,
    QgsPointXY,
    QgsVectorLayer)


def slct2contour(context, feedback, CHID, fds_path, QUANTITY, threshold, t_step, crs, xy_offset):
    """

    Parameters
    ----------
    CHID : str
        run name
    CHID : str
        folder containing run data
    QUANTITY : str
        SCLT quantity
    threshold : float
        threshold value to denote fire arrival
    t_step : float
        time increment between isochrones
    crs: str
        code for coordinate reference system (e.g. EPSG:5070 for NAD83/Conus Albers)
    offset (optional) : QgisPointXY
        [x_o, y_o] adjustment of domain to CRS units, if necessary

    Returns
    -------
    None.

    borrows from /fds/Utilities/Python/scripts

    """

    # parse SMV file
    linesSMV = pfds.zreadlines(fds_path+'/'+CHID+'.smv')

    grids=[]
    SLCTfiles=defaultdict(bool)
    for il in range(len(linesSMV)):
        if ('GRID' in linesSMV[il]):
            gridTRNX, gridTRNY, gridTRNZ = pfds.parseGRID(linesSMV, il)
            grids.append([gridTRNX.copy(),
                          gridTRNY.copy(),
                          gridTRNZ.copy()])

        if ('SLCT' in linesSMV[il]):
            file = '%s.sf'%(linesSMV[il+1][1:].split('.sf')[0])
            feedback.pushInfo(file)
            SLCTfiles[file] = defaultdict(bool)
            SLCTfiles[file]['QUANTITY'] = linesSMV[il+2].strip()
            SLCTfiles[file]['SHORTNAME'] = linesSMV[il+3].strip()
            SLCTfiles[file]['UNITS'] = linesSMV[il+4].strip()
            SLCTfiles[file]['LINETEXT'] = linesSMV[il]
            SLCTfiles[file]['MESH']=int(SLCTfiles[file]['LINETEXT'][5:10])
            SLCTfiles[file]['AGL']=float(SLCTfiles[file]['LINETEXT'][11:20])

    append=False
    maxval=0
    for SLCTfile in SLCTfiles:

        timesSLCF = pfds.readSLCFtimes(fds_path+'/'+SLCTfile)
        times = []
        f = pfds.zopen(fds_path+'/'+SLCTfile)

        qty, sName, uts, iX, eX, iY, eY, iZ, eZ = pfds.readSLCFheader(f)

        # export if correct quantity
        if (qty == QUANTITY):
            (NX, NY, NZ) = (eX-iX, eY-iY, eZ-iZ)
            feedback.pushInfo('Reading from '+SLCTfile+'...')

            shape = (NX+1, NY+1, NZ+1)
            # allocate array to store arrival time
            arrival_data = -1.0*np.ones(shape)
            NT = len(timesSLCF)
            time = np.zeros(NT)
            for i in range(0, NT):
                t, tmp = pfds.readNextTime(f, NX, NY, NZ)
                tmp = np.reshape(tmp, shape, order='F')
                # set arrival time for any points that reach arrival condition
                arrival_data[(tmp >= threshold) & (arrival_data<0)] = t


            # all SCLT should be 2D (x and y)
            arrival_data=np.squeeze(arrival_data)
            # anywhere that the fire arrive gets max val (to make nice contours)
            arrival_data[arrival_data<0]=t

            #create temporary layer containing sample points
            if not append:
                uri='Point?crs='+crs+'&field=id:integer&field=time:double&index=yes'
                pointLayer=QgsVectorLayer(uri, 'fds_sample_points', 'memory')
                pointLayer.startEditing()
                append=True

            pointLayer = _createPointLayer(
                feedback,arrival_data,grids[SLCTfiles[SLCTfile]['MESH']-1],pointLayer,xy_offset)

    f.close()

    pointLayer.commitChanges()
    # for debug, add layer of fds sample points to map
    QgsProject.instance().addMapLayer(pointLayer)
    # layers=QgsProject.instance().layerTreeRoot().findLayerIds()
    # QgsProject.instance().layerTreeRoot().findLayer(pointLayer).setItemVisibilityChecked(False)
    # context.temporaryLayerStore().addMapLayer(pointLayer)
    contourOutput=_createContourLayer(pointLayer.source(),t_step)

    return contourOutput


# get a vector file of points from a 2D numpy array
def _createPointLayer(feedback,data,grid,pointLayer,xy_offset):

    point=QgsFeature()

    feedback.pushInfo('Extracting SCLT data points...')
    total=len(grid[0])*len(grid[1])
    count=0
    for ir in grid[0]:
        for jr in grid[1]:
            # Stop the algorithm if cancel button has been clicked
            if feedback.isCanceled():
                break

            count=count+1
            point.setGeometry(QgsGeometry.fromPointXY(
                QgsPointXY(ir[1]+xy_offset.x(),jr[1]+xy_offset.y())))
            i,j=int(ir[0]),int(jr[0])
            time=data[i,j]
            point.setAttributes([count,time.item()])
            pointLayer.addFeatures([point])

            # Update the progress bar
            feedback.setProgress(int(100*count/total))



    return pointLayer

# get a vector file of points from a 2D numpy array
def _createContourLayer(inputSource,t_step):

    return  processing.run("contourplugin:generatecontours",
            {'InputLayer':inputSource,
            'InputField':'"time"',
            'DuplicatePointTolerance':0,
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
