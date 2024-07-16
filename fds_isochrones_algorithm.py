# -*- coding: utf-8 -*-

"""
/***************************************************************************
 fdsIsochrones
                                 A QGIS plugin
 create isochrones of fire progression from FDS outputs
 Generated by Plugin Builder: http://g-sherman.github.io/Qgis-Plugin-Builder/
                              -------------------
        begin                : 2023-06-27
        copyright            : (C) 2023 by Eric Mueller
        email                : ericvmueller@gmail.com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

__author__ = 'Eric Mueller'
__date__ = '2023-06-27'
__copyright__ = '(C) 2023 by Eric Mueller'

# This will get replaced with a git SHA1 when you do a git archive

__revision__ = '$Format:%H$'

from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import *
import os
from .submodules import processData,fileParse


DEFAULTS = {
    'CHID': '',
    'fds_path': './',
    't_step': '10',
    'threshold': '0',
    'QUANTITY': '0',
    'dateTime': '',
    'origin': 'None',
    'dx_out': '',
    'samplePoints': 'False'
}

quants = {'0':'TIME OF ARRIVAL','1':'LEVEL SET VALUE'}

class fdsIsochronesAlgorithm(QgsProcessingAlgorithm):

    # Constants used to refer to parameters and outputs

    CHID = 'CHID'
    fds_path = 'fds_path'
    QUANTITY = 'QUANTITY'
    threshold = 'threshold'
    t_step = 't_step'
    OUTPUT = 'OUTPUT'
    crs = 'crs'
    origin = 'origin'
    dx_out = 'dx_out'
    dateTime = 'dateTime'
    samplePoints = 'samplePoints'

    def initAlgorithm(self, config):
        """
        define the inputs and output of the algorithm
        """
        project = QgsProject.instance()

        # Define parameter: CHID
        defaultValue, _ = project.readEntry('fds_isochrones', 'CHID', DEFAULTS['CHID'])
        self.addParameter(
            QgsProcessingParameterString(
                name=self.CHID,
                description='FDS case identifier (CHID)',
                multiLine=False,
                defaultValue=defaultValue
            )
        )

        # Define parameter: fds_path
        defaultValue, _ = project.readEntry('fds_isochrones', 'fds_path', DEFAULTS['fds_path'])
        self.addParameter(
            QgsProcessingParameterFile(
                name=self.fds_path,
                description='Working directory',
                behavior=QgsProcessingParameterFile.Folder,
                fileFilter='All files (*.*)',
                defaultValue=defaultValue
            )
        )

        # Define parameter: crs
        defaultValue, _ = project.readEntry('fds_isochrones', 'crs', project.crs().authid())
        self.addParameter(
            QgsProcessingParameterCrs(
                name=self.crs,
                description='CRS for FDS run',
                defaultValue=defaultValue
            )
        )

        # Define parameter: origin
        defaultValue, _ = project.readEntry('fds_isochrones', 'origin', DEFAULTS['origin'])
        print(defaultValue)
        self.addParameter(
            QgsProcessingParameterPoint(
                name=self.origin,
                description='FDS origin in CRS (will be obtained from CHID.out if available)',
                defaultValue=defaultValue,
                optional=True
            )
        )

        # Define parameter: QUANTITY
        defaultValue, _ = project.readEntry('fds_isochrones', 'QUANTITY', DEFAULTS['QUANTITY'])
        self.addParameter(
            QgsProcessingParameterEnum(
                name=self.QUANTITY,
                description='AGL Slice QUANTITY',
                options=['TIME OF ARRIVAL','LEVEL SET VALUE'],
                defaultValue=defaultValue
            )
        )

        # Define parameter: threshold
        defaultValue, _ = project.readEntry('fds_isochrones', 'threshold', DEFAULTS['threshold'])
        self.addParameter(
            QgsProcessingParameterNumber(
                name=self.threshold,
                description='Threshold QUANTITY at fire arrival (e.g. 0 for LEVEL SET VALUE)',
                type=QgsProcessingParameterNumber.Double,
                defaultValue=defaultValue,
                optional=False
            )
        )

        # Define parameter: t_step
        defaultValue, _ = project.readEntry('fds_isochrones', 't_step', DEFAULTS['t_step'])
        self.addParameter(
            QgsProcessingParameterNumber(
                name=self.t_step,
                description='Contour intervals (s)',
                type=QgsProcessingParameterNumber.Double,
                defaultValue=defaultValue,
                optional=False,
                minValue=0.1
            )
        )

        # Define parameter: dx_out
        defaultValue, _ = project.readEntry('fds_isochrones', 'dx_out', DEFAULTS['dx_out'])
        self.addParameter(
            QgsProcessingParameterNumber(
                name=self.dx_out,
                description='Filtered output resolution (m)',
                type=QgsProcessingParameterNumber.Double,
                # defaultValue=defaultValue,
                optional=True,
                minValue=0.01
            )
        )

        # Define parameter: dateTime
        defaultValue, _ = project.readEntry('fds_isochrones', 'dateTime', DEFAULTS['dateTime'])
        self.addParameter(
            QgsProcessingParameterDateTime(
                name=self.dateTime,
                description='Ignition datetime for animation',
                # defaultValue=defaultValue,
                optional=True
            )
        )

        # Define parameter: samplePoints
        defaultValue, _ = project.readEntry('fds_isochrones', 'samplePoints', DEFAULTS['samplePoints'])
        self.addParameter(
            QgsProcessingParameterBoolean(
                name=self.samplePoints,
                description='Save FDS sample points as temporary layer',
                defaultValue=False
            )
        )

        # Add a feature sink in which to store processed features
        self.addParameter(
            QgsProcessingParameterFeatureSink(
                name=self.OUTPUT,
                description='Output layer'
            )
        )

    def processAlgorithm(self, parameters, context, feedback):
        """
        Process the algorithm
        """

        results, outputs, project = {}, {}, QgsProject.instance()

        # Check project crs and save it
        if not project.crs().isValid():
            raise QgsProcessingException(
                f'Project CRS <{project.crs().description()}> is not valid, cannot proceed.'
            )
        project.writeEntry('fds_isochrones', 'project_crs', project.crs().authid())

        # Get parameter: chid
        CHID = self.parameterAsString(parameters, 'CHID', context)
        if not CHID:
            raise QgsProcessingException(self.invalidSourceError(parameters, 'CHID'))
        project.writeEntry('fds_isochrones', 'CHID', CHID)

        # Get parameter: fds_path
        project_path = project.readPath('./')
        if not project_path:
            raise QgsProcessingException(
                'Save the qgis project to disk, cannot proceed.'
            )
        fds_path = self.parameterAsFile(parameters, 'fds_path', context)
        if not fds_path:
            raise QgsProcessingException(
                self.invalidSourceError(parameters, 'fds_path')
            )
        project.writeEntry('fds_isochrones', 'fds_path', fds_path)
        fds_path = os.path.join(project_path, fds_path)  # make abs

        # Get parameter: crs
        crs = self.parameterAsCrs(parameters, 'crs', context)
        if not crs:
            raise QgsProcessingException(
                self.invalidSourceError(crs, 'crs')
            )
        project.writeEntry('fds_isochrones', 'crs', crs.authid())

        # Get parameter: origin
        origin = self.parameterAsPoint(parameters, 'origin', context, crs)
        print('write '+str(origin))
        project.writeEntry('fds_isochrones', 'origin', str(origin))

        # Get parameter: QUANTITY
        QUANTITY = self.parameterAsString(parameters, 'QUANTITY', context)
        project.writeEntry('fds_isochrones', 'QUANTITY', QUANTITY)

        # Get parameter: threshold
        threshold = self.parameterAsDouble(parameters, 'threshold', context)
        project.writeEntry('fds_isochrones', 'threshold', str(threshold))

        # Get parameter: t_step
        t_step = self.parameterAsDouble(parameters, 't_step', context)
        project.writeEntry('fds_isochrones', 't_step', str(t_step))

        # Get parameter: dx_out
        dx_out = self.parameterAsDouble(parameters, 'dx_out', context)
        project.writeEntry('fds_isochrones', 'dx_out', str(dx_out))

        # Get parameter: dateTime
        dateTime = self.parameterAsDateTime(parameters, 'dateTime', context)
        project.writeEntry('fds_isochrones', 'dateTime', str(dateTime))

        # Get parameter: samplePoints
        samplePoints = self.parameterAsBoolean(parameters, 'samplePoints', context)
        project.writeEntry('fds_isochrones', 'samplePoints', str(samplePoints))

        # parse OUT file for geographic information
        origin_lat,origin_lon=fileParse.parseOUT(fds_path+'/'+CHID+'.out')
        # if FDS origin successfully determined, shift the data
        if ((origin_lat>-1e6) and (origin_lon>-1e6)):
            # convert from WGS84 to relevant CRS
            wgs84_crs = QgsCoordinateReferenceSystem("EPSG:4326")
            wgs84_to_fds_tr = QgsCoordinateTransform(
                        wgs84_crs, crs, project)
            xy_offset = QgsPoint(x=origin_lon,y=origin_lat)
            xy_offset.transform(wgs84_to_fds_tr)

        else:
            # xy_offset = QgsPoint(x=0,y=0)
            xy_offset = origin

        contourLayer,maxTime=processData.slct2contour(
            feedback,
            CHID,
            fds_path,
            quants[QUANTITY],
            threshold,
            t_step,
            dx_out,
            crs,
            xy_offset,
            dateTime,
            samplePoints)

        (sink, dest_id) = self.parameterAsSink(parameters, self.OUTPUT,
                context, contourLayer.fields(), contourLayer.wkbType(),
                contourLayer.sourceCrs())

        # Copy over features to the sink (not sure if most effective way to save)
        total = 100.0 / contourLayer.featureCount() if contourLayer.featureCount() else 0
        features = contourLayer.getFeatures()

        for current, feature in enumerate(features):
            # Stop the algorithm if cancel button has been clicked
            if feedback.isCanceled():
                break

            # Add a feature in the sink
            sink.addFeature(feature, QgsFeatureSink.FastInsert)

            # Update the progress bar
            feedback.setProgress(int(current * total))

        # save information for updating symbology
        self.dest_id=dest_id
        self.intervals=round(maxTime/t_step)

        # Return the results of the algorithm
        return {self.OUTPUT: self.dest_id}

    def postProcessAlgorithm(self, context, feedback):
        """
        PostProcessing to define the Symbology
        """
        # identify the output map layer
        layer=QgsProcessingUtils.mapLayerFromString(self.dest_id,context)

        # classify the features according to the time, using fixed t_step interval
        renderer = QgsGraduatedSymbolRenderer('time')
        renderer.setClassificationMethod(QgsClassificationEqualInterval())
        renderer.updateClasses(layer, renderer.mode(), self.intervals)
        # set color ramp to Magma
        default_style = QgsStyle().defaultStyle()
        color_ramp = default_style.colorRamp('Magma')
        renderer.updateColorRamp(color_ramp)
        # increase line width
        renderer.setSymbolSizes(.5,.5)

        # format the labels
        frmt = QgsRendererRangeLabelFormat()
        frmt.setFormat("%1 - %2")
        frmt.setPrecision(-1)
        frmt.setTrimTrailingZeroes(True)
        renderer.setLabelFormat(frmt)

        layer.setRenderer(renderer)
        layer.triggerRepaint()

        return {self.OUTPUT: self.dest_id}

    def name(self):
        """
        Returns the algorithm name, used for identifying the algorithm
        """
        return 'FDS Isochrones'

    def displayName(self):
        """
        Returns the translated algorithm name, which should be used for any
        user-visible display of the algorithm name.
        """
        return self.tr(self.name())

    def group(self):
        """
        Returns the name of the group this algorithm belongs to
        """
        return self.tr(self.groupId())

    def groupId(self):
        """
        Returns the unique ID of the group this algorithm belongs to
        """
        return ''

    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return fdsIsochronesAlgorithm()
