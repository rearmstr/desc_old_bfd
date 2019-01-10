#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2013 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
import numpy as np
import random
import scipy.spatial
from contextlib import contextmanager

import lsst.daf.base as dafBase
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable

from lsst.meas.base.noiseReplacer import NoiseReplacer
from .measureCoadd import MeasureCoaddTask, MeasureCoaddConfig

import lsst.desc.old.bfd as bfd


__all__ = ("MeasurePriorConfig", "MeasurePriorTask")


def matchCatalogs(ra_ref, dec_ref, ra_p, dec_p, dmin):
    ra_ref = np.deg2rad(ra_ref)
    dec_ref = np.deg2rad(dec_ref)
    ra_p = np.deg2rad(ra_p)
    dec_p = np.deg2rad(dec_p)
    posRef = np.dstack([np.sin(dec_ref)*np.cos(ra_ref),
                        np.sin(dec_ref)*np.sin(ra_ref),
                        np.sin(dec_ref)])[0]
    posP = np.dstack([np.sin(dec_p)*np.cos(ra_p),
                      np.sin(dec_p)*np.sin(ra_p),
                      np.sin(dec_p)])[0]
    mytree = scipy.spatial.cKDTree(posRef)
    dist, index = mytree.query(posP)

    # convert to arsec
    dist *= 3600.*(180/np.pi)
    close = dist < dmin
    closeIndices = index[close]
    return close, closeIndices


def matchTreeCatalogs(tree_ref, ra_p, dec_p, dmin):
    ra_p = np.deg2rad(ra_p)
    dec_p = np.deg2rad(dec_p)
    posP = np.dstack([np.sin(dec_p)*np.cos(ra_p),
                      np.sin(dec_p)*np.sin(ra_p),
                      np.sin(dec_p)])[0]
    dist, index = tree_ref.query(posP)

    # convert to arsec
    dist *= 3600*(180./np.pi)

    close = dist < dmin
    closeIndices = index[close]
    return close, closeIndices


def matchXYTreeCatalogs(tree_ref, x, y, dmin):
    posP = np.dstack([x, y])[0]
    dist, index = tree_ref.query(posP)

    close = dist < dmin
    closeIndices = index[close]
    return close, closeIndices


class MeasurePriorConfig(MeasureCoaddConfig):
    snMin = pexConfig.Field(
        dtype=float,
        default=5,
        optional=True,
        doc="Minimun flux S/N"
    )
    snMax = pexConfig.Field(
        dtype=float,
        default=25,
        optional=True,
        doc="Maximum flux S/N"
    )
    fluxMin = pexConfig.Field(
        dtype=float,
        default=None,
        optional=True,
        doc="Minimun flux"
    )
    fluxMax = pexConfig.Field(
        dtype=float,
        default=None,
        optional=True,
        doc="Maximum flux"
    )
    magMin = pexConfig.Field(
        dtype=float,
        default=None,
        optional=True,
        doc="Minimun mag"
    )
    magMax = pexConfig.Field(
        dtype=float,
        default=None,
        optional=True,
        doc="Maximum mag"
    )
    noiseFactor = pexConfig.Field(
        dtype=float,
        default=1,
        optional=True,
        doc="Noise boost factor for kernel smoothing"
    )
    priorSigmaCutoff = pexConfig.Field(
        dtype=float,
        default=5.5,
        optional=True,
        doc="Maximum sigma range when sampling for prior"
    )
    priorSigmaStep = pexConfig.Field(
        dtype=float,
        default=1.,
        optional=True,
        doc="Step size when sampling for prior"
    )
    priorSigmaBuffer = pexConfig.Field(
        dtype=float,
        default=1.,
        optional=True,
        doc="Buffer width of KdTreePrior (in sigma)"
    )
    nSample = pexConfig.Field(
        dtype=int,
        default=30000,
        optional=True,
        doc="Number of templates sampled per target"
    )
    maxXY = pexConfig.Field(
        dtype=float,
        default=4.,
        optional=True,
        doc="Maximum translational displacement in sigma of the nominal covariance matrix"
    )
    sigma = pexConfig.Field(
        dtype=float,
        default=0.8,
        optional=True,
        doc='Sigma used in k-space weight function'
    )
    wIndex = pexConfig.Field(
        dtype=int,
        default=3,
        optional=True,
        doc="index used in k-space weight function"
    )
    centroidName = pexConfig.Field(
        dtype=str,
        default='centroid.sdss',
        optional=True,
        doc="name of centroid to use from the catalog"
    )
    selectionOnly = pexConfig.Field(
        dtype=bool,
        default=False,
        optional=True,
        doc="only do selection in the prior"
    )
    covFile = pexConfig.Field(
        dtype=str,
        default='./test.fits',
        optional=True,
        doc="file that contains the covariance matrices"
    )
    maxVar = pexConfig.Field(
        dtype=float,
        default=0.15,
        optional=True,
        doc="Minimum Variance that will be considered"
    )
    minVar = pexConfig.Field(
        dtype=float,
        default=1e-4,
        optional=True,
        doc="Minimum Variance that will be considered"
    )
    sample = pexConfig.Field(
        dtype=float,
        default=0.2,
        optional=True,
        doc="Only use this fraction of the galaxies"
    )
    invariantCovariance = pexConfig.Field(
        dtype=bool,
        default=False,
        optional=True,
        doc="will all the galaxies have the same covariance"
    )
    reCentroid = pexConfig.Field(
        dtype=bool,
        default=False,
        optional=True,
        doc="recentroid galaxis"
    )
    reCentroidPsf = pexConfig.Field(
        dtype=bool,
        default=False,
        optional=True,
        doc="recentroid galaxis"
    )
    label = pexConfig.Field(
        dtype=str,
        optional=False,
        default='',
        doc="additional label to add to data id"
    )
    matchFile = pexConfig.Field(
        dtype=str,
        optional=True,
        default=None,
        doc="file to match objects to"
    )
    zFile = pexConfig.Field(
        dtype=str,
        optional=True,
        default=None,
        doc="file for redshift information"
    )
    readZFile = pexConfig.Field(
        dtype=bool,
        optional=True,
        default=False,
        doc="read truth z file from butler"
    )
    zPdfFile = pexConfig.Field(
        dtype=str,
        optional=True,
        default=None,
        doc="file for redshift pdf information"
    )
    zMaxCut = pexConfig.Field(
        dtype=float,
        optional=True,
        default=None,
        doc=" maximum redshift cut"
    )
    zMinCut = pexConfig.Field(
        dtype=float,
        optional=True,
        default=None,
        doc=" maximum redshift cut"
    )
    zQualityCut = pexConfig.Field(
        dtype=float,
        optional=True,
        default=None,
        doc="redshift quality cut"
    )
    zField = pexConfig.Field(
        dtype=str,
        optional=True,
        default='frankenz_photoz_mode',
        doc="file for redshift information"
    )
    zRa = pexConfig.Field(
        dtype=str,
        optional=True,
        default='ira',
        doc="ra column from redshift file"
    )
    zDec = pexConfig.Field(
        dtype=str,
        optional=True,
        default='idec',
        doc="ra column from redshift file"
    )
    zId = pexConfig.Field(
        dtype=str,
        optional=True,
        default='object_id',
        doc="ra column from redshift file"
    )
    zQualityField = pexConfig.Field(
        dtype=str,
        optional=True,
        default='frankenz_photoz_zrisk_med',
        doc="field for quality cut"
    )
    zPdfField = pexConfig.Field(
        dtype=str,
        optional=True,
        default='pdf',
        doc="file for redshift information"
    )
    zPdfGridField = pexConfig.Field(
        dtype=str,
        optional=True,
        default='zgrid',
        doc="file for redshift information"
    )
    zBin = pexConfig.Field(
        dtype=int,
        optional=True,
        default=None,
        doc="integer label represnting redshift bin"
    )
    noiseBin = pexConfig.Field(
        dtype=int,
        optional=True,
        default=None,
        doc="integer label represnting noise bin"
    )
    colorBin = pexConfig.Field(
        dtype=int,
        optional=True,
        default=None,
        doc="integer label represnting color bin"
    )
    useLabels = pexConfig.ListField(
        dtype=str,
        default=[],
        optional=True,
        doc="List of labels for which to build the prior. "
    )
    ccFile = pexConfig.Field(
        dtype=str,
        optional=True,
        default=None,
        doc="file for color-color cuts"
    )
    noClobber = pexConfig.Field(
        dtype=bool,
        optional=True,
        default=False,
        doc="check if already exists"
    )
    colorBinFile = pexConfig.Field(
        dtype=str,
        optional=True,
        default=None,
        doc="file for binned color-color cuts"
    )
    useXY = pexConfig.Field(
        dtype=bool,
        default=False,
        optional=True,
        doc="us xy position instead of ra/dec"
    )
    zMatchDist = pexConfig.Field(
        dtype=float,
        optional=True,
        default=0.5,
        doc="redshift quality cut"
    )
    mcSamples = pexConfig.Field(
        dtype=int,
        default=500,
        optional=True,
        doc="number of Monte Carlo samples"
    )
    maxRatio = pexConfig.Field(
        dtype=float,
        default=-1.,
        optional=True,
        doc="apply maximum ratio for selection"
    )
    oldHSC = pexConfig.Field(
        doc="allow ability to run on old HSC data",
        dtype=bool,
        default=False,
    )
    randomSeed = pexConfig.Field(
        dtype=int,
        default=11155,
        optional=True,
        doc="apply maximum ratio for selection"
    )
    selectPeak = pexConfig.Field(
        doc="Select only i-band peaks",
        dtype=bool,
        default=False,
    )
    deblend = pexConfig.Field(
        doc="Was the deblender run?",
        dtype=bool,
        default=True,
    )
    priorType = pexConfig.Field(
        doc="Name of prior type in policy",
        dtype=str,
        default='momentPrior',
    )


class MeasurePriorTask(MeasureCoaddTask):
    ConfigClass = MeasurePriorConfig

    def __init__(self, schema=None, **kwargs):
        """
        """
        MeasureCoaddTask.__init__(self, **kwargs)
        random.seed(self.config.randomSeed)

    def getCovariances(self):
        cat = afwTable.BaseCatalog.readFits(self.config.covFile)

        covList = []

        for rec in cat:
            cov = rec.get('isoCov')
            label = rec.get('label')
            minVariance = rec.get('min')
            maxVariance = rec.get('max')
            cov = np.array(cov.reshape(6, 6), dtype=np.float32)

            covList.append((cov, label, minVariance, maxVariance))

        return covList

    def readInputs(self, dataRef):
        """Return a lsst.pipe.base.Struct containing the Exposure to fit, catalog, measurement
        of the pixel variance and covariance of the matrix as determined by the Psf
        """
        name = self.config.coaddName + "Coadd_calexp"
        if self.config.oldHSC:
            name += "_hsc"
        exposure = dataRef.get(name, immediate=True)

        return pipeBase.Struct(
            sources=dataRef.get(self.dataPrefix + "meas", immediate=True),
            moments=dataRef.get(self.dataPrefix + "moment", immediate=True),
            exposure=exposure,
        )

    def run(self, dataRef):
        """Main driver
        """

        inputs = self.readInputs(dataRef)
        self.log.info("Preparing catalog")
        if self.config.maxObjects is not None:
            first = self.config.firstObject
            last = min(self.config.firstObject + self.config.maxObjects, len(inputs.sources))
            inputs.sources = inputs.sources[first:last]
            inputs.moments = inputs.moments[first:last]
        outCat = self.prepCatalog(inputs)

        self.runMeasure(inputs, outCat, dataRef)

        self.log.info("Writing outputs")
        self.writeOutputs(dataRef, outCat)
        return pipeBase.Struct(outCat=outCat, inputs=inputs)

    def runMeasure(self, inputs, outCat, dataRef):
        # This should tell you how to select objects that apply and also
        # what the covariance matrix is
        covList = self.getCovariances()

        # centroidKey = inputs.sources.schema.find(self.config.centroidName).key
        #noiseKey = outCat.schema.find('noise_variance').key
        noise = 1.

        if self.config.zFile or self.config.readZFile:
            if self.config.readZFile is False:
                import pyfits
                zFile = pyfits.open(self.config.zFile)[1].data
            else:
                zFile = dataRef.get('deepCoadd_truth', immediate=True)
            zRedshift = zFile[self.config.zField]
            zId = zFile[self.config.zId]
            if self.config.zQualityCut is not None:
                zQuality = zFile[self.config.zQualityField]

            if self.config.useXY is False:
                zRa = np.deg2rad(zFile[self.config.zRa])
                zDec = np.deg2rad(zFile[self.config.zDec])

                posRef = np.dstack([np.sin(zDec)*np.cos(zRa), np.sin(zDec)*np.sin(zRa),
                                    np.sin(zDec)])[0]
                zTree = scipy.spatial.cKDTree(posRef)
            else:
                zX = zFile[self.config.zRa]
                zY = zFile[self.config.zDec]

                posRef = np.dstack([zX, zY])[0]
                zTree = scipy.spatial.cKDTree(posRef)

            if self.config.zPdfFile:
                pdfId = pyfits.open(self.config.zPdfFile)[1].data['id']
                pdfData = pyfits.open(self.config.zPdfFile)[1].data[self.config.zPdfField]
                zPdfGrid = pyfits.open(self.config.zPdfFile)[2].data[self.config.zPdfGridField]

                zPdfIndex = (zPdfGrid < self.config.zMaxCut) & (zPdfGrid > self.config.zMinCut)
                self.log.debug('pdf sum in redshift range %0.2f' % np.sum(zPdfIndex))

            if self.config.colorBinFile:

                ccFile = pyfits.open(self.config.colorBinFile)[1].data
                colorBins = ccFile['color_bin']
                assert len(colorBins == len(zId))

        if self.config.ccFile:
            ccFile = np.load(self.config.ccFile)
            ccRa = np.deg2rad(ccFile['ra'])
            ccDec = np.deg2rad(ccFile['dec'])

            posRef = np.dstack([np.sin(ccDec)*np.cos(ccRa), np.sin(ccDec)*np.sin(ccRa),
                                np.sin(ccDec)])[0]
            ccTree = scipy.spatial.cKDTree(posRef)

        for cov, label, varMin, varMax in covList:

            if (label not in self.config.useLabels) and len(self.config.useLabels) > 0:
                self.log.info("Label %s not in %s" % (label, self.config.useLabels))
                continue

            # Build full label
            full_label = label + self.config.label
            if self.config.selectionOnly:
                full_label += '_selection'

            if self.config.noiseBin is not None:
                full_label += '_n%d' % self.config.noiseBin

            if self.config.zBin is not None:
                full_label += '_z%d' % (self.config.zBin)

            if self.config.colorBin is not None:
                full_label += '_c%d' % self.config.colorBin

            # if it already exists and we care don't run anything
            if self.config.noClobber:
                dataRef.dataId['label'] = full_label
                if dataRef.datasetExists(self.dataPrefix + self.config.priorType):
                    self.log.info('Prior already exists %s. skipping' % dataRef.dataId)
                    return

            self.log.info('Processing label %s' % label)
            sigmaFlux = np.sqrt(cov[0, 0])

            minFlux = self.config.snMin*sigmaFlux
            if self.config.fluxMin is not None:
                minFlux = self.config.fluxMin
            elif self.config.magMin is not None:
                # Current assumption is that coadd zeropoint is 27, true for HSC
                minFlux = 10**(-0.4*(self.config.magMin-27))

            maxFlux = self.config.snMax*sigmaFlux
            if self.config.fluxMax is not None:
                maxFlux = self.config.fluxMax
            elif self.config.magMax is not None:
                maxFlux = 10**(-0.4*(self.config.magMax-27))

            momentPrior = bfd.MomentPrior(minFlux, maxFlux,
                                          cov, self.config.invariantCovariance,
                                          self.config.noiseFactor,
                                          self.config.priorSigmaCutoff,
                                          self.config.priorSigmaStep,
                                          self.config.nSample,
                                          self.config.priorSigmaBuffer,
                                          self.config.selectionOnly)
            if self.config.maxRatio > 0:
                momentPrior.setMaxRatio(self.config.maxRatio, self.config.mcSamples)

            momentPrior.setVarianceLimits(varMin, varMax)

            footprints = {measRecord.getId(): (measRecord.getParent(), measRecord.getFootprint())
                          for measRecord in inputs.sources}
            noiseImage = None
            exposureId = None
            noiseReplacer = NoiseReplacer(self.config.noiseReplacer, inputs.exposure, footprints,
                                          noiseImage=noiseImage, log=self.log, exposureId=exposureId)

            ngood = 0
            iteration = 0
            for i, (source, ref, moment) in enumerate(zip(outCat, inputs.sources, inputs.moments)):

                if self.config.sample > 0:
                    if random.uniform(0, 1) > self.config.sample:
                        continue

                with self.noiseContext(inputs, noiseReplacer, ref.getId()) as inputs:
                    if not self.selection(source, ref):
                        self.log.debug("Does not pass selection criteria")
                        continue

                    if moment.get('bfd_flags'):
                        continue

                    # This should happen in the prep catalog stake because we don't want to repeat
                    # I can fix this later...
                    if not self.preMeasure(source, ref, inputs.exposure):
                        continue

                    if (
                        moment.get('bfd_momentsCov')[0] > self.config.maxVar or
                        moment.get('bfd_momentsCov')[0] < self.config.minVar
                    ):
                        self.log.debug(
                            "Does not pass variance cuts %f %0.2f/%0.2f" %
                            (moment.get('bfd_momentsCov')[0], self.config.minVar, self.config.maxVar)
                        )
                        continue

                    sn = moment.get('bfd_moments')[0]/np.sqrt(moment.get('bfd_momentsCov')[0])
                    self.log.debug(
                        "Processing galaxy %d with S/N=%0.2f, flux=%0.2f centroid %0.2f,%0.2f refid %d" %
                        (i, sn, moment.get('bfd_moments')[0], moment.get('bfd_center_x'),
                         moment.get('bfd_center_y'), ref.getId())
                    )

                    # Weight is 1, unless using pdf
                    weight = 1./self.config.sample

                    if self.config.zFile or self.config.readZFile:
                        centroid = ref.getCentroid()
                        sky = inputs.exposure.getWcs().pixelToSky(centroid)
                        if self.config.useXY is False:
                            ra = np.array([sky.toIcrs().getRa().asDegrees()])
                            dec = np.array([sky.toIcrs().getDec().asDegrees()])
                            result = matchTreeCatalogs(zTree, ra, dec, self.config.zMatchDist)
                        else:
                            x = np.array([centroid.getX()])
                            y = np.array([centroid.getY()])
                            result = matchXYTreeCatalogs(zTree, x, y, self.config.zMatchDist)

                        if np.sum(result[0]) < 1:
                            self.log.debug("Failed to match redshift")
                            continue

                        if self.config.zQualityCut is not None:
                            if zQuality[result[1][0]] > self.config.zQualityCut:
                                self.log.debug("Failed redshift quality")
                                continue

                        if (self.config.zMaxCut is not None and self.config.zMinCut is not None and
                                self.config.zPdfFile is None):
                            redshift = zRedshift[result[1][0]]
                            self.log.debug('Redshift %f' % redshift)
                            if (redshift > self.config.zMaxCut) or (redshift < self.config.zMinCut):
                                self.log.debug("Redshift not in range")
                                continue

                        if self.config.zPdfFile is not None:
                            idList = np.where(pdfId == zId[result[1][0]])[0]
                            if len(idList) == 0:
                                self.log.debug("Can't match to pdf")
                                continue

                            weight *= float(np.sum(pdfData[idList[0]][zPdfIndex]))
                            self.log.debug("Using pdf weight %f"%weight)

                        if self.config.colorBinFile is not None:
                            colorBin = colorBins[result[1][0]]
                            if colorBin != self.config.colorBin:
                                continue
                            else:
                                self.log.debug("Selected color bin %d"%colorBin)

                    if self.config.ccFile:
                        centroid = ref.getCentroid()
                        sky = inputs.exposure.getWcs().pixelToSky(centroid)
                        result = matchTreeCatalogs(ccTree, np.array([sky.toIcrs().getRa().asDegrees()]),
                                                   np.array([sky.toIcrs().getDec().asDegrees()]),
                                                   self.config.zMatchDist)

                        if np.sum(result[0]) < 1:
                            self.log.debug('Failed to match color-color')
                            continue

                    if weight < 1e-6:
                        self.log.debug('Weight too low')
                        continue

                    bfd_control = bfd.BfdKMomentControl()
                    bfd_control.sigma = self.config.sigma
                    bfd_control.wIndex = self.config.wIndex
                    bfd_control.maxCentroid = self.config.maxXY
                    bfd_control.ignorePsf = False
                    bfd_control.shift = True
                    bfd_control.reCentroid = self.config.reCentroid
                    bfd_control.reCentroidPsf = self.config.reCentroidPsf

                    try:
                        self.log.debug("Passed cuts")
                        pos = afwGeom.Point2D(moment.get('bfd_center_x'), moment.get('bfd_center_y'))
                        priorGalaxy = bfd.PriorGalaxy(bfd_control)
                        passed = priorGalaxy.addImage(ref, inputs.exposure, pos,
                                                      # This needs to get fixed, it is set to nothing currently
                                                      # source.get(noiseKey),
                                                      noise, True, i)
                        if passed is False:
                            self.log.debug('Failed add image')
                            continue

                        flip = True
                        momentPrior.addPriorGalaxy(priorGalaxy, self.config.maxXY, weight,
                                                   flip, ref.getId())
                        ngood += 1
                    except MemoryError:
                        raise
                    except Exception as err:
                        self.log.warn("Error measuring source %s : %s"
                                      % (source.getId(), err))
                    iteration += 1

            noiseReplacer.end()

            selectionPqr = momentPrior.selectionProbability(cov)
            deselect = selectionPqr.copy()

            # compute not selected probabily 1 - selection
            deselect[0] = 1 - selectionPqr[0]
            for i in range(1, 6):
                deselect[i] *= -1.

            self.log.info('Used %d galaxies in prior' % ngood)
            catalog = momentPrior.getCatalog()
            metadata = dafBase.PropertyList()
            metadata.set('cov', np.array(cov.flatten(), dtype=float))
            metadata.set('selectionPqr', selectionPqr.astype(np.float))
            metadata.set('deselectPqr', deselect.astype(np.float))
            metadata.set('fluxMin', minFlux)
            metadata.set('fluxMax', maxFlux)
            metadata.set('varMin', varMin)
            metadata.set('varMax', varMax)
            if self.config.zFile is not None:
                metadata.set('zFile', self.config.zFile)
                metadata.set('zField', self.config.zField)
            if self.config.zMaxCut is not None:
                metadata.set('zMaxCut', self.config.zMaxCut)
            if self.config.zMinCut is not None:
                metadata.set('zMinCut', self.config.zMinCut)

            if self.config.zPdfFile is not None:
                metadata.set('zPdfFile', self.config.zPdfFile)
            if self.config.colorBinFile is not None:
                metadata.set('colorBinFile', self.config.colorBinFile)

            if self.config.zBin is not None:
                metadata.set('zBin', self.config.zBin)
            if self.config.noiseBin is not None:
                metadata.set('noiseBin', self.config.noiseBin)
            if self.config.colorBin is not None:
                metadata.set('colorBin', self.config.colorBin)

            if self.config.maxRatio > 0:
                metadata.set('maxRatio', self.config.maxRatio)
                metadata.set('mcSamples', self.config.mcSamples)

            metadata.set('noiseFactor', self.config.noiseFactor)
            metadata.set('priorSigmaCutoff', self.config.priorSigmaCutoff)
            metadata.set('priorSigmaStep', self.config.priorSigmaStep)
            metadata.set('priorSigmaBuffer', self.config.priorSigmaBuffer)
            metadata.set('nsample', self.config.nSample)
            metadata.set('selectionOnly', self.config.selectionOnly)
            metadata.set('invariantCovariance', self.config.invariantCovariance)
            metadata.set('maxXY', self.config.maxXY)
            metadata.set('sigma', self.config.sigma)
            metadata.set('wIndex', self.config.wIndex)
            metadata.set('centroidName', self.config.centroidName)
            metadata.set('covFile', self.config.covFile)
            metadata.set('totalWeight', momentPrior.getTotalWeight())

            catalog.getTable().setMetadata(metadata)

            self.log.info(
                'Created %d templates, total weight %f' % (len(catalog), momentPrior.getTotalWeight())
            )
            self.log.info('Selection %s' % selectionPqr)
            dataRef.dataId['label'] = full_label
            dataRef.put(catalog, self.dataPrefix + self.config.priorType)

    def prepCatalog(self, inputs):
        """Prepare the prior and return the output catalog
        """
        outCat = afwTable.SourceCatalog(self.schema)
        srcCat = inputs.sources

        for srcRecord in srcCat:
            outRecord = outCat.addNew()
            outRecord.setId(srcRecord.getId())

        if self.config.calculateVariance:
            x0 = inputs.exposure.getXY0().getX()
            y0 = inputs.exposure.getXY0().getY()
            self.xy0 = afwGeom.Extent2I(-x0, -y0)

        return outCat

    def writeOutputs(self, dataRef, outCat):
        """Write task outputs using the butler.
        """
        return

    def selection(self, source, ref):
        childName = 'deblend_nchild'
        if self.config.oldHSC is False:
            childName = 'deblend_nChild'

        # Don't process blended parent objects
        if self.config.deblend:
            if ref.getParent() == 0 and ref.get(childName) > 0:
                return False
        if ref.getFootprint().getArea() > self.config.maxArea:
            return False
        if ref.getFootprint().getArea() == 0:
            return False
        if self.config.selectPeak:
            if ref.get('merge_peak_i') is False:
                return False

        return True

    @contextmanager
    def noiseContext(self, inputs, noiseReplacer, id):
        """Context manager that applies and removes gain
        """

        noiseReplacer.insertSource(id)
        try:
            yield inputs
        finally:
            noiseReplacer.removeSource(id)
