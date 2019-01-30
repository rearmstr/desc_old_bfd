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
import lsst.daf.base as dafBase
import lsst.pex.config
import lsst.pipe.base
import lsst.afw.image
import lsst.afw.geom
import lsst.afw.table
import numpy as np
from .measureCoadd import MeasureCoaddTask, MeasureCoaddConfig
from lsst.daf.persistence import Butler
import lsst.desc.old.bfd as bfd
import astropy.io.fits as pyfits
import scipy.spatial


__all__ = ("MeasurePqrZConfig", "MeasurePqrZTask")


def matchTreeCatalogs(tree_ref, ra_p, dec_p, dmin):
    ra_p = ra_p*np.pi/180
    dec_p = dec_p*np.pi/180
    posP = np.dstack([np.sin(dec_p)*np.cos(ra_p), np.sin(dec_p)*np.sin(ra_p), np.sin(dec_p)])[0]
    dist, index = tree_ref.query(posP)

    # convert to arsec
    dist *= 3600.*(180/np.pi)

    close = dist < dmin
    closeIndices = index[close]
    return close, closeIndices, dist[close]


def matchXYTreeCatalogs(tree_ref, x, y, dmin):
    posP = np.dstack([x, y])[0]
    dist, index = tree_ref.query(posP)

    close = dist < dmin
    closeIndices = index[close]
    return close, closeIndices, dist[close]


class MeasurePqrZConfig(MeasureCoaddConfig):

    invariantCovariance = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        optional=True,
        doc="will all the galaxies have the same covariance"

    )
    numProc = lsst.pex.config.Field(
        dtype=int,
        default=1,
        optional=True,
        doc="number of processors to use"
    )
    chunk = lsst.pex.config.Field(
        dtype=int,
        default=300,
        optional=True,
        doc="number of processors to use"
    )
    sampleSeed = lsst.pex.config.Field(
        dtype=int,
        default=0,
        optional=True,
        doc="number of processors to use"
    )
    sampleFraction = lsst.pex.config.Field(
        dtype=float,
        default=0.25,
        optional=True,
        doc="number of processors to use"
    )
    priorRerun = lsst.pex.config.Field(
        dtype=str,
        optional=False,
        doc="rerun for the prior"
    )
    priorTracts = lsst.pex.config.ListField(
        dtype=int,
        default=[9813],
        optional=False,
        doc="tract for the prior"
    )
    priorLabel = lsst.pex.config.Field(
        dtype=str,
        default='b0',
        optional=False,
        doc="label for the prior"
    )
    priorFilter = lsst.pex.config.Field(
        dtype=str,
        default='HSC-I',
        optional=False,
        doc="filter for the prior"
    )
    maxPriorFiles = lsst.pex.config.Field(
        dtype=int,
        default=-1,
        optional=False,
        doc="filter for the prior"
    )
    priorPatches = lsst.pex.config.ListField(
        dtype=str,
        default=None,
        optional=True,
        doc="Dictionary for the fraction of the galaxies used for eac"
    )
    checkExists = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        optional=True,
        doc="Dictionary for the fraction of the galaxies used for eac"
    )
    noClobber = lsst.pex.config.Field(
        dtype=bool,
        optional=False,
        default=False,
        doc="check if already exists"
    )
    zFile = lsst.pex.config.Field(
        dtype=str,
        optional=True,
        default=None,
        doc="file for redshift information"
    )
    zField = lsst.pex.config.Field(
        dtype=str,
        optional=False,
        default='mode',
        doc="file for redshift information"
    )
    noClobber = lsst.pex.config.Field(
        dtype=bool,
        optional=False,
        default=False,
        doc="check if already exists"
    )
    ignoreZ = lsst.pex.config.Field(
        dtype=bool,
        optional=False,
        default=False,
        doc="ignore redshift"
    )
    zBins = lsst.pex.config.ListField(
        dtype=float,
        default=[0, 0.55, 0.8, 1.1, 1.6, 3],
        optional=True,
        doc="Dictionary for the fraction of the galaxies used for eac"
    )
    useAllZ = lsst.pex.config.Field(
        dtype=bool,
        optional=False,
        default=True,
        doc="attempt to use all available photo-z estimators to select appropriate prior"
    )
    zType = lsst.pex.config.Field(
        dtype=str,
        optional=False,
        default='frankenz',
        doc="which photo-z estimator if not using all"
    )
    randomizeZ = lsst.pex.config.Field(
        dtype=bool,
        optional=False,
        default=False,
        doc="Randomize redshift assignement"
    )
    zRa = lsst.pex.config.Field(
        dtype=str,
        optional=False,
        default='ira',
        doc="ra column from redshift file"
    )
    zDec = lsst.pex.config.Field(
        dtype=str,
        optional=False,
        default='idec',
        doc="ra column from redshift file"
    )
    zId = lsst.pex.config.Field(
        dtype=str,
        optional=False,
        default='object_id',
        doc="ra column from redshift file"
    )
    useXY = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        optional=True,
        doc="us xy position instead of ra/dec"
    )
    zMatchDist = lsst.pex.config.Field(
        dtype=float,
        optional=False,
        default=0.5,
        doc="redshift quality cut"
    )
    writeEmpty = lsst.pex.config.Field(
        dtype=bool,
        default=True,
        optional=True,
        doc="write empty catalogs"
    )


class MeasurePqrZTask(MeasureCoaddTask):
    ConfigClass = MeasurePqrZConfig

    def __init__(self, schema=None, **kwargs):
        """
        """
        MeasureCoaddTask.__init__(self, **kwargs)

        self.schema = lsst.afw.table.SourceTable.makeMinimalSchema()
        # Should move these into C++?
        self.flagMomKey = self.schema.addField("bfd_flags_moment", doc="flag bad input on moments",
                                               type='Flag')
        self.notSelFluxKey = self.schema.addField("bfd_ns_flux", doc="not selected because of flux",
                                                  type='Flag')
        self.notSelVarKey = self.schema.addField("bfd_ns_var", doc="not selected because of variance",
                                                 type='Flag')
        self.pqrKey = bfd.BfdPqrKey.addFields(self.schema, "bfd")
        self.flagKey = self.schema.find('bfd_flags').key
        self.zBinKey = self.schema.addField("bfd_redshift_bin", doc="Redshift", type=np.int32)
        self.disZKey = self.schema.addField("bfd_redshift_distance", doc="distance to redshift", type=float)
        self.zKey = self.schema.addField("bfd_redshift", doc="redshift used in pqr", type=float)

    def readInputs(self, dataRef):
        """Return a lsst.pipe.base.Struct containing the Exposure to fit, catalog, measurement
        of the pixel variance and covariance of the matrix as determined by the Psf
        """

        return lsst.pipe.base.Struct(
            sources=dataRef.get(self.dataPrefix + "moment", immediate=True,
                                flags=lsst.afw.table.SOURCE_IO_NO_FOOTPRINTS),
        )

    def prepCatalog(self, inputs):
        """Prepare the prior and return the output catalog
        """
        outCat = lsst.afw.table.SourceCatalog(self.schema)
        srcCat = inputs.sources

        for srcRecord in srcCat:
            outRecord = outCat.addNew()
            # outRecord.setId(srcCat.get('id'))

        return outCat

    def prep(self):
        self.prior = bfd.MomentPrior()
        priorFiles = []
        priorButler = Butler(self.config.priorRerun)
        prior_skyMap = priorButler.get('deepCoadd_skyMap')

        for tract in self.config.priorTracts:
            for patchInfo in prior_skyMap[tract]:
                patch = '%d,%d'%patchInfo.getIndex()

                if self.config.priorPatches:
                    if patch not in self.config.priorPatches:
                        continue

                if priorButler.datasetExists('deepCoadd_momentPrior', tract=tract, patch=patch,
                                             filter=self.config.priorFilter, label=self.config.priorLabel):
                    priorFiles.append(priorButler.get('deepCoadd_momentPrior_filename',
                                                      tract=tract, patch=patch,
                                                      filter=self.config.priorFilter,
                                                      label=self.config.priorLabel)[0])

        max_file = len(priorFiles)
        if self.config.maxPriorFiles > 0:
            max_file = self.config.maxPriorFiles

        first = True
        self.zBin = None
        for file in priorFiles[:max_file]:
            if file.find('_parent') > 0:
                self.log.info("Skipping %s, from parent" % file)
                continue
            self.log.info("Adding prior %s" % file)
            try:
                cat = lsst.afw.table.BaseCatalog.readFits(file)
                self.prior.addCatalog(cat, self.config.invariantCovariance,
                                      self.config.sampleFraction, self.config.sampleSeed)
                # Should be same for all prior catalogs
                if first:
                    self.cov = np.array(cat.getTable().getMetadata().getArrayDouble('COV')).reshape(6, 6)
                    self.zBin = cat.getTable().getMetadata().getInt('ZBIN')
                    self.fluxBin = cat.getTable().getMetadata().getInt('NOISEBIN')
                    first = False
            except Exception as e:
                print('Failed to read', e)
                continue

        self.prior.prepare()
        self.fluxMin = self.prior.getFluxMin()
        self.fluxMax = self.prior.getFluxMax()
        self.varMin = self.prior.getVarMin()
        self.varMax = self.prior.getVarMax()
        selectionPqr = self.prior.selectionProbability(self.cov.astype(np.float32))
        deselect = selectionPqr.copy()
        deselect[0] = 1 - selectionPqr[0]
        for i in range(1, 6):
            deselect[i] *= -1.
        self.noSelectPqr = deselect

        self.log.info("Reading redshift file")
        zFile = pyfits.open(self.config.zFile)[1].data
        self.zRedshift = np.zeros(len(zFile))
        self.zId = zFile[self.config.zId]
        if self.config.useXY is False:
            zRa = zFile[self.config.zRa]*np.pi/180
            zDec = zFile[self.config.zDec]*np.pi/180
            if self.config.useAllZ:
                mask = zFile['frankenz_photoz_%s_isnull'%self.config.zField] == False
                self.zRedshift[mask] = zFile['frankenz_photoz_%s'%self.config.zField][mask]
                mask = (zFile['mizuki_photoz_%s_isnull'%self.config.zField] == False) & (self.zRedshift == 0.0)
                self.zRedshift[mask] = zFile['mizuki_photoz_%s'%self.config.zField][mask]
                mask = (zFile['nnpz_photoz_%s_isnull'%self.config.zField] == False) & (self.zRedshift == 0.0)
                self.zRedshift[mask] = zFile['nnpz_photoz_%s'%self.config.zField][mask]
                mask = (zFile['mlz_photoz_%s_isnull'%self.config.zField] == False) & (self.zRedshift == 0.0)
                self.zRedshift[mask] = zFile['mlz_photoz_%s'%self.config.zField][mask]
            else:
                mask = zFile['%s_photoz_%s_isnull'%(self.config.zType, self.config.zField)] == False
                self.zRedshift[mask] = zFile['%s_photoz_%s'%(self.config.zType,self.config.zField)][mask]
            posRef = np.dstack([np.sin(zDec)*np.cos(zRa), np.sin(zDec)*np.sin(zRa), np.sin(zDec)])[0]
            self.log.info('Building redshift treee')
            self.zTree = scipy.spatial.cKDTree(posRef)
        else:
            zX = zFile[self.config.zRa]
            zY = zFile[self.config.zDec]
            posRef = np.dstack([zX, zY])[0]
            self.zTree = scipy.spatial.cKDTree(posRef)
            self.zRedshift = zFile[self.config.zField]

    def run(self, dataRef):
        """Main driver
        """
        self.log.info("Processing %s" % str(dataRef.dataId))

        if self.config.checkExists:
            dataRef.dataId['label'] = self.config.priorLabel
            if dataRef.datasetExists(self.dataPrefix+"pqr"):
                if self.config.noClobber:
                    self.log.info('Pqr already exists %s. skipping' % dataRef.dataId)
                    return
                filename = dataRef.get(self.dataPrefix+"pqr_filename")[0]
                if filename.find('_parent') < 0:
                    self.log.info("Skipping %s, file %s exists" % (str(dataRef.dataId), filename))
                    return

        inputs = self.readInputs(dataRef)
        if self.config.maxObjects is not None:
            first = self.config.firstObject
            last = min(self.config.firstObject + self.config.maxObjects, len(inputs.sources))
            inputs.sources = inputs.sources[first:last]

        outCat = self.runMeasure(inputs.sources, dataRef)
        print('Number of sources', len(outCat))
        if len(outCat) == 0:
            self.log.info("No objects processed")
            if self.config.writeEmpty is False:
                return lsst.pipe.base.Struct(outCat=outCat, inputs=inputs)

        self.writeOutputs(dataRef, outCat)

        return lsst.pipe.base.Struct(outCat=outCat, inputs=inputs)

    def runMeasureMulti(self, args):
        self.runMeasure(self, *args)

    def runMeasure(self, sources, dataRef):

        flags = sources.get('bfd_flags')
        flux = sources.get('bfd_moments')[:, 0]
        noise = sources.get('bfd_momentsCov')[:, 0]
        pqrKey = self.schema.find('bfd_pqr').key

        # Preslection cuts
        pre_sel = flags == False

        # redshift selection
        minimumDist = 0.5

        if self.config.useXY is False:
            ra = np.rad2deg(sources['coord_ra'])
            dec = np.rad2deg(sources['coord_dec'])
            result = matchTreeCatalogs(self.zTree, ra, dec, minimumDist)
        else:
            ra = np.array(sources['bfd.center.x'])
            dec = np.array(sources['bfd.center.y'])
            result = matchXYTreeCatalogs(self.zTree, ra, dec, self.config.zMatchDist)

        redshift = np.zeros(len(ra))
        distance = np.zeros(len(ra))

        redshift[result[0] == False] = -1.
        distance[result[0] == False] = -1.
        redshift[result[0]] = self.zRedshift[result[1]]
        distance[result[0]] = result[2]

        self.log.info("distance length %d,%d,%d"%(len(distance), len(redshift), len(sources)))

        self.log.info("Using redshift bins %s:"%self.config.zBins)
        redshiftBin = np.digitize(redshift, self.config.zBins)
        redshiftBin[redshiftBin == len(self.config.zBins)] = len(self.config.zBins)-1

        if self.config.ignoreZ:
            redshiftBin[:] = self.zBin

        if self.config.randomizeZ:
            np.random.shuffle(redshiftBin)

        noiseBins = np.arange(0.05, 1.25, 0.05)
        noiseBin = np.digitize(noise, noiseBins) - 1
        # Flux selection
        sel = np.logical_and.reduce((pre_sel,
                                        noise > self.varMin,
                                        noise < self.varMax,
                                        flux > self.fluxMin,
                                        flux < self.fluxMax,
                                        redshiftBin == self.zBin
        ))
        self.log.info("Surviving cuts:")
        self.log.info("   presel: %d" % np.sum(pre_sel))
        self.log.info("   noise: %d" % np.sum((noise > self.varMin) & (noise < self.varMax)))
        self.log.info("   flux: %d" % np.sum((flux > self.fluxMin) & (flux < self.fluxMax)))
        self.log.info("   redshift: %d" % np.sum(redshiftBin == self.zBin))
        self.log.info("  total:%d" % np.sum(sel))

        outCat = lsst.afw.table.SourceCatalog(self.schema)

        for ii, (srcRecord, dis, zz, noi) in enumerate(zip(sources[sel], distance[sel], redshift[sel], noiseBin[sel])):
            outRecord = outCat.addNew()
            outRecord.setId(srcRecord.getId())
            outRecord.setRa(srcRecord.getRa())
            outRecord.setDec(srcRecord.getDec())

        self.prior.getPqrCat(sources[sel], outCat, self.config.numProc, self.config.chunk)

        return outCat

    def writeOutputs(self, dataRef, outCat):
        """Write task outputs using the butler.
        """
        dataRef.dataId['label'] = self.config.priorLabel
        dataRef.put(outCat, self.dataPrefix+"pqr")
        return

    def selection(self, source, ref):
        # Don't process blended parent objects
        if ref.getParent() == 0 and ref.get('deblend_nChild') > 0:
            return False
        # For now don't process objects in a blend
        # if ref.getParent()!=0:
        #     return False
        if ref.getFootprint().getArea() > self.config.maxArea:
            return False
        if ref.getFootprint().getArea() == 0:
            return False
        # if ref.get('classification.extendedness') == 0:
        #     return False
        return True
