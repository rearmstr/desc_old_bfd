import sys
import math
import collections

import lsst.daf.base as dafBase
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.afw.cameraGeom as afwCg

from lsst.pipe.base import Struct, ArgumentParser
from lsst.pex.config import Config, Field, ConfigurableField, ListField
from lsst.ctrl.pool.pool  import abortOnError, NODE, Pool, Debugger
from lsst.ctrl.pool.parallel  import BatchPoolTask
from lsst.pipe.drivers.utils import ButlerTaskRunner, getDataRef
#from hsc.meas.tansip.solvetansip import SolveTansipTask
#from .measureCcd import MeasureCcdTask
import numpy as np

Debugger().enabled = True

class SummaryStats:

    def __init__(self, schema):
        self.schema = schema
        self.ccdKey = self.schema.addField('ccd',type=float)
        self.visitKey = self.schema.addField('visit',type=float)
        self.nstarKey = self.schema.addField('nstar',type=float)
        self.medE1Key = self.schema.addField('med_e1',type=float)
        self.stdE1Key = self.schema.addField('std_e1',type=float)
        self.meanE1Key = self.schema.addField('mean_e1',type=float)
        self.madE1Key = self.schema.addField('mad_e1',type=float)
        self.slopeE1XKey = self.schema.addField('slope_e1_x',type=float)
        self.slopeE1YKey = self.schema.addField('slope_e1_y',type=float)
        self.medE2Key = self.schema.addField('med_e2',type=float)
        self.stdE2Key = self.schema.addField('std_e2',type=float)
        self.meanE2Key = self.schema.addField('mean_e2',type=float)
        self.madE2Key = self.schema.addField('mad_e2',type=float)
        self.slopeE2XKey = self.schema.addField('slope_e2_x',type=float)
        self.slopeE2YKey = self.schema.addField('slope_e2_y',type=float)
        self.medSizeKey = self.schema.addField('med_size',type=float)
        self.stdSizeKey = self.schema.addField('std_size',type=float)
        self.meanSizeKey = self.schema.addField('mean_size',type=float)
        self.madSizeKey = self.schema.addField('mad_size',type=float)
        self.slopeSizeXKey = self.schema.addField('slope_size_x',type=float)
        self.slopeSizeYKey = self.schema.addField('slope_size_y',type=float)


    def run(self,dataRef, ccd, visit):

        src = dataRef.get('src')
        #print dataRef.schema.getNames()
        #print src.schema
        mask = src.get('base_PixelFlags_flag_interpolated')==False#np.logical_and(src.get('calib.psf.reserved')==1,
                              #src.get('flags.pixel.interpolated.any')==False)
        sIxxKey = src.schema.find('base_SdssShape_xx').key
        sIyyKey = src.schema.find('base_SdssShape_yy').key
        sIxyKey = src.schema.find('base_SdssShape_xy').key
        mIxxKey = src.schema.find('base_SdssShape_psf_xx').key
        mIyyKey = src.schema.find('base_SdssShape_psf_yy').key
        mIxyKey = src.schema.find('base_SdssShape_psf_xy').key
        sShapeKey = src.schema.find('shape.sdss').key
        mShapeKey = src.schema.find('shape.sdss.psf').key

        stars = src[mask]
        starIxx = src.get(sIxxKey)[mask]
        starIxy = src.get(sIxyKey)[mask]
        starIyy = src.get(sIyyKey)[mask]
        modelIxx = src.get(mIxxKey)[mask]
        modelIxy = src.get(mIxyKey)[mask]
        modelIyy = src.get(mIyyKey)[mask]
        xPos = src.getX()[mask]
        yPos = src.getY()[mask]

        starE1 = (starIxx-starIyy)/(starIxx+starIyy)
        starE2 = (2*starIxy)/(starIxx+starIyy)
        starSize = np.sqrt( 0.5*(starIxx + starIyy))*2.35*0.168 #np.array([a.get(sShapeKey).getDeterminantRadius()*2.35*0.17 for a in stars])
        modelE1 = (modelIxx-modelIyy)/(modelIxx+modelIyy)
        modelE2 = (2*modelIxy)/(modelIxx+modelIyy)
        modelSize = np.sqrt( 0.5*(modelIxx + modelIyy))*2.35*0.168#np.array([a.get(mShapeKey).getDeterminantRadius()*2.35*0.17 for a in stars])

        resE1 = starE1 - modelE1
        resSize = starSize - modelSize
        resE2 = starE2 - modelE2
        medE1 = np.median(resE1)
        medSize = np.median(resSize)
        medE2 = np.median(resE2)
        meanE1 = np.mean(resE1)
        meanSize = np.mean(resSize)
        meanE2 = np.mean(resE2)

        madE1 = sigma_mad(resE1)
        madE2 = sigma_mad(resE2)
        madSize = sigma_mad(resSize)
        stdE1 = np.std(resE1)
        stdE2 = np.std(resE2)
        stdSize = np.std(resSize)
        nStar = len(stars)

        return Struct(nstar=nStar,
                      medE1=medE1,meanE1=meanE1,madE1=madE1,stdE1=stdE1,
                      medE2=medE2,meanE2=meanE2,madE2=madE2,stdE2=stdE2,
                      medSize=medSize,meanSize=meanSize,madSize=madSize,stdSize=stdSize,
                      ccd=ccd,visit=visit)


    def finalize(self, args, visit):
        self.cat=afwTable.BaseCatalog(self.schema)
        for s in args:
            if s is None:
                continue
            rec = self.cat.addNew()
            rec.set(self.nstarKey,s.nstar)
            rec.set(self.visitKey,s.visit)
            rec.set(self.ccdKey,s.ccd)
            rec.set(self.medE1Key,s.medE1)
            rec.set(self.madE1Key,s.madE1)
            rec.set(self.meanE1Key,s.meanE1)
            rec.set(self.stdE1Key,s.stdE1)

            rec.set(self.medE2Key,s.medE2)
            rec.set(self.madE2Key,s.madE2)
            rec.set(self.meanE2Key,s.meanE2)
            rec.set(self.stdE2Key,s.stdE2)

            rec.set(self.medSizeKey,s.medSize)
            rec.set(self.madSizeKey,s.madSize)
            rec.set(self.meanSizeKey,s.meanSize)
            rec.set(self.stdSizeKey,s.stdSize)

        self.cat.writeFits('psf_stat_%d.fits'%visit)



class CreateStarCat:

    def __init__(self, schema):
        self.schema = schema
        self.psfE1Key = self.schema.addField('psf_e1',type=float)
        self.psfE2Key = self.schema.addField('psf_e2',type=float)
        self.psfSizeKey = self.schema.addField('psf_size',type=float)
        self.starE1Key = self.schema.addField('star_e1',type=float)
        self.starE2Key = self.schema.addField('star_e2',type=float)
        self.starSizeKey = self.schema.addField('star_size',type=float)
        self.psfIxxKey = self.schema.addField('psf_ixx',type=float)
        self.psfIyyKey = self.schema.addField('psf_iyy',type=float)
        self.psfIxyKey = self.schema.addField('psf_ixy',type=float)
        self.starIxxKey = self.schema.addField('star_ixx',type=float)
        self.starIyyKey = self.schema.addField('star_iyy',type=float)
        self.starIxyKey = self.schema.addField('star_ixy',type=float)
        self.magKey = self.schema.addField('psf_mag',type=float)
        self.xKey = self.schema.addField('x',type=float)
        self.yKey = self.schema.addField('y',type=float)
        self.fxKey = self.schema.addField('fx',type=float)
        self.fyKey = self.schema.addField('fy',type=float)
        self.visitKey = self.schema.addField('visit',type=np.int32)
        self.ccdKey = self.schema.addField('ccd',type=np.int32)
        self.raKey = self.schema.addField('ra',type=float)
        self.decKey = self.schema.addField('dec',type=float)
        self.resKey = self.schema.addField('reserved',type=np.int32)
        self.isoKey = self.schema.addField('isolated',type=np.int32)


    def run(self,dataRef, ccd, visit, useCalib):

        src = dataRef.get('src')
        # if useCalib:
        #     calexp = dataRef.get('calexp')
        #     calib = calexp.getCalib()
        #     calib.setThrowOnNegativeFlux(0)
        #     ccdTransform = afwCg.cast_Ccd(calexp.getDetector())
        mask = np.logical_and.reduce((
            np.logical_or(src.get('calib_psf_reserved')==1,
                          src.get('calib_psfUsed')==1),
            #src.get('calib.psf.candidate')==1,#),
            src.get('base_PixelFlags_flag_interpolated')==False,
            src.get('base_PixelFlags_flag_saturated')==False
        ))

        sIxxKey = src.schema.find('base_SdssShape_xx').key
        sIyyKey = src.schema.find('base_SdssShape_yy').key
        sIxyKey = src.schema.find('base_SdssShape_xy').key
        mIxxKey = src.schema.find('base_SdssShape_psf_xx').key
        mIyyKey = src.schema.find('base_SdssShape_psf_yy').key
        mIxyKey = src.schema.find('base_SdssShape_psf_xy').key
        magKey = src.schema.find('base_PsfFlux_flux').key
        resKey = src.schema.find('calib_psf_reserved').key

        stars = src[mask]
        starIxx = src.get(sIxxKey)[mask]
        starIxy = src.get(sIxyKey)[mask]
        starIyy = src.get(sIyyKey)[mask]
        modelIxx = src.get(mIxxKey)[mask]
        modelIxy = src.get(mIxyKey)[mask]
        modelIyy = src.get(mIyyKey)[mask]
        resL = src.get(resKey)[mask]
        #resL = [0]*np.sum(mask)
        xPos = src.getX()[mask]
        yPos = src.getY()[mask]
        starE1 = (starIxx-starIyy)/(starIxx+starIyy)
        starE2 = (2*starIxy)/(starIxx+starIyy)
        starSize = np.sqrt( 0.5*(starIxx + starIyy))*2.35*0.168 
        modelE1 = (modelIxx-modelIyy)/(modelIxx+modelIyy)
        modelE2 = (2*modelIxy)/(modelIxx+modelIyy)
        modelSize = np.sqrt( 0.5*(modelIxx + modelIyy))*2.35*0.168 
        raPos = [a.getCoord().getRa().asDegrees() for a in stars]
        decPos = [a.getCoord().getDec().asDegrees() for a in stars]

        #if useCalib:
        #    fp=[ccdTransform.getPositionFromPixel(a.getCentroid()).getMm() for a in stars]
        #    psfMag = calib.getMagnitude(src.get(magKey)[mask])
        #else:
        fp = [a.getCentroid() for a in stars]
        psfMag = -2.5*np.log10(src.get(magKey)[mask])

        xfPos = [a.getX() for a in fp]
        yfPos = [a.getY() for a in fp]
        isos = [(a.getParent()==0 & a.get('deblend_nChild')==0) for a in stars]

        list = []
        for (ra,dec,x,y,fx,fy,me1,me2,msize,mxx,myy,mxy,se1,se2,ssize,sxx,syy,sxy,mag, res,iso) in zip(
            raPos,decPos,xPos,yPos,xfPos,yfPos,modelE1,modelE2,modelSize,modelIxx,modelIyy,modelIxy,
                starE1,starE2,starSize,starIxx,starIyy,starIxy,psfMag, resL,isos):

            list.append(Struct(ra=ra,dec=dec,x=x,y=y,fx=fx,fy=fy,visit=visit,ccd=ccd,
                               me1=me1,me2=me2,msize=msize,mxx=mxx,myy=myy,mxy=mxy,
                               se1=se1,se2=se2,ssize=ssize,sxx=sxx,syy=syy,sxy=sxy,mag=mag,res=res,iso=iso))
        return list

    def finalize(self, args, visit):
        self.cat=afwTable.BaseCatalog(self.schema)
        for cat in args:
            if cat is None: continue
            for s in cat:

                rec = self.cat.addNew()
                rec.set(self.xKey,s.x)
                rec.set(self.yKey,s.y)
                rec.set(self.raKey,s.ra)
                rec.set(self.decKey,s.dec)
                rec.set(self.fxKey,s.fx)
                rec.set(self.fyKey,s.fy)
                rec.set(self.psfE1Key,s.me1)
                rec.set(self.psfE2Key,s.me2)
                rec.set(self.psfSizeKey,s.msize)
                rec.set(self.starE1Key,s.se1)
                rec.set(self.starE2Key,s.se2)
                rec.set(self.starSizeKey,s.ssize)
                rec.set(self.psfIxxKey,s.mxx)
                rec.set(self.psfIyyKey,s.myy)
                rec.set(self.psfIxyKey,s.mxy)
                rec.set(self.starIxxKey,s.sxx)
                rec.set(self.starIyyKey,s.syy)
                rec.set(self.starIxyKey,s.sxy)
                rec.set(self.visitKey,s.visit)
                rec.set(self.ccdKey,s.ccd)
                rec.set(self.magKey,s.mag)
                rec.set(self.resKey,int(s.res))
                rec.set(self.isoKey,int(s.iso))


        self.cat.writeFits('psf_cat_%d.fits'%visit)



def sigma_mad(arr):
    """ Median Absolute Deviation: a "Robust" version of standard deviation.
        Indices variabililty of the sample.
        https://en.wikipedia.org/wiki/Median_absolute_deviation
    """
    arr = np.ma.array(arr).compressed() # should be faster to not use masked arrays.
    med = np.median(arr)
    return 1.4826*np.median(np.abs(arr - med))

class ProcessWlDiagnosticConfig(Config):
    instrument = Field(dtype=str, default="suprimecam", doc="Instrument name, for solvetansip")
    ignoreCcdList = ListField(dtype=int, default=[], doc="List of CCDs to ignore when processing")
    #useCcdList = ListField(dtype=int, default=np.arange(0,104).astype(int), doc="List of CCDs to process")
    useCalib = Field(dtype=bool, default=True, doc="Use calibration zeropoint")

    def setDefaults(self):
        a=1

class ProcessWlDiagnosticTask(BatchPoolTask):
    """Process an entire exposure at once.

    We use MPI to gather the match lists for exposure-wide astrometric and
    photometric solutions.  Note that because of this, different nodes
    see different parts of the code.
    """

    RunnerClass = ButlerTaskRunner
    ConfigClass = ProcessWlDiagnosticConfig
    _DefaultName = "processWlDiagnostic"

    def __init__(self, *args, **kwargs):
        """Constructor.

        All nodes execute this method.
        """
        super(ProcessWlDiagnosticTask, self).__init__(*args, **kwargs)
        self.schema = afwTable.Schema()
        self.runCcd = CreateStarCat(self.schema)

#        self.ccdKey = self.schema.addField('ccd',type=float)
#        self.visitKey = self.schema.addField('visit',type=float)
#        self.nstarKey = self.schema.addField('nstar',type=float)
#        self.medIxxKey = self.schema.addField('med_ixx',type=float)
#        self.stdIxxKey = self.schema.addField('std_ixx',type=float)
#        self.meanIxxKey = self.schema.addField('mean_ixx',type=float)
#        self.madIxxKey = self.schema.addField('mad_ixx',type=float)
#        self.slopeIxxXKey = self.schema.addField('slope_ixx_x',type=float)
#        self.slopeIxxYKey = self.schema.addField('slope_ixx_y',type=float)
#        self.medIyyKey = self.schema.addField('med_iyy',type=float)
#        self.stdIyyKey = self.schema.addField('std_iyy',type=float)
#        self.meanIyyKey = self.schema.addField('mean_iyy',type=float)
#        self.madIyyKey = self.schema.addField('mad_iyy',type=float)
#        self.slopeIyyXKey = self.schema.addField('slope_iyy_x',type=float)
#        self.slopeIyyYKey = self.schema.addField('slope_iyy_y',type=float)
#        self.medIxyKey = self.schema.addField('med_ixy',type=float)
#        self.stdIxyKey = self.schema.addField('std_ixy',type=float)
#        self.meanIxyKey = self.schema.addField('mean_ixy',type=float)
#        self.madIxyKey = self.schema.addField('mad_ixy',type=float)
#        self.slopeIxyXKey = self.schema.addField('slope_ixy_x',type=float)
#        self.slopeIxyYKey = self.schema.addField('slope_ixy_y',type=float)
#        self.cat=afwTable.BaseCatalog(self.schema)
#


    @classmethod
    def batchWallTime(cls, time, parsedCmd, numNodes):
        #numCcds = sum(1 for raft in parsedCmd.butler.get("camera") for ccd in afwCg.cast_Raft(raft))
        #numCycles = int(math.ceil(numCcds/float(numNodes*numProcs)))
        #numExps = len(cls.RunnerClass.getTargetList(parsedCmd))
        return time

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        doBatch = kwargs.pop("doBatch", False)
        parser = ArgumentParser(name="processWlDiagnostic", *args, **kwargs)
        parser.add_id_argument("--id", datasetType="raw", level="visit",
                               help="data ID, e.g. --id visit=12345")
        return parser

    @abortOnError
    def run(self, expRef, butler):
        """Process a single exposure, with scatter-gather-scatter using MPI.

        All nodes execute this method, though the master and slaves have different
        routes through it.  The expRef is only a DummyDataRef on the slaves.
        """
        pool = Pool("processWlDiagnostic")
        pool.cacheClear()
        pool.storeSet(butler=butler)

        dataIdList = dict([(ccdRef.get("ccdExposureId"), ccdRef.dataId)
                           for ccdRef in expRef.subItems("ccd") if ccdRef.datasetExists("raw")])
        dataIdList = collections.OrderedDict(sorted(dataIdList.items()))

        # filter out unwanted CCDs
        # if len(self.config.useCcdList) > 0:
        dataIdListUse = [dataId for dataId in dataIdList.values()
                              if dataId['ccd'] in np.arange(0,104)]
        # else:
        #dataIdListUse = dataIdList.values()
        # Scatter: process CCDs independently
        #self.cat.reserve(len(dataIdListUse))
        for m in dataIdList.values():
            dataId = m
            break
        structList = pool.map(self.process, dataIdListUse)

        numGood = sum(1 for s in structList if s is not None)
        if numGood == 0:
            self.log.warn("All CCDs in exposure failed")
            return
        self.runCcd.finalize(structList, dataId['visit'])
#        self.cat=afwTable.BaseCatalog(self.schema)
#        print 'LEN',len(structList)
#        for s in structList:
#            if s is None:
#                continue
#            rec = self.cat.addNew()
#            rec.set(self.nstarKey,s.nstar)
#            rec.set(self.visitKey,s.visit)
#            rec.set(self.ccdKey,s.ccd)
#            rec.set(self.medIxxKey,s.medIxx)
#            rec.set(self.madIxxKey,s.madIxx)
#            rec.set(self.meanIxxKey,s.meanIxx)
#            rec.set(self.stdIxxKey,s.stdIxx)
#            rec.set(self.slopeIxxXKey,s.slopeIxxX)
#            rec.set(self.slopeIxxYKey,s.slopeIxxY)
#
#            rec.set(self.medIyyKey,s.medIyy)
#            rec.set(self.madIyyKey,s.madIyy)
#            rec.set(self.meanIyyKey,s.meanIyy)
#            rec.set(self.stdIyyKey,s.stdIyy)
#            rec.set(self.slopeIyyXKey,s.slopeIyyX)
#            rec.set(self.slopeIyyYKey,s.slopeIyyY)
#
#            rec.set(self.medIxyKey,s.medIxy)
#            rec.set(self.madIxyKey,s.madIxy)
#            rec.set(self.meanIxyKey,s.meanIxy)
#            rec.set(self.stdIxyKey,s.stdIxy)
#            rec.set(self.slopeIxyXKey,s.slopeIxyX)
#            rec.set(self.slopeIxyYKey,s.slopeIxyY)
#
#        self.cat.writeFits('psf_stat_%d.fits'%dataId['visit'])
#

    def process(self, cache, dataId):
        """Process a single CCD and save the results for a later write.

        Only slaves execute this method.
        """
        cache.result = None
        if dataId["ccd"] in self.config.ignoreCcdList:
            self.log.warn("Ignoring %s: CCD in ignoreCcdList" % (dataId,))
            return None
        dataRef = getDataRef(cache.butler, dataId)
        ccdId = dataRef.get("ccdExposureId")
        with self.logOperation("processing %s (ccdId=%d)" % (dataId, ccdId)):
            try:
                #result = self.runDiagnostic(dataRef,dataId["ccd"],dataId["visit"] )
                result = self.runCcd.run(dataRef,dataId["ccd"],dataId["visit"],self.config.useCalib )
                print('Number of stars',len(result))
            except Exception as e:
                self.log.warn("Failed to process %s: %s\n" % (dataId, e))
                import traceback
                traceback.print_exc()
                return None

            if result is not None:
                cache.result = result
            return result

    def runDiagnostic(self, dataRef,ccd, visit):

        src = dataRef.get('src')
        #mask = np.logical_and(src.get('calib.psf.reserved')==1,
        #                      src.get('flags.pixel.interpolated.any')==False)
        maxk = src.get('flags.pixel.interpolated.any')==False
        sIxxKey = src.schema.find('shape.sdss.xx').key
        sIyyKey = src.schema.find('shape.sdss.yy').key
        sIxyKey = src.schema.find('shape.sdss.xy').key
        mIxxKey = src.schema.find('shape.sdss.psf.xx').key
        mIyyKey = src.schema.find('shape.sdss.psf.yy').key
        mIxyKey = src.schema.find('shape.sdss.psf.xy').key

        stars = src[mask]
        starIxx = src.get(sIxxKey)[mask]
        starIxy = src.get(sIxyKey)[mask]
        starIyy = src.get(sIyyKey)[mask]
        modelIxx = src.get(mIxxKey)[mask]
        modelIxy = src.get(mIxyKey)[mask]
        modelIyy = src.get(mIyyKey)[mask]
        xPos = src.getX()[mask]
        yPos = src.getY()[mask]

        resIxx = starIxx - modelIxx
        resIxy = starIxy - modelIxy
        resIyy = starIyy - modelIyy
        medIxx = np.median(resIxx)
        medIxy = np.median(resIxy)
        medIyy = np.median(resIyy)
        meanIxx = np.mean(resIxx)
        meanIxy = np.mean(resIxy)
        meanIyy = np.mean(resIyy)

        madIxx = sigma_mad(resIxx)
        madIyy = sigma_mad(resIyy)
        madIxy = sigma_mad(resIxy)
        stdIxx = np.std(resIxx)
        stdIyy = np.std(resIyy)
        stdIxy = np.std(resIxy)
        nStar = len(stars)

        slopeIxxX,cIxxX=np.polyfit(xPos,resIxx,1)
        slopeIyyX,cIyyX=np.polyfit(xPos,resIyy,1)
        slopeIxyX,cIxyX=np.polyfit(xPos,resIxy,1)

        slopeIxxY,cIxxY=np.polyfit(yPos,resIxx,1)
        slopeIyyY,cIyyY=np.polyfit(yPos,resIyy,1)
        slopeIxyY,cIxyY=np.polyfit(yPos,resIxy,1)

#        rec.set(self.nstarKey,nStar) 
#        rec.set(self.medIxxKey,medIxx) 
#        rec.set(self.madIxxKey,madIxx) 
#        rec.set(self.meanIxxKey,meanIxx) 
#        rec.set(self.stdIxxKey,stdIxx) 
#        rec.set(self.slopeIxxXKey,slopeIxxX) 
#        rec.set(self.slopeIxxYKey,slopeIxxY) 
#
#        rec.set(self.medIyyKey,medIyy) 
#        rec.set(self.madIyyKey,madIyy) 
#        rec.set(self.meanIyyKey,meanIyy) 
#        rec.set(self.stdIyyKey,stdIyy) 
#        rec.set(self.slopeIyyXKey,slopeIyyX) 
#        rec.set(self.slopeIyyYKey,slopeIyyY) 
#
#        rec.set(self.medIxyKey,medIxy) 
#        rec.set(self.madIxyKey,madIxy) 
#        rec.set(self.meanIxyKey,meanIxy) 
#        rec.set(self.stdIxyKey,stdIxy) 
#        rec.set(self.slopeIxyXKey,slopeIxyX) 
#        rec.set(self.slopeIxyYKey,slopeIxyY) 
#


        return Struct(nstar=nStar,
                      medIxx=medIxx,meanIxx=meanIxx,madIxx=madIxx,stdIxx=stdIxx,
                      slopeIxxX=slopeIxxX,slopeIxxY=slopeIxxY,
                      medIyy=medIyy,meanIyy=meanIyy,madIyy=madIyy,stdIyy=stdIyy,
                      slopeIyyX=slopeIyyX,slopeIyyY=slopeIyyY,
                      medIxy=medIxy,meanIxy=meanIxy,madIxy=madIxy,stdIxy=stdIxy,
                      slopeIxyX=slopeIxyX,slopeIxyY=slopeIxyY,ccd=ccd,visit=visit)
                      

    def _getConfigName(self):
        return None
    def _getEupsVersionsName(self):
        return None
    def _getMetadataName(self):
        return None



