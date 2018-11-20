
# from lsst.pipe.base import ArgumentParser, ButlerInitializedTaskRunner, ConfigDatasetType
# from lsst.pipe.tasks.processCcd import ProcessCcdTask
# from lsst.pex.config import Config, Field, ConfigurableField, ListField
# from lsst.ctrl.pool.parallel import BatchParallelTask, BatchTaskRunner


# class ProcessBfdCoaddConfig(Config):
#     processCcd = ConfigurableField(
#         target=ProcessCcdTask, doc="CCD processing task")
#     ignoreCcdList = ListField(dtype=int, default=[],
#                               doc="List of CCDs to ignore when processing")
#     ccdKey = Field(dtype=str, default="ccd",
#                    doc="DataId key corresponding to a single sensor")


# class ProcessBfdCoaddTaskRunner(BatchTaskRunner, ButlerInitializedTaskRunner):
#     """Run batches, and initialize Task using a butler"""
#     pass


# class ProcessBfdCoaddTask(BatchParallelTask):
#     """Process CCDs in parallel
#     """
#     ConfigClass = ProcessBfdCoaddConfig
#     _DefaultName = "processBfdCoadd"
#     RunnerClass = ProcessBfdCoaddTaskRunner

#     def __init__(self, butler=None, psfRefObjLoader=None, astromRefObjLoader=None, photoRefObjLoader=None,
#                  *args, **kwargs):
#         """!
#         Constructor
#         The psfRefObjLoader, astromRefObjLoader, photoRefObjLoader should
#         be an instance of LoadReferenceObjectsTasks that supplies an external
#         reference catalog. They may be None if the butler argument is
#         provided or the particular reference catalog is not required.
#         @param[in] butler  The butler is passed to the refObjLoader constructor in case it is
#             needed.  Ignored if the refObjLoader argument provides a loader directly.
#         @param[in] psfRefObjLoader  Reference catalog loader for PSF determination.
#         @param[in] astromRefObjLoader  Reference catalog loader for astrometric calibration.
#         @param[in] photoRefObjLoader Reference catalog loader for photometric calibration.
#         @param[in,out] kwargs  other keyword arguments for lsst.ctrl.pool.BatchParallelTask
#         """
#         BatchParallelTask.__init__(self, *args, **kwargs)
#         self.ignoreCcds = set(self.config.ignoreCcdList)
#         self.makeSubtask("processCcd", butler=butler, psfRefObjLoader=psfRefObjLoader,
#                          astromRefObjLoader=astromRefObjLoader, photoRefObjLoader=photoRefObjLoader)

#     @classmethod
#     def _makeArgumentParser(cls, *args, **kwargs):
#         kwargs.pop("doBatch", False)
#         parser = ArgumentParser(name="processBfdCoadd", *args, **kwargs)
#         parser.add_id_argument("--id",
#                                datasetType=ConfigDatasetType(
#                                    name="processCcd.isr.datasetType"),
#                                level="sensor",
#                                help="data ID, e.g. --id visit=12345 ccd=67")
#         return parser

#     def run(self, sensorRef):
#         """Process a single CCD, with scatter-gather-scatter using MPI.
#         """
#         if sensorRef.dataId[self.config.ccdKey] in self.ignoreCcds:
#             self.log.warn("Ignoring %s: CCD in ignoreCcdList" %
#                           (sensorRef.dataId))
#             return None

#         with self.logOperation("processing %s" % (sensorRef.dataId,)):
#             return self.processCcd.runDataRef(sensorRef)

import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.pipe.base as pipeBase
from lsst.pipe.base import Struct, ArgumentParser, TaskRunner
from lsst.ctrl.pool.parallel import BatchPoolTask #,BatchTaskRunner, BatchParallelTask
#from lsst.pipe.drivers.coaddDriver import CoaddDriverTaskRunner
#from lsst.pipe.drivers.multiBandDriver import MultiBandDriverTaskRunner
#from lsst.pipe.tasks.coaddBase import CoaddTaskRunner
from lsst.pipe.drivers.utils import TractDataIdContainer
from lsst.ctrl.pool.pool import Pool, abortOnError
from lsst.pex.config import Config, Field, ConfigurableField

from .measureCoadd import MeasureCoaddTask
from .processBfdPatch import ProcessBfdPatchTask



class ProcessBfdCoaddConfig(Config):
    processPatch = ConfigurableField(target=ProcessBfdPatchTask, doc="Patch processing task")
    coaddName = Field(dtype=str, default="deep", doc="Name of coadd")


def unpickle(factory, args, kwargs):
    """Unpickle something by calling a factory"""
    return factory(*args, **kwargs)


class BfdCoaddDriverTaskRunner(TaskRunner):
    """TaskRunner for running MultiBandTask
    This is similar to the lsst.pipe.base.ButlerInitializedTaskRunner,
    except that we have a list of data references instead of a single
    data reference being passed to the Task.run, and we pass the results
    of the '--reuse-outputs-from' command option to the Task constructor.
    """

    def __init__(self, TaskClass, parsedCmd, doReturnResults=False):
        TaskRunner.__init__(self, TaskClass, parsedCmd, doReturnResults)

    def makeTask(self, parsedCmd=None, args=None):
        """A variant of the base version that passes a butler argument to the task's constructor
        parsedCmd or args must be specified.
        """
        if parsedCmd is not None:
            butler = parsedCmd.butler
        elif args is not None:
            dataRefList, kwargs = args
            butler = dataRefList[0].butlerSubset.butler
        else:
            raise RuntimeError("parsedCmd or args must be specified")
        return self.TaskClass(config=self.config, log=self.log, butler=butler)



class ProcessBfdCoaddTask(BatchPoolTask):
    ConfigClass = ProcessBfdCoaddConfig
    _DefaultName = "processBfdCoadd"
    RunnerClass = BfdCoaddDriverTaskRunner

    def __init__(self, butler=None, *args, **kwargs):
        BatchPoolTask.__init__(self, **kwargs)

        self.makeSubtask('processPatch')
        self.bulter = butler

    def __reduce__(self):
        """Pickler"""
        return unpickle, (self.__class__, [], dict(config=self.config, name=self._name,
                                                   parentTask=self._parentTask, log=self.log,
                                                   ))
    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        kwargs.pop("doBatch", False)
        parser = ArgumentParser(name=cls._DefaultName, *args, **kwargs)
        parser.add_id_argument("--id", "deepCoadd", help="data ID, e.g. --id tract=12345 patch=1,2",
                               ContainerClass=TractDataIdContainer)
        return parser

    @classmethod
    def batchWallTime(cls, time, parsedCmd, numCores):
        return time

    @abortOnError
    def run(self, patchRefList):
        """!Determine which tracts are non-empty before processing
        @param tractPatchRefList: List of tracts and patches to include in the coaddition
        @param butler: butler reference object
        @param selectIdList: List of data Ids (i.e. visit, ccd) to consider when making the coadd
        @return list of references to sel.runTract function evaluation for each tractPatchRefList member
        """
        pool = Pool("all")

        patches = {}
        tract = None
        for patchRef in patchRefList:
            dataId = patchRef.dataId
            if tract is None:
                tract = dataId["tract"]
            else:
                assert tract == dataId["tract"]

            patch = dataId["patch"]
            if patch not in patches:
                patches[patch] = []
            patches[patch].append(patchRef)
        #import pdb;pdb.set_trace()
        
        pool.map(self.runPatch, patches.values())#[patch for patch in patches.values()])
        #with self.logOperation("processing %s" % (sensorRef.dataId,)):
        #    return [self.map(self.measure, patchRefList) for patchRefList in tractPatchRefList]

    def runPatch(self, cache, dataRefList):
        
        for dataRef in dataRefList:
            with self.logOperation("processing %s" % (dataRef.dataId,)):
                self.processPatch.run(dataRef)


    def _getConfigName(self):
        return None
    def _getEupsVersionsName(self):
        return None
    def _getMetadataName(self):
        return None

# class DriverTaskRunner(pipeBase.TaskRunner):

#     # def __init__(self, TaskClass, parsedCmd, doReturnResults=False):
#     #     CoaddTaskRunner.__init__(self, TaskClass, parsedCmd, doReturnResults)

#     # def makeTask(self, parsedCmd=None, args=None):
#     #     return self.TaskClass(config=self.config, log=self.log)

#     @staticmethod
#     def getTargetList(parsedCmd, **kwargs):
#         """!Get bare butler into Task
#         @param parsedCmd results of parsing command input
#         """
#         kwargs["butler"] = parsedCmd.butler
#         return [(parsedCmd.id.refList, kwargs), ]



# class ProcessBfdCoaddTask(BatchPoolTask):
#     ConfigClass = ProcessBfdCoaddConfig
#     _DefaultName = "processBfdCoadd"
#     RunnerClass = DriverTaskRunner

#     def __init__(self, **kwargs):
#         BatchPoolTask.__init__(self, **kwargs)

#         self.makeSubtask('measure')

#     def __reduce__(self):
#         """Pickler"""
#         return unpickle, (self.__class__, [], dict(config=self.config, name=self._name,
#                                                    parentTask=self._parentTask, log=self.log,
#                                                    ))
#     @classmethod
#     def _makeArgumentParser(cls, *args, **kwargs):
#         kwargs.pop("doBatch", False)
#         parser = ArgumentParser(name=cls._DefaultName, *args, **kwargs)
#         parser.add_id_argument("--id", "deepCoadd", help="data ID, e.g. --id tract=12345 patch=1,2",
#                                ContainerClass=TractDataIdContainer)
#         return parser

#     @classmethod
#     def batchWallTime(cls, time, parsedCmd, numCores):
#         return time

# #    @abortOnError
#     def run(self, tractPatchRefList, butler):
#         """!Determine which tracts are non-empty before processing
#         @param tractPatchRefList: List of tracts and patches to include in the coaddition
#         @param butler: butler reference object
#         @param selectIdList: List of data Ids (i.e. visit, ccd) to consider when making the coadd
#         @return list of references to sel.runTract function evaluation for each tractPatchRefList member
#         """
#         pool = Pool("tracts")
#         pool.storeSet(butler=butler, skymap=butler.get(self.config.coaddName + "Coadd_skyMap"))
#         print(tractPatchRefList)
        
#         for patchRefList in tractPatchRefList:
#              pool.map(self.runBfd, patchRefList)

#         #return [self.runTract(patchRefList, butler) for patchRefList in tractPatchRefList]

#     def runTract(self, patchRefList, butler):

#         pool = Pool("patches")
#         pool.map(self.runBfd, patchRefList)

#     def runBfd(self, cache, patchRef):

#         with self.logOperation("processing %s " % (patchRef.dataId)):
#             try:
#                 result = self.measure.run(patchRef)
#             except Exception as e:
#                 self.log.warn("Failed to process %s: %s\n" % (patchRef.dataId, e))
#                 import traceback
#                 traceback.print_exc()
#                 return None

#     def _getConfigName(self):
#         return None
#     def _getEupsVersionsName(self):
#         return None
#     def _getMetadataName(self):
#         return None


