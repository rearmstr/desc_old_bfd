
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
from lsst.pipe.base import Struct, ArgumentParser
from lsst.ctrl.pool.parallel import BatchPoolTask
from lsst.pipe.drivers.coaddDriver import CoaddDriverTaskRunner
from lsst.pipe.tasks.coaddBase import CoaddTaskRunner
from lsst.pipe.drivers.utils import getDataRef, TractDataIdContainer
from lsst.ctrl.pool.pool import Pool, abortOnError, NODE
from lsst.pex.config import Config, Field, ConfigurableField, ListField
from .measureCoadd import MeasureCoaddTask
from .processBfdPatch import PatchRunner

# class ProcessBfdCoaddConfig(Config):
#     processPatch = ConfigurableField(target=ProcessBfdPatchTask, doc="Patch processing task")
#     coaddName = Field(dtype=str, default="deep", doc="Name of coadd")


# def unpickle(factory, args, kwargs):
#     """Unpickle something by calling a factory"""
#     return factory(*args, **kwargs)

# class ProcessBfdCoaddTask(BatchParallelTask):
#     """Process Patches in parallel
#     """
#     ConfigClass = ProcessBfdCoaddConfig
#     _DefaultName = "processBfdCoadd"
#     RunnerClass = BatchTaskRunner

#     def __init__(self, butler=None, *args, **kwargs):
#         """!
#         Constructor
#         """
#         BatchParallelTask.__init__(self, *args, **kwargs)
#         self.makeSubtask("processPatch")

#     @classmethod
#     def _makeArgumentParser(cls, *args, **kwargs):
#         kwargs.pop("doBatch", False)
#         parser = ArgumentParser(name="processBfdCoadd", *args, **kwargs)
#         parser.add_id_argument("--id", "deepCoadd", help="data ID, e.g. --id tract=12345 patch=1,2",
#                                 ContainerClass=CoaddDataIdContainer)
#         return parser

#     def run(self, patchRef):
#         """Process a patch, with scatter-gather-scatter using MPI.
#         """
#         with self.logOperation("processing %s" % (patchRef.dataId,)):
#             self.processPatch.run(patchRef)

#     def _getConfigName(self):
#         return None
#     def _getEupsVersionsName(self):
#         return None
#     def _getMetadataName(self):
#         return None

from lsst.pipe.base import ArgumentParser, ButlerInitializedTaskRunner, ConfigDatasetType
from lsst.pipe.tasks.processCcd import ProcessCcdTask
from lsst.pex.config import Config, Field, ConfigurableField, ListField
from lsst.ctrl.pool.parallel import BatchParallelTask, BatchTaskRunner


class ProcessBfdCoaddConfig(Config):
    processCcd = ConfigurableField(
        target=ProcessCcdTask, doc="CCD processing task")
    ignoreCcdList = ListField(dtype=int, default=[],
                              doc="List of CCDs to ignore when processing")
    ccdKey = Field(dtype=str, default="ccd",
                   doc="DataId key corresponding to a single sensor")


class SingleFrameTaskRunner(BatchTaskRunner, ButlerInitializedTaskRunner):
    """Run batches, and initialize Task using a butler"""
    pass


class ProcessBfdCoaddTask(BatchParallelTask):
    """Process CCDs in parallel
    """
    ConfigClass = ProcessBfdCoaddConfig
    _DefaultName = "processBfdCoadd"
    RunnerClass = SingleFrameTaskRunner

    def __init__(self, butler=None, psfRefObjLoader=None, astromRefObjLoader=None, photoRefObjLoader=None,
                 *args, **kwargs):
        """!
        Constructor
        The psfRefObjLoader, astromRefObjLoader, photoRefObjLoader should
        be an instance of LoadReferenceObjectsTasks that supplies an external
        reference catalog. They may be None if the butler argument is
        provided or the particular reference catalog is not required.
        @param[in] butler  The butler is passed to the refObjLoader constructor in case it is
            needed.  Ignored if the refObjLoader argument provides a loader directly.
        @param[in] psfRefObjLoader  Reference catalog loader for PSF determination.
        @param[in] astromRefObjLoader  Reference catalog loader for astrometric calibration.
        @param[in] photoRefObjLoader Reference catalog loader for photometric calibration.
        @param[in,out] kwargs  other keyword arguments for lsst.ctrl.pool.BatchParallelTask
        """
        BatchParallelTask.__init__(self, *args, **kwargs)
        self.ignoreCcds = set(self.config.ignoreCcdList)
        self.makeSubtask("processCcd", butler=butler, psfRefObjLoader=psfRefObjLoader,
                         astromRefObjLoader=astromRefObjLoader, photoRefObjLoader=photoRefObjLoader)

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        kwargs.pop("doBatch", False)
        parser = ArgumentParser(name="processBfdCoadd", *args, **kwargs)
        parser.add_id_argument("--id",
                               datasetType=ConfigDatasetType(
                                   name="processCcd.isr.datasetType"),
                               level="sensor",
                               help="data ID, e.g. --id visit=12345 ccd=67")
        return parser

    def runDataRef(self, sensorRef):
        """Process a single CCD, with scatter-gather-scatter using MPI.
        """
        if sensorRef.dataId[self.config.ccdKey] in self.ignoreCcds:
            self.log.warn("Ignoring %s: CCD in ignoreCcdList" %
                          (sensorRef.dataId))
            return None

        with self.logOperation("processing %s" % (sensorRef.dataId,)):
            return self.processCcd.runDataRef(sensorRef)
>>>>>>> Stashed changes

    def _getConfigName(self):
        return None
    def _getEupsVersionsName(self):
        return None
    def _getMetadataName(self):
        return None
