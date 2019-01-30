
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
from lsst.ctrl.pool.parallel import BatchParallelTask, BatchTaskRunner
from lsst.coadd.utils.coaddDataIdContainer import CoaddDataIdContainer

from .momentSummary import MomentSummaryTask
from astropy.table import Table

class ProcessBfdMomentPriorConfig(Config):
    momentSummary = ConfigurableField(target=MomentSummaryTask, doc="Patch processing task")
    coaddName = Field(dtype=str, default="deep", doc="Name of coadd")


def unpickle(factory, args, kwargs):
    """Unpickle something by calling a factory"""
    return factory(*args, **kwargs)

class ProcessBfdMomentPriorTask(BatchParallelTask):
    """Process Patches in parallel
    """
    ConfigClass = ProcessBfdMomentPriorConfig
    _DefaultName = "processBfdMomentPrior"
    RunnerClass = BatchTaskRunner

    def __init__(self, butler=None, *args, **kwargs):
        """!
        Constructor
        """
        BatchParallelTask.__init__(self, *args, **kwargs)
        self.makeSubtask("momentSummary")

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        kwargs.pop("doBatch", False)
        parser = ArgumentParser(name="processBfdCoadd", *args, **kwargs)
        parser.add_id_argument("--id", "deepCoadd", help="data ID, e.g. --id tract=12345 patch=1,2",
                               ContainerClass=CoaddDataIdContainer)
        return parser

    def run(self, patchRef):
        """Process a patch, with scatter-gather-scatter using MPI.
        """
        priorRefList = []

        for inoise in range(0, 24):
            for iflux in range(1, 3):
                for zbin in range(1, 6):
                    label = 'b%d_n%d_z%d' % (inoise, iflux, zbin)
                    try:
                        uri = patchRef.get('deepCoadd_momentPrior_filename', label=label)[0]
                    except:
                        self.log.info('No data for %s %s', patchRef.dataId, label)
                        continue
                    priorRefList.append(uri)

        for inoise in range(0, 24):
                for zbin in range(1, 6):
                    label = 'b%d_selection_z%d' % (inoise, zbin)
                    try:
                        uri = patchRef.get('deepCoadd_momentPrior_filename', label=label)[0]
                    except:
                        self.log.info('No data for %s %s', patchRef.dataId, label)
                        continue
                    priorRefList.append(uri)

        with self.logOperation("processing %s" % (patchRef.dataId,)):
            self.momentSummary.run(priorRefList)

    def _getConfigName(self):
        return None
    def _getEupsVersionsName(self):
        return None
    def _getMetadataName(self):
        return None
