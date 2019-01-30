import lsst.pipe.base as pipeBase
from lsst.pipe.base import ArgumentParser
from lsst.pex.config import Config, Field, ConfigurableField
from lsst.pipe.tasks.coaddBase import CoaddDataIdContainer, CoaddTaskRunner

from .measureCoadd import MeasureCoaddTask


class ProcessBfdPatchSingleConfig(Config):
    measure = ConfigurableField(target=MeasureCoaddTask, doc="Patch Processing")
    coaddName = Field(dtype=str, default="deep", doc="Name of coadd")


class PatchRunner(CoaddTaskRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        """Get bare butler into Task"""
        kwargs["butler"] = parsedCmd.butler
        return [(parsedCmd.id.refList, kwargs), ]


class ProcessBfdPatchSingleTask(pipeBase.CmdLineTask):
    """Process an tract at once.
    """

    ConfigClass = ProcessBfdPatchSingleConfig
    _DefaultName = "processBfdPatchSingle"
    RunnerClass = PatchRunner

    def __init__(self, *args, **kwargs):
        """Constructor.

        All nodes execute this method.
        """
        super(ProcessBfdPatchSingleTask, self).__init__(*args, **kwargs)
        self.makeSubtask('measure')

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        parser = ArgumentParser(name="processBfdPatchSingle", *args, **kwargs)
        parser.add_id_argument("--id", "deepCoadd", help="data ID, e.g. --id tract=12345 patch=1,2",
                               ContainerClass=CoaddDataIdContainer)
        return parser

    def run(self, tractPatchRefList, butler):
        """
        """
        try:
            self.measure.prep()
        except Exception as e:
            self.log.warn("Failed to prepare: %s\n" % (e))
            return

        #for tract in tractPatchRefList:
        for patchRef in tractPatchRefList:
            result = self.process(patchRef)

    def process(self, dataRef):
        """Process a single Patch
        """
        try:
            result = self.measure.run(dataRef)
        except Exception as e:
            self.log.warn("Failed to process %s: %s\n" % (dataRef.dataId, e))
            import traceback
            traceback.print_exc()
            return None

        return result

    def write(self, cache, struct, focusMd=None):
        """Write the outputs.
        """

    def _getConfigName(self):
        return None

    def _getEupsVersionsName(self):
        return None

    def _getMetadataName(self):
        return None
