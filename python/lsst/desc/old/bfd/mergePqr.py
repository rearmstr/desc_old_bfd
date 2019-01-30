import numpy

from lsst.pipe.base import CmdLineTask, Struct, TaskRunner, ArgumentParser
from lsst.pex.config import Config, Field, ListField
from lsst.pipe.tasks.coaddBase import CoaddDataIdContainer
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.afw.detection as afwDetect
from lsst.daf.base import PropertyList
import re

class MergePqrRunner(TaskRunner):
    def makeTask(self, parsedCmd=None, args=None):
        """Provide a butler to the Task constructor"""
        if parsedCmd is not None:
            butler = parsedCmd.butler
        elif args is not None:
            dataRefList, kwargs = args
            butler = dataRefList[0].getButler()
        else:
            raise RuntimeError("Neither parsedCmd or args specified")
        return self.TaskClass(config=self.config, log=self.log, butler=butler)

    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        """Provide a list of patch references for each patch

        The patch references within the list will have different filters.
        """
        refList = {} # Will index this as refList[tract][patch][label] = ref
        for ref in parsedCmd.id.refList:
            tract = ref.dataId["tract"]
            patch = ref.dataId["patch"]
            filter = ref.dataId["filter"]
            label = ref.dataId["label"]
            if not tract in refList:
                refList[tract] = {}
            if not patch in refList[tract]:
                refList[tract][patch] = {}
            if label in refList[tract][patch]:
                raise RuntimeError("Multiple versions of %s" % (ref.dataId,))
                #refList[tract][patch][filter] = {}
            #if filter in refList[tract][patch][filter]:
            #    raise RuntimeError("Multiple versions of %s" % (ref.dataId,))
            refList[tract][patch][label] = ref
        #import pdb;pdb.set_trace()
        return [(list(ref.values()), kwargs) for tract,patch in refList.items() for label,ref in patch.items()]


class MergePqrConfig(Config):
    labelList = ListField(dtype=str, default=[],
                          doc="List of labels to merge into one pqr file")
    coaddName = Field(dtype=str, default="deep", doc="Name of coadd")
    #selectFile = Field(dtype=str, default="/tigress/rea3/bfd_work/priorSelection.fits", doc="Name of coadd")
    #selectFileBin1 = Field(dtype=str, default=None, doc="Name of coadd")
    #selectFileBin2 = Field(dtype=str, default=None, doc="Name of coadd")


class MergePqrTask(CmdLineTask):
    """A base class for merging source catalogs

    Merging detections (MergeDetectionsTask) and merging measurements
    (MergeMeasurementsTask) are currently so similar that it makes
    sense to re-use the code, in the form of this abstract base class.

    Sub-classes should set the following class variables:
    * _DefaultName: name of Task
    * inputDataset: name of dataset to read
    * outputDataset: name of dataset to write
    * getSchemaCatalogs to the output of _makeGetSchemaCatalogs(outputDataset)

    In addition, sub-classes must implement the mergeCatalogs method.
    """
    _DefaultName = 'MergePqr'
    ConfigClass = MergePqrConfig
    RunnerClass = MergePqrRunner

    @classmethod
    def _makeArgumentParser(cls):
        """Create a suitable ArgumentParser

        We will use the ArgumentParser to get a provide a list of data
        references for patches; the RunnerClass will sort them into lists
        of data references for the same patch
        """
        parser = ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "deepCoadd_pqr",
                               ContainerClass=CoaddDataIdContainer,
                               help="data ID, e.g. --id tract=12345 patch=1,2 filter=g label=l1^l2^l3")
        return parser

    def __init__(self, butler=None, schema=None, **kwargs):
        """Initialize the task.

        Keyword arguments (in addition to those forwarded to CmdLineTask.__init__):
         - schema: the schema of the detection catalogs used as input to this one
         - butler: a butler used to read the input schema from disk, if schema is None

        Derived classes should use the getInputSchema() method to handle the additional
        arguments and retreive the actual input schema.
        """
        CmdLineTask.__init__(self, **kwargs)

        # self.select = {}
        # self.deselect = {}
        # selectFile = afwTable.BaseCatalog.readFits(self.config.selectFile)
        # for rec in selectFile:
        #     self.select[rec.get('label')] = rec.get('selection')
        #     self.deselect[rec.get('label')] = rec.get('deselection')

        # if self.config.selectFileBin1:
        #     selectFile1 = afwTable.BaseCatalog.readFits(self.config.selectFileBin1)
        #     for rec in selectFile1:
        #         self.select[rec.get('label')+'_0'] = rec.get('selection')
        #         self.deselect[rec.get('label')+'_0'] = rec.get('deselection')

        # if self.config.selectFileBin2:
        #     selectFile2 = afwTable.BaseCatalog.readFits(self.config.selectFileBin2)
        #     for rec in selectFile2:
        #         self.select[rec.get('label')+'_1'] = rec.get('selection')
        #         self.deselect[rec.get('label')+'_1'] = rec.get('deselection')

        # How many bins of S/N are there?
        self.snBins = 2

    def run(self, patchRefList):
        """Merge coadd sources from multiple bands

        patchRefList: list of patch data reference for each filter
        """
        catalogs = dict(self.readCatalog(patchRef) for patchRef in patchRefList)
        #import pdb;pdb.set_trace()
        mergedCatalog = self.mergeCatalogs(catalogs, patchRefList[0])
        self.write(patchRefList[0], mergedCatalog)

    def readCatalog(self, patchRef):
        """Read input catalog

        We read the input dataset provided by the 'inputDataset'
        class variable.
        """
        labelName = patchRef.dataId["label"]
        catalog = patchRef.get(self.config.coaddName + "Coadd_pqr", immediate=True)
        self.log.info("Read %d pqr sources for label %s: %s" % (len(catalog), labelName, patchRef.dataId))
        return labelName, catalog

    def mergeCatalogs(self, catalogs, patchRef):
        """Merge multiple catalogs

        catalogs: dict mapping filter name to source catalog

        Returns: merged catalog
        """
        cat = catalogs[list(catalogs.keys())[0]]
        schema = afwTable.Schema()

        mapper = afwTable.SchemaMapper(cat.schema)
        mapper.addMinimalSchema(afwTable.SourceTable.makeMinimalSchema(),True)
        for item in cat.schema:
            mapper.addMapping(item.key, item.field.getName())

        outSchema = mapper.getOutputSchema()
        labelKey = outSchema.addField('prior_label',type=str,doc='label of prior',size=15)
        nsKey = {}
        nsPqrKey = {}
        selPqrKey = {}
        if self.snBins > 1:
            for bin in range(self.snBins):
                nsPqrKey[bin] = outSchema.addField('bfd.ns.pqr.%d' % bin,type='ArrayF',
                                                   doc='non selection term for S/N bin %d'%bin,size=6)

                selPqrKey[bin] = outSchema.addField('bfd.sel.pqr.%d'%bin, type='ArrayF',
                                                    doc='selection term for S/N bin %d'%bin,size=6)
                nsKey[bin] = outSchema.addField('bfd.ns.flux.%d'%bin,type='Flag',doc='non selection S/N flag bin %d'%bin)
        selKey  = outSchema.addField('bfd.sel', type='ArrayF',
                                     doc='selection term for this object',size=6)

        catalog = afwTable.SourceCatalog(outSchema)

        flagKey = cat.schema.find('bfd.flags').key
        fluxKey = cat.schema.find('bfd.ns.flux').key
        pqrKey = cat.schema.find('bfd.pqr').key
        labelNames = catalogs.keys()
        self.log.info('combining %s'% str(labelNames))

        for i,records in enumerate(zip(*catalogs.values())):
            new_record = catalog.addNew()
            new_record.set(flagKey,True)
            if self.snBins > 1:
                for bin in range(self.snBins):
                    new_record.set(nsKey[bin],False)
            goodMeasure = False
            flags = numpy.array([((record.get(flagKey)==False)&(record.get(fluxKey)==False)) for record in records])
            if numpy.sum(flags)>1:
                self.log.warn("Found more than one measurement for this record, this shouldn't happen")
                continue

            index = numpy.where(flags==True)[0]
            if len(index) == 1:
                new_record.set(labelKey,labelNames[index[0]])
                new_record.assign(records[index[0]],mapper)
                name = labelNames[index[0]]
                sname = (re.search('(b\d*)',name)).groups()[0]
                new_record.set(flagKey, False)
                #new_record.set(selKey, self.select[sname])

            # Identify objects that are in one S/N bin but not others
            ns_flags = numpy.array([ ((record.get(flagKey)==False)&(record.get(fluxKey)==True)) for record in records])
            index = numpy.where(ns_flags==True)[0]

            if self.snBins > 1:
                if numpy.sum(ns_flags) == 1:
                    name = labelNames[index[0]]
                    sname = (re.search('(b\d*)',name)).groups()[0]
                    if len(name)<4:
                        snBin=0
                        selBin=1
                    else:
                        snBin=1
                        selBin=0
                    new_record.set(nsPqrKey[snBin], records[index[0]].get(pqrKey))
                    #new_record.set(selPqrKey[selBin], self.select[sname+'_'+str(selBin)])
                    new_record.set(nsKey[snBin], True)
            # Identify objects that are not in the flux selection region of all S/N bins
            if numpy.sum(ns_flags) > 1:
                names = numpy.array(labelNames)[ns_flags]
                lnames = [len(a) for a in names]
                ns_label = numpy.array(names)[lnames == numpy.min(lnames)][0]
                #new_record.set(pqrKey, self.deselect[ns_label])
                #new_record.set(selKey, self.select[ns_label])
                new_record.set(fluxKey, True)
                new_record.set(flagKey, False)
                new_record.set(labelKey, ns_label)


        return catalog

    def write(self, patchRef, catalog):
        """Write the output

        We write as the dataset provided by the 'outputDataset'
        class variable.
        """
        print(patchRef.dataId)
        patchRef.put(catalog, "deepCoadd_mergedPqr")
        #catalog.writeFits('test.fits')
        # since the filter isn't actually part of the data ID for the dataset we're saving,
        # it's confusing to see it in the log message, even if the butler simply ignores it.
        #mergeDataId = patchRef.dataId.copy()
        #del mergeDataId["filter"]
        #self.log.info("Wrote merged catalog: %s" % (mergeDataId,))

    def writeMetadata(self, dataRefList):
        """No metadata to write, and not sure how to write it for a list of dataRefs"""
        return None

    def _getConfigName(self):
        return None
    def _getEupsVersionsName(self):
        return None
    def _getMetadataName(self):
        return None
