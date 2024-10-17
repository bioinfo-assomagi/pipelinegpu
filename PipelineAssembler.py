from Pipeline import Pipeline

from Pipes.RawReadsPipe import RawReadsPipe
from Pipes.QualityReadsPipe import QualityReadsPipe
from Pipes.InputPipe import Setup

from Pipes.DiagnosysCorePipe import ReadFastQFilesPipe
#from Pipes.DiagnosysCorePipe import ProcessFastQFilesPipe
from Pipes.DiagnosysCorePipe import UnzipFastQFilesPipe
from Pipes.DiagnosysCorePipe import ProcessFastQFilesPipe2
from Pipes.DiagnosysCorePipe import PreAlignmentPipe
from Pipes.DiagnosysCorePipe import SetSamplesPipe
from Pipes.DiagnosysCorePipe import ResyncDBPipe
from Pipes.DiagnosysCorePipe import ProcessPipe
from Pipes.DiagnosysCorePipe import CoverageWrapperPipe
from Pipes.ParallelWrapper import ParallelWrapper

from Pipes.EndPipe import EndPipe
from Pipes.StartPipe import StartPipe
from Pipes.Fq2BamPipe import Fq2BamPipe
from Pipes.InterStage import SampleListFam
from Pipes.CoveragePipe import CoveragePipe2
from Pipes.VariantCallPipe import VariantCallPipe
#from Pipes.AnnotationPipe import AnnotationPipe
from Pipes.VariantFilterWrapperPipe import VariantFilterWrapperPipe


class PipelineAssembler():

    def __init__(self):
        pass

    @staticmethod
    def factory(pipeline_type):
        if pipeline_type == "resource":
            return Pipeline(Setup()).assemblePipe(ResyncDBPipe()).assemblePipe(ReadFastQFilesPipe()).assemblePipe(EndPipe())
        elif pipeline_type == "process":
            return Pipeline(StartPipe()).assemblePipe(UnzipFastQFilesPipe()).assemblePipe(ProcessFastQFilesPipe2()).assemblePipe(SetSamplesPipe()).assemblePipe(PreAlignmentPipe()).assemblePipe(EndPipe())
        elif pipeline_type == "interstage":
            return Pipeline(SampleListFam())
        elif pipeline_type == "fq2bam":
            return Pipeline(SampleListFam()).assemblePipe(Fq2BamPipe()).assemblePipe(EndPipe())
        elif pipeline_type == "coverage":
            return Pipeline(CoveragePipe2()).assemblePipe(EndPipe())
        elif pipeline_type == "variantcall":
            #return Pipeline(VariantCallPipe()).assemblePipe(AnnotationPipe()).assemblePipe(EndPipe())
            pass
        elif pipeline_type == "test":
            return Pipeline(Setup()).assemblePipe(ResyncDBPipe()).assemblePipe(ReadFastQFilesPipe()).assemblePipe(ProcessPipe()).assemblePipe(SampleListFam()).assemblePipe(Fq2BamPipe()).assemblePipe(CoverageWrapperPipe()).assemblePipe(VariantCallPipe()).assemblePipe(VariantFilterWrapperPipe()).assemblePipe(EndPipe())
        elif pipeline_type == "vcf_first":
            return Pipeline(Setup()).assemblePipe(ResyncDBPipe()).assemblePipe(ReadFastQFilesPipe()).assemblePipe(ProcessPipe()).assemblePipe(SampleListFam()).assemblePipe(Fq2BamPipe()).assemblePipe(VariantCallPipe()).assemblePipe(CoverageWrapperPipe()).assemblePipe(VariantFilterWrapperPipe()).assemblePipe(EndPipe())
        elif pipeline_type == "bamstart":
            return Pipeline(Setup()).assemblePipe(ResyncDBPipe()).assemblePipe(CoverageWrapperPipe()).assemblePipe(VariantCallPipe()).assemblePipe(VariantFilterWrapperPipe()).assemblePipe(EndPipe())

    def independent_pipeline(pipeline_type):
        if pipeline_type == "coverage":
            return Pipeline(Setup()).assemblePipe(CoveragePipe2())
        elif pipeline_type == "fq2bam":
            return Pipeline(Setup()).assemblePipe(SampleListFam()).assemblePipe(Fq2BamPipe()).assemblePipe(EndPipe())

    def GetParentPipeline():
        pass
        