

"""
Module Description:
-------------------

This module provides the `PipelineAssembler` class, which is responsible for assembling various types
of pipelines based on the provided pipeline type. It utilizes multiple `Pipe` subclasses to construct
complex processing pipelines for tasks such as data setup, processing, variant calling, and annotation.

Available pipeline types include:

- **resource**: Sets up resources and synchronizes the database.
- **process**: Processes raw data files and prepares samples.
- **interstage**: Handles inter-stage data transfer.
- **fq2bam**: Converts FASTQ files to BAM format.
- **coverage**: Computes coverage statistics.
- **variantcall**: (Not implemented) Intended for variant calling processes.
- **test**: A comprehensive test pipeline that includes multiple processing stages.
- **vcf_first**: Processes VCF files with parallel annotation.
- **bamstart**: Starts processing from BAM files.

The module also provides a factory method to create independent pipelines for specific tasks.
"""

from Pipeline import Pipeline
from Pipes.InputPipe import Setup

from Pipes.DiagnosysCorePipe import ReadFastQFilesPipe

from Pipes.DiagnosysCorePipe import ResyncDBPipe
from Pipes.DiagnosysCorePipe import ProcessPipe
from Pipes.DiagnosysCorePipe import CoverageWrapperPipe
from Pipes.CoverageStatisticsWrapperPipe import CoverageStatisticsWrapperPipe
from Pipes.VariantFilterPipe import VariantFilterPipe
from Pipes.AnnotationPipe import AnnotationPipe


from Pipes.ParallelWrapper import ParallelWrapper

from Pipes.EndPipe import EndPipe
from Pipes.Fq2BamPipe import Fq2BamPipe
from Pipes.InterStage import SampleListFam
from Pipes.VariantCallPipe import VariantCallPipe
from Pipes.VariantFilterWrapperPipe import VariantFilterWrapperPipe
from Pipes.IndelWrapperPipe import IndelWrapperPipe
from Pipes.MultipleDiseaseChildSampleHandlerPipe import MultipleDiseaseSampleHandlerPipe


class PipelineAssembler():

    """
    The `PipelineAssembler` class provides factory methods to assemble and retrieve different types
    of pipelines based on the provided pipeline type. It utilizes various `Pipe` subclasses to build
    complex processing pipelines for tasks such as data processing, variant calling, and annotation.

    **Usage Example:**

    .. code-block:: python

        assembler = PipelineAssembler()
        pipeline = assembler.factory('process')
        result = pipeline.start(data='Initial data')

    **Available Pipeline Types:**

    - **resource**: Sets up resources and synchronizes the database.
    - **process**: Processes raw data files and prepares samples.
    - **interstage**: Handles inter-stage data transfer (TBR).
    - **fq2bam**: Converts FASTQ files to BAM format.
    - **coverage**: Computes coverage statistics.
    - **test**: A comprehensive test pipeline that includes multiple processing stages.
    - **vcf_first**: Processes VCF files with parallel annotation.
    - **bamstart**: Starts processing from BAM files.
    """

    def __init__(self):
        pass

    @staticmethod
    def factory(pipeline_type):
        """
         Factory method to assemble and return a pipeline based on the provided `pipeline_type`.

        :param pipeline_type: The type of pipeline to assemble.
        :type pipeline_type: str
        :return: An assembled `Pipeline` instance corresponding to the specified type.
        :rtype: Pipeline
        :raises ValueError: If an unknown pipeline type is provided.
        """
        
        if pipeline_type == "resource":
            pass
        elif pipeline_type == "test":
            return Pipeline(Setup()).assemblePipe(ResyncDBPipe()).assemblePipe(ReadFastQFilesPipe()).assemblePipe(ProcessPipe()).assemblePipe(SampleListFam()).assemblePipe(Fq2BamPipe()).assemblePipe(CoverageWrapperPipe()).assemblePipe(VariantCallPipe()).assemblePipe(CoverageStatisticsWrapperPipe()).assemblePipe(VariantFilterWrapperPipe()).assemblePipe(EndPipe())
        elif pipeline_type == "vcf_first":
            # NOTE: right now all of the following are here due to testing
            #return Pipeline(Setup()).assemblePipe(ResyncDBPipe()).assemblePipe(ReadFastQFilesPipe()).assemblePipe(ProcessPipe()).assemblePipe(SampleListFam()).assemblePipe(Fq2BamPipe()).assemblePipe(VariantCallPipe()).assemblePipe(CoverageWrapperPipe()).assemblePipe(VariantFilterWrapperPipe()).assemblePipe(CoverageStatisticsWrapperPipe()).assemblePipe(IndelWrapperPipe()).assemblePipe(EndPipe())
            #return Pipeline(Setup()).assemblePipe(VariantFilterWrapperPipe()).assemblePipe(EndPipe())
            #return Pipeline(Setup()).assemblePipe(ParallelWrapper(Pipeline(VariantFilterPipe()).assemblePipe(AnnotationPipe()))).assemblePipe(EndPipe())
            #return Pipeline(Setup()).assemblePipe(CoverageStatisticsWrapperPipe()).assemblePipe(EndPipe())
            #return Pipeline(Setup()).assemblePipe(IndelWrapperPipe()).assemblePipe(EndPipe())
            return Pipeline(Setup()).assemblePipe(ResyncDBPipe()).assemblePipe(ReadFastQFilesPipe()).assemblePipe(ProcessPipe()).assemblePipe(Fq2BamPipe()).assemblePipe(VariantCallPipe()).assemblePipe(MultipleDiseaseSampleHandlerPipe()).assemblePipe(CoverageWrapperPipe()).assemblePipe(CoverageStatisticsWrapperPipe()).assemblePipe(IndelWrapperPipe()).assemblePipe(VariantFilterWrapperPipe()).assemblePipe(EndPipe())
        
        elif pipeline_type == "bamstart":
            return Pipeline(Setup()).assemblePipe(ResyncDBPipe()).assemblePipe(CoverageWrapperPipe()).assemblePipe(VariantCallPipe()).assemblePipe(VariantFilterWrapperPipe()).assemblePipe(EndPipe())
        