from Pipes.Pipe import Pipe
from Pipes.ParallelPipe import ParallelPipe
import os
import csv

import config
import sys
import utils
from os.path import join
import glob
import re
import shutil
import tools
from DBContext import DBContext
import pandas as pd
import numpy as np

import dir_tree
from Entities.Sample import Sample

# TODO: if files not found inside samples list, try reading them from principal_directory
class VariantCallPipe():
    """ Class responsible for running varinat calling (pb_gatk and pb_deepvariant). 
    TODO: create a wrapper for this class. I.e. the wrapper should be responsible for handling the samples (receiving them from the
    directories, and supplying them to the this class). The VariantCallPipe should act on one sample at time, so that the convetion
    of each Pipe being responsible for one sample at a time is satisfied."""

    def __init__(self):
        super().__init__()

    def process(self, **kwargs):
        # self.principal_directory = kwargs.pop("principal_directory", None)
        self.principal_directory = dir_tree.principal_directory.path
        sample_jsons = glob.glob(join(dir_tree.principal_directory.sample_data.path, "*.json"))
        samples = [Sample.fromJSON(json_file) for json_file in sample_jsons]
        self.panel = kwargs.pop("panel", None)

        print("TESTING VARIANT CALL PIPE: {}".format(samples))
        self.HaplotypeCaller(samples)
        self.DeepVariant(samples)

        kwargs.update(
            {"principal_directory": self.principal_directory, "samples": samples, "panel": self.panel})
        return kwargs
    

    def HaplotypeCaller(self, samples):

        docker_input_parabricks = os.path.join(config.DOCKER_WORKDIR, 'bam')

        for sample in samples:
            sample_name = str(sample.name)
            try: 
                bam_filename = sample.bam.split('/')[-1]
            except Exception as e:
                raise("In variant calling pipe, BAM file coulnd not be found for sample {}".format(sample_name))

            command = ' '.join(['docker', 'run', '--gpus', 'all', '--rm',
                                '--volume', "{}/:{}".format(config.REF, config.DOCKER_REFDIR),
                                '--volume', "{}/:{}".format(self.principal_directory, config.DOCKER_WORKDIR),
                                '--volume',
                                "{}/:{}".format(os.path.join(self.principal_directory, "temp"), config.DOCKER_OUTPUTDIR),
                                "{}".format(config.PARABRICKS_VERSION),
                                'pbrun', 'haplotypecaller',
                                '--ref', "{}/{}".format(config.DOCKER_REFDIR, config.REF_GENOME_NAME),
                                "--in-bam", os.path.join(docker_input_parabricks, bam_filename),
                                '--out-variants', "{}/{}_pb_gatk.vcf".format(config.DOCKER_OUTPUTDIR, sample_name)])
            
            if os.system(command) != 0:
                raise Exception("Haplotypecaller failed!")

            sample.vcf_path_haplotypecaller = "{}/{}_pb_gatk.vcf".format(os.path.join(self.principal_directory, "temp"), sample_name)
            sample.saveJSON()
            
    def DeepVariant(self, samples):
        docker_input_parabricks = os.path.join(config.DOCKER_WORKDIR, 'bam')

        for sample in samples:
            sample_name = str(sample.name)
            try: 
                bam_filename = sample.bam.split('/')[-1]
            except Exception as e:
                raise("In variant calling pipe, BAM file coulnd not be found for sample {}".format(sample_name))

            command = ' '.join(['docker', 'run', '--gpus', 'all', '--rm', 
                                '--volume', "{}/:{}".format(config.REF, config.DOCKER_REFDIR), 
                                '--volume', "{}/:{}".format(self.principal_directory, config.DOCKER_WORKDIR), 
                                '--volume', "{}/:{}".format(os.path.join(self.principal_directory, "temp"), config.DOCKER_OUTPUTDIR), 
                                "{}".format(config.PARABRICKS_VERSION_DEEPVARIANT), 
                            'pbrun', 'deepvariant', 
                            '--ref', "{}/{}".format(config.DOCKER_REFDIR, config.REF_GENOME_NAME), 
                            "--in-bam", os.path.join(docker_input_parabricks, bam_filename),
                            '--out-variants', "{}/{}_pb_deepvariant.vcf".format(config.DOCKER_OUTPUTDIR, sample_name)])

            if os.system(command) != 0:
                raise Exception("Deepvariant failed!")
            
            sample.vcf_path_deepvariant = "{}/{}_pb_deepvariant.vcf".format(os.path.join(self.principal_directory, "temp"), sample_name)
            sample.saveJSON()
        
