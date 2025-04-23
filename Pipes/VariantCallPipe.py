from Pipes.Pipe import Pipe
from Pipes.ParallelPipe import ParallelPipe
import os
import csv
import subprocess
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
    """ Class responsible for running variant calling (pb_gatk and pb_deepvariant). 
    TODO: create a wrapper for this class. I.e. the wrapper should be responsible for handling the samples (receiving them from the
    directories, and supplying them to the this class). The VariantCallPipe should act on one sample at time, so that the convetion
    of each Pipe being responsible for one sample at a time is satisfied."""

    def __init__(self):
        super().__init__()

    def process(self, **kwargs):
        print("PROGRESS_FLAG:50% - Running Variant Calling (GATK HaplotypeCaller and DeepVariant)...", flush=True)
        self.principal_directory = dir_tree.principal_directory.path

        # carica tutti i sample
        sample_jsons = glob.glob(os.path.join(self.principal_directory, 'sample_data', '*.json'))
        samples = [Sample.fromJSON(jf) for jf in sample_jsons]
        self.panel = kwargs.get("panel", None)
        
        # 1) chiamata varianti
        self.HaplotypeCaller(samples)
        self.DeepVariant(samples)

        # 2) split & merge SNP/INDEL
        self._split_variants(samples)

        # aggiorna kwargs
        kwargs.update({
            'principal_directory': self.principal_directory,
            'samples': samples
        })
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
                                '--haplotypecaller-options', '"-A StrandBiasBySample -A DepthPerAlleleBySample"',
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
        
    def _split_variants(self, samples):
        """
        Splits raw VCFs into SNP and INDEL files for each caller, without merging.
        """
        temp_dir = os.path.join(self.principal_directory, 'temp')
        os.makedirs(temp_dir, exist_ok=True)
        ref_fa = os.path.join(config.REF, config.REF_GENOME_NAME)

        for sample in samples:
            base = sample.name
            hap_vcf = sample.vcf_path_haplotypecaller
            dv_vcf = sample.vcf_path_deepvariant

            # HaplotypeCaller SNP/INDEL
            hap_snp = os.path.join(temp_dir, f"{base}_hap_snp.vcf")
            hap_indel = os.path.join(temp_dir, f"{base}_hap_indel.vcf")
            subprocess.run([
                config.GATK, 'SelectVariants', '-R', ref_fa, '-V', hap_vcf,
                '--select-type-to-include', 'SNP', '-O', hap_snp
            ], check=True)
            subprocess.run([
                config.GATK, 'SelectVariants', '-R', ref_fa, '-V', hap_vcf,
                '--select-type-to-include', 'INDEL', '-O', hap_indel
            ], check=True)

            # # DeepVariant SNP/INDEL
            # dv_snp = os.path.join(temp_dir, f"{base}_dv_snp.vcf")
            # dv_indel = os.path.join(temp_dir, f"{base}_dv_indel.vcf")
            # subprocess.run([
            #     config.GATK, 'SelectVariants', '-R', ref_fa, '-V', dv_vcf,
            #     '--select-type-to-include', 'SNP', '-O', dv_snp
            # ], check=True)
            # subprocess.run([
            #     config.GATK, 'SelectVariants', '-R', ref_fa, '-V', dv_vcf,
            #     '--select-type-to-include', 'INDEL', '-O', dv_indel
            # ], check=True)

            # Update sample with split paths
            sample.vcf_hap_snp = hap_snp
            sample.vcf_hap_indel = hap_indel
            # sample.vcf_dv_snp = dv_snp
            # sample.vcf_dv_indel = dv_indel
            sample.saveJSON()