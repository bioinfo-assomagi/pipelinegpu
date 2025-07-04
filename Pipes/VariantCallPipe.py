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


        self.vcf_quality_filter(samples, self.principal_directory)
        self.strand_bias_filter(samples)

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
                                '--out-variants', "{}/{}_pb_gatk_raw.vcf".format(config.DOCKER_OUTPUTDIR, sample_name)])
            
            if os.system(command) != 0:
                raise Exception("Haplotypecaller failed!")

            sample.vcf_path_haplotypecaller_raw = "{}/{}_pb_gatk_raw.vcf".format(os.path.join(self.principal_directory, "temp"), sample_name)
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
            hap_vcf = sample.vcf_path_haplotypecaller_raw
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


    def vcf_quality_filter(self, samples, principal_directory):
         for sample in samples:
            # sample = str(sample_obj.name)
            temp_dir = os.path.join(principal_directory, 'temp')
            os.makedirs(temp_dir, exist_ok=True)
            ref_fa = os.path.join(config.REF, config.REF_GENOME_NAME)

            # Input VCFs from HaplotypeCaller split
            snp_vcf = getattr(sample, 'vcf_hap_snp', None)
            indel_vcf = getattr(sample, 'vcf_hap_indel', None)
            if not snp_vcf or not os.path.exists(snp_vcf):
                raise FileNotFoundError(f"SNP VCF not found for sample {sample.name}: {snp_vcf}")
            if not indel_vcf or not os.path.exists(indel_vcf):
                raise FileNotFoundError(f"INDEL VCF not found for sample {sample.name}: {indel_vcf}")

            # Filtration for SNP
            snp_filtered = os.path.join(temp_dir, f"{sample.name}_hap_snp.filtered.vcf")
            cmd_filt_snp = [
                config.GATK, 'VariantFiltration',
                '-R', ref_fa,
                '-V', snp_vcf,
                '-O', snp_filtered,
                '--filter-name', 'LowQD', '--filter-expression', 'QD < 2.0',
                '--filter-name', 'LowQUAL', '--filter-expression', 'QUAL < 30.0',
                '--filter-name', 'HighSOR', '--filter-expression', 'SOR > 3.0',
                '--filter-name', 'HighFS', '--filter-expression', 'FS > 60.0',
                '--filter-name', 'LowMQ', '--filter-expression', 'MQ < 40.0',
                '--filter-name', 'LowMQRankSum', '--filter-expression', 'MQRankSum < -12.5',
                '--filter-name', 'LowReadPosRankSum', '--filter-expression', 'ReadPosRankSum < -8.0'
            ]
            subprocess.run(cmd_filt_snp, check=True)
            sample.vcf_hap_snp_filtered = snp_filtered

            # Filtration for INDEL
            indel_filtered = os.path.join(temp_dir, f"{sample.name}_hap_indel.filtered.vcf")
            cmd_filt_indel = [
                config.GATK, 'VariantFiltration',
                '-R', ref_fa,
                '-V', indel_vcf,
                '-O', indel_filtered,
                '--filter-name', 'LowQD', '--filter-expression', 'QD < 2.0',
                '--filter-name', 'LowQUAL', '--filter-expression', 'QUAL < 30.0',
                '--filter-name', 'HighFS', '--filter-expression', 'FS > 200.0',
                '--filter-name', 'LowReadPosRankSum', '--filter-expression', 'ReadPosRankSum < -20.0'
            ]
            subprocess.run(cmd_filt_indel, check=True)
            sample.vcf_hap_indel_filtered = indel_filtered

            # Merge filtered SNP and INDEL into single VCF
            vcf_dir = os.path.join(principal_directory, 'vcf')
            os.makedirs(vcf_dir, exist_ok=True)
            merged_vcf = os.path.join(vcf_dir, f"{sample.name}_qualityfilter_gatk.vcf")
            cmd_merge = [
                config.GATK, 'MergeVcfs',
                '-I', snp_filtered,
                '-I', indel_filtered,
                '-O', merged_vcf
            ]
            subprocess.run(cmd_merge, check=True)
            sample.vcf_quality_filtered = merged_vcf
            sample.saveJSON()


    def strand_bias_filter(self, samples):
        """ THis is a prefilter of the outputted vcfs from the gatk_filter pipe, removing variants with strand bias.
        Actually, the HighSOR variants are written into a seperate vcf which will be later used for exclusive join with the main
        vcf file. """
        for sample in samples:
            vcf_quality_filtered = sample.vcf_quality_filtered
            df = pd.read_csv(vcf_quality_filtered, sep='\t', comment='#', names=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', str(sample.name)])
            df_filtered = df[~df["FILTER"].str.contains('HighSOR', na=False)]
            sample.vcf_path_haplotypecaller = "{}/{}_pb_gatk.vcf".format(os.path.join(self.principal_directory, "temp"), sample.name)

            
            df_filtered.to_csv(sample.vcf_path_haplotypecaller, sep='\t', index=False)
            sample.saveJSON()
        