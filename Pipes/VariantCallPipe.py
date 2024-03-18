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

class VariantCallPipe(ParallelPipe):
    path = os.getcwd()
    geno37 = join('/home/magi/', 'dataset/GENOME/37/all_chr37.fa')
    geno38 = join('/home/magi/', 'dataset/GENOME/38/all_chr38.fa')
    dbsnp144_37 = join('/home/magi/', 'dataset/dbsnp144/37/common_all_20150605_2.vcf')
    dbsnp144_38 = join('/home/magi/', 'dataset/dbsnp144/38/common_all_20150603_2.vcf')
    indel_37 = join('/home/magi/', 'dataset/dbsnp144/37/common_all_20150605_indel.vcf')
    indel_38 = join('/home/magi/', 'dataset/dbsnp144/38/common_all_20150603_indel.vcf')
    clinvar = join('/home/magi/', 'dataset/dbsnp144/38/clinvar_20140929_2.vcf')
    clinvar_indel = join('/home/magi/', 'dataset/dbsnp144/38/clinvar_20140929_indel.vcf')
    BUCHIARTIFICIALI = join('/home/magi/', 'PROJECT/diagnosys/bin/BUCHIARTIFICIALI.txt')

    def __init__(self):
        super().__init__()

    def process(self, **kwargs):
        self.principal_directory = kwargs.pop("principal_directory", None)
        samples: list = kwargs.pop("samples", None)

        print("TESTING VARIANT CALL PIPE: {}".format(samples))
        self.ParabricksHaplotypeCaller(samples)
        self.DeepVariant(samples)
        self.Filter(samples)

        kwargs.update(
            {"principal_directory": self.principal_directory, "samples": samples})
        return kwargs
    

    def ParabricksHaplotypeCaller(self, samples):

        docker_input_parabricks = os.path.join(config.DOCKER_WORKDIR, 'bam')

        for sample in samples:
            sample_name = str(sample['name'])
            bam_filename = sample['bam'].split('/')[-1]

            os.system(' '.join(['docker', 'run', '--gpus', 'all', '--rm',
                                '--volume', "{}/:{}".format(config.REF, config.DOCKER_REFDIR),
                                '--volume', "{}/:{}".format(self.principal_directory, config.DOCKER_WORKDIR),
                                '--volume',
                                "{}/:{}".format(os.path.join(self.principal_directory, "temp"), config.DOCKER_OUTPUTDIR),
                                "{}".format(config.PARABRICKS_VERSION),
                                'pbrun', 'haplotypecaller',
                                '--ref', "{}/Homo_sapiens_assembly38.fasta".format(config.DOCKER_REFDIR),
                                "--in-bam", os.path.join(docker_input_parabricks, bam_filename),
                                '--out-variants', "{}/{}_pb_gatk.vcf".format(config.DOCKER_OUTPUTDIR, sample_name)]))

            sample['vcf_path_haplotypecaller'] = "{}/{}_pb_gatk.vcf".format(os.path.join(self.principal_directory, "temp"), sample_name)
            
    def DeepVariant(self, samples):
        docker_input_parabricks = os.path.join(config.DOCKER_WORKDIR, 'bam')

        for sample in samples:
            sample_name = str(sample['name'])
            bam_filename = sample['bam'].split('/')[-1]

            os.system(' '.join(['docker', 'run', '--gpus', 'all', '--rm', 
                                '--volume', "{}/:{}".format(config.REF, config.DOCKER_REFDIR), 
                                '--volume', "{}/:{}".format(self.principal_directory, config.DOCKER_WORKDIR), 
                                '--volume', "{}/:{}".format(os.path.join(self.principal_directory, "temp"), config.DOCKER_OUTPUTDIR), 
                                "{}".format(config.PARABRICKS_VERSION_DEEPVARIANT), 
                            'pbrun', 'deepvariant', 
                            '--ref', "{}/Homo_sapiens_assembly38.fasta".format(config.DOCKER_REFDIR), 
                            "--in-bam", os.path.join(docker_input_parabricks, bam_filename),
                            '--out-variants', "{}/{}_pb_deepvariant.vcf".format(config.DOCKER_OUTPUTDIR, sample_name)]))

            sample['vcf_path_deepvariant'] = "{}/{}_pb_deepvariant.vcf".format(os.path.join(self.principal_directory, "temp"), sample_name)
        
    def Filter(self, samples):
        for sample in samples:
           
            filter_intronic_vcf_haplotypecaller, filter_cds_vcf_haplotypecaller = self.vcf_vertical_filter(sample, 'haplotypecaller')
            filter_intronic_vcf_deepvariant, filter_cds_vcf_deepvariant = self.vcf_vertical_filter(sample, 'deepvariant')
            
            filter_cds_vcf_haplotypecaller = self.buchiartificiali_filter(filter_cds_vcf_haplotypecaller)
            filter_cds_vcf_deepvariant = self.buchiartificiali_filter(filter_cds_vcf_deepvariant)
            
                  
            filter_intronic_vcf_haplotypecaller.to_csv(join(self.principal_directory, 'vcf/', sample['name']) + '_unfied_all.vcf', sep='\t', index=False)
            filter_cds_vcf_haplotypecaller.to_csv(join(self.principal_directory, 'vcf/', sample['name']) + '_unfied_only_CDS.vcf', sep='\t', index=False)
            
            filter_intronic_vcf_deepvariant.to_csv(join(self.principal_directory, 'vcf/', sample['name']) + '_samt_all.vcf', sep='\t', index=False)
            filter_cds_vcf_deepvariant.to_csv(join(self.principal_directory, 'vcf/', sample['name']) + '_samt_only_CDS.vcf', sep='\t', index=False)
            
            # TODO: add them to samples
                                     
    def buchiartificiali_filter(self, filtered_cds_vcf):
        buchiartificiali = pd.read_csv(self.BUCHIARTIFICIALI, sep='\t', header=0)
        
        for index, row in buchiartificiali.iterrows():
            chrom_number = row['#CHROM']
            start = row['START']
            end = row['END']
            mask = ~((filtered_cds_vcf["#CHROM"] == chrom_number) & (filtered_cds_vcf['POS'] >= start) & (filtered_cds_vcf['POS'] <= end))
            filtered_cds_vcf = filtered_cds_vcf[mask]
        
        return filtered_cds_vcf           
            
    def vcf_vertical_filter(self, sample, vcf_type):
        vertical = sample['vertical']
        if type(vertical) is str:
            vertical = pd.read_csv(vertical, sep='\t')
        #verticalX = sample['verticalx']
        
        if vcf_type == 'haplotypecaller':
            vcf = pd.read_csv(sample['vcf_path_haplotypecaller'], sep='\t', comment='#', names=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', str(sample['name'])])
        elif vcf_type == 'deepvariant':
            vcf = pd.read_csv(sample['vcf_path_deepvariant'], sep='\t', comment='#', names=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', str(sample['name'])])
        else:
            raise Exception("VCF_VERTICAL_FILTER: Invalid vcf_type option.")    
        
        filter_intronic = pd.merge(vcf, vertical, on=['#CHROM', 'POS'], how='left')
        filter_cds = pd.merge(vcf, vertical, on=['#CHROM', 'POS'], how='inner')
        
        return filter_intronic, filter_cds



        # # Instead of generating the bed and vertical, we read them, as they should already be present in their respective folders (generated by CoveragePipe)
        # bed_ = join(folder_pheno, 'bed_' + sample_name)
        # bed = pd.read_csv(bed_, sep='\t')
        # bedx_ = join(folder_pheno, 'bedX_' + sample_name)
        # bedx = pd.read_csv(bedx_, sep='\t')

        # vertical_ = join(folder_pheno, 'vertical_' + sample_name)
        # vertical = pd.read_csv(vertical_, sep='\t')
        # verticalx_ = join(folder_pheno, 'verticalX_' + sample_name)
        # verticalX = pd.read_csv(verticalx_, sep='\t')

        # if not vertical.empty:
        #     buchiartificiali = pd.read_csv(self.BUCHIARTIFICIALI,sep='\t',header=0)
        #     if 'OCULAR' in panel: _vertical_macro=join(path,'bin','VERTICAL','vertical_ONA20509.2020')
        #     elif 'VASCULAR' in panel:_vertical_macro=join(path,'bin','VERTICAL','vertical_VASCNA20509.2020')
        #     elif 'NEUROLOGY' in panel: _vertical_macro=join(path,'bin','VERTICAL','vertical_NEURNA20509.2020')
        #     elif 'INTEGRACARDIOSTANCHEZZA' in panel: _vertical_macro=join(path,'bin','VERTICAL','vertical_INTEGRACARDIOSTANCHEZZA.2023')
        #     elif 'MIXED' in panel: _vertical_macro=join(path,'bin','VERTICAL','vertical_MIX120509.2020')
        #     elif 'LYMPHOBESITY' in panel: _vertical_macro=join(path,'bin','VERTICAL','vertical_LIMPHNA20509.2020')
        #     elif 'INFERTILIT' in panel: _vertical_macro=join(path,'bin','VERTICAL','vertical_INA20509.2020')
        #     elif 'CHERATOCONO' in panel: _vertical_macro=join(path,'bin','VERTICAL','vertical_CHERATOCONO.2020')
        #     elif 'PCDH19' in panel: _vertical_macro=join(path,'bin','VERTICAL','vertical_PCDH19')
        #     elif 'GENEOBNEW' in panel: _vertical_macro=join(path,'bin','VERTICAL','vertical_GENEOBNEW')
        #     elif 'CUSTOM' in panel: _vertical_macro=join(folder_pheno,'vertical_' + sample_name)
        #     elif 'CANCER' in panel:
        #         print ('CANCERRRRRRRRRRRRRRRRRRRRRRRR',join(path,'bin','VERTICAL','vertical_CANC20509.2022'))
        #         _vertical_macro=join(path,'bin','VERTICAL','vertical_CANC20509.2022')
        #     else: _vertical_macro = ''  #'' #print('VERTICAL not found')

        #     try: vertical_macro =  pd.read_csv(_vertical_macro, sep='\t')
        #     except: print ('ERRORE!!!') #verti
        #     sample = vcf_from_samtools(sample, quality, genome, vertical, verticalX, principal_directory, phenotype, vertical_macro, buchiartificiali)
        #     name = join(principal_directory,'temp/',sample+'_final_disease')
        #     GATK_unfied(sample, genome, gatk, quality, threads, principal_directory, dbsnp, verticalX, buchiartificiali)
        #     #plot(args,folder_name,sample)