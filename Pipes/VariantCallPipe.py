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

# TODO: if files not found inside samples list, try reading them from principal_directory
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
        self.panel = kwargs.pop("panel", None)

        print("TESTING VARIANT CALL PIPE: {}".format(samples))
        self.ParabricksHaplotypeCaller(samples)
        self.DeepVariant(samples)
        self.Filter(samples)
        self.Merge(samples)

        kwargs.update(
            {"principal_directory": self.principal_directory, "samples": samples, "panel": self.panel})
        return kwargs
    

    def ParabricksHaplotypeCaller(self, samples):

        docker_input_parabricks = os.path.join(config.DOCKER_WORKDIR, 'bam')

        for sample in samples:
            sample_name = str(sample['name'])
            bam_filename = sample['bam'].split('/')[-1]

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

            sample['vcf_path_haplotypecaller'] = "{}/{}_pb_gatk.vcf".format(os.path.join(self.principal_directory, "temp"), sample_name)
            
    def DeepVariant(self, samples):
        docker_input_parabricks = os.path.join(config.DOCKER_WORKDIR, 'bam')

        for sample in samples:
            sample_name = str(sample['name'])
            bam_filename = sample['bam'].split('/')[-1]

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
            
            sample['vcf_path_deepvariant'] = "{}/{}_pb_deepvariant.vcf".format(os.path.join(self.principal_directory, "temp"), sample_name)
        
    """ Filter exones. Note, the dataframes named intronic contain all regions, exones + introns. """
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

            sample['intronic_vcf_haplotypecaller'] = filter_intronic_vcf_haplotypecaller
            sample['cds_vcf_haplotypecaller'] = filter_cds_vcf_haplotypecaller
            sample['intronic_vcf_deepvariant'] = filter_intronic_vcf_deepvariant
            sample['cds_vcf_deepvariant'] = filter_cds_vcf_deepvariant
            
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
    
    def Merge(self, samples):
        
        temp_directory_path = os.path.join(self.principal_directory, "temp/")
        vcf_directory_path = os.path.join(self.principal_directory, "vcf/")
        annotation_directory_path = os.path.join(self.principal_directory, "annotation/")

        for sample in samples:
            
            sample_name = str(sample["name"])

            cds_path = os.path.join(temp_directory_path, sample_name + "_temp_CDS_annot.csv")
            all_path = os.path.join(vcf_directory_path, sample_name + "_temp_all_annot.csv")

            cds_samtools = sample['cds_vcf_deepvariant']
            cds_gatk = sample['cds_vcf_haplotypecaller']
            all_samtools = sample['intronic_vcf_deepvariant']
            all_gatk = sample['intronic_vcf_haplotypecaller']

            # Just in case samples will contain the paths to the vcf instead of the dataframe
            if type(cds_samtools) is str:
                cds_samtools = pd.read_csv(cds_samtools, sep='\t')
            
            if type(cds_gatk) is str:
                cds_gatk = pd.read_csv(cds_gatk, sep='\t')
            
            if type(all_samtools) is str:
                all_samtools = pd.read_csv(all_samtools, sep='\t')
            
            if type(all_gatk) is str:
                all_gatk = pd.read_csv(all_gatk, sep='\t')

            all_samtools['refseq'] = all_samtools['refseq'].astype(str)
            all_samtools['hgmd'] = all_samtools['hgmd'].astype(str)
            all_samtools['GENE'] = all_samtools['GENE'].astype(str)
            all_gatk['refseq'] = all_gatk['refseq'].astype(str)
            all_gatk['hgmd'] = all_gatk['hgmd'].astype(str)
            all_gatk['GENE'] = all_gatk['GENE'].astype(str)

            CDS = pd.merge(cds_samtools, cds_gatk, on=['#CHROM', 'POS', 'ID', 'GENE', 'exone', 'length', 'strand',
								'refseq', 'hgmd'], how='outer', suffixes=('_samt', '_gatk'))
            
            ALL = pd.merge(all_samtools, all_gatk, on=['#CHROM', 'POS', 'ID', 'GENE', 'exone', 'length', 'strand',
								'refseq', 'hgmd'], how='outer', suffixes=('_samt', '_gatk'))
            
            CDS['QUAL_samt'].fillna(-999,inplace=True)
            CDS['QUAL'] = CDS['QUAL_samt'].astype(int)
            CDS['INFO'] = CDS['INFO_samt']
            CDS['FILTER'] = CDS['FILTER_gatk']
            CDS.loc[CDS['FILTER'].isnull(),'FILTER'] = '.'
            CDS['REF'] = CDS['REF_gatk']
            CDS.loc[CDS['REF'].isnull(),'REF'] = CDS['REF_samt']
            CDS['ALT'] = CDS['ALT_gatk']
            CDS.loc[CDS['ALT'].isnull(),'ALT'] = CDS['ALT_samt']

            CDS['ALT'] = np.where(CDS['ALT'].str.split(',').str.get(0) == '*',
					CDS['ALT'].str.split(',').str.get(1),
					CDS['ALT'].str.split(',').str.get(0))
            CDS = CDS[CDS['ALT']!='*']
            mask1 = CDS['INFO'].isnull()
            CDS.loc[mask1,'INFO'] = CDS['INFO_gatk']
            CDS['FORMAT'] = 'GT:PL:DP'

            CDS.drop(['QUAL_samt','INFO_samt','FORMAT_samt'], axis=1,inplace=True)

            ALL['QUAL_samt'].fillna(-999,inplace=True)
            ALL['QUAL'] = ALL['QUAL_samt'].astype(int)
            ALL['FILTER'] = ALL['FILTER_gatk']
            ALL.loc[ALL['FILTER'].isnull(),'FILTER'] = '.'
            ALL['INFO'] = ALL['INFO_samt']
            ALL['REF'] = ALL['REF_gatk']
            ALL.loc[ALL['REF'].isnull(),'REF'] = ALL['REF_samt']
            ALL['ALT'] = ALL['ALT_gatk']
            ALL.loc[ALL['ALT'].isnull(),'ALT'] = ALL['ALT_samt']
            ALL['ALT'] = np.where(ALL['ALT'].str.split(',').str.get(0) == '*',
                                ALL['ALT'].str.split(',').str.get(1),
                                ALL['ALT'].str.split(',').str.get(0))
            ALL = ALL[ALL['ALT']!='*']
            
            maskalt = ALL['ALT'] != '*'
            ALL.loc[maskalt,'ALT'] = ALL['ALT']
            mask1 = ALL['INFO'].isnull()
            ALL.loc[mask1,'INFO'] = ALL['INFO_gatk']
            ALL['FORMAT'] = 'GT:PL:DP'
            ALL.drop(['QUAL_samt','INFO_samt','FORMAT_samt'], axis=1,inplace=True)
            CDS.drop(['QUAL_gatk','INFO_gatk','FORMAT_gatk'], axis=1,inplace=True)
            ALL.drop(['QUAL_gatk','INFO_gatk','FORMAT_gatk'], axis=1,inplace=True)
            ALL['ALT'] = ALL['ALT'].str.split(',').str.get(0)
            CDS['ALT'] = CDS['ALT'].str.split(',').str.get(0)
            CDS['POS'] = CDS['POS'].astype(int)
            ALL['POS'] = ALL['POS'].astype(int)
            ALL['DEPTH'] = ALL['INFO'].str.split('DP=').str.get(1).str.split(';').str.get(0)

            samtools = sample_name + '_samt'
            gatk = sample_name + '_gatk'

            if ((self.panel == 'ALLGENE') or (self.panel == 'trusightone') or (self.panel == 'CANCER')):
                ALL = ALL[( pd.to_numeric(ALL['DEPTH'])>=10)]
                # ALL = ALL[( pd.to_numeric(ALL['DEPTH'])>=20)]
                ALL = ALL[ pd.to_numeric(ALL['QUAL'])>=18]
            else:
                ALL = ALL[( pd.to_numeric(ALL['DEPTH'])>=1) | ( pd.to_numeric(ALL['DEPTH'])==1)]
                #ALL = ALL[ALL['QUAL']>=18]

            CDS = CDS[['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT', samtools, gatk,
					'GENE','exone','length','strand','refseq','hgmd']]
            ALL = ALL[['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT', samtools, gatk,
                        'GENE','exone','length','strand','refseq','hgmd']]
            
            CDS.iloc[:,:].fillna('.',inplace=True)
            ALL.iloc[:,:].fillna('.',inplace=True)
            CDS[samtools].fillna('.',inplace=True)
            ALL[samtools].fillna('.',inplace=True)
            CDS[gatk].fillna('.',inplace=True)
            ALL[gatk].fillna('.',inplace=True)
            CDS.to_csv(cds_path, sep='\t', header=True, index=False,encoding='utf-8')
            ALL.to_csv(all_path, sep='\t', header=True, index=False,encoding='utf-8')
            vcf_annotate_all = join(annotation_directory_path, sample_name + '_all_annot.vcf')
            vcf_annotate_CDS = join(annotation_directory_path, sample_name + '_CDS_annot.vcf')
            vcf_annotate_CDS2 = join(annotation_directory_path, sample_name + '_CDS_annot_dbNSFP.vcf')

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