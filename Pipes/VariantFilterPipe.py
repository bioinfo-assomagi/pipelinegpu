from Pipes.Pipe import Pipe
from Pipes.ParallelPipe import ParallelPipe
import os


import config

import utils
from os.path import join
import glob
import pandas as pd
import numpy as np

import dir_tree
from Entities.Sample import Sample

# TODO: if files not found inside samples list, try reading them from principal_directory
class VariantFilterPipe(ParallelPipe):

    BUCHIARTIFICIALI = join('/home/magi/', 'PROJECT/diagnosys/bin/BUCHIARTIFICIALI.txt')

    def __init__(self):
        super().__init__()

    def process(self, **kwargs):
        self.thread_id = os.getpid()

        self.sample = kwargs.pop("sample")
        self.panel = kwargs.pop("panel", None)
        self.genome_type = kwargs.pop("genome", "geno38")

        self.thread_print("TESTING VARIANT Filter PIPE: {}".format(self.sample.name))

        self.Filter()
        self.Merge()
        self.VEP()

        self.sample.saveJSON()
        
        kwargs.update(
            {"sample": self.sample, "panel": self.panel, "genome": self.genome_type})
        
        return kwargs
      
    
    def Filter(self):
        """ Filter the vcf files with the regions present in the BED file (exones/CDS).
            Additionally it removes variants that fall in buchi regions (regions specified in the buchiartificiali file.)
            Note, the dataframes named intronic contain all regions, exones + introns. 

            Returns:
                pd.DataFrame: fitered variants 
        """

            
        if not hasattr(self.sample, "vcf_path_haplotypecaller"):
            raise Exception("No haplotypecaller output registered for this sample!")
        elif not  os.path.exists(self.sample.vcf_path_haplotypecaller):
            raise Exception("The registered haplotypecaller output file doesn't exist in the project directory!")
        
        if not hasattr(self.sample, "vcf_path_deepvariant"):
            raise Exception("No deepvariant output registered for this sample!")
        elif not  os.path.exists(self.sample.vcf_path_deepvariant):
            raise Exception("The registered deepvariant output file doesn't exist in the project directory!")
        
        filter_intronic_vcf_haplotypecaller, filter_cds_vcf_haplotypecaller = self.vcf_vertical_filter('haplotypecaller')
        filter_intronic_vcf_deepvariant, filter_cds_vcf_deepvariant = self.vcf_vertical_filter('deepvariant')
        
        filter_cds_vcf_haplotypecaller = self.buchiartificiali_filter(filter_cds_vcf_haplotypecaller)
        filter_cds_vcf_deepvariant = self.buchiartificiali_filter(filter_cds_vcf_deepvariant)

        # Remove 0/0 from deepvariant results
        filter_cds_vcf_deepvariant['genotype'] = filter_cds_vcf_deepvariant[str(self.sample.name)].str.split(':').str[0]
        filter_intronic_vcf_deepvariant['genotype'] = filter_intronic_vcf_deepvariant[str(self.sample.name)].str.split(':').str[0]
        self.thread_print("Removing low quality vcfs from deepvariant ... ")
        filter_cds_vcf_deepvariant = filter_cds_vcf_deepvariant[(~filter_cds_vcf_deepvariant['genotype'].isin(['./.', '0/0'])) & (filter_cds_vcf_deepvariant['QUAL'] != 0.0)]
        filter_intronic_vcf_deepvariant = filter_intronic_vcf_deepvariant[(~filter_intronic_vcf_deepvariant['genotype'].isin(['./.', '0/0'])) & (filter_intronic_vcf_deepvariant['QUAL'] != 0.0)]
        
        filter_cds_vcf_deepvariant = filter_cds_vcf_deepvariant.drop('genotype', axis = 1)
        filter_intronic_vcf_deepvariant = filter_intronic_vcf_deepvariant.drop('genotype', axis = 1)


        vcf_paths = utils.define_vcf_paths(self.sample.name, dir_tree)

        filter_intronic_vcf_haplotypecaller.to_csv(vcf_paths["haplotypecaller_intronic_vcf_path"], sep='\t', index=False)
        filter_cds_vcf_haplotypecaller.to_csv(vcf_paths["haplotypecaller_cds_vcf_path"], sep='\t', index=False)
        
        filter_intronic_vcf_deepvariant.to_csv(vcf_paths["deepvariant_intronic_vcf_path"], sep='\t', index=False)
        filter_cds_vcf_deepvariant.to_csv(vcf_paths["deepvariant_cds_vcf_path"], sep='\t', index=False)

        # TODO: add them to samples
        self.sample.intronic_vcf_haplotypecaller = vcf_paths["haplotypecaller_intronic_vcf_path"]
        self.sample.cds_vcf_haplotypecaller = vcf_paths["haplotypecaller_cds_vcf_path"]
        self.sample.intronic_vcf_deepvariant = vcf_paths["deepvariant_intronic_vcf_path"]
        self.sample.cds_vcf_deepvariant = vcf_paths["deepvariant_cds_vcf_path"]
        
      
                                     
    def buchiartificiali_filter(self, filtered_cds_vcf):
        """ Remove from the VCF file, variants that fall in positions that are present in the buchiartificiali file.

                Args:
                    filtered_cds_vcf (pandas.DataFrame): dataframe containing the variants filtered by the BED file

                Returns:
                    pandas.DataFrame: dataframe with variants that fall outside the regions present in buchiartificiali. Doesn't write to any file.
        """
        buchiartificiali = pd.read_csv(self.BUCHIARTIFICIALI, sep='\t', header=0)
        
        for index, row in buchiartificiali.iterrows():
            chrom_number = row['#CHROM']
            start = row['START']
            end = row['END']
            mask = ~((filtered_cds_vcf["#CHROM"] == chrom_number) & (filtered_cds_vcf['POS'] >= start) & (filtered_cds_vcf['POS'] <= end))
            filtered_cds_vcf = filtered_cds_vcf[mask]
        
        return filtered_cds_vcf           
            
    # def intersect(self, sample, vcf_type):
    #     import pybedtools
    #     a = pybedtools.BedTool(sample.vcf)
    #     b = pybedtools.BedTool(sample.bed)
    #     intersection = a.intersect(b, u=True)
        

    def vcf_vertical_filter(self, vcf_type):
        """ Intersect withthe BED and BEDX files, which are transformed to verticals, so pandas.merge can be used.
        
                Returns:
                    pandas.Dataframe: the filtered dataframe. Doesn't write to any file.
        """
        sample = self.sample

        vertical = sample.verticalX
        if type(vertical) is str:
            vertical = pd.read_csv(vertical, sep='\t')
        #verticalX = sample['verticalx']
        
        if vcf_type == 'haplotypecaller':
            vcf = pd.read_csv(sample.vcf_path_haplotypecaller, sep='\t', comment='#', names=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', str(sample.name)])
        elif vcf_type == 'deepvariant':
            vcf = pd.read_csv(sample.vcf_path_deepvariant, sep='\t', comment='#', names=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', str(sample.name)])
        else:
            raise Exception("VCF_VERTICAL_FILTER: Invalid vcf_type option.")    
        
        filter_intronic = pd.merge(vcf, vertical, on=['#CHROM', 'POS'], how='left')
        filter_cds = pd.merge(vcf, vertical, on=['#CHROM', 'POS'], how='inner')
        
        return filter_intronic, filter_cds
    
    def Merge(self):
        
        temp_directory_path = dir_tree.principal_directory.temp.path
        vcf_directory_path = dir_tree.principal_directory.vcf.path
        annotation_directory_path = dir_tree.principal_directory.annotation.path

        sample = self.sample
            
        sample_name = str(sample.name)

        cds_path = os.path.join(temp_directory_path, sample_name + "_temp_CDS_annot.csv")
        all_path = os.path.join(vcf_directory_path, sample_name + "_temp_all_annot.csv")

        cds_samtools = sample.cds_vcf_deepvariant
        cds_gatk = sample.cds_vcf_haplotypecaller
        all_samtools = sample.intronic_vcf_deepvariant
        all_gatk = sample.intronic_vcf_haplotypecaller

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
        
        CDS['QUAL_gatk'].fillna(-999,inplace=True)
        CDS['QUAL'] = CDS['QUAL_gatk'].astype(int)
        CDS['INFO'] = CDS['INFO_gatk']
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
        #CDS.loc[mask1,'INFO'] = CDS['INFO_gatk']
        
        CDS.loc[mask1, 'INFO'] = CDS.loc[mask1].apply(
            lambda row: ';'.join(
                [f"{key}={value}" for key, value in zip(row['FORMAT_samt'].split(':'), row["{}_samt".format(sample_name)].split(':'))]
                                ), axis=1)

        CDS['FORMAT'] = 'GT:PL:DP'


        ALL['QUAL_gatk'].fillna(-999,inplace=True)
        ALL['QUAL'] = ALL['QUAL_gatk'].astype(int)
        ALL['INFO'] = ALL['INFO_gatk']
        ALL['FILTER'] = ALL['FILTER_gatk']
        ALL.loc[ALL['FILTER'].isnull(),'FILTER'] = '.'
        ALL['REF'] = ALL['REF_gatk']
        ALL.loc[ALL['REF'].isnull(),'REF'] = ALL['REF_samt']
        ALL['ALT'] = ALL['ALT_gatk']
        ALL.loc[ALL['ALT'].isnull(),'ALT'] = ALL['ALT_samt']
        ALL['ALT'] = np.where(ALL['ALT'].str.split(',').str.get(0) == '*',
                            ALL['ALT'].str.split(',').str.get(1),
                            ALL['ALT'].str.split(',').str.get(0))
        ALL = ALL[ALL['ALT'] != '*']

        mask1 = ALL['INFO'].isnull()
        #ALL.loc[mask1,'INFO'] = ALL['INFO_gatk']
        ALL.loc[mask1, 'INFO'] = ALL.loc[mask1].apply(
             lambda row: ';'.join(
                [f"{key}={value}" for key, value in zip(row['FORMAT_samt'].split(':'), row["{}_samt".format(sample_name)].split(':'))]
                                  ), axis=1)

        ALL['FORMAT'] = 'GT:PL:DP'
        ALL.drop(['QUAL_samt','INFO_samt','FORMAT_samt'], axis=1,inplace=True)
        ALL.drop(['QUAL_gatk','INFO_gatk','FORMAT_gatk'], axis=1,inplace=True)
        CDS.drop(['QUAL_samt','INFO_samt','FORMAT_samt'], axis=1,inplace=True)
        CDS.drop(['QUAL_gatk','INFO_gatk','FORMAT_gatk'], axis=1,inplace=True)
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
        sample.cds_path = cds_path
        sample.all_path = all_path
        sample.vcf_annot_CDS = vcf_annotate_CDS


    def VEP(self):
        annotation_dir = dir_tree.principal_directory.annotation.path 
        

        if self.genome_type == 'geno38':
            os.system(' '.join([config.VEPS,'--cache --refseq --offline --use_given_ref --assembly GRCh38 --fasta', config.GENO38,
				'--fork 8 --force_overwrite --vcf --format vcf --everything --af_1kg',
				'--buffer_size 500 --force --xref_refseq --exclude_predicted --use_transcript_ref',
				'--plugin dbNSFP,%s,'
				'Eigen-PC-raw_coding_rankscore,CADD_raw_rankscore,DANN_rankscore,MetaSVM_rankscore,'
				'MetaLR_rankscore,FATHMM_converted_rankscore,MutationTaster_converted_rankscore,'
				'MutationAssessor_rankscore,Polyphen2_HDIV_rankscore,Polyphen2_HVAR_rankscore,'
				'SIFT_converted_rankscore,LRT_converted_rankscore,MutPred_rankscore,'
				'PROVEAN_converted_rankscore,VEST4_rankscore,REVEL_rankscore,'
				'SIFT_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,LRT_pred,MutationTaster_pred,MutationAssessor_pred,FATHMM_pred,'
				'PROVEAN_pred,MetaSVM_pred,MetaLR_pred,M-CAP_pred,fathmm-MKL_coding_pred,'
				'GERP++_RS_rankscore,phyloP100way_vertebrate_rankscore,phyloP30way_mammalian_rankscore,'
				'phyloP30way_mammalian_rankscore,phastCons100way_vertebrate_rankscore,'
				'phastCons30way_mammalian_rankscore,SiPhy_29way_logOdds_rankscore,'
				'clinvar_clnsig,clinvar_review,Interpro_domain,'
				'gnomAD_exomes_POPMAX_AF,gnomAD_exomes_POPMAX_nhomalt,gnomAD_exomes_controls_AF,'
				'gnomAD_genomes_POPMAX_AF,gnomAD_genomes_POPMAX_nhomalt,'
				'gnomAD_exomes_controls_AC,'
				'clinvar_review,clinvar_MedGen_id,clinvar_OMIM_id,clinvar_Orphanet_id,gnomAD_exomes_controls_nhomalt,',
				'--plugin dbscSNV,%s,GRCh38,ada_score,rf_score',
				'-i',self.sample.cds_path,'-o', self.sample.vcf_annot_CDS]) % (config.dbNSFP38_gz, config.dbscSNV11_gz))