from Pipes.Pipe import Pipe
import os
import pandas as pd
import config
import glob
import dir_tree
import numpy as np
from Bio.Seq import Seq

""" Pipe that will modify the results of VEP. Responsible for creating _CDS_annot.csv"""
class AnnotationPipe(Pipe):

    eredita = None
    hgmd = None
    appris = None
    eccezioni = None

    def __init__(self):
        pass

    def process(self, **kwargs):
        self.principal_directory = kwargs.pop("principal_directory", None)
        #self.samples = kwargs.pop("samples", None)
        self.sample = kwargs.pop("sample")
        self.genome_type = kwargs.pop("genome", None)
        self.panel = kwargs.pop("panel", None)

        self.input_vcf = os.path.join(self.principal_directory, "vcf/")
        self.output_annotation = os.path.join(self.principal_directory, "annotation/")
        self.output_final = os.path.join(self.principal_directory, "final/")
        self.folder_temp = os.path.join(self.principal_directory, "temp/")
        self.folder_COV = os.path.join(self.principal_directory, "coverage/")
        self.input_phenotype = os.path.join(self.principal_directory, "pheno/phenotype")

        # Signature of this pipe TODO: ...
        self.ext = "_CDS_annot.csv"

        self.run()

        self.sample.saveJSON()

        kwargs.update({"principal_directory": self.principal_directory, "sample": self.sample, "genome_type": self.genome_type, "panel": self.panel})
        return kwargs

    def filter_vep_refseq(self):
        APPRIS = pd.read_csv(config.APPRIS, sep="\t")
        ECCEZIONI = pd.read_csv(config.ECCEZIONI, sep="\t")

        APPRIS['nm_refseq'] = APPRIS['refseq'].str.split(".").str.get(0)
        ECCEZIONI['nm_refseq'] = ECCEZIONI['refseq'].str.split(".").str.get(0)

        folder_COV = dir_tree.principal_directory.coverage.path
        sample_name = str(self.sample.name)
        coverage = pd.read_csv(os.path.join(folder_COV, sample_name, sample_name + "_all"), sep='\t', header=0)
        if self.genome_type == "geno38":
            eredita = pd.read_csv(config.EREDITA38, sep="\t", header=0)
        elif self.genome_type == "geno37":
            eredita = pd.read_csv(config.EREDITA37, sep="\t", header=0)

        coverage = coverage[['#CHROM','POS','C%','G%','T%','A%','ins%','del%','sum']]
        CDS = pd.read_csv(self.sample.vcf_annot_CDS, sep="\t", header=None, comment='#', names=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 
                                                                                                sample_name + '_samt', sample_name + '_gatk', 'GENE', 'exone', 'length', 'strand', 'refseq', 'hgmd'])

        cds_coverage = pd.merge(CDS, coverage, on=["#CHROM", "POS"], how='left')
        cds_coverage_eredita = pd.merge(cds_coverage, eredita, on=['GENE'], how='left')

        result = cds_coverage_eredita
        # HERE, HERE IS WHERE YOU SET INFO1 AND INFO2

        info_split = result['INFO'].str.split('CSQ=')
        result["INFO1"] = info_split.str.get(0)
        result["INFO2"] = info_split.str.get(1).str.split(",")

        # Explode results to unique NM_s
        exploded_results = result.explode("INFO2")
        exploded_results["nm_refseq"] = exploded_results["INFO2"].str.split("|").str.get(6).str.split(".").str.get(0)

        # Change strategy ... second ....
        # Merge exploded results with APPRIS and ECCEZIONI
        exploded_results_appris = pd.merge(exploded_results, APPRIS, how="inner", on=["GENE", "nm_refseq"], suffixes=("_result", "_appris")) # shouldn't expect suffixes (except refseq_appris), since APPRIS doesn't have any in common column names except GENE and nm_refseq
        exploded_results_eccezioni = pd.merge(exploded_results, ECCEZIONI, how="inner", on=["GENE", "nm_refseq"], suffixes=("_result", "_eccezioni")) # only refseq column name is expected to be in common

        # Now merge those two together, outer join, and you will get single row per variant, with the corresponding refseqs that were found
        exploded_results_appris_eccezioni = pd.merge(exploded_results_appris, exploded_results_eccezioni, on=["#CHROM", "POS"], how="outer", suffixes=("_appris", "_eccezioni"))

        # Create the INFO column in the "processed results" - careful, there are no null here since there already is an INFO in exploded_results
        exploded_results_appris_eccezioni["INFO_processed"] = exploded_results_appris_eccezioni["INFO1_appris"] + "CSQ=" + exploded_results_appris_eccezioni["INFO2_appris"]
        exploded_results_appris_eccezioni["INFO_processed"].fillna(exploded_results_appris_eccezioni["INFO1_eccezioni"] + "CSQ=" + exploded_results_appris_eccezioni["INFO2_eccezioni"], inplace=True)

        # When joining with the original results, the info column of the processed will have _processed
        final_result2 = pd.merge(result, exploded_results_appris_eccezioni, how="left", on=["#CHROM", "POS"], suffixes=("", "_processed"))
        final_result2["INFO_processed"].fillna(final_result2["INFO"].str.split(",").str.get(0), inplace=True) # what if there is no "|," i.e. only one refseq
        final_result2["INFO"] = final_result2["INFO_processed"]

        cleaned_result = final_result2[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
            sample_name + '_samt', sample_name + '_gatk', 'GENE', 'exone', 'length', 'strand',
            'refseq', 'hgmd', 'C%', 'G%', 'T%', 'A%', 'ins%', 'del%', 'sum',
            'INHERITANCE', 'VERBOSE']]
        

        return cleaned_result

    def adapt_annotation(self, result):
        name = str(self.sample.name)
        gatk = name + "_gatk"
        samtools = name + "_samt"
        folder_final = dir_tree.principal_directory.final.path
        
        cols = ['sample','HGVS','CHROM','POS','GENE','ALTGENE','ID','types','QUAL','DEPTH','mapquality','INHERITANCE',
		'samtools_geno','gatk_geno','unbalance','consequence','impact','clin_sign','strand',
		'refseq','hgmd','num_exon','num_intron','HGVS_c','HGVS_p','CDS_position',
		'amino_acids','codons','variation_cosmic','sift','polyphen','GMAF','ExAC_MAF','MAX_MAF','Adj_MAF','AFR_MAF','AMR_MAF','EAS_MAF',
		'EUR_MAF','SAS_MAF','AA_MAF','EA_MAF','pubmed',
		'hgmd_mutation','hgmd_function','hgmd_phenotype','hgmd_pubmed',

		'CADD_rankscore','DANN_rankscore','EigenPC_rankscore','FATHMM_rankscore','FATHMM_pred','GERP_rankscore','Interpro_domain',
		'LRT_rankscore','LRT_pred','MCAP_pred','MetaLR_pred','MetaLR_rankscore','MetaSVM_pred','MetaSVM_rankscore',
		'MutPred_rankscore','MutationAssessor_pred','MutationAssessor_rankscore','MutationTaster_rankscore','MutationTaster_pred',
		'PROVEAN_rankscore','PROVEAN_pred','Polyphen2HDIV_pred','Polyphen2HDIV_rankscore','Polyphen2HVAR_pred',
		'Polyphen2HVAR_rankscore','REVEL_rankscore','SIFT_rankscore','SIFT_pred','SiPhy29way_rankscore',
		'VEST3_rankscore','clinvar_clnsig','fathmmMKL_pred',
		'phastCons100way_rankscore','phastCons20way_rankscore','phyloP100way_rankscore','phyloP20way_rankscore',

		'ada_score','rf_score',

		'clinvar_MedGen_id','clinvar_OMIM_id','clinvar_Orphanet_id','clinvar_clnsig','clinvar_review',
		'gnomAD_exomes_POPMAX_AF','gnomAD_exomes_POPMAX_nhomalt','gnomAD_exomes_controls_AF',
		'gnomAD_exomes_controls_nhomalt','gnomAD_genomes_controls_nhomalt','gnomAD_genomes_POPMAX_AF','gnomAD_genomes_POPMAX_nhomalt',
		'gnomAD_genomes_controls_POPMAX_AF','gnomAD_genomes_controls_POPMAX_nhomalt','gnomAD_exomes_controls_AC',
		'gnomAD_genomes_controls_AC'
		]

        if self.genome_type == "geno38":
            HGMD = pd.read_csv(config.HGMD38, spe="\t", header=None, names=["#CHROM", "START", "END", "hgmd"], encoding="latin")
        elif self.genome_type == "geno37":
            HGMD = pd.read_csv(config.HGMD37, spe="\t", header=None, names=["#CHROM", "START", "END", "hgmd"])

        if len(result) != 0:
            for index,row in result.iterrows():
                try:
                    # if ((row['G%'] != 0.0) & (row['ALT'] == 'G')):
                    # 	result.loc[index,'unbalance'] = 'G='+unicode(row['G%'])
                    # elif ((row['C%'] != 0.0) & (row['ALT'] == 'C')):
                    # 	result.loc[index,'unbalance'] = 'C='+unicode(row['C%'])
                    # elif ((row['T%'] != 0.0) & (row['ALT'] == 'T')):
                    # 	result.loc[index,'unbalance'] = 'T='+unicode(row['T%'])
                    # elif ((row['A%'] != 0.0) & (row['ALT'] == 'A')):
                    # 	result.loc[index,'unbalance'] = 'A='+unicode(row['A%'])
                    # if ((row['ins%'] != 0.0)):
                    # 	result.loc[index,'unbalance'] = 'ins='+unicode(row['ins%'])
                    # if ((row['del%'] != 0.0)):
                    # 	result.loc[index,'unbalance'] = 'del='+unicode(row['del%'])
                    if ((row['G%'] != 0.0) & (row['ALT'] == 'G')):
                        result.loc[index,'unbalance'] = 'G='+str(row['G%'])
                    elif ((row['C%'] != 0.0) & (row['ALT'] == 'C')):
                        result.loc[index,'unbalance'] = 'C='+str(row['C%'])
                    elif ((row['T%'] != 0.0) & (row['ALT'] == 'T')):
                        result.loc[index,'unbalance'] = 'T='+str(row['T%'])
                    elif ((row['A%'] != 0.0) & (row['ALT'] == 'A')):
                        result.loc[index,'unbalance'] = 'A='+str(row['A%'])
                    if ((row['ins%'] != 0.0)):
                        result.loc[index,'unbalance'] = 'ins='+str(row['ins%'])
                    if ((row['del%'] != 0.0)):
                        result.loc[index,'unbalance'] = 'del='+str(row['del%'])


                    result.loc[index,'DEPTH2'] = row['sum']
                except KeyError:
                    print ('TROVA VARIANTE TERMINATO!!!')
    ######################################################################################################
            result['samtools'] = result[samtools].str.split(':').str.get(0)
            result['gatk'] = result[gatk].str.split(':').str.get(0)

            # result['DEPTH'] = np.where(result['INFO'].str.split(';').str.get(0) == 'INDEL',
            # 			result['INFO'].str.split(';').str.get(3).str.split('=').str.get(1),
            # 			result['INFO'].str.split(';').str.get(0).str.split('=').str.get(1))
            result['DEPTH'] = result['INFO'].str.split('DP=').str.get(1).str.split(';').str.get(0)

            result['mapquality'] = result['INFO'].str.split('CSQ=').str.get(0).str.split('MQ=').str.get(1).str.split(';').str.get(0)

            result['consequence'] = result['INFO'].str.split('|').str.get(1)
            result['impact'] = result['INFO'].str.split('|').str.get(2)
            result['symbol'] = result['INFO'].str.split('|').str.get(3)
            #result['transcript_id'] = result['INFO'].str.split('|').str.get(6)
            result['num_exon'] = result['INFO'].str.split('|').str.get(8)
            result['num_intron'] = result['INFO'].str.split('|').str.get(9)
            result['HGVS_c'] = result['INFO'].str.split('|').str.get(10)
            result['HGVS_p'] = result['INFO'].str.split('|').str.get(11)
            result['cDNA_position'] = result['INFO'].str.split('|').str.get(12)
            result['CDS_position'] = result['INFO'].str.split('|').str.get(13)
            result['protein_position'] = result['INFO'].str.split('|').str.get(14)
            result['amino_acids'] = result['INFO'].str.split('|').str.get(15)
            result['codons'] = result['INFO'].str.split('|').str.get(16)

            result['variation'] = result['INFO'].str.split('|').str.get(17).str.split('&').str.get(0)
            try: result['variation2'] = result['INFO'].str.split('|').str.get(17).str.split('&').str.get(1)
            except: result['variation2'] = ''

            result['variation'].fillna('',inplace=True)
            result['variation2'].fillna('',inplace=True)
            try:
                result['ID'] = np.where(result['variation'].str.contains('rs'),result['variation'],
                                            np.where(result['variation2'].str.contains('rs'),result['variation2'],np.nan))
            except:
                print ('ERROR IN NP.WHERE!!!!')
                result['ID'] = np.nan
            result['cosmic'] = result['INFO'].str.split('|').str.get(17).str.split('&').str.get(1)
            result['cosmic'].fillna('',inplace=True)

            try:
                result['variation_cosmic'] = np.where(result['cosmic'].str.contains('COSM'),result['cosmic'],np.nan)
            except:
                print ('ERROR COSMIC IN NP.WHERE!!!!')
                result['variation_cosmic'] = np.nan
            result['distance'] = result['INFO'].str.split('|').str.get(18)
            result['strand'] = result['INFO'].str.split('|').str.get(19)
            result['hgnc_id'] = result['INFO'].str.split('|').str.get(23)
            result['canonical'] = result['INFO'].str.split('|').str.get(24)

            result['gene_pheno'] = result['INFO'].str.split('|').str.get(38)
            result['sift'] = result['INFO'].str.split('|').str.get(39)
            result['polyphen'] = result['INFO'].str.split('|').str.get(40)
            result['domain'] = result['INFO'].str.split('|').str.get(41)
            result['GMAF'] = result['INFO'].str.split('|').str.get(44)
            result['AFR_MAF'] = result['INFO'].str.split('|').str.get(51)
            result['AMR_MAF'] = result['INFO'].str.split('|').str.get(52)
            result['EAS_MAF'] = result['INFO'].str.split('|').str.get(54)
            result['EUR_MAF'] = result['INFO'].str.split('|').str.get(56)
            result['SAS_MAF'] = result['INFO'].str.split('|').str.get(58)
            result['AA_MAF'] = result['INFO'].str.split('|').str.get(48)
            result['EA_MAF'] = result['INFO'].str.split('|').str.get(59)
            #result['ExAC_MAF'] = result['INFO'].str.split('|').str.get(47)
            #result['Adj_MAF'] = result['INFO'].str.split('|').str.get(48)
            result['ExAC_MAF'] = result['INFO'].str.split('|').str.get(50) #GnomadAF
            result['MAX_MAF'] = result['INFO'].str.split('|').str.get(70)  #MAX_MAF
            result['Adj_MAF'] = result['INFO'].str.split('|').str.get(71)  #MAX_AF_POPS
            result['clin_sign'] = result['INFO'].str.split('|').str.get(72)
            result['somatic'] = result['INFO'].str.split('|').str.get(73)
            result['PHENO'] = result['INFO'].str.split('|').str.get(74)
            result['pubmed'] = result['INFO'].str.split('|').str.get(75)
            result['pubmed'] = result['pubmed'].str.split('&').str.get(0)+' '+result['pubmed'].str.split('&').str.get(1)+' '+result['pubmed'].str.split('&').str.get(2)
            try:result['CADD_rankscore'] = result['INFO'].str.split('|').str.get(81).str.split(',').str.get(0)
            except: result['CADD_rankscore'] = result['INFO'].str.split('|').str.get(81)
            try: result['DANN_rankscore'] = result['INFO'].str.split('|').str.get(82).str.split(',').str.get(0)
            except: result['DANN_rankscore'] = result['INFO'].str.split('|').str.get(82)
            try: result['EigenPC_rankscore'] = result['INFO'].str.split('|').str.get(83).str.split(',').str.get(0)
            except: result['EigenPC_rankscore'] = result['INFO'].str.split('|').str.get(83)
            try: result['FATHMM_rankscore'] = result['INFO'].str.split('|').str.get(84).str.split(',').str.get(0)
            except: result['FATHMM_rankscore'] = result['INFO'].str.split('|').str.get(84)
            result['FATHMM_pred'] = result['INFO'].str.split('|').str.get(85)
            try: result['GERP_rankscore'] = result['INFO'].str.split('|').str.get(86).str.split(',').str.get(0)
            except: result['GERP_rankscore'] = result['INFO'].str.split('|').str.get(86)
            result['Interpro_domain'] = result['INFO'].str.split('|').str.get(87)
            try: result['LRT_rankscore'] = result['INFO'].str.split('|').str.get(88).str.split(',').str.get(0)
            except: result['LRT_rankscore'] = result['INFO'].str.split('|').str.get(88)
            result['LRT_pred'] = result['INFO'].str.split('|').str.get(89)
            result['MCAP_pred'] = result['INFO'].str.split('|').str.get(90)
            result['MetaLR_pred'] = result['INFO'].str.split('|').str.get(91)
            try: result['MetaLR_rankscore'] = result['INFO'].str.split('|').str.get(92).str.split(',').str.get(0)
            except:result['MetaLR_rankscore'] = result['INFO'].str.split('|').str.get(92)
            result['MetaSVM_pred'] = result['INFO'].str.split('|').str.get(93)
            try: result['MetaSVM_rankscore'] = result['INFO'].str.split('|').str.get(94).str.split(',').str.get(0)
            except:result['MetaSVM_rankscore'] = result['INFO'].str.split('|').str.get(94)
            try: result['MutPred_rankscore'] = result['INFO'].str.split('|').str.get(95).str.split(',').str.get(0)
            except:result['MutPred_rankscore'] = result['INFO'].str.split('|').str.get(95)
            result['MutationAssessor_pred'] = result['INFO'].str.split('|').str.get(96)
            try: result['MutationAssessor_rankscore'] = result['INFO'].str.split('|').str.get(97).str.split(',').str.get(0)
            except: result['MutationAssessor_rankscore'] = result['INFO'].str.split('|').str.get(97)
            try: result['MutationTaster_rankscore'] = result['INFO'].str.split('|').str.get(98).str.split(',').str.get(0)
            except:result['MutationTaster_rankscore'] = result['INFO'].str.split('|').str.get(98)
            result['MutationTaster_pred'] = result['INFO'].str.split('|').str.get(99)
            try: result['PROVEAN_rankscore'] = result['INFO'].str.split('|').str.get(100).str.split(',').str.get(0)
            except:result['PROVEAN_rankscore'] = result['INFO'].str.split('|').str.get(100)
            result['PROVEAN_pred'] = result['INFO'].str.split('|').str.get(101)
            result['Polyphen2HDIV_pred'] = result['INFO'].str.split('|').str.get(102)
            try: result['Polyphen2HDIV_rankscore'] = result['INFO'].str.split('|').str.get(103).str.split(',').str.get(0)
            except:result['Polyphen2HDIV_rankscore'] = result['INFO'].str.split('|').str.get(103)
            result['Polyphen2HVAR_pred'] = result['INFO'].str.split('|').str.get(104)
            try: result['Polyphen2HVAR_rankscore'] = result['INFO'].str.split('|').str.get(105).str.split(',').str.get(0)
            except:result['Polyphen2HVAR_rankscore'] = result['INFO'].str.split('|').str.get(105)
            try: result['REVEL_rankscore'] = result['INFO'].str.split('|').str.get(106).str.split(',').str.get(0)
            except:result['REVEL_rankscore'] = result['INFO'].str.split('|').str.get(106)
            try: result['SIFT_rankscore'] = result['INFO'].str.split('|').str.get(107).str.split(',').str.get(0)
            except:result['SIFT_rankscore'] = result['INFO'].str.split('|').str.get(107)
            result['SIFT_pred'] = result['INFO'].str.split('|').str.get(108)
            try: result['SiPhy29way_rankscore'] = result['INFO'].str.split('|').str.get(109).str.split(',').str.get(0)
            except:result['SiPhy29way_rankscore'] = result['INFO'].str.split('|').str.get(109)
            try: result['VEST3_rankscore'] = result['INFO'].str.split('|').str.get(110).str.split(',').str.get(0)
            except:result['VEST3_rankscore'] = result['INFO'].str.split('|').str.get(110)
            try: result['clinvar_MedGen_id'] = result['INFO'].str.split('|').str.get(111).str.split(',').str.get(0)
            except:result['clinvar_MedGen_id'] = result['INFO'].str.split('|').str.get(111)
            try: result['clinvar_OMIM_id'] = result['INFO'].str.split('|').str.get(112).str.split(',').str.get(0)
            except:result['clinvar_OMIM_id'] = result['INFO'].str.split('|').str.get(112)
            try: result['clinvar_Orphanet_id'] = result['INFO'].str.split('|').str.get(113).str.split(',').str.get(0)
            except:result['clinvar_Orphanet_id'] = result['INFO'].str.split('|').str.get(113)
            try: result['clinvar_clnsig'] = result['INFO'].str.split('|').str.get(114).str.split(',').str.get(0)
            except:result['clinvar_clnsig'] = result['INFO'].str.split('|').str.get(114)
            try: result['clinvar_review'] = result['INFO'].str.split('|').str.get(115).str.split(',').str.get(0)
            except:result['clinvar_review'] = result['INFO'].str.split('|').str.get(115)
            try: result['fathmmMKL_pred'] = result['INFO'].str.split('|').str.get(116).str.split(',').str.get(0)
            except: result['fathmmMKL_pred'] = result['INFO'].str.split('|').str.get(116)
            try: result['gnomAD_exomes_POPMAX_AF'] = result['INFO'].str.split('|').str.get(117).str.split(',').str.get(0)
            except:result['gnomAD_exomes_POPMAX_AF'] = result['INFO'].str.split('|').str.get(117)
            try: result['gnomAD_exomes_POPMAX_nhomalt'] = result['INFO'].str.split('|').str.get(118).str.split(',').str.get(0)
            except:result['gnomAD_exomes_POPMAX_nhomalt'] = result['INFO'].str.split('|').str.get(118)
            try: result['gnomAD_exomes_controls_AC'] = result['INFO'].str.split('|').str.get(119).str.split(',').str.get(0)
            except:result['gnomAD_exomes_controls_AC'] = result['INFO'].str.split('|').str.get(119)
            try: result['gnomAD_exomes_controls_AF'] = result['INFO'].str.split('|').str.get(120).str.split(',').str.get(0)
            except:result['gnomAD_exomes_controls_AF'] = result['INFO'].str.split('|').str.get(120)
            try: result['gnomAD_exomes_controls_nhomalt'] = result['INFO'].str.split('|').str.get(121).str.split(',').str.get(0)
            except:result['gnomAD_exomes_controls_nhomalt'] = result['INFO'].str.split('|').str.get(121)
            try: result['gnomAD_genomes_POPMAX_AF'] = result['INFO'].str.split('|').str.get(122).str.split(',').str.get(0)
            except:result['gnomAD_genomes_POPMAX_AF'] = result['INFO'].str.split('|').str.get(122)
            try: result['gnomAD_genomes_POPMAX_nhomalt'] = result['INFO'].str.split('|').str.get(123).str.split(',').str.get(0)
            except:result['gnomAD_genomes_POPMAX_nhomalt'] = result['INFO'].str.split('|').str.get(123)
            try: result['gnomAD_genomes_controls_AC'] = -999 #result['INFO'].str.split('\|').str.get(115).str.split(',').str.get(0)
            except:result['gnomAD_genomes_controls_AC'] =  -999 #result['INFO'].str.split('\|').str.get(115)
            try: result['gnomAD_genomes_controls_POPMAX_AF'] =  -999 #result['INFO'].str.split('\|').str.get(116).str.split(',').str.get(0)
            except:result['gnomAD_genomes_controls_POPMAX_AF'] =  -999 #result['INFO'].str.split('\|').str.get(116)
            try: result['gnomAD_genomes_controls_POPMAX_nhomalt'] =  -999 #result['INFO'].str.split('\|').str.get(117).str.split(',').str.get(0)
            except:result['gnomAD_genomes_controls_POPMAX_nhomalt'] =  -999 #result['INFO'].str.split('\|').str.get(117)
            try: result['gnomAD_genomes_controls_nhomalt'] =  -999 #result['INFO'].str.split('\|').str.get(118).str.split(',').str.get(0)
            except:result['gnomAD_genomes_controls_nhomalt'] =  -999#result['INFO'].str.split('\|').str.get(118)
            try: result['phastCons100way_rankscore'] = result['INFO'].str.split('\|').str.get(124).str.split(',').str.get(0)
            except:result['phastCons100way_rankscore'] = result['INFO'].str.split('\|').str.get(124)
            try: result['phastCons20way_rankscore'] = result['INFO'].str.split('\|').str.get(125).str.split(',').str.get(0)
            except:result['phastCons20way_rankscore'] = result['INFO'].str.split('\|').str.get(125)
            try: result['phyloP100way_rankscore'] = result['INFO'].str.split('\|').str.get(126).str.split(',').str.get(0)
            except:result['phyloP100way_rankscore'] = result['INFO'].str.split('\|').str.get(126)
            try: result['phyloP20way_rankscore'] = result['INFO'].str.split('\|').str.get(127).str.split(',').str.get(0)
            except:result['phyloP20way_rankscore'] = result['INFO'].str.split('\|').str.get(127)
            try: result['ada_score'] = result['INFO'].str.split('\|').str.get(128).str.split(',').str.get(0)
            except:result['ada_score'] = result['INFO'].str.split('\|').str.get(128)
            try: result['rf_score'] = result['INFO'].str.split('\|').str.get(129).str.split(',').str.get(0)
            except:result['rf_score'] = result['INFO'].str.split('\|').str.get(129)

            mask1a = result['samtools'] == '0/1'
            mask1a2 = result['samtools'] =='1/2'
            mask1b = result['gatk'] == '0/1'
            mask2a = result['samtools'] == '0/0'
            mask2b = result['gatk'] == '0/0'
            mask3a = result['samtools'] == '1/1'
            mask3b = result['gatk'] == '1/1'
            result.loc[mask1a, 'samtools_geno'] = 'het'
            result.loc[mask1a2, 'samtools_geno'] = 'het'
            result.loc[mask1b, 'gatk_geno'] = 'het'
            result.loc[mask2a, 'samtools_geno'] = 'homo_wild'
            result.loc[mask2b, 'gatk_geno'] = 'homo_wild'
            result.loc[mask3a, 'samtools_geno'] = 'homo'
            result.loc[mask3b, 'gatk_geno'] = 'homo'

            HGMD['POS'] = HGMD['END']
            result2 = pd.merge(result,HGMD,on=['#CHROM','POS'],how='left')
            result2.drop_duplicates(subset=['#CHROM','POS'],inplace=True)
            result2['hgmd'] = result2['hgmd_x']
            try:
                result2['hgmd_mutation'] = result2['hgmd_y'].str.split('|').str.get(2)
                result2['hgmd_function'] = result2['hgmd_y'].str.split('|').str.get(3)
                result2['hgmd_phenotype'] = result2['hgmd_y'].str.split('|').str.get(4)
                result2['hgmd_pubmed'] = result2['hgmd_y'].str.split('|').str.get(5)
            except AttributeError:
                result2['hgmd_mutation'] = ''
                result2['hgmd_function'] = ''
                result2['hgmd_phenotype'] = ''
                result2['hgmd_pubmed'] = ''
                print ('AttributeError!!!')

            result2.drop('hgmd_y',axis=1,inplace=True)
            result2.drop('hgmd_x',axis=1,inplace=True)
            result2['length'] = result2['length'].map('{:,.0f}'.format)

            for index, row in result2.iterrows():
                count_ref = len(row['REF'])
                count_alt = len(row['ALT'])
                result2.loc[index,'count_ref'] = int(count_ref)
                result2.loc[index,'count_alt'] = int(count_alt)
            result2['START'] = np.where((result2['count_ref']>1),
                        result2['POS']+1,result2['POS'])

            result2['END'] = np.where((result2['count_ref'] > result2['count_alt']),
                        result2['POS']+(result2['count_ref']-1),
                        np.where((result2['count_ref'] < result2['count_alt']),
                        result2['START']+((result2['count_alt']-1)-(result2['count_ref']-1)),
                        result2['POS']))

            result2['count_ref'] = result2['count_ref'].astype(int)
            result2['count_alt'] = result2['count_alt'].astype(int)
            for index, row in result2.iterrows():
                if (row['count_ref'] > row['count_alt'])&(row['count_ref']>1):
                    row['REF_2'] = row['REF'][row['count_alt']:]
                    row['ALT_2'] = '-'
                    result2.loc[index,'REF_2'] = row['REF_2']
                    result2.loc[index,'ALT_2'] = row['ALT_2']
                    result2.loc[index,'types'] = 'DELETION'
                elif (row['count_ref'] < row['count_alt'])&(row['count_alt']>1):
                    row['REF_2'] = '-'
                    row['ALT_2'] = row['ALT'][row['count_ref']:]
                    result2.loc[index,'REF_2'] = row['REF_2']
                    result2.loc[index,'ALT_2'] = row['ALT_2']
                    result2.loc[index,'types'] = 'INSERTION'
                elif (row['count_ref'] == row['count_alt'])&(row['count_ref'] == 1):
                    row['REF_2'] = row['REF']
                    row['ALT_2'] = row['ALT']
                    result2.loc[index,'REF_2'] = row['REF_2']
                    result2.loc[index,'ALT_2'] = row['ALT_2']
                    result2.loc[index,'types'] = 'SVN'

            result2['START'] = result2['START'].astype(int)
            result2['END'] = result2['END'].astype(int)
            result2['HGVS'] = result2['#CHROM']+':'+result2['START'].astype(str)+'-'+result2['END'].astype(str)+':'+result2['REF_2']+'/'+result2['ALT_2']

            result2['HGVS_c'].fillna('',inplace=True)
            result2['HGVS_p'].fillna('',inplace=True)


            for index, row in result2.iterrows():
                if row['strand'] == str(-1):
                    try: seq_x = Seq(row['REF_2']) #, IUPAC.unambiguous_dna)
                    except: seq_x = Seq(row['REF']) #, IUPAC.unambiguous_dna)
                    try: seq_y = Seq(row['ALT_2']) #, IUPAC.unambiguous_dna)
                    except: seq_x = Seq(row['ALT']) #, IUPAC.unambiguous_dna)
                    seq_ref = seq_x.reverse_complement()
                    seq_alt = seq_y.reverse_complement()
                    row['REF_2'] = seq_ref
                    row['ALT_2'] = seq_alt
                #try:
                row['HGVS_c'] = row['HGVS_c'].replace('N>',str(row['REF_2']+'>'))
                DEL = 'del'+'N'*(len(row['REF_2']))
                row['HGVS_c'] = row['HGVS_c'].replace(DEL,'del'+str(row['REF_2']))
                #except:  #row['HGVS_c'] = ''
                result2.loc[index,'HGVS_c'] = row['HGVS_c']

            mask1 = result2['HGVS_p'].str.contains('p.%3D')
            result2.loc[mask1,'HGVS_p'] = '(p.%3D)'

            result2['HGVS_c'] = result2['HGVS_c'].replace(r'NN*','',regex=True)
            result2['HGVS_c'] = result2['HGVS_c'].replace(r'delN+','',regex=True)
            result2['HGVS_c'] = result2['HGVS_c'].replace(r'del-','',regex=True)
            result2['HGVS_c'] = result2['HGVS_c'].replace(r'^M_','NM_',regex=True)
            result2['HGVS_c'] = result2['HGVS_c'].replace(r'^R_','NR_',regex=True)

            result2['refseq'].fillna(result2['HGVS_c'].str.split('.').str.get(0),inplace=True)
            result2['variation_cosmic'].fillna('RISN',inplace=True)

            mask1 = result2['variation_cosmic'].str.contains('RISN')
            result2.loc[mask1,'variation_cosmic'] = ''
            mask2 = result2['variation_cosmic'].str.contains('rs')
            result2.loc[mask2,'variation_cosmic'] = ''
            mask3 = result2['variation_cosmic'].str.contains('RPGR')
            result2.loc[mask3,'variation_cosmic'] = ''

            result2['hgmd_mutation'].replace('null;null','',inplace=True)
            result2['CHROM'] = result2['#CHROM']
            result2['ALTGENE'] = result2['GENE']
            result2['GENE'] = result2['symbol']

            #result2['GENE'] = result2['GENE'].fillna(result2['symbol'])
            #result2['GENE'] = result2['symbol']

            result2['sample'] = str(name)
            #result2b = result2[result2['QUAL'] >= 18]
            result2b = result2
            result2b['hgmd_mutation'].fillna('',inplace=True)
            result2b.sort_values(by=['#CHROM','POS'],ascending=[True,True],inplace=True)

            result2b['DEPTH'] = result2b['DEPTH'].str.split(',').str.get(0)
            result2b['DEPTH2'].fillna(result2b['DEPTH'],inplace=True)
            result2b['DEPTH2'].fillna(-999,inplace=True)
            result2b['DEPTH'] = result2b['DEPTH2'].astype(int)

            result2b['mapquality'].fillna(-999,inplace=True)
            result2b['mapquality2'] = result2b['mapquality'].astype(float)
            result2b['mapquality'] = result2b['mapquality2'].astype(int)
            result2b['GENE'].replace('DFNB31','WHRN',inplace=True)
            result3 = result2b[cols]
            print ('Annotation CDS len:', len(result3),'->',str(name))

        else:
            result3 = pd.DataFrame(columns=cols)
            print ('Annotation CDS len:', len(result3),'->',str(name))

        return result3

    def run(self):
        refseq_filter = self.filter_vep_refseq()
        adapt_annot = self.adapt_annotation(refseq_filter)

        # write it to final/ dir
        final_dir = dir_tree.principal_directory.final.path
        out_filepath = os.join(final_dir, str(self.sample.name) + "_CDS_annot.csv")
        adapt_annot.to_csv(adapt_annot, out_filepath)



    



