from Pipes.Pipe import Pipe
import os
import pandas as pd
import config
import glob
import dir_tree
import numpy as np
from Bio.Seq import Seq

""" Pipe that will modify the results of VEP. """
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

        self.run()

        self.sample.saveJSON()

        kwargs.update({"principal_directory": self.principal_directory, "samples": self.samples, "genome_type": self.genome_type, "panel": self.panel})
        return kwargs

    def final_annotation(self):
        folder_COV = dir_tree.principal_directory.coverage.path
        sample_name = str(self.sample.name)
        coverage = pd.read_csv(os.path.join(folder_COV, sample_name, sample_name + "_all"), sep='\t', header=0)
        eredita = pd.read_csv(config.EREDITA37, sep="\t", header=0)
        coverage = coverage[['#CHROM','POS','C%','G%','T%','A%','ins%','del%','sum']]
        CDS = pd.read_csv(self.sample.vcf_annot_CDS, sep="\t", header=None, comment='#', names=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 
                                                                                                sample_name + '_samt', sample_name + '_gatk', 'GENE', 'exone', 'length', 'strand', 'refseq', 'hgmd'])

        cds_coverage = pd.merge(CDS, coverage, on=["#CHROM", "POS"], how='left')
        cds_coverage_eredita = pd.merge(cds_coverage, eredita, on=['GENE'], how='left')

        result = cds_coverage_eredita
        info_split = result['INFO'].str.split('CSQ=') #  LHS of CSQ contains VariantCalling INFO; RHS of CSQ contains VEP data;
        result["INFO1"] = info_split.str.get(0)
        result["INFO2"] = info_split.str.get(1).str.split(",")


        #TODO:
        result3['PRINCIPAL'] = None

        REFSEQPRINCIPAL = APPRIS['refseq'].str.split('.').str.get(0)+'.'
        REFSEQPRINCIPALECCEZIONI = ECCEZIONI['refseq'].str.split('.').str.get(0)+'.'
        ECCEZIONI.dropna(subset=['refseq'],inplace=True)
        REFSEQPRINCIPALECCEZIONI.dropna(inplace=True)
        for index,row in result3.iterrows():
            for x in row.INFO2:
                for refseq in REFSEQPRINCIPALECCEZIONI:
                    if refseq in x:
                        #print (str(refseq))
                        if (ECCEZIONI['GENE'][ECCEZIONI['refseq'].str.split('.').str.get(0)+'.' == refseq]==row['GENE']).item():
                            result3.loc[index,'PRINCIPAL'] = str(row['INFOFIRST'])+'CSQ='+x
                        else: pass
                for refseq in REFSEQPRINCIPAL:
                    if refseq in x:
                        if (APPRIS['GENE'][APPRIS['refseq'].str.split('.').str.get(0)+'.' == refseq]==row['GENE']).item():
                            result3.loc[index,'PRINCIPAL'] = str(row['INFOFIRST'])+'CSQ='+x
                        else: pass

        result3['PRINCIPAL2'] = np.where(result3['PRINCIPAL'].isnull(),
                        result3['INFO3'],np.nan)

        result3['INFO'] = result3['PRINCIPAL']
        result3['INFO'].fillna(result3['PRINCIPAL2'],inplace=True)
        result3.drop('PRINCIPAL',axis=1,inplace=True)
        result3.drop('INFOFIRST',axis=1,inplace=True)
        result3.drop('INFO2',axis=1,inplace=True)
        result3.drop('INFO3',axis=1,inplace=True)
        result3.drop('PRINCIPAL2',axis=1,inplace=True)
        


    def run(self):
        pass
        #self.call_vep()

        #self.annotate_vcfs()

        # self.phenotype = pd.read_csv(
        #     self.input_phenotype, sep="\t", header=0, encoding="utf-8"
        # )

        # if self.genome_type == "geno37":
        #     # phenotype = pd.read_csv(self.input_phenotype, sep='\t', header=0)
        #     self.eredita = pd.read_csv(config.EREDITA37, sep="\t", header=0)
        #     self.hgmd = pd.read_csv(
        #         config.HGMD37,
        #         sep="\t",
        #         header=None,
        #         names=["#CHROM", "START", "END", "hgmd"],
        #     )

        # elif self.genome_type == "geno38":
        #     self.eredita = pd.read_csv(config.EREDITA38, sep="\t", header=0)
        #     self.hgmd = pd.read_csv(
        #         config.HGMD38,
        #         sep="\t",
        #         header=None,
        #         names=["#CHROM", "START", "END", "hgmd"],
        #         encoding="latin",
        #     )
        #     self.appris = pd.read_csv(config.APPRIS, sep="\t", header=0)
        #     self.eccezioni = pd.read_csv(config.ECCEZIONI, sep="\t", header=0)


    



