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

        kwargs.update({"principal_directory": self.principal_directory, "sample": self.sample, "genome_type": self.genome_type, "panel": self.panel})
        return kwargs

    def final_annotation(self):
        APPRIS = pd.read_csv(config.APPRIS, sep="\t")
        ECCEZIONI = pd.read_csv(config.ECCEZIONI, sep="\t")

        APPRIS['nm_refseq'] = APPRIS['refseq'].str.split(".").str.get(0)
        ECCEZIONI['nm_refseq'] = ECCEZIONI['refseq'].str.split(".").str.get(0)

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

    def run(self):
        final_annot = self.final_annotation()
        final_annot.to_csv("annotationpipe_result.csv", sep="\t")

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


    



