from Pipes.Pipe import Pipe
import os
import pandas as pd
import config
import glob
import dir_tree
import numpy as np
from Bio.Seq import Seq


class AnnotationPipe(Pipe):

    eredita = None
    hgmd = None
    appris = None
    eccezioni = None

    def __init__(self):
        pass

    def process(self, **kwargs):
        self.principal_directory = kwargs.pop("principal_directory", None)
        self.samples = kwargs.pop("samples", None)
        self.genome_type = kwargs.pop("genome", None)
        self.panel = kwargs.pop("panel", None)

        self.input_vcf = os.path.join(self.principal_directory, "vcf/")
        self.output_annotation = os.path.join(self.principal_directory, "annotation/")
        self.output_final = os.path.join(self.principal_directory, "final/")
        self.folder_temp = os.path.join(self.principal_directory, "temp/")
        self.folder_COV = os.path.join(self.principal_directory, "coverage/")
        self.input_phenotype = os.path.join(self.principal_directory, "pheno/phenotype")

        self.run()

        kwargs.update({"principal_directory": self.principal_directory, "samples": self.samples, "genome_type": self.genome_type, "panel": self.panel})
        return kwargs

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


    



