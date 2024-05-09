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

import CDSHandler
import cutCDS_jurgen

genome_dict = {
    "geno37": join("/home/magi/", "dataset/GENOME/37/all_chr37.fa"),
    "geno38": join("/home/magi/", "dataset/GENOME/38/all_chr38.fa"),
}

##########################################################################################
dbsnp144_37 = join("/home/magi/", "dataset/dbsnp144/37/common_all_20150605_2.vcf")
dbsnp144_38 = join("/home/magi/", "dataset/dbsnp144/38/common_all_20150603_2.vcf")
indel_37 = join("/home/magi/", "dataset/dbsnp144/37/common_all_20150605_indel.vcf")
indel_38 = join("/home/magi/", "dataset/dbsnp144/38/common_all_20150603_indel.vcf")
clinvar = join("/home/magi/", "dataset/dbsnp144/38/clinvar_20140929_2.vcf")
clinvar_indel = join("/home/magi/", "dataset/dbsnp144/38/clinvar_20140929_indel.vcf")
#############################################################################################
BUCHIARTIFICIALI = join("/home/magi/", "PROJECT/diagnosys/bin/BUCHIARTIFICIALI.txt")
##############################################################################################


class CoveragePipe2(ParallelPipe):
    def __init__(self) -> None:
        super().__init__()

    def process(self, **kwargs):
        self.thread_id = kwargs.pop("thread_id", None)
        principal_directory = kwargs.pop("principal_directory", None)
        panel = kwargs.pop("panel", None)
        folder_pheno = join(principal_directory, "pheno")
        # phenotype = kwargs.pop("phenotype")
        genome_type = kwargs.pop("genome", None)
        dest = kwargs.pop("dest", None)
        path = kwargs.pop("path")
        input_phenotype = join(principal_directory, "pheno/phenotype")
        phenotype = pd.read_csv(input_phenotype, sep="\t", header=0, dtype=str)

        dbsnp = dbsnp144_38 if genome_type == "geno38" else dbsnp144_37
        genome = genome_dict[genome_type]
        queue = kwargs.pop("queue", None)

        sample = kwargs.pop("sample", None)
        sample_name = str(sample["name"])
        genelist = list(phenotype["gene"][phenotype["sample"] == sample_name])
        sospetto = phenotype["malattia"][phenotype["sample"] == sample_name].unique()[0]

        vertical = pd.DataFrame()
        verticalX = pd.DataFrame()
        if genome_type == "geno38":
            vertical, verticalX, BED = self.cutCDS(genelist, sample, folder_pheno, dest)
        elif genome_type == "geno37":
            vertical = self.cutCDS37(genelist, sample, folder_pheno, dest)

        if sospetto == "CONNETTIVOPATIE":
            print("connettivopatie")
            bed_ = join(folder_pheno, "bed_" + sample_name)
            bed = pd.read_csv(bed_, sep="\t")
            region = pd.DataFrame(
                {
                    "#CHROM": pd.Series(
                        [
                            "chr9",
                            "chr9",
                            "chr9",
                            "chr9",
                            "chr9",
                            "chr9",
                            "chr9",
                            "chr9",
                            "chr9",
                            "chr9",
                            "chr9",
                            "chr9",
                            "chr9",
                            "chr9",
                            "chr9",
                        ]
                    ),
                    "START": pd.Series(
                        [
                            86584322,
                            86585077,
                            86585812,
                            86585652,
                            86586188,
                            86586587,
                            86586797,
                            86587759,
                            86588201,
                            86588817,
                            86589432,
                            86590377,
                            86591910,
                            86592604,
                            86593110,
                        ]
                    ),
                    "END": pd.Series(
                        [
                            86584355,
                            86585246,
                            86585827,
                            86585734,
                            86586271,
                            86586641,
                            86587104,
                            86587887,
                            86588314,
                            86588888,
                            86589504,
                            86590420,
                            86591966,
                            86592701,
                            86593167,
                        ]
                    ),
                    "GENE": pd.Series(
                        [
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                        ]
                    ),
                    "length": pd.Series(
                        [34, 170, 16, 83, 84, 55, 308, 129, 114, 72, 73, 44, 57, 98, 58]
                    ),
                    "exone": pd.Series([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
                    "strand": pd.Series([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]),
                    "refseq": pd.Series(
                        [
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                        ]
                    ),
                    "hgmd": pd.Series(
                        [
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                        ]
                    ),
                }
            )
            bed = pd.concat([bed, region])
            bed.to_csv(bed_, sep="\t", index=False)
            bedx_ = join(folder_pheno, "bedX_" + sample_name)
            bedx = pd.read_csv(bedx_, sep="\t")
            vertical, verticalX = cutCDS_jurgen.create_vertical(
                bed, bedx, sample_name, folder_pheno
            )
            print(vertical)
            vertical_ = join(folder_pheno, "vertical_" + sample_name)
            vertical = pd.read_csv(vertical_, sep="\t")
            verticalx_ = join(folder_pheno, "verticalX_" + sample_name)
            verticalX = pd.read_csv(verticalx_, sep="\t")
        elif sospetto == "CROMATINOPATIE":
            print("CROMATINOPATIE, in ")
            bed_ = join(folder_pheno, "bed_" + sample_name)
            bed = pd.read_csv(bed_, sep="\t")
            region = pd.DataFrame(
                {
                    "#CHROM": pd.Series(
                        [
                            "chr11",
                            "chr12",
                            "chr15",
                            "chr15",
                            "chr15",
                            "chr15",
                            "chr16",
                            "chr16",
                            "chr22",
                            "chr3",
                            "chr3",
                            "chr3",
                            "chr5",
                            "chr5",
                            "chr7",
                            "chr7",
                            "chr7",
                            "chr7",
                            "chr7",
                            "chr7",
                            "chr3",
                            "chr22",
                            "chr22",
                            "chr22",
                            "chr22",
                            "chr22",
                            "chr22",
                            "chrX",
                        ]
                    ),
                    "START": pd.Series(
                        [
                            119077253,
                            112910723,
                            66995991,
                            66995821,
                            67001018,
                            96875310,
                            30134513,
                            55515763,
                            19770411,
                            20202572,
                            20202683,
                            8591535,
                            78280735,
                            78280904,
                            140624418,
                            150693467,
                            150700214,
                            150700331,
                            150700720,
                            150700214,
                            20216072,
                            19748403,
                            19753887,
                            19754000,
                            19754265,
                            19754305,
                            19444417,
                            136649004,
                        ]
                    ),
                    "END": pd.Series(
                        [
                            119077347,
                            112910816,
                            66996065,
                            66995870,
                            67001101,
                            96875443,
                            30134555,
                            55515815,
                            19770570,
                            20202743,
                            20202743,
                            8591605,
                            78280783,
                            78280957,
                            140624528,
                            150693561,
                            150700516,
                            150700413,
                            150700808,
                            150700494,
                            20216167,
                            19748627,
                            19754144,
                            19754071,
                            19754311,
                            19754311,
                            19444466,
                            136649122,
                        ]
                    ),
                    "GENE": pd.Series(
                        [
                            "CBL",
                            "PTPN11",
                            "int",
                            "SMAD6",
                            "int",
                            "NR2F2",
                            "MAPK3",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "ARSB",
                            "BRAF",
                            "NOS3",
                            "int",
                            "int",
                            "int",
                            "int",
                            "SGOL1",
                            "TBX1",
                            "TBX1",
                            "int",
                            "int",
                            "int",
                            "UFD1L",
                            "ZIC3",
                        ]
                    ),
                    "length": pd.Series(
                        [
                            95,
                            94,
                            75,
                            50,
                            84,
                            134,
                            43,
                            53,
                            160,
                            172,
                            61,
                            71,
                            49,
                            54,
                            111,
                            95,
                            303,
                            83,
                            89,
                            281,
                            96,
                            225,
                            258,
                            72,
                            47,
                            7,
                            50,
                            119,
                        ]
                    ),
                    "exone": pd.Series(
                        [
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                        ]
                    ),
                    "strand": pd.Series(
                        [
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                            1,
                        ]
                    ),
                    "refseq": pd.Series(
                        [
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                            "N",
                        ]
                    ),
                    "hgmd": pd.Series(
                        [
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                            "int",
                        ]
                    ),
                }
            )

            c = pd.concat([bed, region])
            bed.to_csv(bed_, sep="\t", index=False)
            bedx_ = join(folder_pheno, "bedX_" + sample_name)
            bedx = pd.read_csv(bedx_, sep="\t")
            vertical, verticalX = cutCDS_jurgen.create_vertical(
                bed, bedx, sample_name, folder_pheno
            )
            vertical_ = join(folder_pheno, "vertical_" + sample_name)
            vertical = pd.read_csv(vertical_, sep="\t")
            verticalx_ = join(folder_pheno, "verticalX_" + sample_name)
            verticalX = pd.read_csv(verticalx_, sep="\t")

        if not vertical.empty:
            buchiartificiali = pd.read_csv(BUCHIARTIFICIALI, sep="\t", header=0)
            vertical_macro_path = utils.get_vertical_macro(panel)
            if vertical_macro_path is None:
                vertical_macro = None
                print("No vertical macro for this panel!")
            else:
                vertical_macro = pd.read_csv(vertical_macro_path, sep="\t")

            count_unbalance(
                genome,
                vertical,
                verticalX,
                principal_directory,
                sample,
                phenotype,
                vertical_macro,
            )

            cols = [
                "#CHROM",
                "POS",
                "C%",
                "G%",
                "T%",
                "A%",
                "ins%",
                "del%",
                "sum",
                "DEPTH",
                "GENE",
                "exone",
                "length",
                "strand",
                "refseq",
                "hgmd",
                "sample",
            ]
            folder_to_count = join(principal_directory, "temp/", "to_count/")

            folder_to_count_file = join(
                principal_directory, "temp/", "to_count/", sample_name + "*_disease"
            )
            glob_ = glob.glob(folder_to_count_file)
            glob_ = str(glob_)[1:-1].replace(",", " ")
            glob_ = str(glob_)[1:-1].replace("'", "")
            os.system(
                " ".join(
                    [
                        "cat",
                        glob_,
                        ">",
                        join(
                            principal_directory, "temp/", sample_name + "_final_disease"
                        ),
                    ]
                )
            )
            name = join(principal_directory, "temp/", sample_name + "_final_disease")
            count = pd.read_csv(name, sep="\t", header=None, names=cols)
            count.fillna(0, inplace=True)

            folder_to_count_sex = join(
                principal_directory, "temp/", "to_count/", sample_name + "*_SEX"
            )
            glob_sex = glob.glob(folder_to_count_sex)
            glob_sex = str(glob_sex)[1:-1].replace(",", " ")
            glob_sex = str(glob_sex)[1:-1].replace("'", "")
            os.system(
                " ".join(
                    [
                        "cat",
                        glob_sex,
                        ">",
                        join(principal_directory, "temp/", sample_name + "_final_sex"),
                    ]
                )
            )
            namesex = join(principal_directory, "temp/", sample_name + "_final_sex")
            count_sex = pd.read_csv(namesex, sep="\t", header=None, names=cols)
            count_sex.fillna(0, inplace=True)

            filter_disease_files(
                BED,
                vertical,
                principal_directory,
                phenotype,
                count,
                sample_name,
                buchiartificiali,
            )

            file_to_count = join(
                principal_directory, "temp/", sample_name + "_to_count"
            )
            a = pd.read_csv(
                file_to_count,
                sep="\t",
                header=None,
                quoting=csv.QUOTE_NONE,
                encoding="utf-8",
                low_memory=False,
                on_bad_lines="skip",
                names=["CHROM", "POS", "info", "DEPTH", "CALL", "quality"],
                chunksize=40 * 100024,
            )
            try:
                filter_macroarea_files(
                    vertical_macro, principal_directory, sample_name, buchiartificiali
                )
            except:
                filter_macroarea_files(
                    vertical, principal_directory, sample_name, buchiartificiali
                )
            print("outplot")
            print(
                join(
                    principal_directory,
                    "coverage/",
                    sample_name,
                    sample_name + "_final_sex",
                )
            )
            count_sex.to_csv(
                join(
                    principal_directory,
                    "coverage/",
                    sample_name,
                    sample_name + "_final_sex",
                ),
                sep="\t",
                index=False,
            )
            print("to_Csv")

        # os.system(' '.join(['rm -r',join(folder_to_count)])) # TODO: fix this

        sample["vertical"] = vertical
        sample["verticalx"] = verticalX
        sample["bed"] = BED
        if queue is not None:
            queue.put(sample)

        self.thread_print("Coverage Analysis finished, check bed and vertical files.")
        kwargs.update(
            {
                "thread_id": self.thread_id,
                "principal_directory": principal_directory,
                "panel": panel,
                "genome": genome_type,
                "sample": sample,
                "dest": dest,
                "path": path,
                "queue": queue,
            }
        )

        return kwargs

    def cutCDS(self, genelist, sample, folder_pheno, server_id):
        sample_name = sample["name"]
        vertical, verticalX, BED = cutCDS_jurgen.cutCDS(
            genelist, sample_name, folder_pheno, server_id
        )
        # Additional processing ... TODO

        return vertical, verticalX, BED

    def cutCDS37(self, genelist, sample, folder_pheno, server_id):
        pass


def count_unbalance(
    genome, vertical, verticalX, folder, sample, phenotype, vertical_macro
):

    sample_name = str(sample["name"])
    bam = sample["bam"]

    mpileup_out = join(folder, "temp/", sample_name + "_to_count")
    file_vcf = join(folder, "temp/", sample_name + "_samt.vcf")
    folder_to_count = join(folder, "temp/", "to_count/")
    folder_to_macroarea = join(folder, "temp/", "to_macroarea/")


    # OPTIONAL to replace by parabricks
    os.system(
        " ".join(
            [
                "samtools",
                "mpileup",
                "-Q 0 -q 0 -d10000000 -L 100000 -A",
                bam,
                ">",
                mpileup_out,
            ]
        )
    )

    print("load count!!!")
    a = pd.read_csv(
        mpileup_out,
        sep="\t",
        header=None,
        quoting=csv.QUOTE_NONE,
        encoding="utf-8",
        low_memory=False,
        on_bad_lines="skip",
        names=["CHROM", "POS", "info", "DEPTH", "CALL", "quality"],
        chunksize=40 * 100024,
    )

    i = 0
    for chunk in a:
        i += 1
        chunk["sample"] = sample_name
        count = count_coverage(folder, sample_name, chunk)
        name = join(folder_to_count, sample_name + "_" + str(i))
        count_ = ""
        count_sex = ""
        count_ = filter_for_panel(vertical, folder, name, sample_name)
        count_sex = filter_for_SEX(verticalX, folder, name, sample_name)
        if vertical_macro is None:
            filter_for_macroarea(verticalX, folder, sample_name, i)
        else:
            filter_for_macroarea(vertical_macro, folder, sample_name, i)
        print("filter_disease!!!")
        filter_disease(folder, name, phenotype, count_, sample_name)
        print(i)

    print("load vcf!!!")
    ########################################################################
    print("------------------------")
    return str(sample_name)


def count_coverage(folder, sample, TO_COUNT):
    TO_COUNT = TO_COUNT[["CHROM", "POS", "DEPTH", "CALL"]]

    TO_COUNT["."] = TO_COUNT["CALL"].str.count(r"[.]")
    TO_COUNT[","] = TO_COUNT["CALL"].str.count(r"[,]")
    TO_COUNT["*"] = TO_COUNT["CALL"].str.count(r"[*]")
    TO_COUNT["$"] = TO_COUNT["CALL"].str.count(r"[$]")
    TO_COUNT["C"] = TO_COUNT["CALL"].str.count(r"[cC]")
    TO_COUNT["G"] = TO_COUNT["CALL"].str.count(r"[gG]")
    TO_COUNT["T"] = TO_COUNT["CALL"].str.count(r"[tT]")
    TO_COUNT["A"] = TO_COUNT["CALL"].str.count(r"[aA]")
    TO_COUNT["ins"] = TO_COUNT["CALL"].str.count(r"\+[0-9]+[ACGTNacgtn]+")
    TO_COUNT["del"] = TO_COUNT["CALL"].str.count(r"\-[0-9]+[ACGTNacgtn]+")
    TO_COUNT["other"] = TO_COUNT["CALL"].str.count(r"[*><$^]")

    TO_COUNT["sum"] = (
        TO_COUNT["."]
        + TO_COUNT[","]
        + TO_COUNT["."]
        + TO_COUNT["C"]
        + TO_COUNT["T"]
        + TO_COUNT["A"]
        + TO_COUNT["G"]
        + TO_COUNT["ins"]
        + TO_COUNT["del"]
    )

    TO_COUNT["#CHROM"] = TO_COUNT["CHROM"]
    TO_COUNT = TO_COUNT[
        ["#CHROM", "POS", "DEPTH", "C", "G", "T", "A", "ins", "del", "sum"]
    ]

    TO_COUNT["DEPTH"].fillna(0, inplace=True)
    TO_COUNT["sum"].fillna(0, inplace=True)

    TO_COUNT["sum"] = TO_COUNT["sum"].astype(int)
    TO_COUNT["POS"] = TO_COUNT["POS"].astype(int)
    TO_COUNT["DEPTH"] = TO_COUNT["DEPTH"].astype(int)

    TO_COUNT["selection"] = 1
    TO_COUNT_final = TO_COUNT[
        ["#CHROM", "POS", "DEPTH", "C", "G", "T", "A", "ins", "del", "sum", "selection"]
    ]
    print("COUNT OK!!!")
    print("all count:", len(TO_COUNT_final))

    name_filter_intermedio = sample + "_unbalance_inter.txt"
    file_filtered_intermedio = join(folder, "temp/", name_filter_intermedio)

    TO_COUNT_final.to_csv(file_filtered_intermedio, sep="\t", index=False)
    return


def filter_for_panel(vertical, folder, name, sample):
    name_filter_intermedio = sample + "_unbalance_inter.txt"
    file_filtered_intermedio = join(folder, "temp/", name_filter_intermedio)
    COUNT = pd.read_csv(file_filtered_intermedio, sep="\t", header=0)
    vertical["filt"] = 1
    name_filter = sample + "_unbalance_filter.txt"
    file_filtered = join(folder, "temp/", name_filter)
    print("-------FILTERING-------")
    FILTER = pd.merge(COUNT, vertical, on=["#CHROM", "POS"], how="left")
    FILTER["DEPTH"].fillna(0, inplace=True)
    FILTER["sum"].fillna(0, inplace=True)
    FILTER["selection"].fillna(0, inplace=True)
    FILTER["filt"].fillna(0, inplace=True)
    FILTER = FILTER[FILTER["filt"] == 1]
    FILTER["C%"] = FILTER["C"].astype(float) / FILTER["sum"].astype(float)
    FILTER["G%"] = FILTER["G"].astype(float) / FILTER["sum"].astype(float)
    FILTER["T%"] = FILTER["T"].astype(float) / FILTER["sum"].astype(float)
    FILTER["A%"] = FILTER["A"].astype(float) / FILTER["sum"].astype(float)
    FILTER["ins%"] = FILTER["ins"].astype(float) / FILTER["sum"].astype(float)
    FILTER["del%"] = FILTER["del"].astype(float) / FILTER["sum"].astype(float)
    FILTER["C%"] = FILTER["C%"].map("{:,.3f}".format)
    FILTER["G%"] = FILTER["G%"].map("{:,.3f}".format)
    FILTER["T%"] = FILTER["T%"].map("{:,.3f}".format)
    FILTER["A%"] = FILTER["A%"].map("{:,.3f}".format)
    FILTER["ins%"] = FILTER["ins%"].map("{:,.3f}".format)
    FILTER["del%"] = FILTER["del%"].map("{:,.3f}".format)
    FILTER.drop_duplicates(subset=["#CHROM", "POS", "DEPTH"], keep="last")
    FILTER = FILTER[
        [
            "#CHROM",
            "POS",
            "C%",
            "G%",
            "T%",
            "A%",
            "ins%",
            "del%",
            "sum",
            "DEPTH",
            "GENE",
            "exone",
            "length",
            "strand",
            "refseq",
            "hgmd",
            "selection",
            "filt",
        ]
    ].sort_values(by=["#CHROM", "POS"], ascending=[True, True])
    FILTER["sum"] = FILTER["sum"].astype(int)
    FILTER["POS"] = FILTER["POS"].astype(int)
    FILTER["DEPTH"] = FILTER["DEPTH"].astype(int)
    FILTER["selection"] = FILTER["selection"].astype(int)
    FILTER["filt"] = FILTER["filt"].astype(int)
    FILTER["exone"] = FILTER["exone"].astype(int)
    FILTER["length"] = FILTER["length"].astype(int)
    FILTER["strand"] = FILTER["strand"].astype(int)
    FILTER2 = FILTER[FILTER["filt"] == 1]
    print("len FILTER:", len(FILTER2))
    print("FILTER OK!!!")
    return FILTER2


def filter_for_SEX(verticalx, folder, name, sample):
    name_filter_intermedio = sample + "_unbalance_inter.txt"
    file_filtered_intermedio = join(folder, "temp/", name_filter_intermedio)
    COUNT = pd.read_csv(file_filtered_intermedio, sep="\t", header=0)
    vertical = verticalx[verticalx["GENE"].isin(["AMELX", "AMELY", "SRY"])]
    vertical["filt"] = 1
    name_filter = sample + "_unbalance_filter.txt"
    file_filtered = join(folder, "temp/", name_filter)
    print("-------FILTERING-------")
    FILTER = pd.merge(COUNT, vertical, on=["#CHROM", "POS"], how="left")
    FILTER["DEPTH"].fillna(0, inplace=True)
    FILTER["sum"].fillna(0, inplace=True)
    FILTER["selection"].fillna(0, inplace=True)
    FILTER["filt"].fillna(0, inplace=True)
    FILTER = FILTER[FILTER["DEPTH"] >= 5]
    FILTER = FILTER[FILTER["filt"] == 1]
    FILTER["C%"] = FILTER["C"].astype(float) / FILTER["sum"].astype(float)
    FILTER["G%"] = FILTER["G"].astype(float) / FILTER["sum"].astype(float)
    FILTER["T%"] = FILTER["T"].astype(float) / FILTER["sum"].astype(float)
    FILTER["A%"] = FILTER["A"].astype(float) / FILTER["sum"].astype(float)
    FILTER["ins%"] = FILTER["ins"].astype(float) / FILTER["sum"].astype(float)
    FILTER["del%"] = FILTER["del"].astype(float) / FILTER["sum"].astype(float)
    FILTER["C%"] = FILTER["C%"].map("{:,.1f}".format)
    FILTER["G%"] = FILTER["G%"].map("{:,.1f}".format)
    FILTER["T%"] = FILTER["T%"].map("{:,.1f}".format)
    FILTER["A%"] = FILTER["A%"].map("{:,.1f}".format)
    FILTER["ins%"] = FILTER["ins%"].map("{:,.1f}".format)
    FILTER["del%"] = FILTER["del%"].map("{:,.1f}".format)
    FILTER.drop_duplicates(subset=["#CHROM", "POS", "DEPTH"], keep="last")
    FILTER.drop("selection", axis=1, inplace=True)
    FILTER.drop("filt", axis=1, inplace=True)
    FILTER["sum"] = FILTER["sum"].astype(int)
    FILTER["POS"] = FILTER["POS"].astype(int)
    FILTER["DEPTH"] = FILTER["DEPTH"].astype(int)
    FILTER["exone"] = FILTER["exone"].astype(int)
    FILTER["length"] = FILTER["length"].astype(int)
    FILTER["strand"] = FILTER["strand"].astype(int)
    FILTER_SEX = FILTER[
        [
            "#CHROM",
            "POS",
            "C%",
            "G%",
            "T%",
            "A%",
            "ins%",
            "del%",
            "sum",
            "DEPTH",
            "GENE",
            "exone",
            "length",
            "strand",
            "refseq",
            "hgmd",
        ]
    ].sort_values(by=["#CHROM", "POS"], ascending=[True, True])

    FILTER_SEX["sample"] = str(sample)
    print("len FILTER SEX:", len(FILTER_SEX))
    print("FILTER SEX OK!!!")
    FILTER_SEX.to_csv(name + "_SEX", sep="\t", header=False, index=False)
    return FILTER_SEX


def filter_for_macroarea(vertical, folder, sample, i):
    print("filter_for_macroarea")
    name_filter_intermedio = sample + "_unbalance_inter.txt"
    file_filtered_intermedio = join(folder, "temp/", name_filter_intermedio)
    COUNT = pd.read_csv(file_filtered_intermedio, sep="\t", header=0)

    vertical["filt"] = 1
    print("-------FILTERING-------")
    FILTER = pd.merge(COUNT, vertical, on=["#CHROM", "POS"], how="left")

    FILTER["DEPTH"].fillna(0, inplace=True)
    FILTER["sum"].fillna(0, inplace=True)
    FILTER["selection"].fillna(0, inplace=True)
    FILTER["filt"].fillna(0, inplace=True)

    FILTER = FILTER[FILTER["filt"] == 1]

    FILTER["C%"] = FILTER["C"].astype(float) / FILTER["sum"].astype(float)
    FILTER["G%"] = FILTER["G"].astype(float) / FILTER["sum"].astype(float)
    FILTER["T%"] = FILTER["T"].astype(float) / FILTER["sum"].astype(float)
    FILTER["A%"] = FILTER["A"].astype(float) / FILTER["sum"].astype(float)
    FILTER["ins%"] = FILTER["ins"].astype(float) / FILTER["sum"].astype(float)
    FILTER["del%"] = FILTER["del"].astype(float) / FILTER["sum"].astype(float)

    FILTER["C%"] = FILTER["C%"].map("{:,.3f}".format)
    FILTER["G%"] = FILTER["G%"].map("{:,.3f}".format)
    FILTER["T%"] = FILTER["T%"].map("{:,.3f}".format)
    FILTER["A%"] = FILTER["A%"].map("{:,.3f}".format)
    FILTER["ins%"] = FILTER["ins%"].map("{:,.3f}".format)
    FILTER["del%"] = FILTER["del%"].map("{:,.3f}".format)
    FILTER.drop_duplicates(subset=["#CHROM", "POS", "DEPTH"], keep="last")
    FILTER = FILTER[
        [
            "#CHROM",
            "POS",
            "C%",
            "G%",
            "T%",
            "A%",
            "ins%",
            "del%",
            "sum",
            "DEPTH",
            "GENE",
            "exone",
            "length",
            "strand",
            "refseq",
            "hgmd",
            "selection",
            "filt",
        ]
    ].sort_values(by=["#CHROM", "POS"], ascending=[True, True])

    FILTER["sum"] = FILTER["sum"].astype(int)
    FILTER["POS"] = FILTER["POS"].astype(int)
    FILTER["DEPTH"] = FILTER["DEPTH"].astype(int)
    FILTER["selection"] = FILTER["selection"].astype(int)
    FILTER["filt"] = FILTER["filt"].astype(int)
    FILTER["exone"] = FILTER["exone"].astype(int)
    FILTER["length"] = FILTER["length"].astype(int)
    FILTER["strand"] = FILTER["strand"].astype(int)
    FILTER2 = FILTER[FILTER["filt"] == 1]

    print("len FILTER:", len(FILTER2))
    print("FILTER OK!!!")

    folder_to_macroarea = join(folder, "temp/", "to_macroarea")
    FILTER2.to_csv(
        join(folder_to_macroarea, sample + "_" + str(i) + "_macroarea"),
        sep="\t",
        index=False,
    )
    return FILTER2


def filter_disease(folder, name, phenotype, COUNT, sample):
    folder_coverage = join(folder, "coverage/", sample)
    result = join(folder_coverage, sample + "_all")
    phenotype_ = phenotype[phenotype["sample"].astype(str) == str(sample)]
    a = phenotype_[["malattia", "gene"]].drop_duplicates()
    x1 = pd.DataFrame(a["gene"])
    b = COUNT  # [COUNT['GENE'].isin(x1['gene'])]
    b.drop("selection", axis=1, inplace=True)
    b.drop("filt", axis=1, inplace=True)
    b.to_csv(name + "_disease", sep="\t", header=False, index=False)
    return


def filter_disease_files(
    BED, vertical, folder, phenotype, COUNT, sample, buchiartificiali
):
    print("start final filter!!!")
    COUNT = COUNT[
        ["#CHROM", "POS", "C%", "G%", "T%", "A%", "ins%", "del%", "sum", "DEPTH"]
    ].sort_values(by=["#CHROM", "POS"], ascending=[True, True])
    phenotype_ = phenotype[phenotype["sample"].astype(str) == str(sample)]
    a = phenotype_[["malattia", "gene"]].drop_duplicates()
    x1 = pd.DataFrame(BED["GENE"].drop_duplicates())
    # x1 = pd.DataFrame(a['gene'])
    x2 = pd.DataFrame({"GENE": pd.Series(["CTD-3074O7.11", "TM4SF2"])})
    x = x1._append(x2)
    FILTER = vertical[vertical["GENE"].isin(x["GENE"])]

    b = pd.merge(COUNT, FILTER, on=["#CHROM", "POS"], how="right")
    for index, row in buchiartificiali.iterrows():
        x = row["#CHROM"]
        y = row["START"]
        z = row["END"]
        mask1 = (b["#CHROM"] == x) & (b["POS"] >= y) & (b["POS"] <= z)
        b.loc[mask1, "DEPTH"] = 0

    b["DEPTH"].fillna(0, inplace=True)
    b["sum"].fillna(0, inplace=True)
    try:
        b["filt"].fillna(0, inplace=True)
    except:
        b["filt"] = 0
    b["sum"] = b["sum"].astype(int)
    b["POS"] = b["POS"].astype(int)
    b["DEPTH"] = b["DEPTH"].astype(int)
    b["filt"] = b["filt"].astype(int)
    b["exone"] = b["exone"].astype(int)
    b["length"] = b["length"].astype(int)
    b["strand"] = b["strand"].astype(int)

    COUNT = b
    COUNT["sample"] = sample
    COUNT.fillna(0, inplace=True)
    COUNT.drop(["filt"], axis=1, inplace=True)
    print(len(COUNT))
    print("Sequence: ", len(COUNT))

    folder_coverage = join(folder, "coverage/", sample)
    if not os.path.exists(folder_coverage):
        os.makedirs(folder_coverage)

    result = join(folder_coverage, sample + "_all")
    result_buchi = join(folder_coverage, sample + "_buchi")
    result_non_buchi = join(folder_coverage, sample + "_marked")
    result_not_marked = join(folder_coverage, sample + "_only_0")

    cov = len(COUNT[COUNT["DEPTH"] >= 10])
    buchi = len(COUNT[COUNT["DEPTH"] < 10])
    # cov = len(COUNT[COUNT['DEPTH'] >= 20])
    # buchi = len(COUNT[COUNT['DEPTH'] < 20])
    tot = cov + buchi
    COUNT["sample"] = str(sample)

    # buchi = COUNT[COUNT['DEPTH'] < 20]
    # marked = COUNT[COUNT['DEPTH'] >= 20]
    buchi = COUNT[COUNT["DEPTH"] < 10]
    marked = COUNT[COUNT["DEPTH"] >= 10]
    not_marked = COUNT[COUNT["DEPTH"] == 0]

    buchi.to_csv(result_buchi, sep="\t", index=False)
    marked.to_csv(result_non_buchi, sep="\t", index=False)
    not_marked.to_csv(result_not_marked, sep="\t", index=False)
    COUNT.to_csv(result, sep="\t", index=False)

    print("Len COV > 10: ", cov)
    print("Len BUCHI: ", len(buchi))
    try:
        print(
            "% Sequence COVERED: ", "{:,.1f}".format(float(cov) / float(tot) * 100), "%"
        )
    except ZeroDivisionError:
        print(0)

    print("Len disease: ", len(COUNT))
    return


def filter_macroarea_files(vertical_macro, folder, sample, buchiartificiali):
    print("start filter_macroarea_files!!!")
    folder_to_macroarea = join(folder, "temp/", "to_macroarea")
    path_macro = join(folder_to_macroarea, sample + "_*_macroarea")
    # print(path_)
    files_macro = glob.glob(path_macro)
    macro_count = pd.DataFrame(
        columns=[
            "#CHROM",
            "POS",
            "C%",
            "G%",
            "T%",
            "A%",
            "ins%",
            "del%",
            "sum",
            "DEPTH",
            "GENE",
            "exone",
            "length",
            "strand",
            "refseq",
            "hgmd",
            "selection",
            "filt",
        ]
    )
    for file_macro in files_macro:
        # concat files
        print(file_macro)
        macro = pd.read_csv(file_macro, sep="\t")
        macro_count = macro_count._append(macro, ignore_index=True)
    print(macro_count)
    if True:
        macro_count = macro_count.sort_values(
            by=["#CHROM", "POS"], ascending=[True, True]
        )
        print(macro_count)
        b = macro_count
        for index, row in buchiartificiali.iterrows():
            x = row["#CHROM"]
            y = row["START"]
            z = row["END"]
            mask1 = (b["#CHROM"] == x) & (b["POS"] >= y) & (b["POS"] <= z)
            b.loc[mask1, "DEPTH"] = 0

        b["DEPTH"].fillna(0, inplace=True)
        b["sum"].fillna(0, inplace=True)
        b["filt"].fillna(0, inplace=True)
        b["sum"] = b["sum"].astype(int)
        b["POS"] = b["POS"].astype(int)
        b["DEPTH"] = b["DEPTH"].astype(int)
        b["filt"] = b["filt"].astype(int)
        b["exone"] = b["exone"].astype(int)
        b["length"] = b["length"].astype(int)
        b["strand"] = b["strand"].astype(int)

        COUNT = b
        COUNT["sample"] = sample
        COUNT.fillna(0, inplace=True)
        print(len(COUNT))
        print("Sequence: ", len(COUNT))

        folder_coverage = join(folder, "coverage/", sample)
        if not os.path.exists(folder_coverage):
            os.makedirs(folder_coverage)

        result = join(folder_coverage, sample + "_all_macro")
        result_buchi = join(folder_coverage, sample + "_buchi_macro")
        result_non_buchi = join(folder_coverage, sample + "_marked_macro")
        result_not_marked = join(folder_coverage, sample + "_only_0_macro")

        cov = len(COUNT[COUNT["DEPTH"] >= 10])
        # cov = len(COUNT[COUNT['DEPTH'] >= 20])
        buchi = len(COUNT[COUNT["DEPTH"] < 10])
        # buchi = len(COUNT[COUNT['DEPTH'] < 20])
        tot = cov + buchi
        COUNT["sample"] = str(sample)

        buchi = COUNT[COUNT["DEPTH"] < 10]
        marked = COUNT[COUNT["DEPTH"] >= 10]
        # buchi = COUNT[COUNT['DEPTH'] < 10]
        # marked = COUNT[COUNT['DEPTH'] >= 10]
        not_marked = COUNT[COUNT["DEPTH"] == 0]

        buchi.to_csv(result_buchi, sep="\t", index=False)
        marked.to_csv(result_non_buchi, sep="\t", index=False)
        not_marked.to_csv(result_not_marked, sep="\t", index=False)
        COUNT.to_csv(result, sep="\t", index=False)

        # print ('Len COV > 10: ',cov)
        print("Len BUCHI: ", len(buchi))
        try:
            print(
                "% Sequence COVERED: ",
                "{:,.1f}".format(float(cov) / float(tot) * 100),
                "%",
            )
        except ZeroDivisionError:
            print(0)

        print("Len disease: ", len(COUNT))
    return
