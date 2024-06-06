
from Pipe import Pipe
import utils
import os
import pandas as pd
import csv

import dir_tree

cols = ["#CHROM", "POS", "C%", "G%", "T%", "A%", "ins%", "del%", "sum", "DEPTH", "GENE", "exone", "length", "strand", "refseq", "hgmd", "sample"]

class CoveragePipeTest(Pipe):

    def __init__(self) -> None:
        super().__init__()

    def process(self, **kwargs):
        self.panel = kwargs.pop("panel")
        self.sample = kwargs.pop("sample")
        self.load_prereq()
        #sample_name = kwargs.pop("sample")
        # sample_data_dir = dir_tree.sample_data.path
        # sample_json_filename = "{}.json".format(sample_name)
        # sample_json = os.path.join(sample_data_dir, sample_json_filename)
        # self.sample_data = Sample.fromJSON(sample_json(sample_json))

        self.count_unbalance()
        count_disease, count_sex = self.merge_chunks()
        

        kwargs.update({"panel": self.panel, "sample": self.sample})
        return kwargs

    def filter_disease_files(BED, vertical, folder, phenotype, COUNT, sample, buchiartificiali):
        COUNT = COUNT[["#CHROM", "POS", "C%", "G%", "T%", "A%", "ins%", "del%", "sum", "DEPTH"]].sort_values(by=["#CHROM", "POS"], ascending=[True, True])
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

        folder_coverage = os.path.join(folder, "coverage/", sample)
        if not os.path.exists(folder_coverage):
            os.makedirs(folder_coverage)

        result = os.path.join(folder_coverage, sample + "_all")
        result_buchi = os.path.join(folder_coverage, sample + "_buchi")
        result_non_buchi = os.path.join(folder_coverage, sample + "_marked")
        result_not_marked = os.path.join(folder_coverage, sample + "_only_0")

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


    def load_prereq(self):
        self.vertical_df, self.verticalX_df, self.vertical_macro_df = self.get_vertical_dataframes()
        self.buchiartificali = utils.get_buchiartificiali()

    def merge_chunks(self):
        #TODO: add full paths to the *_disease and *_SEX
        final_disease_filename = os.path.join(dir_tree.temp.path, "{}_final_disease".format(str(self.sample["name"])))
        os.system(" ".join(["cat", "*_disease", ">", final_disease_filename]))
        count_disease = pd.read_csv(final_disease_filename, sep="\t", header=None, names=cols)
        count_disease.fillna(0, inplace=True)

        final_sex_filename = os.path.join(dir_tree.temp.path, "{}_final_sex".format(str(self.sample["name"])))
        os.system(" ".join(["cat", "*_SEX", ">", final_sex_filename]))
        count_sex = pd.read_csv(final_sex_filename, sep="\t", header=None, names=cols)
        count_sex.fillna(0, inplace=True)

        return count_disease, count_sex


    def count_unbalance(self):
        bam = self.get_bam_filename()
        sample_name = str(self.sample["name"])
        mpileup_out = os.path.join(dir_tree.temp_dir, sample_name + "_to_count")
        folder_to_count = dir_tree.temp_dir.to_count_dir
        folder_to_macroarea = dir_tree.temp_dir.to_macroarea_dir
        #file_vcf = join(folder, "temp/", sample_name + "_samt.vcf")

        os.system(" ".join(["samtools", "mpileup", "-Q 0 -q 0 -d10000000 -L 100000 -A", bam, ">", mpileup_out]))

        mpileup_out_chunks = pd.read_csv(
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
        
        for index, chunk in enumerate(mpileup_out_chunks):
            chunk["sample"] = sample_name
            chunk.rename(columns={"CHROM": "#CHROM"}, inplace=True)

            # INNER JOIN WITH VERTICAL
            chunky_vertical = pd.merge(chunk, self.vertical_df, on=["#CHROM", "POS"], how='inner')
            chunky_verticalX = pd.merge(chunk, self.verticalX_df, on=["#CHROM", "POS"], how='inner')

            coverage_vertical = self.count_coverage(chunky_vertical)
            coverage_verticalX = self.count_coverage(chunky_verticalX)

            coverage_vertical.to_csv(os.path.join(folder_to_count, sample_name + "_" + str(index)) + "_disease", sep="\t", header=False, index=False)
            coverage_verticalX.to_csv(os.path.join(folder_to_count, sample_name + "_" + str(index)) + "_SEX", sep="\t", header=False, index=False)
            
            if self.vertical_macro_df is not None:
                chunky_vertical_macro = pd.merge(chunk, self.vertical_macro_df, on=["#CHROM", "POS"], how='inner')
                coverage_vertical_macro = self.count_coverage(chunky_vertical_macro)
                coverage_vertical_macro.to_csv(folder_to_macroarea, sample_name + "_" + str(index) + "_macroarea", sep="\t", index=False,)

    def count_coverage(self, chunk):
        #coverage = chunk[["CHROM", "POS", "DEPTH", "CALL"]]
    
        # col 5 of samtools mpileup output, each row represents the bases in the reads,
        # covering that position in the reference genome
        read_base_patterns = {".": r"[.]", 
                            ",": r"[*]", 
                            "$": r"[$]", 
                            "C": r"[cC]", 
                            "G": r"[gG]", 
                            "T": r"[tT]", 
                            "A": r"[aA]", 
                            "ins": r"\+[0-9]+[ACGTNacgtn]+", 
                            "del": r"\-[0-9]+[ACGTNacgtn]+", 
                            "other": r"[*><$^]"
                            }
        
        data = []
        for key, value in read_base_patterns.items():
            data.append(chunk["CALL"].str.count(value).rename(key))
            
        coverage_counts = pd.concat(data, axis=1)
        coverage_counts['sum'] = coverage_counts[[".", ",", "C", "G", "T", "A", "ins", "del"]].sum(axis = 1)
        df = pd.concat([chunk[["#CHROM", "POS", "DEPTH", "CALL"]], coverage_counts], axis = 1)
    
        df = df[["#CHROM", "POS", "DEPTH", "C", "G", "T", "A", "ins", "del", "sum"]]
        
        df["DEPTH"].fillna(0, inplace=True)
        df["sum"].fillna(0, inplace=True)
        df["sum"] = df["sum"].astype(int)
        df["POS"] = df["POS"].astype(int)
        df["DEPTH"] = df["DEPTH"].astype(int)
        df["selection"] = 1
        
        df["C%"] = df["C"].astype(float) / df["sum"].astype(float)
        df["G%"] = df["G"].astype(float) / df["sum"].astype(float)
        df["T%"] = df["T"].astype(float) / df["sum"].astype(float)
        df["A%"] = df["A"].astype(float) / df["sum"].astype(float)
        df["ins%"] = df["ins"].astype(float) / df["sum"].astype(float)
        df["del%"] = df["del"].astype(float) / df["sum"].astype(float)
        df["C%"] = df["C%"].map("{:,.3f}".format)
        df["G%"] = df["G%"].map("{:,.3f}".format)
        df["T%"] = df["T%"].map("{:,.3f}".format)
        df["A%"] = df["A%"].map("{:,.3f}".format)
        df["ins%"] = df["ins%"].map("{:,.3f}".format)
        df["del%"] = df["del%"].map("{:,.3f}".format)

        return df

    def get_vertical_dataframes(self):
        vertical_macro_df = None
        
        vertical_path, verticalX_path, vertical_macro_path = self.get_vertical_files()
        if vertical_path is None:
            raise Exception("No vertical files for sample {}".format(str(self.sample["name"])))
        if vertical_macro_path is not None:
            vertical_macro_df = pd.read_csv(vertical_macro_path, sep="\t")
        else:
            print("No vertical_macro for panel {}".format(self.panel))
        
        vertical_df = pd.read_csv(vertical_path, sep="\t")
        verticalX_df = pd.read_csv(verticalX_path, sep="\t")

        return vertical_df, verticalX_df, vertical_macro_df

    def get_vertical_files(self):
        """ Check the pheno/ directory if the vertical and verticalX are present. Then, check if virtual_macro is also present. See config.cfg for the
        directories of vertical_macro files."""
        sample = self.sample
        _vertical_path = None
        _verticalX_path = None
        _vertical_macro_path = None
        pheno_dir = dir_tree.pheno_dir # ideally, we want to store each directory globally during runtime

        _vertical_macro_path = utils.get_vertical_macro(self.panel)

        if "vertical" not in sample.keys() or sample["vertical"] is None:
            vertical_filename = "vertical_{}".format(str(sample["name"]))
            verticalX_filename = "verticalX_{}".format(str(sample["name"]))
            if os.path.isfile(os.path.join(pheno_dir, vertical_filename)):
                _vertical_path = os.path.join(pheno_dir, vertical_filename)
            if os.path.isfile(os.path.join(pheno_dir, verticalX_filename)):
                _verticalX_path = os.path.join(pheno_dir, verticalX_filename)

            return _vertical_path, _verticalX_path, _vertical_macro_path

        return sample["vertical"], sample["verticalX"], _vertical_macro_path

    def get_bam_filename(self):
        sample = self.sample
        bam_dir = dir_tree.bam_dir

        if "bam" not in sample.keys() or sample["bam"] is None:
            # Search in the directory tree if we can find bam files for our sample
            bam_filename = "{}_final.bam".format(str(sample["name"]))
            if os.path.isfile(os.path.join(bam_dir, bam_filename)):
                return os.path.join(bam_dir, bam_filename)
            
        return sample["bam"]