import pandas as pd
import os
import numpy as np
import csv

def create_vertical(bed_df):
    bed_df['POS'] = [[x for x in range(bed_df.at[i, 'START'], bed_df.at[i, 'END'] + 1)] for i in bed_df.index]
    vertical = bed_df.drop(['START', 'END'], axis=1).explode('POS').reset_index()
    return vertical.drop_duplicates(vertical.drop_duplicates(subset=['#CHROM','POS'],inplace=True))

def count_coverage_old(TO_COUNT):
    TO_COUNT = TO_COUNT[["#CHROM", "POS", "DEPTH", "CALL"]]

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
    
    return TO_COUNT_final

def count_coverage(chunk):
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
    df = pd.concat([chunk[["#CHROM", "POS", "DEPTH", "CALL", "sample"]], coverage_counts], axis = 1)
  
    df = df[["#CHROM", "POS", "DEPTH", "C", "G", "T", "A", "ins", "del", "sum", "sample"]]
    
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

def filter_for_panel(vertical, coverage):
    vertical["filt"] = 1
    print("-------FILTERING-------")
    FILTER = pd.merge(coverage, vertical, on=["#CHROM", "POS"], how="left")
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

bam = '/home/magi/PROJECT/diagnosys/RESULT/13_Feb_2024_OCULARE/bam/E675.2023_final.bam'
mpileup_out = '/home/magi/PROJECT/diagnosys/bin_jurgen/E675.2023_to_count'

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

vertical = pd.read_csv("/home/magi/PROJECT/diagnosys/RESULT/13_Feb_2024_OCULARE/pheno/vertical_E675.2023", sep="\t")
print(vertical)
chunky = None
chunky_vertical = None
sample_name = "E675.2023"
for index, chunk in enumerate(a):
    
    chunk["sample"] = sample_name
    chunk.rename(columns={"CHROM": "#CHROM"}, inplace=True)
    # INNER JOIN WITH VERTICAL
    chunky_vertical = pd.merge(chunk, vertical, on=["#CHROM", "POS"], how='inner')
    chunky = chunk
    coverage = count_coverage(chunk)
    coverage_old = count_coverage_old(chunk)
    coverage_filt = count_coverage(chunky_vertical)
    #coverage.to_csv("E530.2019_coverage_counts_{}".format(index), sep="\t", index=False)
    if len(chunky_vertical) != 0:
        break

coverage_old_filtered = pd.merge(coverage_old, vertical, on=["#CHROM", "POS"], how='inner')
coverage_old_filtered.to_csv("coverage_old_filtered_inner_join.csv", sep="\t")

coverage_old_filter_for_panel = filter_for_panel(vertical, coverage_old)
coverage_old_filter_for_panel.to_csv("coverage_old_filter_for_panel.csv", sep="\t")

chunky_vertical.to_csv("chunky_vertical.csv", sep="\t")
coverage_filt.to_csv("coverage_filtered.csv", sep="\t")
print(chunky_vertical)