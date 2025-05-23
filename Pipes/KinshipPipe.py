from Pipes.ParallelPipe import ParallelPipe
import os
import pandas as pd
import numpy as np
import glob
from os.path import join, isfile
import dir_tree
import utils
from typing import Optional

class KinshipPipe(ParallelPipe):
    """Single-sample kinship analysis (paternity / maternity test).

    This Pipe is the OO/pipe translation of the legacy *bin/kinship_test.py*.
    One instance handles **one** proband sample. If parents are present in the
    *sample_list_FAM.csv* file, their VCFs are fetched automatically.

    Output (written under the *phase/* directory of the project):
        <sample>_all_annot.csv     # filtered, panel-restricted VCF for every family member
        <sample>_family.csv        # merged trio before incompatibility filter
        <sample>_family_cleaned.vcf# Mendelian-consistent variants only
    """

    def __init__(self):
        super().__init__()

    # ------------------------------------------------------------------ #
    # Public entry‐point (Pipeline contract)
    # ------------------------------------------------------------------ #
    def process(self, **kwargs):
        # multiprocessing diagnostics helpers
        self.thread_id = os.getpid()

        # mandatory kwargs provided by previous pipes / wrapper
        self.sample = kwargs.pop("sample")            # Entities.Sample object
        self.panel = kwargs.pop("panel", None)        # e.g. CANCER, OCULAR…
        self.dest   = kwargs.pop("dest")              # server id (b/r/s/z)
        self.genome = kwargs.pop("genome", "geno38") # geno37 | geno38
        self.lock   = kwargs.pop("lock", None)        # optional multiprocessing lock

        self.run()

        # push enriched kwargs downstream
        kwargs.update({"sample": self.sample,
                       "panel" : self.panel,
                       "dest"  : self.dest,
                       "genome": self.genome})
        return kwargs

    # ------------------------------------------------------------------ #
    # Main logic
    # ------------------------------------------------------------------ #
    def run(self):
        # Resolve project directories from dir_tree helper
        folder_root   = dir_tree.principal_directory.path
        folder_vcf    = dir_tree.principal_directory.vcf.path
        self.folder_phase = dir_tree.principal_directory.phase.path

        os.makedirs(self.folder_phase, exist_ok=True)

        # Vertical panel (gene list)
        vertical_path = utils.get_vertical_macro(self.panel)
        self.vertical = (pd.read_csv(vertical_path, sep="\t")
                         if vertical_path and isfile(vertical_path)
                         else pd.DataFrame(columns=["#CHROM", "POS"]))

        # Family information
        fam_csv = os.path.join(folder_root, "sample_list_FAM.csv")
        if not isfile(fam_csv):
            self.thread_print(f"No sample_list_FAM.csv – skip kinship for {self.sample.name}")
            return
        sample_family = pd.read_csv(fam_csv, sep="\t")
        sample_family = sample_family.astype({"Sample Id": str, "Riferimento": str})

        # Identify trio (proband, father, mother)
        prob_id = str(self.sample.name)
        family_df = sample_family[sample_family["Riferimento"] == prob_id]
        if family_df.empty:
            self.thread_print(f"{prob_id} not found in family list – skip.")
            return

        try:
            father_id = family_df[family_df["Familiarita"] == "PADRE"].iloc[0, 0]
            mother_id = family_df[family_df["Familiarita"] == "MADRE"].iloc[0, 0]
        except IndexError:
            self.thread_print(f"Parents not found for {prob_id} – skip.")
            return

        # ------------------------------------------------------------------ #
        # 1. Generate <sample>_all_annot.csv for every family member
        # ------------------------------------------------------------------ #
        for fam_id in family_df["Sample Id"].astype(str):
            unified_p = join(folder_vcf, f"{fam_id}_unfied_all.vcf")
            samt_p    = join(folder_vcf, f"{fam_id}_samt_all.vcf")
            if not (isfile(unified_p) and isfile(samt_p)):
                self.thread_print(f"VCFs missing for {fam_id} – skipped.")
                continue
            unified_df = pd.read_csv(unified_p, sep="\t", header=0)
            samt_df    = pd.read_csv(samt_p   , sep="\t", header=0)
            self.filter_vcf(samt_df, unified_df, fam_id)

        # ------------------------------------------------------------------ #
        # 2. Perform kinship checks on the trio
        # ------------------------------------------------------------------ #
        self.pat_mat_test(prob_id, mother_id, father_id)

    # ------------------------------------------------------------------ #
    # Helper functions (adapted from kinship_test.py)
    # ------------------------------------------------------------------ #
    def filter_vcf(self, samt_df: pd.DataFrame, unified_df: pd.DataFrame, sample_id: str):
        """Merge samtools & haplotypecaller VCFs, restrict to vertical panel,
        hard-filter on QUAL/DP, normalise GT field, write *_all_annot.csv."""

        # In our implementation we will change the order, because we will not give priority
        # to the samt_df, but to the unified_df. So instead, of going trhough each line, 
        # and replacing _x with _y, we will change the order of the suffixes, from ("_x", "_y")
        # to ("_y", "_x")
        merged = pd.merge(samt_df, unified_df,
                          on=["#CHROM", "POS", "ID", "GENE", "exone", "length",
                              "strand", "refseq", "hgmd"], how="outer", suffixes=("_y", "_x"))

        # Consolidate genotype and core VCF fields
        merged[sample_id] = merged[f"{sample_id}_x"].combine_first(merged[f"{sample_id}_y"])
        merged["QUAL"]   = merged["QUAL_x"].fillna(-999).astype(int)
        merged["FILTER"] = merged["FILTER_x"].fillna(".")
        merged["INFO"]   = merged["INFO_x"].combine_first(merged["INFO_y"])
        merged["REF"]    = merged["REF_x"].combine_first(merged["REF_y"])
        merged["ALT"]    = merged["ALT_x"].combine_first(merged["ALT_y"])

        # First ALT allele, drop symbolic, remove rows with only "*"
        merged["ALT"] = np.where(merged["ALT"].str.split(",").str.get(0) == "*",
                                 merged["ALT"].str.split(",").str.get(1),
                                 merged["ALT"].str.split(",").str.get(0))
        merged = merged[merged["ALT"] != "*"]

        # Depth & quality filters and INDEL exclusion
        merged["POS"]   = merged["POS"].astype(int)
        merged         = merged[merged["INFO"].str.split(";").str.get(0) != "INDEL"]
        merged["DEPTH"] = (merged["INFO"].str.split("DP=").str.get(1)
                                               .str.split(";").str.get(0).fillna(-999).astype(int))
        merged         = merged[(merged["QUAL"] >= 30) & (merged["DEPTH"] >= 20)]

        # Column selection & genotype normalisation
        merged = merged[["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                         "INFO", "FORMAT", sample_id, "GENE", "exone", "refseq"]]
        merged[sample_id] = (merged[sample_id].fillna("./.:")
                                            .str.split(":").str.get(0))
        # Map non-canonical allele codes to ./.
        for bad in ["2", "3", "4"]:
            merged[sample_id] = (merged[sample_id]
                                 .str.replace(fr"{bad}/", "./", regex=True)
                                 .str.replace(fr"/{bad}", "/.", regex=True))

        # Intersect with vertical panel if provided
        if not self.vertical.empty:
            merged = pd.merge(merged, self.vertical[["#CHROM", "POS"]], on=["#CHROM", "POS"], how="inner")

        out_csv = join(self.folder_phase, f"{sample_id}_all_annot.csv")
        merged.to_csv(out_csv, sep="\t", index=False)
        self.thread_print(f"Wrote {out_csv}")

    def pat_mat_test(self, proband: str, mother: str, father: str):
        """Merge *_all_annot.csv of the trio and remove Mendelian-incompatible
        variants. Result is written to <proband>_family_cleaned.vcf"""

        # Helper function to read the csv files complying to the output format of filter_vcf (ending in _all_annot.csv)
        def read_csv(sid):
            return pd.read_csv(join(self.folder_phase, f"{sid}_all_annot.csv"), sep="\t")

        try:
            prob_df  = read_csv(proband)
            fath_df  = read_csv(father)
            moth_df  = read_csv(mother)
        except FileNotFoundError:
            self.thread_print("Annotated CSV missing – cannot perform kinship test.")
            return

        # Build trio table (similar to legacy script but vectorised)
        fam = (prob_df.merge(fath_df, how="outer", on=["#CHROM", "POS", "ID", "REF", "ALT"], suffixes=("", "_f"))
                      .merge(moth_df, how="outer", on=["#CHROM", "POS", "ID", "REF", "ALT"], suffixes=("", "_m")))

        # Unify annotations columns
        for col in ["QUAL", "FILTER", "INFO", "FORMAT", "GENE", "refseq", "exone"]:
            alt_cols = [f"{col}_f", f"{col}_m"]
            base = fam[col].combine_first(fam[alt_cols[0]]).combine_first(fam[alt_cols[1]])
            fam[col] = base
            fam.drop(columns=[c for c in alt_cols if c in fam], inplace=True)

        # Fill missing GT with 0/0
        for sid in [proband, father, mother]:
            fam[sid] = fam[sid].fillna("0/0")

        fam = fam[["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                   father, proband, mother, "GENE", "refseq", "exone"]]

        trio_csv = join(self.folder_phase, f"{proband}_family.csv")
        fam.sort_values(["#CHROM", "POS"]).to_csv(trio_csv, sep="\t", index=False)

        # ------------------------------------------------------------------ #
        # Mendelian incompatibility filter (vectorised replacement for the
        # giant list of conditions of the legacy script)
        # ------------------------------------------------------------------ #
        incompatible = [
            ("1/1", "1/1", "1/1"), ("1/1", "1/1", "0/0"),
            ("0/0", "1/1", "1/1"), ("0/0", "1/1", "0/1"),
            ("0/0", "1/1", "1/0"), ("0/1", "1/1", "0/0"),
            ("1/0", "1/1", "0/0"), ("0/0", "1/1", "0/0"),
            ("1/1", "0/1", "1/1"), ("1/1", "1/0", "1/1"),
            ("0/1", "0/0", "0/0"), ("1/0", "0/0", "0/0"),
            ("0/0", "0/0", "0/1"), ("0/0", "0/0", "1/0"),
            ("1/1", "0/0", "1/1"), ("0/1", "0/0", "0/1"),
            ("0/1", "0/0", "1/0"), ("1/0", "0/0", "0/1"),
            ("1/0", "0/0", "1/0"), ("1/1", "0/0", "0/1"),
            ("1/1", "0/0", "1/0"), ("0/1", "0/0", "1/1")
        ]
        mask = pd.Series(False, index=fam.index)
        for f, p, m in incompatible:
            mask |= (fam[father] == f) & (fam[proband] == p) & (fam[mother] == m)
        fam_clean = fam[~mask]

        out_vcf = join(self.folder_phase, f"{proband}_family_cleaned.vcf")
        fam_clean.to_csv(out_vcf, sep="\t", index=False)
        self.thread_print(f"Wrote {out_vcf}")

        self._write_parentage_stats(fam_clean, proband, mother, father)


    # ------------------------------------------------------------------ #
    # Private methods
    # ------------------------------------------------------------------ #
    def _write_parentage_stats(self, fam: pd.DataFrame, proband: str, mother: str, father: str):
        """
        Derive the same statistics as the legacy compatibilita_pat / compatibilita_mat
        lines and append them to <coverage>/<proband>/<proband>_stat_cov.csv.
        """

        # Autosomals only
        fam_a = fam[~fam["#CHROM"].isin(["chrX","chrY"])]

        # Homo variants counts
        v_fat = fam_a[(fam_a[father] == "1/1") & (fam_a[mother] == "0/0")]
        v_mot = fam_a[(fam_a[mother] == "1/1") & (fam_a[father] == "0/0")]
        
        wrong_fat = v_fat[v_fat[proband] == "0/0"]
        wrong_mot = v_mot[v_mot[proband] == "0/0"]

        # Helper to calculate percentage (kinship_score)
        def pct(wrong, total) -> str:
            if total.shape[0] == 0:
                return "0.0% (0/0)"
            
            good = 100.0 - (wrong.shape[0] / total.shape[0]) * 100.0
            return f"{round(good, 2)}% ({wrong.shape[0]}/{total.shape[0]})"
        
    
        # ------------------------------------------------------------------#
        # Update the stat_cov file
        # ------------------------------------------------------------------#
        cov_dir = join(dir_tree.principal_directory.coverage.path, proband)
        os.makedirs(cov_dir, exist_ok=True)
        cov_csv = join(cov_dir, f"{proband}_stat_cov.csv")

        cov_df: pd.DataFrame
        if isfile(cov_csv):
            cov_df = pd.read_csv(cov_csv, sep="\t")
        else:
            cov_df = pd.DataFrame(index=[0])

        cov_df["compatibilita_pat"] = compat_pat
        cov_df["compatibilita_mat"] = compat_mat
        cov_df.to_csv(cov_csv, sep="\t", index=False)

        self.thread_print(
            f"Parentage stats → {compat_pat} | {compat_mat}  →  {cov_csv}"
        )



    # ------------------------------------------------------------------ #
    # TODO: add the evaluate_contamination method
    # ------------------------------------------------------------------ #

    # ------------------------------------------------------------------ #
    # Convenience wrappers around utils helpers
    # ------------------------------------------------------------------ #
    def thread_print(self, msg):
        utils.thread_print(getattr(self, "thread_id", None), msg)
