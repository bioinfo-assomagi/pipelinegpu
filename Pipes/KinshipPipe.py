from Pipes.ParallelPipe import ParallelPipe
import os
import pandas as pd
import numpy as np
import glob
from os.path import join, isfile
import dir_tree
import utils
from typing import Optional
from typing import Union
from pathlib import Path
import logging
import sys

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
        self._logger = logging.getLogger(__name__)
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
        sample = str(self.sample.name)
        # Resolve project directories from dir_tree helper
        folder_root   = dir_tree.principal_directory.path
        folder_vcf    = dir_tree.principal_directory.vcf.path
        self.folder_phase = dir_tree.principal_directory.phase.path
        folder_cov = dir_tree.principal_directory.coverage.path

        cov_dir = Path(folder_cov) / sample
        vcf_dir = Path(folder_vcf)
        stat_cov_path = cov_dir / f"{sample}_stat_cov.csv"
        vcf_unified_path = vcf_dir / f"{sample}_unfied_all.vcf"

        os.makedirs(self.folder_phase, exist_ok=True)

        # Vertical panel (gene list)
        vertical_path = utils.get_vertical_macro(self.panel)
        self.vertical = (pd.read_csv(vertical_path, sep="\t")
                         if vertical_path and isfile(vertical_path)
                         else pd.DataFrame(columns=["#CHROM", "POS"]))
        # self._logger.debug(f"Vertical panel: {self.vertical}")
        

        # ---------------------------------------------------------------------------------------------------------------------------- # 
        # 0. Contamination check - this will run for the input sample - no bothering with searching for other members of the family
        # ---------------------------------------------------------------------------------------------------------------------------- #

        try:
            stat_cov = pd.read_csv(stat_cov_path, sep="\t")
        except FileNotFoundError as err:
            self._logger.error("({} thread): Coverage file not found  → {}".format(self.sample.name,  err.filename))
            return
        stat_cov = (
            stat_cov.assign(
                compatibilita_pat="-999",
                compatibilita_mat="-999",
                probcontaminazione="-999",
            )
        )

        stat_cov['cov_medio'] = stat_cov['cov_medio'].replace(',','',regex=True)

        self._logger.debug("({} thread): stat_cov prepared".format(self.sample.name))
        
        try:
            samt = pd.read_csv(vcf_unified_path, sep="\t")
            prob_cont = self.evaluate_contamination(
                samt, sample_col=sample, folder_cov=folder_cov
            )
            stat_cov["probcontaminazione"] = f"{prob_cont:.2f}"
        except FileNotFoundError:
            self._logger.error("({} thread): Unified VCF not found  → {}".format(self.sample.name,  vcf_unified_path))
        except Exception as exc:                     
            self._logger.error("({} thread): Contamination eval failed → {}".format(self.sample.name, exc))        

        # write back both to pipeline folder and Django export
        stat_cov.to_csv(stat_cov_path, sep="\t", index=False)
        #(path_django / f"{sample}_stat_cov.csv").write_text(stat_cov.to_csv(sep="\t", index=False))
        self._logger.debug("(%s thread):Updated stat_cov for %s; the next methods will work on the updated stat_cov found in %s", self.sample.name, self.sample.name, stat_cov_path)    



        # ---------------------------------------------------------------------------------------------------------------------------- # 
        # 0.5. Family information
        # ---------------------------------------------------------------------------------------------------------------------------- #
        fam_csv = os.path.join(folder_root, "sample_list_FAM.csv")
        if not isfile(fam_csv):
            self._logger.warning(f"No sample_list_FAM.csv – skip kinship for {self.sample.name}")
            return
        sample_family = pd.read_csv(fam_csv, sep="\t")
        sample_family = sample_family.astype({"Sample Id": str, "Riferimento": str})

        # Identify trio (proband, father, mother)
        # What happens here is that if the sample is not the proband, it will not be analysed further to asses paternity, maternity, or undergo any type of filtering
        # Only probands pass here. 
        # So basically, if the pipeline input is a list of samples belonging to the same family; with one proband, only the proband will pass here; and undergo filtering
        # and paternity matenrity tests; additionally, the siblings will undergo patternity and maternity tests, but no filtering.

        # If there are two probands, they will be executed in parallel, called from the WrapperPipe (VariantFilterPipe)
        # and TODO: the pat_mat_test of siblings must be considered cautiosly
        prob_id = str(self.sample.name)
        family_df = sample_family[sample_family["Riferimento"] == prob_id]
        if family_df.empty:
            self._logger.warning("({} thread): {} not found as riferimento in family list.".format(self.sample.name, prob_id))
            self._logger.info("[{} thread] Checking if we are dealing with a sibling ... ".format(self.sample.name))

            # check if we are dealing with a sibling
            current_sample_familarity_data = sample_family[sample_family["Sample Id"] == prob_id]
            current_sample_familarita = current_sample_familarity_data["Familiarita"].iloc[0]
            current_sample_affeto = current_sample_familarity_data["Affeto"].iloc[0]
            if not ((current_sample_familarita == "FRATELLO" or current_sample_familarita == "SORELLA") and current_sample_affeto == "SI"):
                self._logger.info("({} thread): {} is not a sibling, neither a proband, or not affected - skip kinship test.".format(self.sample.name, prob_id))
                return
            else:
                # get the nucleo of the family (XXX:nucleo or proband? - the riferimento)
                current_sample_riferimento = current_sample_familarity_data["Riferimento"].iloc[0]

                self._logger.info("({} thread): {} is a sibling of {} - compute pat_mat_test (kinship test), then return.".format(self.sample.name, self.sample.name, current_sample_riferimento))
                # get the father and the mother of the current_sample by looking at the father and mother that has the same riferimento as the current_sample_riferimento
                # if the current sample is a sibling to the nucleo (riferimento) than they are supposed to share the same parents

                try:
                    current_sample_father = sample_family[(sample_family["Riferimento"] == current_sample_riferimento) & (sample_family["Familiarita"] == "PADRE")]["Sample Id"].iloc[0]
                    current_sample_mother = sample_family[(sample_family["Riferimento"] == current_sample_riferimento) & (sample_family["Familiarita"] == "MADRE")]["Sample Id"].iloc[0]
                except IndexError:
                    self._logger.warning("({} thread): Father or Mother not found for {} – skip.".format(self.sample.name, prob_id))
                    return
                self.pat_mat_test(prob_id, current_sample_mother, current_sample_father)

            return

        try:
            father_id = family_df[family_df["Familiarita"] == "PADRE"].iloc[0, 0]
            mother_id = family_df[family_df["Familiarita"] == "MADRE"].iloc[0, 0]
        except IndexError:
            self._logger.warning("({} thread): Parents not found for {} – skip.".format(self.sample.name, prob_id))
            return


       
        # ------------------------------------------------------------------ #
        # 1. Generate <sample>_all_annot.csv for every family member
        # ------------------------------------------------------------------ #
        for fam_id in family_df["Sample Id"].astype(str):
            unified_p = join(folder_vcf, f"{fam_id}_unfied_all.vcf")
            samt_p    = join(folder_vcf, f"{fam_id}_samt_all.vcf")
            #self._logger.debug("({} thread): Processing {}".format(self.sample.name, fam_id))
            if not (isfile(unified_p) and isfile(samt_p)):
                self._logger.warning("({} thread): VCFs missing for {} – skipped.".format(self.sample.name, fam_id))
                continue
            unified_df = pd.read_csv(unified_p, sep="\t", header=0)
            samt_df    = pd.read_csv(samt_p   , sep="\t", header=0)
            self.filter_vcf(samt_df, unified_df, fam_id, self.vertical, output_dir=self.folder_phase)

        # ------------------------------------------------------------------ #
        # 2. Perform kinship checks on the trio
        # ------------------------------------------------------------------ #
        # We should perform it for all children TODO: i.e. brothers and sisters of the probando?


        self.pat_mat_test(prob_id, mother_id, father_id)

        
        

    # ------------------------------------------------------------------ #
    # Functions adapted from kinship_test.py in the /diagnosys/bin dir
    # ------------------------------------------------------------------ #
    def filter_vcf(self,               
        samt_sample: pd.DataFrame,
        unified_sample: pd.DataFrame,
        sample: str,
        vertical: pd.DataFrame,
        output_dir: Union[str, Path, None] = None,
        qual_threshold: int = 30,
        depth_threshold: int = 20,
    ) -> pd.DataFrame:
        """
        Merge SAMtools and UnifiedGenotyper calls, tidy the record and perform
        simple quality filters.

        Parameters
        ----------
        samt_sample, unified_sample
            VCF–like tables produced by SAMtools and UnifiedGenotyper.  They must
            contain *at least* the columns used as merge keys plus the per-sample
            genotype column named in *sample*.
        sample
            Name of the genotype column to keep (usually the sample identifier).
        vertical
            Extra annotations keyed by ``#CHROM`` + ``POS`` (e.g., target-region
            specification).
        output_dir
            If given, write the resulting table as ``<sample>_all_annot.csv`` to
            that folder (directory is created if necessary).
        qual_threshold, depth_threshold
            Minimum values for the ``QUAL`` and ``DP`` (depth) fields.

        Returns
        -------
        pd.DataFrame
            The cleaned and filtered data.
        """
        # ------------------------------------------------------------------ merge
        merge_keys = [
            "#CHROM",
            "POS",
            "ID",
            "GENE",
            "exone",
            "length",
            "strand",
            "refseq",
            "hgmd",
        ]
        df = pd.merge(
            samt_sample,
            unified_sample,
            how="outer",
            on=merge_keys,
            suffixes=("_sam", "_uni"),
        )
        
        # -------------------------------------------------------- value “coalesce”
        def coalesce(col_left: str, col_right: str, dest: str, dtype = None):
            """Take value from left column, fall back to right column."""
            series = df[col_left].where(df[col_left].notna(), df[col_right])
            df[dest] = series.astype(dtype) if dtype else series

        coalesce(sample + "_sam", sample + "_uni", sample)
        coalesce("QUAL_sam", "QUAL_uni", "QUAL", dtype=int)
        coalesce("FILTER_uni", "FILTER_uni", "FILTER")            # just copy
        coalesce("INFO_sam", "INFO_uni", "INFO")
        coalesce("REF_sam", "REF_uni", "REF")
        coalesce("ALT_sam", "ALT_uni", "ALT")

        # -------------------------------------------------------- ALT field tidy-up
        alt_split = df["ALT"].str.split(",", n=1, expand=True)
        df["ALT"] = np.where(alt_split[0] == "*", alt_split[1], alt_split[0])
        df = df.query("ALT != '*'")

        df["POS"] = df["POS"].astype("Int64")

        # ------------------------------ drop INDEL records, create DEPTH, FILTER QC
        df = df[df["INFO"].str.split(";", n=1).str[0] != "INDEL"]

        df["DEPTH"] = (
            df["INFO"]
            .str.extract(r"DP=(\d+)", expand=False)
            .astype("Int64")
            .fillna(-999)
        )

        df = df.query("QUAL >= @qual_threshold and DEPTH >= @depth_threshold")

        # ------------------------------------------------------------ genotype fix
        genotype = df[sample].fillna("./.:").str.split(":", n=1).str[0]

        # Replace secondary alleles (2, 3, 4 …) with missing
        # 2/0 -> ./0, 0/3 -> 0/. , 3/. -> ./.
        genotype = genotype.replace(
            {
                r"([23-9])/\.": "./.",
                r"\./([23-9])": "./.",
                r"([23-9])/([01])": "./\\2",
                r"([01])/([23-9])": "\\1/.",
            },
            regex=True,
        )
        df[sample] = genotype
        # --------------------------------------------------- final column ordering
        base_cols = [
            "#CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO",
            "FORMAT",
            sample,
            "GENE",
            "exone",
            "refseq",
        ]
        df["FORMAT"] = "GT"
        df = df[base_cols].fillna(".")

        # --------------------------------------------- merge with vertical table
        df = df.merge(vertical, how="inner", on=["#CHROM", "POS"], suffixes=("", "_y"))

        for col in ("GENE", "refseq", "exone"):
            # ensure both sides are object dtype
            left = df[col].astype(object).values  
            right = df[f"{col}_y"].astype(object).values

            # 1) Wherever left is "." or missing, take the right value
            combined = np.where((left == ".") | pd.isna(left), right, left)

            # 2) Wherever still missing, restore literal "."
            combined = np.where(pd.isna(combined), ".", combined)

            df[col] = combined

        # drop the helper columns
        df.drop(columns=[f"{c}_y" for c in ("GENE", "refseq", "exone")], inplace=True)

        # ------------------------------------------------------------- write CSV
        if output_dir is not None:
            out_dir = Path(output_dir).expanduser().resolve()
            out_dir.mkdir(parents=True, exist_ok=True)
            out_path = out_dir / f"{sample}_all_annot.csv"
            df.to_csv(out_path, sep="\t", index=False)
            self._logger.debug("({} thread) [KinshipPipe.filter_vcf]: Wrote \033[94m {} \033[0m".format(self.sample.name, out_path))

        return df 

    def pat_mat_test(self, proband: str, mother: str, father: str
    ):
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
            self._logger.warning("({} thread): Annotated CSV missing – cannot perform kinship test.".format(self.sample.name))
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

        
        # Mendelian incompatibility filter (vectorised replacement for the giant list of conditions of the legacy script)      
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

        self._logger.debug("({} thread) [KinshipPipe.pat_mat_test]: Wrote \033[94m {} \033[0m".format(self.sample.name, out_vcf))
        
        self._update_parental_compatibility(fam_clean, str(father), str(mother), str(proband), dir_tree.principal_directory.path)

    def _update_parental_compatibility(self,
        fam_clean: pd.DataFrame,
        father_col: str,
        mother_col: str,
        proband_col: str,
        project_root: Union[str, Path],
    ) -> pd.DataFrame:
        """
        Compute paternal / maternal genotype–compatibility rates and update the
        existing *stat_cov* file for *proband_col*.

        Parameters
        ----------
        fam_clean : pd.DataFrame
            VCF-like dataframe (already cleaned) that contains genotype columns
            for *father_col*, *mother_col* and *proband_col*.
        father_col, mother_col, proband_col : str
            Column names for father, mother and proband.
        project_root : str or pathlib.Path
            Project folder that contains::
                coverage/<proband>/<proband>_stat_cov.csv

        Returns
        -------
        pd.DataFrame
            The updated statistics table (also saved to disk).
        """
        
        # -------------------------------------------------------------- autosomes
        autosomal = fam_clean.loc[~fam_clean["#CHROM"].isin({"chrX", "chrY"})]

        # ------------------------------------------------ parental segregation sets
        pat_vars = autosomal.query(
            f"`{father_col}` == '1/1' and `{mother_col}` == '0/0'"
        )
        mat_vars = autosomal.query(
            f"`{father_col}` == '0/0' and `{mother_col}` == '1/1'"
        )

        
        # print(fam_clean[(fam_clean[father_col] == '1/1') & (fam_clean[mother_col] == '0/0')])


        wrong_pat = pat_vars.query(f"`{proband_col}` == '0/0'")
        wrong_mat = mat_vars.query(f"`{proband_col}` == '0/0'")
        
        def _fmt_compat(wrong: int, total: int) -> str:
            """Return 'xx.xx% (wrong/total)' or 'N/A (0/0)' when total is 0."""
            if total == 0:
                return "N/A (0/0)"
            pct = round(100.0 * (1 - wrong / total), 2)
            return f"{pct}% ({wrong}/{total})"

        # -------------------------------------------------------- load --> update
        proband_id = proband_col           # keep the external name
        stats_path = (
            Path(project_root)
            / "coverage"
            / proband_id
            / f"{proband_id}_stat_cov.csv"
        )

        if not stats_path.is_file():
            raise FileNotFoundError(stats_path)

        stat_cov = pd.read_csv(stats_path, sep="\t")

        stat_cov["compatibilita_pat"] = _fmt_compat(len(wrong_pat), len(pat_vars))
        stat_cov["compatibilita_mat"] = _fmt_compat(len(wrong_mat), len(mat_vars))

        stat_cov.to_csv(stats_path, sep="\t", index=False)
        self._logger.info("[{} thread] [KinshipPipe._update_parental_compatibility]: Compatibility stats updated → \033[94m {} \033[0m".format(self.sample.name, stats_path))

        return stat_cov

    # ------------------------------------------------------------------ #
    # TODO: add the evaluate_contamination method
    # ------------------------------------------------------------------ #
    def evaluate_contamination(self, samtall: pd.DataFrame, sample_col: str, folder_cov: Union[str, Path], qual_min: int = 18,
        depth_min: int = 25,
        hmean: float = 0.471707,
        hstd: float = 0.067168,
        low_thresh: float = 0.340057441652,
        high_thresh: float = 0.603356263732
    ):
        """
        Estimate sample-level DNA contamination.

        Parameters
        ----------
        samtall : pd.DataFrame
            VCF-like DataFrame with at least: ``POS``, ``REF``, ``ALT``, ``QUAL``,
            ``INFO`` and the genotype column named in *sample_col*.
        sample_col : str
            Column containing *GT:DP…* strings for the sample under analysis.
        folder_cov : str | pathlib.Path
            Directory that holds ``<sample>/<sample>_all_macro`` (coverage stats).
        qual_min, depth_min : int, default (18, 25)
            Minimum quality and depth thresholds for heterozygous sites.
        hmean, hstd : float
            Mean and std-dev used for Z-scoring the unbalance ratio.
        low_thresh, high_thresh : float
            Legacy absolute thresholds for extreme unbalance.

        Returns
        -------
        float
            Estimated contamination probability **as a percentage** (0‒100),
            rounded to two decimals.  If no suitable records are found, 0.0 is
            returned.
        """

        # hmean: float = 0.471707
        # hstd: float = 0.067168
        # low_thresh: float = 0.340057441652
        # high_thresh: float = 0.603356263732

        # preprocess

        df = samtall.copy()                       
        df["samtools"] = df[sample_col].str.split(":", n=1).str[0]

        geno_map = {"0/1": "het", "1/2": "het", "0/0": "homo_wild", "1/1": "homo"}
        df["samtools_geno"] = df["samtools"].map(geno_map).fillna("other")

        df["count_ref"] = df["REF"].str.len().astype(int)
        df["count_alt"] = df["ALT"].str.len().astype(int)

        df["START"] = np.where(df["count_ref"] > 1, df["POS"] + 1, df["POS"])

        ref_gt_alt = df["count_ref"] > df["count_alt"]
        alt_gt_ref = df["count_alt"] > df["count_ref"]
        df["END"] = np.select(
            [ref_gt_alt, alt_gt_ref],
            [
                df["POS"] + (df["count_ref"] - 1),          # deletion
                df["START"] + (df["count_alt"] - df["count_ref"]),  # insertion
            ],
            default=df["POS"],                              # single-base SVN
        )

        # -------------------------------------------------------------- type
        df["types"] = "SVN"
        del_mask = ref_gt_alt & (df["count_ref"] > 1)
        ins_mask = alt_gt_ref & (df["count_alt"] > 1)


        # helpers to slice variable-length strings row-wise
        df.loc[del_mask, "REF_2"] = df.loc[del_mask].apply(
            lambda r: r["REF"][r["count_alt"] :], axis=1
        )


        
        df["ALT_2"]  = pd.Series(dtype="object")
        df["types"]  = pd.Series(dtype="object")

        df.loc[del_mask, ["ALT_2", "types"]] = "-", "DELETION"

        df.loc[ins_mask, "ALT_2"] = df.loc[ins_mask].apply(
            lambda r: r["ALT"][r["count_ref"] :], axis=1
        )
        df.loc[ins_mask, ["REF_2", "types"]] = "-", "INSERTION"

        svn_mask = ~(del_mask | ins_mask)
        df.loc[svn_mask, ["REF_2", "ALT_2"]] = df.loc[svn_mask, ["REF", "ALT"]].values

        # -------------------------------------------------------------- depth
        df["DEPTH"] = (
            df["INFO"]
            .str.extract(r"DP=(\d+)", expand=False)
            .astype("Int64")
            .fillna(-999)
        )

        # -------------------------------------------------------------- keep only good heterozygous calls
        df = df.query(
            "samtools_geno == 'het' and QUAL >= @qual_min and DEPTH > @depth_min"
        )

        # -------------------------------------------------------- coverage file
        cov_path = Path(folder_cov) / sample_col / f"{sample_col}_all_macro"
        if not cov_path.exists():
            raise FileNotFoundError(cov_path)
        self._logger.debug(cov_path)
        cov = pd.read_csv(cov_path, sep="\t")
        merged = df.merge(cov, on=["#CHROM", "POS"], how="inner")
        
        # preferred (x) columns, TODO: fallback to (_y) where missing
        merged = merged.rename(columns={'DEPTH_x':'DEPTH','GENE_x': 'GENE','refseq_x': 'refseq','exone_x': 'exone', 'length_x': 'length', 'strand_x': 'strand', 'hgmd_x': 'hgmd'})

        # optional debug dump
        debug_out = cov_path.with_name(f"{sample_col}_all_macrotest")
        merged.to_csv(debug_out, sep="\t", index=False)
        self._logger.debug("(%s thread): Debug macro written → %s", self.sample.name, debug_out)

        # ----------------------------------------------------- unbalance string
        merged["unbalance"] = pd.Series(index=merged.index, dtype="object")

        # base-specific imbalance
        for b in "GCTA":
            pct = f"{b}%"
            mask = (merged["ALT"] == b) & (merged[pct] != 0)
            merged.loc[mask, "unbalance"] = (
                b + "=" + merged.loc[mask, pct].astype(str)
            )

        # indel imbalance
        for tag in ("ins", "del"):
            pct = f"{tag}%"
            if pct in merged.columns:
                mask = merged[pct] != 0
                merged.loc[mask, "unbalance"] = f"{tag}=" + merged.loc[mask, pct].astype(
                    str
                )

        # ----------------------------------------------------- contamination
        data = (
            merged.query("types == 'SVN' and unbalance != 'del=nan'")
            .assign(UNBAL=lambda d: d["unbalance"].str.split("=").str[1].astype(float))
        )

        if data.empty:
            return 0.0

        data["zscore"] = (data["UNBAL"] - hmean) / hstd
        #self._logger.debug(data["zscore"])
        

        # “extreme” by absolute threshold
        extremes = ((data["UNBAL"] < low_thresh) | (data["UNBAL"] > high_thresh)).sum()
        perc_extremes = 100.0 * extremes / len(data)

        # “extreme” by z-score
        extremes_z = ((data["zscore"] < -2) | (data["zscore"] > 2)).sum()
        perc_z = 100.0 * extremes_z / len(data)

        # The old pipeline used the z-score percentage
        probability = round(perc_z, 2)
        self._logger.debug(
            "(%s thread) [KinshipPipe.evaluate_contamination]: Contamination for %s → %.2f%% (|UNBAL|: %.2f %%, z-score: %.2f %%)",
            self.sample.name,
            sample_col,
            probability,
            perc_extremes,
            perc_z,
        )
        return probability



    # ------------------------------------------------------------------ #
    # Convenience wrappers around utils helpers
    # ------------------------------------------------------------------ #
