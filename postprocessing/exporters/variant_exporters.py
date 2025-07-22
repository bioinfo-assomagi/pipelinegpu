
import pandas as pd
import psycopg2
from ..file_table_header_map import pheno_predict_map, other_annot_map, pindel_map
from psycopg2.extras import execute_values
from ...db.stored_procedures import reset_for_new_run as reset_for_new_run_proc, upsert_preselection as upsert_preselection_proc, insert_sanger_allvariation as insert_sanger_allvariation_proc
import csv
import math

import glob
import logging
import sys
from ..loaders.acept_loaders import load_sample_status
from ...LogFormatter import SampleAlignedFormatter
from pathlib import Path
from ..pp_utils import sample_log_justifier
import os
import datetime

# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------

RESULTS_PATH = '/home/magi/PROJECT/VIRTUAL/...'

conn = psycopg2.connect(
    dbname=os.getenv("dbname"), 
    user=os.getenv("user"), 
    password=os.getenv("password"), 
    host=os.getenv("host"), 
    port=os.getenv("port")
    )


# builds a private handler for this module (variant_exporters) only
logger = logging.getLogger(__name__)    # e.g. bin_jurgen.postprocessing.exporters.variant_exporters

_priv_handler = logging.StreamHandler(sys.stdout)
_priv_handler.setFormatter(SampleAlignedFormatter())

logger.addHandler(_priv_handler)        # attach new handler
logger.setLevel(logging.DEBUG)          
logger.propagate = False                # <-   **key line**: don't let it reach the root handler


# ------------------------------
# HELPERS
# ------------------------------

def as_text(v):
    """Return None for NaN/None, else a str."""
    if pd.isna(v):
        return None
    return str(v)

def none_if_nan(x):
    """Convert NaN / NaT → None so psycopg2 will send NULL."""
    if (x is None
        or (isinstance(x, float) and math.isnan(x))
        or (isinstance(x, str) and x.lower() in ("nan", "nat"))):
        return None
    return x


# ----------------------------------------
# MAIN FUNCTIONS
# ----------------------------------------

# The following functions will run in order.


def upsert_variant(table_name: str, variant_df_filepath: Path, col_map: dict):
    sample_id = variant_df_filepath.name.split("_")[0]
    referto = load_sample_status(sample_id, conn)
    if referto != 'Lavorazione':
        logger.error(f"{sample_log_justifier(sample_id)}: Not Lavorazione: upsert_variant {table_name} aborted!")
        return


    with conn.cursor() as cur:
        delete_sql = """ DELETE FROM {table_name} WHERE sample_id = '{sample_id}' """.format(table_name=table_name, sample_id=sample_id)
        cur.execute(delete_sql)
        #logger.info(f"  – removed {cur.rowcount} old rows for sample {sample_id} in table {table_name}")

    
    TABLE_COLS = list(col_map.values()) 
    HEADER_TO_IDX = {h: i for i, h in enumerate(col_map.keys())}
    
    insert_sql = """
        INSERT INTO {table_name} ({cols})
        VALUES ({placeholders})
        ON CONFLICT (sample_id, hgvs) DO UPDATE
        SET {updates}
        """.format(
            table_name   = table_name,
            cols         = ", ".join(TABLE_COLS),
            placeholders = ", ".join(["%s"] * len(TABLE_COLS)),
            updates      = ", ".join([
                f"{col} = EXCLUDED.{col}"
                for col in TABLE_COLS
                if col not in ("sample_id", "hgvs")   # don’t update PK
            ])
        )

    with conn.cursor() as cur, open(variant_df_filepath) as fp:
        reader = csv.DictReader(fp, delimiter="\t")
        cnt = 0
        for raw in reader:
            row = [ raw[k] for k in col_map.keys() ]   # reorder
            cur.execute(insert_sql, row)
            cnt += 1

    conn.commit()

    logger.info(f"{sample_log_justifier(sample_id)}: upsert_variant {table_name} completed!")
 


def reset_for_new_run(df_filepath: str):
    sample_id = df_filepath.name.split("_")[0]
    referto = load_sample_status(sample_id, conn)
    if referto != 'Lavorazione':
        logger.error(f"{sample_log_justifier(sample_id)}: Not Lavorazione: reset_for_new_run aborted!")
        return

    reset_for_new_run_proc(sample_id, conn)
    
    logger.info(f"{sample_log_justifier(sample_id)}: reset_for_new_run completed!")
    


def upsert_preselection(df_filepath: str):
    sample_id = df_filepath.name.split("_")[0]
    referto = load_sample_status(sample_id, conn)
    if referto != 'Lavorazione':
        logger.error(f"{sample_log_justifier(sample_id)}: Not Lavorazione: upsert_preselection aborted!")
        return
    
    df = pd.read_csv(df_filepath, sep="\t")
    for row in df.iterrows():
        upsert_preselection_proc(row, conn)        
    
    logger.info(f"{sample_log_justifier(sample_id)}: upsert_preselection completed!")

def upsert_sanger_allvariation(df_filepath : str):
    sample_id = df_filepath.name.split("_")[0]
    referto = load_sample_status(sample_id, conn)
    if referto != 'Lavorazione':
        logger.error(f"{sample_log_justifier(sample_id)}: Not Lavorazione: upsert_sanger_allvariation aborted!")
        return False
    df = pd.read_csv(df_filepath, sep ="\t", header=0)
    for i, row in df.iterrows():
        insert_sanger_allvariation_proc(row, conn)
        
    logger.info(f"{sample_log_justifier(sample_id)}: upsert_sanger_allvariation completed!")
    return True


def update_sample_timestamp(sample_id):
    with conn.cursor() as cur:
        cur.execute("UPDATE acept_sample SET analisi_bioinfo = %s WHERE sample = %s", (datetime.datetime.now(), sample_id))
    conn.commit()
    logger.info(f"{sample_log_justifier(sample_id)}: update_sample_timestamp completed!")

def bulk_upsert_genestat(df_filepath : Path):
    """
    Delete old rows for *sample_name* and bulk-insert *df* into genestat_cov.
    """
    sample_name = df_filepath.name.split("_")[0]
    sample_id = sample_name
    TABLE='ngs_genestatcov_new'

    referto = load_sample_status(sample_id, conn)
    if referto != 'Lavorazione':
        logger.error(f"{sample_log_justifier(sample_id)}: Not Lavorazione: bulk_upsert_genestat aborted!")
        return False


    with conn.cursor() as cur:

        cur.execute(f"DELETE FROM {TABLE} WHERE sample_id = %s", (sample_name,))
        #logger.info(f"  – removed {cur.rowcount} old rows for sample {sample_name} in table {TABLE}")

        cols = [
            'sample', 'GENE', 'strand', 'refseq', 'hgmd', 'exone', 'CHROM',
            'Exon_Start', 'Exon_End', 'start_buco', 'end_buco',
            'len_gene', 'len_exon', 'buco_exone', 'note', 'buco_gene',
            '%for_exon', '%for_gene', 'mean', 'variante', #'linea'
        ]

        table_cols = [
            'sample_id', '\"GENE\"', 'strand', 'refseq', 'hgmd', 'exone', 'chrm',
            'exone_start', 'exone_end', 'start_buco', 'end_buco',
            'len_gene', 'len_exone', 'buco_exone', 'note', 'buco_gene',
            'perc_exone', 'perc_gene', 'cov_mean', 'variante', #'linea'
        ]

        df = pd.read_csv(df_filepath, sep ="\t", header=0)
        values = [
            tuple(none_if_nan(row[c]) for c in cols)
            for _, row in df[cols].iterrows()
        ]

        col_list = ", ".join(table_cols)
        execute_values(
            cur,
            f"INSERT INTO {TABLE} ({col_list}) VALUES %s",
            values
        )
    conn.commit()
    #logger.info(f"  – inserted {len(values)} rows for {sample_name} in table {TABLE}\n")
    logger.info(f"{sample_log_justifier(sample_id)}: bulk_upsert_genestat completed!")
    return True

def bulk_upsert_statcov(df_filepath : Path):

    sample_name = df_filepath.name.split("_")[0]
    sample_id = sample_name
    STATCOV_TABLE='ngs_statcov_new'

    referto = load_sample_status(sample_id, conn)
    if referto != 'Lavorazione':
        logger.error(f"{sample_log_justifier(sample_id)}: Not Lavorazione: bulk_upsert_statcov aborted!")
        return False

    df = pd.read_csv(df_filepath, sep="\t", dtype={"sample": str})

    # normalise sample id
    df["sample"] = (
        df["sample"].str.replace(r"\.202$", ".2020", regex=True).astype(str)
    )

    # add somatic columns if they’re missing
    som_cols = ["cutoff_somatic", "count_somatic",
                "perc_somatic", "cov_medio_somatic"]
    for c in som_cols:
        if c not in df.columns:
            df[c] = -999

    """DELETE old rows for sample, then bulk-insert df into STATCOV_TABLE."""
    cols = [
        "sample", "cutoff", "count", "perc", "cov_medio", "sesso",
        "compatibilita_pat", "cutoff_somatic", "compatibilita_mat",
        "count_somatic", "probcontaminazione", "perc_somatic",
        "cov_medio_somatic",
    ]

    table_cols = [
        "sample_id", "cutoff", "count", "perc", "cov_medio", "sesso",
        "compatibilita_pat", "cutoff_somatic", "compatibilita_mat",
        "count_somatic", "probcontaminazione", "perc_somatic",
        "cov_medio_somatic",
    ]
    
    tuples = [tuple(none_if_nan(r[c]) for c in cols) for _, r in df[cols].iterrows()]
    cols_sql = ", ".join(table_cols)
    with conn.cursor() as cur:
        cur.execute(f"DELETE FROM {STATCOV_TABLE} WHERE sample_id = %s", (sample_name,))
        logger.info(f"{sample_log_justifier(sample_id)}: deleted %s old rows for sample %s", cur.rowcount, sample_name)

        execute_values(cur,
            f"INSERT INTO {STATCOV_TABLE} ({cols_sql}) VALUES %s", tuples
        )
    conn.commit()
    logger.info(f"{sample_log_justifier(sample_id)}: bulk_upsert_statcov completed!")
    return True


def main(STAGING_DIR: Path):
    
    file_list_variant = list(STAGING_DIR.glob('annot/*'))
    file_list_indel = list(STAGING_DIR.glob('INDEL/*'))
    file_list_genestat = list(STAGING_DIR.glob('coverage/*'))

    # splitting just for log clarity
    print(''.join('-' for _ in range(220)))
    #print('------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------')
    # Let's group files by samples first
    files_by_sample = {}
    for file in file_list_variant:
        sample_id = file.name.split("_")[0]
        if sample_id not in files_by_sample:
            files_by_sample[sample_id] = []
        files_by_sample[sample_id].append(file)
    
    for file in file_list_indel:
        sample_id = file.name.split("_")[0]
        if sample_id not in files_by_sample:
            files_by_sample[sample_id] = []
        files_by_sample[sample_id].append(file)
    
    for file in file_list_genestat:
        sample_id = file.name.split("_")[0]
        if sample_id not in files_by_sample:
            files_by_sample[sample_id] = []
        files_by_sample[sample_id].append(file)


    for sample_id, files in files_by_sample.items():
        # again, redundant loops just for log clarity and order of calls
        for file in files:
            if file.name.endswith("_pheno_predict.csv"):
                upsert_variant("ngs_ngsvariants_new", file, pheno_predict_map)

        for file in files:
            if file.name.endswith("_other_annot.csv"):
                upsert_variant("ngs_ngsotherannotations_new", file, other_annot_map)
        
        for file in files:
            # print(file)
            if file.name.endswith("_final_indel.csv"):
                upsert_variant("ngs_indelvariants_new", file, pindel_map)

        for file in files:
            if file.name.endswith("_variant_selection.csv"):
                reset_for_new_run(file)
                upsert_preselection(file)
                sanger_status = upsert_sanger_allvariation(file)
                if sanger_status:
                    update_sample_timestamp(sample_id)
                
                

        for file in files:
            if file.name.endswith("_gene_cov.csv"):
                bulk_upsert_genestat(file)

        for file in files:
            if file.name.endswith("_stat_cov.csv"):
                bulk_upsert_statcov(file)
                
        print(''.join('-' for _ in range(220)))
        #print('------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------')
       



def run():
    pass

def test(STAGING_DIR: Path):

    file_list_variant = list(STAGING_DIR.glob('annot/*'))
    for file in file_list_variant:
        try:
            if file.name.endswith("_pheno_predict.csv"):
                print("ngs_ngsvariants_new", file, "pheno_predict_map")
            elif file.name.endswith("_other_annot.csv"):
                print("ngs_ngsotherannotations_new", file, "other_annot_map")
        except Exception as e:
            print(str(e))

    # upsert_variant("ngs_ngsvariants_new", "/home/magi/PROJECT/diagnosys/bin_jurgen/tests/135.2025/final/135.2025_pheno_predict.csv", pheno_predict_map)
    # upsert_variant("ngs_ngsotherannotations_new", "/home/magi/PROJECT/diagnosys/bin_jurgen/tests/135.2025/final/135.2025_other_annot.csv", other_annot_map)
    
    file_list_indel = list(STAGING_DIR.glob('INDEL/*'))
    for file in file_list_indel:
        if file.name.endswith("_final_indel.csv"):
            print("ngs_indelvariants_new", file, "pindel_map")

    for file in file_list_variant:
        if file.name.endswith("_variant_selection.csv"):
            print("reset_for_new_run", file)
            print("upsert_preselection", file)
            print("upsert_sanger_allvariation", file)

    file_list_genestat = list(STAGING_DIR.glob('coverage/*'))
    for file in file_list_genestat:
        if file.name.endswith("_gene_cov.csv"):
            print("bulk_upsert_genestat", file)

if __name__ == "__main__":
    main()