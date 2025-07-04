import csv
import numpy as np
import os, sys
from os import system
import psycopg2
from psycopg2.extras import execute_values
import math
from ..loaders.acept_loaders import load_sample_status

import pandas as pd
import glob
from pathlib import Path
import logging
from ...LogFormatter import SampleAlignedFormatter
import datetime
from ..pp_utils import sample_log_justifier

today_pre_date = datetime.datetime.now()
from datetime import datetime
today_date = datetime.date(datetime.now())


# -----------------------------------------------------------------------------
# LOGGER SETUP
# -----------------------------------------------------------------------------
# builds a private handler for this module (variant_exporters) only
logger = logging.getLogger(__name__)    # e.g. bin_jurgen.postprocessing.exporters.variant_exporters

_priv_handler = logging.StreamHandler(sys.stdout)
_priv_handler.setFormatter(SampleAlignedFormatter())

logger.addHandler(_priv_handler)        # attach new handler
logger.setLevel(logging.INFO)          
logger.propagate = False                # <-   **key line**: don't let it reach the root handler

#path_checklist = '/home/bioinfo/VIRTUAL38/apimagi_prod/NGS_RESULT/CHECKLIST/'
#file_list_checklist = glob.glob(path_checklist+'/*')


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



# -----------------------------------------------------------------------------
# CONFIGS and CONSTANTS
# -----------------------------------------------------------------------------

db_config_checklist = {
    "dbname": os.getenv("dbname"), 
    "user": os.getenv("user"), 
    "password": os.getenv("password"), 
    "host": os.getenv("host"), 
    "port": os.getenv("port")
}

from .constants import dict_BA1, dict_vusstatus_it, diction_cause, columns_criteri_cause_final, dict_choice_interpretation


# -----------------------------------------------------------------------------
# FUNCTIONS
# -----------------------------------------------------------------------------

def upsert_checklist(file: Path):


    sample_id = file.name.split('_')[0]
    referto = load_sample_status(sample_id, psycopg2.connect(**db_config_checklist)) 
    if referto != 'Lavorazione':
        logger.error(f"{sample_log_justifier(sample_id)}: Not Lavorazione: upsert_checklist aborted!")
        return

    df = pd.read_csv(file, sep='\t', dtype='object')
    df['BA1_intensita_final'] = (
        df['BA1_intensita_final'].replace('Stand-Alone', 'Stand Alone')
    )
    df['choice_interpretation'] = df['choice_interpretation'].replace(dict_choice_interpretation)
    df['BA1_intensita_final'] = df['BA1_intensita_final'].replace(dict_BA1)
    df['vusstatus'] = df['vusstatus'].replace(dict_vusstatus_it)

    # for col in columns_criteri_cause_final:
    #     # If the fill value is a string, cast the column to object first
    #     if isinstance(diction_cause[col], str):
    #         df[col] = df[col].astype('object')

    df[columns_criteri_cause_final] = df[columns_criteri_cause_final].fillna(diction_cause)

    #print(df.iloc[1:5,1:5])
    #logger.debug("Upserting checklist_unit for sample {}...".format(file.name.__str__()))
    _upsert_checklist_unit(df, db_config_checklist)
    logger.info(f"{sample_log_justifier(sample_id)}: upsert_checklist completed!")
    

def _upsert_checklist_unit(datareader, db_config):
    """
    Upsert all rows from your pandas DataFrame into the checklist table.
    - `datareader`: pd.DataFrame, cleaned and ready
    - `db_config`: dict with keys 'dbname','user','password','host','port'
    """

    # column names present in the file - should be the same as the ones in the db, if not: match them
    cols = [
        'sample_id','HGVS','BA1_final','GENE','PMPOT_final',
        'BA1_intensita_final','BA1_cause_final','BS1_final','BS1_intensita_final','BS1_cause_final',
        'BS2_final','BS2_intensita_final','BS2_cause_final','BS3_final','BS3_intensita_final',
        'BS3_cause_final','BS4_final','BS4_intensita_final','BS4_cause_final','BP1_final',
        'BP1_intensita_final','BP1_cause_final','BP2_final','BP2_intensita_final','BP2_cause_final',
        'BP3_final','BP3_intensita_final','BP3_cause_final','BP4_final','BP4_intensita_final',
        'BP4_cause_final','BP5_final','BP5_intensita_final','BP5_cause_final','BP6_final',
        'BP6_intensita_final','BP6_cause_final','BP7_final','BP7_intensita_final','BP7_cause_final',
        'PVS1_final','PVS1_intensita_final','PVS1_cause_final','PS1_final','PS1_intensita_final',
        'PS1_cause_final','PS2_final','PS2_intensita_final','PS2_cause_final','PS3_final',
        'PS3_intensita_final','PS3_cause_final','PS4_final','PS4_intensita_final','PS4_cause_final',
        'PM1_final','PM1_intensita_final','PM1_cause_final','PM2_final','PM2_intensita_final',
        'PM2_cause_final','PM3_final','PM3_intensita_final','PM3_cause_final','PM4_final',
        'PM4_intensita_final','PM4_cause_final','PM5_final','PM5_intensita_final','PM5_cause_final',
        'PM6_final','PM6_intensita_final','PM6_cause_final','PP1_final','PP1_intensita_final',
        'PP1_cause_final','PP2_final','PP2_intensita_final','PP2_cause_final','PP3_final',
        'PP3_intensita_final','PP3_cause_final','PP4_final','PP4_intensita_final','PP4_cause_final',
        'PP5_final','PP5_intensita_final','PP5_cause_final',
        'choice_interpretation','cadd_score','ada_score','rf_score','revel_score',
        'var_on_gene','allinheritance','zigosita','importanza','consequence',
        'choice_interpretation_final','vusstatus','sample_date'
    ]

    db_table_cols = [
        'sample_id','hgvs','\"BA1_final\"','gene','\"PMPOT_final\"',
        '\"BA1_intensita_final\"','\"BA1_cause_final\"','\"BS1_final\"','\"BS1_intensita_final\"','\"BS1_cause_final\"',
        '\"BS2_final\"','\"BS2_intensita_final\"','\"BS2_cause_final\"','\"BS3_final\"','\"BS3_intensita_final\"',
        '\"BS3_cause_final\"','\"BS4_final\"','\"BS4_intensita_final\"','\"BS4_cause_final\"','\"BP1_final\"',
        '\"BP1_intensita_final\"','\"BP1_cause_final\"','\"BP2_final\"','\"BP2_intensita_final\"','\"BP2_cause_final\"',
        '\"BP3_final\"','\"BP3_intensita_final\"','\"BP3_cause_final\"','\"BP4_final\"','\"BP4_intensita_final\"',
        '\"BP4_cause_final\"','\"BP5_final\"','\"BP5_intensita_final\"','\"BP5_cause_final\"','\"BP6_final\"',
        '\"BP6_intensita_final\"','\"BP6_cause_final\"','\"BP7_final\"','\"BP7_intensita_final\"','\"BP7_cause_final\"',
        '\"PVS1_final\"','\"PVS1_intensita_final\"','\"PVS1_cause_final\"','\"PS1_final\"','\"PS1_intensita_final\"',
        '\"PS1_cause_final\"','\"PS2_final\"','\"PS2_intensita_final\"','\"PS2_cause_final\"','\"PS3_final\"',
        '\"PS3_intensita_final\"','\"PS3_cause_final\"','\"PS4_final\"','\"PS4_intensita_final\"','\"PS4_cause_final\"',
        '\"PM1_final\"','\"PM1_intensita_final\"','\"PM1_cause_final\"','\"PM2_final\"','\"PM2_intensita_final\"',
        '\"PM2_cause_final\"','\"PM3_final\"','\"PM3_intensita_final\"','\"PM3_cause_final\"','\"PM4_final\"',
        '\"PM4_intensita_final\"','\"PM4_cause_final\"','\"PM5_final\"','\"PM5_intensita_final\"','\"PM5_cause_final\"',
        '\"PM6_final\"','\"PM6_intensita_final\"','\"PM6_cause_final\"','\"PP1_final\"','\"PP1_intensita_final\"',
        '\"PP1_cause_final\"','\"PP2_final\"','\"PP2_intensita_final\"','\"PP2_cause_final\"','\"PP3_final\"',
        '\"PP3_intensita_final\"','\"PP3_cause_final\"','\"PP4_final\"','\"PP4_intensita_final\"','\"PP4_cause_final\"',
        '\"PP5_final\"','\"PP5_intensita_final\"','\"PP5_cause_final\"',
        'choice_interpretation','cadd_score','ada_score','rf_score','revel_score',
        'var_on_gene','allinheritance','zigosita','importanza','consequence',
        'choice_interpretation_final','vusstatus','sample_date'
    ]

    today_date = datetime.date(datetime.now())

    # building the INSERT statement with placeholders and an ON CONFLICT upsert.
    col_list_sql = ', '.join(c for c in db_table_cols)  
    #placeholder_sql = ', '.join(['%s'] * len(cols))

    # Exclude the conflict keys from the update set
    update_cols   = [c for c in db_table_cols if c not in ('sample_id', 'hgvs')]
    update_sql    = ', '.join(f"{c} = EXCLUDED.{c}" for c in update_cols)

    insert_sql = f"""
    INSERT INTO checklist_resultinterpretation_new ({col_list_sql})
    VALUES %s
    ON CONFLICT (sample_id, hgvs)
      DO UPDATE SET {update_sql};
    """

    # gather rows
    data_tuples = [
        tuple(row[c] if c != 'sample_date' else today_date for c in cols)
        for _, row in datareader.iterrows()
    ]

    logger.debug("Length of data_tuples = {}".format(len(data_tuples)))
    logger.debug("Just before executing the insert statement...")
    # Bulk execution (execute_values usage)
    conn = psycopg2.connect(**db_config)
    try:
        with conn:
            with conn.cursor() as cur:
                # optional: template= for explicit casting, see below
                execute_values(cur, insert_sql, data_tuples, page_size=500)
        logger.debug("Upserted %s rows into checklist.", len(data_tuples))
    finally:
        conn.close()



    

def main(STAGING_DIR: Path):
    
    file_list_checklist = list(STAGING_DIR.glob('CHECKLIST/*'))
    for file in file_list_checklist:
        if file.name.endswith("interpretation_varsome.csv"):
            #logger.debug("Upserting checklist for sample {}...".format(file.name))
            upsert_checklist(file)
            print(''.join('-' for _ in range(200)))
    


if __name__ == "__main__":
    #main(Path("/home/magi/VIRTUAL/NGS_RESULTS_TEST/"))
    pass