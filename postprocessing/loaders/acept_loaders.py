import pandas as pd
from sqlalchemy import create_engine
import logging
import os
from dotenv import load_dotenv
from pathlib import Path
from ...LogFormatter import SampleAlignedFormatter
import psycopg2
import sys

logger = logging.getLogger(__name__)


# __file__ is the location of this script; adjust if your .env lives elsewhere
env_path = Path(__file__).parent / "ppcfg.env"
load_dotenv(dotenv_path=env_path)

# builds a private handler for this module (variant_exporters) only
logger = logging.getLogger(__name__)    # e.g. bin_jurgen.postprocessing.exporters.variant_exporters

_priv_handler = logging.StreamHandler(sys.stdout)
_priv_handler.setFormatter(SampleAlignedFormatter())

logger.addHandler(_priv_handler)        # attach new handler
logger.setLevel(logging.DEBUG)          
logger.propagate = False                # <-   **key line**: don't let it reach the root handler

#conn = create_engine(os.getenv("DB_CONNECTION_STRING"))

conn = psycopg2.connect(
    dbname=os.getenv("dbname"), 
    user=os.getenv("user"), 
    password=os.getenv("password"), 
    host=os.getenv("host"), 
    port=os.getenv("port")
    )

def load_sample_status(sample_id, conn=conn):
    # check if sample is in 'Lavorazione', otherwise abort
    with conn.cursor() as cur:
        # ─── 0) verify the sample exists and is in “Lavorazione” ─────────────
        cur.execute(
            """
            SELECT referto
              FROM acept_sample
             WHERE sample = %s
            """,
            (sample_id,)
        )

        result = cur.fetchone()
        if result is None:
            logger.error(f"[Sample {sample_id}]: Sample not found in acept_sample!")
            return None
        referto, = result
        
    return referto