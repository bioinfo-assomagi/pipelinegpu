import pandas as pd
import psycopg2
from psycopg2.extras import RealDictCursor 
from sqlalchemy import create_engine
import os


def load_campioni_accettazione():

    # conn = psycopg2.connect(
    #         dbname="postgres", 
    #         user="bioinfo", 
    #         password="password_sicura", 
    #         host="192.168.1.143", 
    #         port="5432"
    #     )

    engine = create_engine(
        "postgresql+psycopg2://bioinfo:password_sicura@192.168.1.143:5432/postgres"
    )

    # engine = create_engine(
    #     create_engine(os.getenv("DB_CONNECTION_STRING"))
    # )

    SQL = """
            SELECT
                a.sample_id::text                                AS sample_id,          -- cast once, not in pandas
                a.riferimento as riferimento_id,
                a.familiarita,
                a.familiarita_relativa,
                a.campionebiologico1 as campionebiologico,
                a.cognome,
                a.nome,
                a.sesso,
                a.nazione,
                a.comunenascita as luogonascita,
                a.datanascita,
                a.codice_fiscale,

                u_gen.username                                    AS operatoregenetista,
                u_san.username                                    AS operatoresanger,

                c.clinico                                         AS clinico,
                c.ospedale                                        AS istituto,

                a.panelid,
                a.malattia_text as fenotipo,
                a.operatoregenetista_id      -- keep the raw id if you still need it
                    AS operatoregenetista_id,
                a.affetto,
                a.panel AS pannello,
                a.select_nucleofamiliare,
                a.select_nucleosenzaprobando::bool,
                a.select_nucleoprobtwist,
                a.somatico,

                s.strutturainviante                               AS strutturainviante
            FROM
                access_accettazione                     a
            LEFT JOIN auth_user                  u_gen  ON u_gen.id = a.operatoregenetista_id::integer
            LEFT JOIN auth_user                  u_san  ON u_san.id = a.operatoresanger_id::integer
            LEFT JOIN access_clinico                    c      ON c.id      = a.clinicoaccettazione_id::integer
            LEFT JOIN access_strutturainviante         s      ON s.id      = a.strutturainvianteaccettazione_id::integer
            ORDER BY a.sample_id;
        """

    campioni_accetazione_df = pd.read_sql(SQL, engine)
    campioni_accetazione_df.columns = (
        col.lower() for col in campioni_accetazione_df.columns
    )

    engine.dispose()

    return campioni_accetazione_df




def load_acc():

    """
    sample_id | timestamp_accettazione            | select_nucleosenzaprobando | affetto_flag | riferimento | familiarita | familiarita_relativa | nome       | cognome
    285.2025  | 2025-02-27 11:04:25.928517 +00:00 | false                      | NO           | 284.2025    | MADRE       | MADRE                | ANTONINA   | LANDO
    286.2025  | 2025-02-27 11:19:02.217091 +00:00 | false                      | SI           | 286.2025    | PROBANDO    | PROBANDO             | MARIA      | CIRIANNI
    287.2025  | 2025-02-27 11:27:59.194686 +00:00 | false                      | SI           | 287.2025    | PROBANDO    | PROBANDO             | MARCO      | MARCON
    290.2025  | 2025-03-07 11:10:26.921312 +00:00 | true                       | NO           | 290.2025    | SORELLA     | SORELLA              | GEMMA      | CICCIOLI
    291.2025  | 2025-03-07 11:23:17.051104 +00:00 | true                       | NO           | 291.2025    | SORELLA     | SORELLA              | MARINA     | PORTULANO
    292.2025  | 2025-03-07 11:51:54.468829 +00:00 | true                       | NO           | 537.2022    | MADRE       | MADRE                | GIULIA     | DESSI
    295.2025  | 2025-03-17 13:18:18.595987 +00:00 | false                      | SI           | 295.2025    | PROBANDO    | PROBANDO             | ANNA       | TRANI
    296.2025  | 2025-03-18 09:20:41.040149 +00:00 | false                      | SI           | 296.2025    | PROBANDO    | PROBANDO             | SIMONETTA  | COMINCILIO
    
    """

    engine = create_engine(
        "postgresql+psycopg2://bioinfo:password_sicura@192.168.1.143:5432/postgres"
    )

    # pre-load every field we will need from ACCETTAZIONE **once**
    ACC_SQL = """
    SELECT sample_id::text,
        timestamp_accettazione as data_creazione,
        select_nucleosenzaprobando::bool,
        affetto                       AS affetto_flag,          -- avoid name clash
        riferimento,
        familiarita,
        familiarita_relativa,
        nome,
        cognome,
        datanascita,
        sesso,
        panel,
        select_nucleofamiliare
    FROM   access_accettazione;
    """
    acc_df = pd.read_sql(ACC_SQL, engine)
    acc_df.columns = (
        col.lower() for col in acc_df.columns
    )
    engine.dispose()

    return acc_df



def load_fam_acc():

    """
    # +-------------+---------------------------------------------------------------------------------------+
    # | riferimento | famiglianucleo                                                                        |
    # +-------------+---------------------------------------------------------------------------------------+
    # |   100.2023  | 101.2023: MADRE: PIPERNI: LORETA: 1959-05-22: F: NO,102.2023: PADRE: STORNELLI:...    |
    # |  1004.2024  | 1005.2024: MADRE: SPERANZA: MARIA CARMELA: 1962-09-12: F: NO                          |
    # |  1006.2023  | 1056.2023: MADRE: GEAN: RACHELE: 1957-01-16: F: NO,1132.2023: SORELLA: REBUZZI:...    |
    # |  1006.2024  | 1127.2024: MADRE: PUOTI: LILIANA: 1967-11-03: F: SI                                   |
    # |  1007.2023  | 1008.2023: FIGLIA: BARBARESI: GIORGIA: 1987-01-11: F: SI                              |
    # |  1008.2024  | 1009.2024: SORELLA: MUSCATELLO: ANNA: 1984-11-11: F: NO                               |
    # |  1009.2023  | 1009.2023: SORELLA: PELLICCIA: GEMMA: 1970-11-06: F: NO                               |
    # +-------------+---------------------------------------------------------------------------------------+

    """


    engine = create_engine(
        "postgresql+psycopg2://bioinfo:password_sicura@192.168.1.143:5432/postgres"
    )
    
    FAM_SQL = """
        SELECT riferimento as riferimento_id,
            STRING_AGG(
                sample_id        || ':' ||
                COALESCE(NULLIF(familiarita_relativa,''), familiarita) || ':' ||
                cognome          || ':' ||
                nome             || ':' ||
                datanascita::text|| ':' ||
                sesso            || ':' ||
                COALESCE(NULLIF(affetto,''),'NO')
            , ',')              AS famiglianucleo
        FROM   access_accettazione
        WHERE  select_nucleofamiliare::boolean
        AND NOT (familiarita = 'PROBANDO' AND riferimento = sample_id)
        GROUP  BY riferimento;
    """
    fam_df = pd.read_sql(FAM_SQL, engine)          # one row per “riferimento”
    engine.dispose()

    return fam_df