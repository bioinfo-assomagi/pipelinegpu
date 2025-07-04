import pandas as pd
import psycopg2
from psycopg2.extras import RealDictCursor 
from psycopg2.extras import execute_values
from psycopg2 import DatabaseError
import logging
import sys
from ...LogFormatter import SampleAlignedFormatter

from pathlib import Path
import math

from ..loaders.access_loaders import load_campioni_accettazione, load_acc, load_fam_acc
from ..loaders.acept_loaders import load_sample_status

import glob


# builds a private handler for this module (variant_exporters) only
logger = logging.getLogger(__name__)    # e.g. bin_jurgen.postprocessing.exporters.variant_exporters

_priv_handler = logging.StreamHandler(sys.stdout)
_priv_handler.setFormatter(SampleAlignedFormatter())

logger.addHandler(_priv_handler)        # attach new handler
logger.setLevel(logging.DEBUG)          
logger.propagate = False                # <-   **key line**: don't let it reach the root handler

principal_directory = "dummy"
#path_refertazioneautomatica = '/home/bioinfo/VIRTUAL/MAGIS/NGS_RESULT/REFERTAZIONE/'
#path_refertazioneautomatica = '/home/magi/PROJECT/diagnosys/bin_jurgen/tests/'
#file_list_refertazioneautomatica = glob.glob(path_refertazioneautomatica+'/*')


METODI_FEB2023 = 'L\'estrazione e purificazione di DNA genomico viene eseguita in modalita\' automatica \
                    (ZINEXTS MagPurix 12 System) o manuale con kit specifico a seconda della matrice biologica fornita \
                    (sangue periferico, saliva o altro tessuto). La preparazione di librerie di DNA genomico, la cattura \
                    e l\'arricchimento delle regioni target [regioni codificanti, giunzioni introne-esone dei geni inclusi \
                    nel pannello e regioni non codificanti specifiche (vedi sezione "geni nel sospetto")] avviene mediante \
                    Twist Custom Panel EF Workflow (Twist Bioscience). Il sequenziamento massivo in parallelo delle regioni \
                    target viene eseguito su piattaforma Illumina (sequencing by synthesis, SBS) in modalita\' paired-end \
                    con lunghezza delle reads di 150bp. L\'eventuale sequenziamento Sanger su seconda estrazione di DNA per \
                    la conferma di varianti che non soddisfano i parametri di qualita\' validati internamente viene effettuato \
                    in conformita\' con la procedura interna (PMID:31599443). Durante la lavorazione vengono \
                    eseguiti controlli interni di qualita\' per l\'esclusione di contaminazione e scambio dei campioni biologici \
                    (PMID:33170178).'

ANALISI_FEB2023 = 'I dati grezzi di sequenziamento (Fastq file) vengono filtrati in base alla qualita\' di ciascuna read \
                    (Phred-Score) Le reads di buona qualita\' vengono allineate alla sequenza di riferimento umana GRCh38 (Hg38) \
                    utilizzando il software BWA. Dalle reads allineate e ordinate con software SAMTOOLS vengono rimossi i \
                    duplicati e quindi riallineate e ricalibrate localmente intorno agli indel e alle varianti puntiformi \
                    note utilizzando il software GATK. Il VCF (Variant Call Format) viene generato da due distinti software: \
                    SAMTOOLS (Bcftools) e GATK (HaplotypeCaller). Vengono selezionate le varianti nelle regioni codificanti \
                    fino a +/- 15 bp e in regioni non codificanti specifiche. Tali varianti vengono annotate con il software \
                    VEP utilizzando dbSNP di RefSeq ed il database APPRIS per determinare il trascritto principale del gene \
                    interessato (PMID: 31378919). Le varianti vengono filtrate per frequenza allelica (Minor Allele Frequency \
                    <3%, con alcune eccezioni definite internamente) e interpretate mediante un algoritmo bioinformatico \
                    (PMID:34946832) che prevede l\'assegnazione di criteri di interpretazione sulla base delle informazioni \
                    fornite dalla piattaforma VarSome tramite l\'ambiente Stable-API con alcune modifiche: (1) utilizzo del \
                    tool AutoPVS1 (PMID:32442321) per l\' assegnazione del criterio PVS1 per varianti con predetta perdita \
                    di funzione; (2) assegnazione dei criteri PP3 e BP4 sulla base dei predittori REVEL (PMID: 27666373), \
                    CADD (PMID: 24487276) e di altri 15 predittori funzionali per le varianti missenso; dei predittori \
                    ADA Boost splicing e Random Forest (RF) splicing per varianti introniche in siti non canonici di splicing, \
                    missenso e sinonime in siti di regolazione dello splicing); (3) ripristino della forza di default \
                    (Moderate) del criterio PM2 in caso di diversa assegnazione. Sono previste modifiche ai criteri assegnati \
                    tramite VarSome sulla base di procedure validate internamente secondo le linee guida ACMG (PMID: 25741868) \
                    e in accordo con le raccomandazioni ACGS (Ellard S. et al. 2020). \
                    L\'analisi delle CNVs viene eseguita utilizzando il tool CoNVaDING (PMID: 26864275) secondo quanto descritto \
                    in PMID:34890029. L\'interpretazione delle CNVs viene eseguita sulla base di informazioni  fornite dall\'ambiente \
                    Stable-API di VarSome. Sono previste modifiche in accordo con le raccomandazioni ACMG e ClinGen (PMID: 31690835).'

SPECIFICHE_OCT2023 = 'In accordo con le raccomandazioni ACGS (Ellard S. et al. 2020), le VUS del tipo SNVs/small indels \
                    vengono distinte in tre classi (calde, tiepide, fredde) sulla base della loro possibile rilevanza \
                    clinica (PMID:37628650). Non vengono riportate: SNVs/small indels interpretate come B/LB/VUS fredde/VUS tiepide \
                    (non in grado di assumere lo status di VUS calda in seguito all\'eventuale attribuzione di criteri \
                    di patogenicita\' potenziali relativi a storia familiare, stato allelico o specificita\' del fenotipo) \
                    e VUS calde in geni associati a patologie a trasmissione recessiva o digenica in assenza di seconda \
                    variante di rilevanza clinica nello stesso gene o nel gene coinvolto in digenia. Non vengono riportate \
                    CNV interpretate come B/LB; CNV interpretate come VUS in geni associati a patologie a trasmissione \
                    autosomica recessiva o digenica in assenza di altra variante di rilevanza clinica nello stesso gene \
                    o nel gene coinvolto in digenia. L\'interpretazione delle varianti avviene sulla base dei dati \
                    disponibili alla data dell\'analisi e pertanto puo\' cambiare nel tempo in seguito all\'acquisizione \
                    di nuove conoscenze e all\'aggiornamento dei dati di consultazione. Acronimi: AR=autosomico recessiva; \
                    AD=autosomico dominante; XLR=X linked recessiva; XLD=X linked dominante; YL=Y linked; DIG=digenica; \
                    Mu=multifattoriale; SNVs=Varianti a singolo nucleotide; Small Indels= piccole inserzioni/delezioni; \
                    CNVs=Variazioni del numero di copie; P=Patogenetica; LP=Probabilmente Patogenetica; VUS=Variante di \
                    Significato Incerto; LB=Probabilmente Benigna; B=Benigna.'

SPECIFICHE_NOPROB_OCT2023 = 'La non contestuale analisi del probando potrebbe limitare l\'utilita\' diagnostica del test. Pertanto \
                    si fa presente che nei familiari riferiti clinicamente non affetti, inclusi nel nucleo familiare, viene \
                    eseguita esclusivamente la ricerca delle SNVs/small indels/CNVs individuate precedentemente nel probando \
                    affetto per le quali e\' stata richiesta analisi di segregazione e delle eventuali varianti aggiuntive \
                    individuate in altri membri consanguinei affetti del nucleo familiare. La segregazione di tali varianti \
                    aggiuntive non puo\' essere eseguita nel probando. \
                    Vengono riportate anche eventuali varianti P/LP identificate nel probando di riferimento, ma non causative \
                    del fenotipo clinico (stato di portatore). In accordo con le raccomandazioni ACGS (Ellard S. et al. 2020), \
                    le VUS del tipo SNVs/small indels vengono ulteriormente distinte in tre classi (calde, tiepide, fredde) \
                    sulla base della loro possibile rilevanza clinica (PMID:37628650). Non vengono riportate nella tabella \
                    dei risultati varianti assenti nel soggetto analizzato oppure varianti precedentemente identificate \
                    nel probando riclassificate come B/LB/VUS fredde. L\'interpretazione delle varianti avviene sulla base \
                    dei dati disponibili alla data dell\'analisi e pertanto puo\' cambiare nel tempo in seguito \
                    all\'acquisizione di nuove conoscenze e all\'aggiornamento dei dati di consultazione. Acronimi: \
                    AR=autosomico recessiva; AD=autosomico dominante; XLR=X linked recessiva; XLD=X linked dominante; \
                    YL=Y linked; DIG=digenica; Mu=multifattoriale; SNVs=Varianti a singolo nucleotide; Small Indels=inserzioni \
                    e delezioni; CNVs=Variazioni del numero di copie; P=Patogenetica; LP=Probabilmente Patogenetica; \
                    VUS=Variante di Significato Incerto; LB=Probabilmente Benigna; B=Benigna.'

LIMITI_OCT2023 = 'Il test non permette la rilevazione di: varianti al di fuori della sequenza codificante +/- 15bp \
                se non diversamente indicato (vedere sezione "Geni analizzati"), SNVs e small indels in regioni \
                con coverage inferiori a 10X, ripetizioni di triplette, mosaicismi della linea germinale, \
                delezioni mono- e biesoniche in eterozigosi e duplicazioni fino a 5 esoni nei cromosomi autosomici, \
                delezioni in eterozigosi e duplicazioni nei cromosomi sessuali ed altri riarrangiamenti strutturali \
                complessi. ll test non garantisce la rilevazione di: small indels superiori ai 20 bp, delezioni mono- e \
                biesoniche in omozigosi e in emizigosi sul cromosoma X. Il presente test non permette di definire \
                l\'esatta estensione ed i punti di rottura delle eventuali CNVs rilevate. Per le duplicazioni non e\' \
                possibile stabilire l\'esatto numero di copie e la loro localizzazione nel genoma. Pertanto, se ne \
                consiglia la conferma mediante una metodica molecolare quantitativa specifica, possibilmente coinvolgendo \
                gli esoni fiancheggianti.'

TWIST_99586838 = '99.7% [CI (99.35-99.89)]'
TWIST_97871383 = '99.51% [CI (99.00-99.8)]'
TWIST_91957849 = '99.83% [CI (99.08-99.99)]'
TWIST_97673193 = '99.56% [CI (99.1-99.82)]'
TWIST_96246486 = '99.79% [CI (99.25-99.97)]'
TWIST_97436393 = '98.88% [CI (96.0-99.86)]'
TWIST_91441162 = '99.81% [CI (98.93-99.99)]'
TWIST_98071606 = '99.72% [CI (99.43-99.89)]'


def test(STAGING_DIR: Path):

    refertazione_dir = STAGING_DIR / "REFERTAZIONE"
    file_list_refertazioneautomatica = refertazione_dir.glob("*")
    file_list_refertazioneautomatica = list(file_list_refertazioneautomatica)
    #print(file_list_refertazioneautomatica)

def main(STAGING_DIR: Path):

    refertazione_dir = STAGING_DIR / "REFERTAZIONE"
    file_list_refertazioneautomatica = refertazione_dir.glob("*")
    file_list_refertazioneautomatica = list(file_list_refertazioneautomatica)

    campioni_accettazione_df = load_campioni_accettazione()
    acc_df = load_acc()
    fam_df = load_fam_acc()

    for path in file_list_refertazioneautomatica:
        p = Path(path)
        parts = p.name.split('_', 2)
        
        # process only “…/refertazioneit_*.tsv”
        if parts[0] != "refertazioneit":
            continue

        logger.info("Processing file: {}".format(path.name))

        campioni = pd.read_csv(path,
                            sep='\t',
                            dtype={'sample_id': str})
        campioni['varianti_principali'] = campioni['varianti_principali'].astype(str)
        #print(campioni_accettazione_df)


        campioni = campioni.merge(
            acc_df[['sample_id', 'data_creazione']],
            on='sample_id', how='left'
        )

        # ── 3. keep only samples that exist in both tables
        samples = (campioni_accettazione_df        # ← produced in the previous step
                .merge(campioni, on='sample_id', how='inner'))

        #print(samples)
        # ── 4. add constant columns
        samples = samples.assign(
            gestione                    = 'DIAGNOSTICA',
            referto                     = 'Preanalisi',
            analisi_eseguite            = ANALISI_FEB2023,
            metodi_utilizzati           = METODI_FEB2023,
            limiti_ngs                  = LIMITI_OCT2023,
            specifiche_refertazione     = SPECIFICHE_OCT2023
        )

        # ── 5. tidy strings / NaNs
        samples['conclusioni_secondarie'] = samples['conclusioni_secondarie'].fillna('')
        samples[['operatoregenetista', 'operatoresanger']] = (
            samples[['operatoregenetista', 'operatoresanger']].fillna('')
        )

        samples['luogonascita'] = (
            samples['luogonascita'].replace('', pd.NA)
                                    .fillna(samples['nazione'])
        )

        # ── 6. nr_varianti_principali  (vectorised)
        samples['nr_varianti_principali'] = (
            samples['varianti_principali']
                .fillna('')
                .apply(lambda s: 0 if s in ('', ',') else len(s.split(',')) - 1)
        )

        # ── 7. sensibilità analitica  (dict-based map)
        sens_map = {
            'INFERTILITY':   TWIST_98071606,
            'MIXED1':        TWIST_98071606,
            'OCULARE':       TWIST_91441162,
            'CANCER':        TWIST_97436393,
            'VASCULAR':      TWIST_98071606,
            'NEUROLOGY':     TWIST_98071606,
            'LYMPHOBESITY':  TWIST_98071606
        }
        samples['sensibilita_analitica'] = (
            samples['pannello'].map(sens_map).fillna('-999')
        )

        # ── 8. override “specifiche_refertazione”
        #     (join the two boolean flags already cached in acc_df, then vectorise)
        flag_cols = acc_df[['sample_id',
                            'affetto_flag']]
        samples = samples.merge(flag_cols, on='sample_id', how='left')

        no_prob   = samples['select_nucleosenzaprobando']
        affetto   = samples['affetto_flag'] == 'SI'

        samples.loc[no_prob & ~affetto,
                    'specifiche_refertazione'] = SPECIFICHE_NOPROB_OCT2023


        # ── append the pre-built “famiglianucleo” text
        samples = samples.merge(fam_df, on='riferimento_id', how='left')

        # ── optional clean-up of empty variant strings
        samples.loc[samples['varianti_principali'].str.len() == 1,
                    'varianti_principali'] = ''
        samples.loc[samples['altre_varianti'].str.len() == 1,
                    'altre_varianti']      = ''


        samples.to_csv('/home/magi/PROJECT/diagnosys/bin_jurgen/tests/refertazion_exporters_sample.csv', sep="\t", index=False)
        
        """ Here is the place where the Lavorazione check will be done for each sample; if that sample is not in Lavorazione, the push_to_db function will not be called. """
        #print(samples)
        for idx, sample in samples.iterrows():
            # row is a Series; wrap it in a list to make a 1-row DataFrame
            sample_df = pd.DataFrame([sample])
            #print(sample_df)

            referto = load_sample_status(str(sample['sample_id']))
            if referto == 'Lavorazione':
                push_to_db(sample_df)
            else:
                logger.error(f"[Sample {sample['sample_id']}]: Not Lavorazione: push_to_db aborted!")
            
       


def push_to_db(sample_df):
    """ Note that the input is a DataFrame with a single row. 
    
    The function will update the refertazione_refertazione table for the given sample. And the plural row
    references are just dummy statements. 
    In general, we expect one referto for one sample.
    """

    UPSERT_SQL = None

    
    # -------------------------------------------------------- #
    # Mapping of the columns of the DB and the sample_df       |
    #                                                          |
    # Sample_df           -> Refertazione DB table             |
    #                                                          |
    # varianti_principali -> mutazioni_trovate                 |
    # altre_varianti      -> mutazioni_trovate_secondary       |
    # affetto_flag        -> affetto                           |
    # -------------------------------------------------------- #

    table_cols_non_affetto = [
            'sample_id', 'gestione', 'referto',
            'riferimento_id', 'familiarita', 'campionebiologico',
            'cognome', 'nome', 'sesso', 'nazione', 'luogonascita',
            'datanascita', 'codice_fiscale', 'clinico',
            'strutturainviante', 'pannello', 'fenotipo', 'geni_in_fenotipo',
            'tecnologia', 'coverage_medio', 'coverage_10', 'coverage_25', 'lista_buchi',
            'mutazioni_trovate', 'mutazioni_trovate_secondary',
            'conclusioni_screening',
            '\"AFFETTO\"', 'data_creazione',
            '\"operatoreGENETISTA\"', 'istituto', 'familiarita_relativa',
            'sensibilita_analitica', 'analisi_eseguite', 'metodi_utilizzati',
            'limiti_ngs', 'versionapi',
            'select_nucleofamiliare', 'select_nucleosenzaprobando',
            'select_nucleoprobtwist', 'somatico', '\"operatoreSANGER\"',
            'famiglianucleo', 'specifiche_refertazione'
        ]
       
    table_cols_affetto = [
            'sample_id', 'gestione', 'referto',
            'riferimento_id', 'familiarita', 'campionebiologico',
            'cognome', 'nome', 'sesso', 'nazione', 'luogonascita',
            'datanascita', 'codice_fiscale', 'clinico',
            'strutturainviante', 'pannello', 'fenotipo',
            'tecnologia', 'coverage_medio', 'coverage_10', 'lista_buchi',
            'mutazioni_trovate', 'mutazioni_trovate_secondary',
            'conclusioni_screening',
            '\"AFFETTO\"', 'data_creazione',
            '\"operatoreGENETISTA\"', 'istituto', 'familiarita_relativa',
            'sensibilita_analitica', 'analisi_eseguite', 'metodi_utilizzati',
            'limiti_ngs', 'versionapi',
            'select_nucleofamiliare', 'select_nucleosenzaprobando',
            'select_nucleoprobtwist', 'somatico', '\"operatoreSANGER\"',
            'famiglianucleo', 'specifiche_refertazione'
        ]


    # Make the column names of samples_df match the column names of the DB table
    sample_df.rename(columns={"APIversion": "versionapi"}, inplace=True)   
    sample_df.rename(columns={"varianti_principali": "mutazioni_trovate"}, inplace=True)   
    sample_df.rename(columns={"altre_varianti": "mutazioni_trovate_secondary"}, inplace=True)   
    
    # dummy plural rows
    def df_to_rows(df, cols):
        rows = []
        for _, r in df.iterrows():
            rows.append(tuple(none_if_nan(r[col]) for col in cols))
        
        #print(len(rows))
        return rows

    conn = psycopg2.connect(
            dbname="postgres", 
            user="bioinfo", 
            password="password_sicura", 
            host="192.168.1.143", 
            port="5432"
        )

    controllo = "Riuscito"

    # ----------------------------- #
    # Update access_accettazione    #
    # ----------------------------- #
    # ids = [str(s) for s in sample_df['sample_id']]

    # dummy plural ids
    ids = []
    # "dummy" loops
    for s in sample_df['sample_id']:
        ids.append(str(s))
        

    try:
        with conn.cursor() as cur:
            execute_values(
                cur,
                """
                UPDATE access_accettazione AS a
                SET    tipo_referto = 'Preanalisi'
                FROM  (VALUES %s) AS v(sample_id)
                WHERE  a.sample_id::text = v.sample_id;
                """,
                [(sid,) for sid in ids]
            )
        conn.commit()
        logger.info(f"[Sample {ids}]: push_to_db successeded!")

    except DatabaseError as e:
        raise(e)

    
    
    # ------------------------------------------------ #
    # Handle refertazione_refertazione table           #
    # ------------------------------------------------ #

    affetto_samples_df = sample_df.loc[
        (sample_df['affetto'] != 'NO')
    ].copy()
    
    if affetto_samples_df.empty:
        
        non_affetto = sample_df.loc[
            (sample_df['affetto'] != 'SI')
        ].copy()

        if not non_affetto.empty:

            # ---------------------------------------------- #
            # Upsert refertazione - for NON AFFETTO          #
            # ---------------------------------------------- #

            non_affetto = non_affetto.reindex(columns=table_cols_non_affetto)
            rows = df_to_rows(non_affetto, table_cols_non_affetto)
            quoted_cols = ', '.join(table_cols_non_affetto)

            UPSERT_SQL = f"""
            INSERT INTO refertazione_refertazione_new ({quoted_cols})
            VALUES %s
            ON CONFLICT (sample_id) DO UPDATE SET
                conclusioni_screening   = EXCLUDED.conclusioni_screening,
                sensibilita_analitica   = EXCLUDED.sensibilita_analitica,
                analisi_eseguite        = EXCLUDED.analisi_eseguite,
                metodi_utilizzati       = EXCLUDED.metodi_utilizzati,
                limiti_ngs              = EXCLUDED.limiti_ngs,
                versionapi              = EXCLUDED.versionapi,
                specifiche_refertazione = EXCLUDED.specifiche_refertazione;
            """  

    else:

        # ---------------------------------------------- #
        # Upsert refertazione - for AFFETTO              #
        # ---------------------------------------------- #
            
        affetto_samples_df = affetto_samples_df.reindex(columns=table_cols_affetto)
        rows = df_to_rows(affetto_samples_df, table_cols_affetto)
        quoted_cols = ', '.join(table_cols_affetto)

        UPSERT_SQL = f"""
        INSERT INTO refertazione_refertazione_new ({quoted_cols})
        VALUES %s
        ON CONFLICT (sample_id) DO UPDATE SET
            conclusioni_screening   = EXCLUDED.conclusioni_screening,
            sensibilita_analitica   = EXCLUDED.sensibilita_analitica,
            analisi_eseguite        = EXCLUDED.analisi_eseguite,
            metodi_utilizzati       = EXCLUDED.metodi_utilizzati,
            limiti_ngs              = EXCLUDED.limiti_ngs,
            versionapi              = EXCLUDED.versionapi,
            specifiche_refertazione = EXCLUDED.specifiche_refertazione;
        """
    
    if UPSERT_SQL:
        #print(UPSERT_SQL)
        with conn.cursor() as cur:
            execute_values(cur, UPSERT_SQL, rows)
        conn.commit()
    else:
        sys.exit("No rows to upsert")



def upsert_refertazione(samples_df : pd.DataFrame):
    pass


def none_if_nan(x):
    """Convert pandas/NumPy NaN → None so they go to PostgreSQL as NULL."""
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return None
    return x


if __name__ == "__main__":
    main()