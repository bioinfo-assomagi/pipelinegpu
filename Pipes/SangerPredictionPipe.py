#!/usr/bin/env python3


# ora posso fare:
import pandas as pd
import numpy as np
from joblib import load
import os
import json
from Pipes.Pipe import Pipe
from Entities.Sample import Sample


class SangerPredictionPipe(Pipe):
    """
    Pipeline step to predict Sanger confirmation, estimate contamination,
    decide MAF and write two output tables.
    """

    def __init__(self):
        # Carica il modello (hard‑coded path)
        model_path = "/home/alessandro/PROJECT/SKLEARN/logistic_model.pkl"
        self.model = load(model_path)

        # Colonne per i due output
        self.cols1 = [
            'sample','HGVS','CHROM','POS','GENE','ALTGENE','ID','types','QUAL','DEPTH','mapquality',
            'INHERITANCE','samtools_geno','gatk_geno','unbalance','consequence','impact','clin_sign',
            'strand','refseq','hgmd','num_exon','num_intron','HGVS_c','HGVS_p','CDS_position',
            'amino_acids','codons','variation_cosmic','sift','polyphen','decisionmaf','decisionINFO',
            'GMAF','ExAC_MAF','MAX_MAF','Adj_MAF','AFR_MAF','AMR_MAF','EAS_MAF','EUR_MAF','SAS_MAF',
            'AA_MAF','EA_MAF','pubmed','hgmd_mutation','hgmd_function','hgmd_phenotype','hgmd_pubmed',
            'probsanger','probcontaminazione','clinvar_MedGen_id','clinvar_OMIM_id',
            'clinvar_Orphanet_id','clinvar_clnsig','clinvar_review','gnomAD_exomes_POPMAX_AF',
            'gnomAD_exomes_POPMAX_nhomalt','gnomAD_exomes_controls_AF','gnomAD_exomes_controls_nhomalt',
            'gnomAD_genomes_controls_nhomalt','gnomAD_genomes_POPMAX_AF','gnomAD_genomes_POPMAX_nhomalt',
            'gnomAD_genomes_controls_POPMAX_AF','gnomAD_genomes_controls_POPMAX_nhomalt',
            'gnomAD_exomes_controls_AC','gnomAD_genomes_controls_AC'
        ]

        self.cols2 = [
            'sample','HGVS','CHROM','POS','GENE','ID','CADD_rankscore','DANN_rankscore',
            'EigenPC_rankscore','FATHMM_rankscore','FATHMM_pred','GERP_rankscore','Interpro_domain',
            'LRT_rankscore','LRT_pred','MCAP_pred','MetaLR_pred','MetaLR_rankscore','MetaSVM_pred',
            'MetaSVM_rankscore','MutPred_rankscore','MutationAssessor_pred',
            'MutationAssessor_rankscore','MutationTaster_rankscore','MutationTaster_pred',
            'PROVEAN_rankscore','PROVEAN_pred','Polyphen2HDIV_pred','Polyphen2HDIV_rankscore',
            'Polyphen2HVAR_pred','Polyphen2HVAR_rankscore','REVEL_rankscore','SIFT_rankscore',
            'SIFT_pred','SiPhy29way_rankscore','VEST3_rankscore','clinvar_clnsig','fathmmMKL_pred',
            'phastCons100way_rankscore','phastCons20way_rankscore','phyloP100way_rankscore',
            'phyloP20way_rankscore','ada_score','rf_score'
        ]

    def add_prediction(self, df: pd.DataFrame, p_threshold: float = 0.9) -> pd.DataFrame:
        df = df.copy()
        df['QUAL'] = pd.to_numeric(df['QUAL'], errors='coerce')
        df['DEPTH'] = pd.to_numeric(df['DEPTH'], errors='coerce')

        mask = (df['types'] == 'SVN') & (df['QUAL'] > 0)
        feats = df.loc[mask, ['QUAL', 'DEPTH']].assign(
            qual=lambda x: np.log(x['QUAL']),
            depth=lambda x: x['DEPTH']
        )[['qual', 'depth']]

        df['probsanger'] = 0
        if not feats.empty:
            probs = self.model.predict_proba(feats.values)[:, 1]
            df.loc[mask, 'probsanger'] = (probs >= p_threshold).astype(int)

        return df

    def evaluate_contamination(
        self,
        df: pd.DataFrame,
        hmean: float = 0.471707,
        hstd: float = 0.067168
    ) -> pd.DataFrame:
        df = df.copy()
        kval = pd.to_numeric(df['QUAL'], errors='coerce')
        dval = pd.to_numeric(df['DEPTH'], errors='coerce')
        mask = (
            (df['gatk_geno'] == 'het') &
            (df['types'] == 'SVN') &
            (kval >= 60) &
            (dval > 20) &
            df['unbalance'].notna() &
            (df['unbalance'] != 'del=nan')
        )

        unbal = df.loc[mask, 'unbalance'].str.split('=', n=1).str[1].astype(float)
        z = (unbal - hmean) / hstd
        percz = (z.abs().gt(2).sum() / len(z) * 100) if len(z) else 0.0
        df['probcontaminazione'] = f"{percz:.2f}"

        return df

    def make_maf_decisional(self, df: pd.DataFrame) -> pd.DataFrame:
        df = df.copy()
        exclude_pops = ['gnomAD_ASJ', 'gnomAD_FIN', 'gnomAD_OTHER', 'gnomAD_OTH']
        df['Adj_MAF2'] = (
            df['Adj_MAF'].fillna('unknown').astype(str)
              .str.split('&', n=1).str[0]
              .replace({pop: 'exclude' for pop in exclude_pops})
        )
        df['MAX_MAF2'] = np.where(
            df['Adj_MAF2'].str.contains('gnomAD', na=False),
            df['MAX_MAF'].replace(1, np.nan),
            np.nan
        )

        df['decisionmaf'] = (
            df['gnomAD_exomes_POPMAX_AF']
            .combine_first(df['MAX_MAF2'])
            .combine_first(df['ExAC_MAF'])
        )

        conds = [
            df['gnomAD_exomes_POPMAX_AF'].notna(),
            df['decisionmaf'].eq(df['MAX_MAF2'])
        ]
        df['decisionINFO'] = np.select(conds, ['POPMAX', 'POP'], default='ALL')

        return df

    def process(self, **kwargs):
        self.principal_directory = kwargs.pop("principal_directory")
        self.sample = kwargs.pop("sample")

        final_dir = os.path.join(self.principal_directory, "final")
        if not os.path.isdir(final_dir):
            raise FileNotFoundError(f"Directory finale non trovata: {final_dir}")
        if not os.listdir(final_dir):
            raise RuntimeError(f"Directory finale esiste ma è vuota: {final_dir}")

        sample_name = self.sample.name
        sample_base = sample_name.split('_')[0]

        # Carica manualmente il JSON dei metadati
        # 1) Se Sample ha attributo json_path (o simile), usalo, altrimenti ricomponi il percorso
        json_path = getattr(self.sample, 'json_path', None)
        if not json_path:
            # assume che il JSON si chiami <sample.name>.json nella cartella sample_data
            json_path = os.path.join(
                self.principal_directory,
                'sample_data',
                f"{self.sample.name}.json"
            )
        with open(json_path, 'r', encoding='utf-8') as jf:
            meta = json.load(jf)
        
        pheno_path = meta.get("pheno_annot")
        if not pheno_path:
            raise ValueError(f"Metadata 'pheno_annot' mancante per {sample_name}")

        df = pd.read_csv(pheno_path, sep="\t", encoding="utf-8")

        df_pred = self.add_prediction(df)
        df_cont = self.evaluate_contamination(df_pred)
        samplefinalnew = self.make_maf_decisional(df_cont)

        # CSV 1: pheno_predict
        valid_cols1 = [c for c in self.cols1 if c in samplefinalnew.columns]
        df_pheno = samplefinalnew.reindex(columns=valid_cols1)
        df_pheno['sample'] = sample_base
        out1 = os.path.join(final_dir, f"{sample_base}_pheno_predict.csv")
        df_pheno.to_csv(out1, sep="\t", index=False, encoding="utf-8")
        print(f"Wrote {len(df_pheno)} rows to {out1}")

        # CSV 2: other_annot
        exclude_consequences = {
            'synonymous_variant','5_prime_UTR_variant','upstream_gene_variant',
            'downstream_gene_variant','3_prime_UTR_variant','intergenic_variant',
            'non_coding_transcript_exon_variant&non_coding_transcript_variant',
            'intron_variant&non_coding_transcript_variant','regulatory_region_variant',
            'splice_region_variant&intron_variant',
            'splice_region_variant&intron_variant&non_coding_transcript_variant',
            'splice_region_variant&non_coding_transcript_exon_variant&non_coding_transcript_variant',
            'splice_region_variant&synonymous_variant'
        }
        mask = ~samplefinalnew['consequence'].isin(exclude_consequences)
        valid_cols2 = [c for c in self.cols2 if c in samplefinalnew.columns]
        df_other = (
            samplefinalnew.loc[mask]
            .reindex(columns=valid_cols2)
            .fillna(-999)
        )
        df_other['sample'] = sample_base
        out2 = os.path.join(final_dir, f"{sample_base}_other_annot.csv")
        df_other.to_csv(out2, sep="\t", index=False, encoding="utf-8")
        print(f"Wrote {len(df_other)} rows to {out2}")


        self.sample.pheno_predict = out1
        self.sample.other_annot  = out2
        sample_json_dir = os.path.join(self.principal_directory, 'sample_data')
        self.sample.set_filepath(sample_json_dir)
        self.sample.saveJSON()

        # 4) prosegui con il resto
        kwargs.update({
            "principal_directory": self.principal_directory,
            "sample": self.sample
        })
        return kwargs






