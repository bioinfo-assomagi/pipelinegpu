import os
import sys
import datetime
import pandas as pd
import dir_tree
import glob
import numpy as np
import re
from __future__ import annotations
from pathlib import Path


from Pipes import Pipe


iconfig = {
    "exception": "/home/magi/PROJECT/diagnosys/bin/SELEZIONEAUTOMATICA/eccezioni.csv",
    "regioniproblematiche": "/home/magi/PROJECT/diagnosys/bin/SELEZIONEAUTOMATICA/regioniproblematiche.txt",
    "gene_varsome": "/home/magi/PROJECT/diagnosys/bin/SELEZIONEAUTOMATICA/GENEVAR/",
	"cols1" : ['sample','CHROM','POS','GENE','HGVS','ID','consequence','unbalance','decisionmaf','decisionINFO','definitivemaf',
			'MAF<20','comments','Description','MAF_selection_AD','MAF_selection_AR',
			'MAF_selection','predictors_decision','predictors_B','predictors_D','controls_nhomalt','allphenotype',
			'allinheritance','inheritance','prevalence','ada_score','rf_score','HGVS_c','HGVS_p','final_inheritance','STATO','zigosita'],
    "traduttore" : {'D':'D','P':'D','T':'B','N':'B','H':'D','B':'B','A':'D','L':'B','M':'D','-999.0':'NULL','-999':'NULL','.':'NULL','U':'NULL'},
    "traduttore2" : {'D':'D','P':'B','T':'B','N':'B','H':'D','B':'B','A':'D','L':'B','M':'D','-999.0':'NULL','-999':'NULL','.':'NULL','U':'NULL'},
	"PRED" : ['DANN','EigenPC','FATHMM',
            'LRT','MCAP','MetaLR',
            'MetaSVM','MutPred','MutationAssessor',
            'MutationTaster','PROVEAN','Polyphen2HDIV',
            'Polyphen2HVAR','SIFT','VEST3'],
	
    "prevalenza_path" : '/home/magi/PROJECT/diagnosys/bin/SELEZIONEAUTOMATICA/Prevalenza.csv',
	"prevalzenza_freq_maf_path" : '/home/magi/PROJECT/diagnosys/bin/SELEZIONEAUTOMATICA/Freq_maf.csv'

}


# --------------------------------------------------------------------------- #
# 1. I/O helpers – read CSV / TSV and return DataFrames
# --------------------------------------------------------------------------- #
def io_load_prevalenza(p: Path) -> pd.DataFrame:
    df = pd.read_csv(p, sep=",")
    return df.assign(**{"Sospetto diagnostico": lambda d: d["Sospetto diagnostico"].str.upper()}) \
             .fillna("NA")

def io_load_prev_maf(p: Path) -> pd.DataFrame:
    return pd.read_csv(p, sep="\t").fillna("NA")

def io_load_sample_ds(p: Path) -> pd.DataFrame:
    df = pd.read_csv(p, sep="\t")
    return (df.assign(sample=lambda d: d["sample"].astype(str)   # assure str dtype
                               .str.replace(r"\.202$", ".2020", regex=True)))

def io_load_problematic_region(p: Path) -> pd.DataFrame:
	regioniproblematiche = pd.read_csv(p, sep='\t',header=0)[['tipo','hgvs','categoria']]
	regioniproblematiche['HGVS'] = regioniproblematiche['hgvs']
	
def io_load_exception(p: Path) -> pd.DataFrame:
	df = pd.read_csv(p, sep='\t',header=0)

# --------------------------------------------------------------------------- #
# 2. Pure transformation helpers – never touch disk, never mutate caller
# --------------------------------------------------------------------------- #
def tf_load_sample_data(pheno: pd.DataFrame, other: pd.DataFrame) -> pd.DataFrame:
    merged = (pd.merge(pheno, other, on=["sample", "HGVS", "CHROM", "POS"],
                       how="outer", indicator=True)
                .assign(ID=lambda d: d["ID_x"].combine_first(d["ID_y"]),
                        GENE=lambda d: d["GENE_x"].combine_first(d["GENE_y"]))
                .drop(columns=["ID_x", "ID_y", "GENE_x", "GENE_y", "_merge"])
                .assign(GENE=lambda d: d["ALTGENE"],
                        decisionmaf=lambda d: d["decisionmaf"].astype(float)))
    return merged



folder_name = dir_tree.principal_directory.path
folder_coverage = dir_tree.principal_directory.coverage.path
folder_final = dir_tree.principal_directory.final.path
folder_control = dir_tree.principal_directory.control.path
folder_pheno = dir_tree.principal_directory.pheno.path
folder_temp = dir_tree.principal_directory.temp.path


def io_select_omim(panel):
    if 'OCULAR' in panel: sospetti_omim='/home/magi/PROJECT/diagnosys/bin/SELEZIONEAUTOMATICA/OMIM/OMIM_ALL.csv'
    elif 'CANCER' in panel: sospetti_omim='/home/magi/PROJECT/diagnosys/bin/SELEZIONEAUTOMATICA/OMIM/OMIM_ALL.csv'
    elif 'VASCULAR' in panel: sospetti_omim='/home/magi/PROJECT/diagnosys/bin/SELEZIONEAUTOMATICA/OMIM/OMIM_ALL.csv'
    elif 'NEUROLOGY' in panel: sospetti_omim='/home/magi/PROJECT/diagnosys/bin/SELEZIONEAUTOMATICA/OMIM/OMIM_ALL.csv'
    elif 'MIXED' in panel: sospetti_omim='/home/magi/PROJECT/diagnosys/bin/SELEZIONEAUTOMATICA/OMIM/OMIM_ALL.csv'
    elif 'INTEGRACARDIOSTANCHEZZA' in panel: sospetti_omim='/home/magi/PROJECT/diagnosys/bin/SELEZIONEAUTOMATICA/OMIM/OMIM_ALL.csv'
    elif 'LYMPHOBESITY' in panel: sospetti_omim='/home/magi/PROJECT/diagnosys/bin/SELEZIONEAUTOMATICA/OMIM/OMIM_ALL.csv'
    elif 'INFERTILITA' in panel: sospetti_omim='/home/magi/PROJECT/diagnosys/bin/SELEZIONEAUTOMATICA/OMIM/OMIM_ALL.csv'
    else: sospetti_omim = '/home/magi/PROJECT/diagnosys/bin/SELEZIONEAUTOMATICA/OMIM/OMIM_ALL.csv'

    return sospetti_omim

def io_load_omim(self, omim_path):
        SOSPETTI=pd.read_csv(omim_path, sep='\t', header=0)
        SOSPETTI.fillna('Unkwnown', inplace=True)
        new=SOSPETTI["GENE (geneMimNumber)"].str.split(" ", n = 1, expand = True)
        # making separate first name column from new data frame
        SOSPETTI["GENE"]= new[0]
        # making separate last name column from new data frame
        SOSPETTI["geneMimNumber"]= new[1]
        SOSPETTI.drop(columns =["GENE (geneMimNumber)"], inplace = True)

def load_SAMPLE_data(pheno_predict, other_annot):
	SAMPLE=pd.merge(pheno_predict, other_annot, on=['sample', 'HGVS', 'CHROM', 'POS'], how='outer', indicator=True)

	SAMPLE.loc[:,'ID'] = SAMPLE['ID_x']
	SAMPLE.loc[:,'GENE'] = SAMPLE['GENE_x']
	SAMPLE['ID'].fillna(SAMPLE['ID_y'])
	SAMPLE['GENE'].fillna(SAMPLE['GENE_y'])

	SAMPLE.drop('ID_x', axis=1, inplace=True)
	SAMPLE.drop('ID_y', axis=1, inplace=True)
	SAMPLE.drop('GENE_x', axis=1, inplace=True)
	SAMPLE.drop('GENE_y', axis=1, inplace=True)

	SAMPLE.loc[:,'GENE'] = SAMPLE['ALTGENE']
	#SAMPLE=SAMPLE[SAMPLE['_merge']=='both']
	SAMPLE.drop('_merge', axis=1, inplace=True)
	# name=SAMPLE_to_analyze[0].split('/')[-1]
	# name=name.split('.')[0]+'.'+name.split('.')[1]
	# name=name.split('_')[0]
	# ds=SAMPLE_DS[SAMPLE_DS['Sample']==name]['ds']
	SAMPLE.loc[:,'decisionmaf'] = SAMPLE['decisionmaf'].astype(float)
	return SAMPLE


def invert_MAF(SAMPLE, cols1, t=0.9):
	if len(SAMPLE) >=1:
		SAMPLE['MAX_MAF'] = SAMPLE['decisionmaf']
		SAMPLE.loc[(SAMPLE['MAX_MAF']>=t) & (SAMPLE['MAX_MAF']!=1),'comments']='MAF>90%: inversion applied'
		SAMPLE.loc[(SAMPLE['MAX_MAF']>=t) & (SAMPLE['MAX_MAF']!=1),'MAX_MAF']=1-SAMPLE.loc[SAMPLE['MAX_MAF']>=t,'MAX_MAF']
		#SAMPLE = SAMPLE[(SAMPLE['comments']!='MAF>90%: inversion applied')] #& (SAMPLE['unbalance'].str.split('=').str.get(1).astype(float)<=0.95)]
	else: SAMPLE = pd.DataFrame(columns=cols1)
	return SAMPLE


def tf_exclude_problemregion(df: pd.DataFrame, regions: pd.DataFrame) -> pd.DataFrame:
    """
    Exclude or flag variants overlapping known problematic regions.

    Parameters
    ----------
    df : pd.DataFrame
        Variant table with an 'HGVS' column of the form 'chrN:start-end:ref/alt'.
    region_df : pd.DataFrame
        Table with columns ['tipo','hgvs','categoria'] specifying
        genomic intervals and their types/categories.

    Returns
    -------
    pd.DataFrame
        A copy of `df` where:
          - PROBLEMATICA and NON REFERTABILE regions are marked as NOT SELECTED
          - OMOLOGIA regions are labeled as OMOLOGIA
          - FUNCTIONAL IMPACT regions containing 'patogenetica' select variants
            (MAF<20 = 1, Description='SELECTED')
          - FUNCTIONAL IMPACT regions containing 'benigna' label variants BENIGNA
    """
	
    df = df.copy()
    df['fromFILE'] = 1
	
    regs = regions.copy()
    regs[['CHROM','coords']] = regs['hgvs'].str.split(':', expand=True)
    regs[['start','end']] = regs['coords'].str.split('-', expand=True).astype(int)
    regs['CHROM'] = regs['CHROM'].str.replace('chr','', case=False).str.upper()
	
    # parse sample HGVS into numeric position
    samp = df.copy()
    samp[['CHROM','coords']] = samp['HGVS'].str.split(':', expand=True, n=1)
    samp[['pos_str','rest']] = samp['coords'].str.split('-', expand=True, n=1)
    samp['pos'] = pd.to_numeric(samp['pos_str'], errors='coerce')
    samp['CHROM'] = samp['CHROM'].str.replace('chr','', case=False).str.upper()
	
    # find variants located in problematic regions using vectorized masks, and apply the labels
    for tipo, group in regs.groupby('tipo'):
		# iterate through the group of regions of that typology
        for _, row in group.iterrows():
            mask = (
                (samp['CHROM'] == row['CHROM']) &
                (samp['pos'] >= row['start']) &
                (samp['pos'] <= row['end'])
            )
            if tipo == 'PROBLEMATICA':
                df.loc[mask, ['Description','fromFILE','comments']] = ['NOT SELECTED', 0, 'Reg. Problematica from FILE']
            elif tipo == 'NON REFERTABILE':
                df.loc[mask, ['Description','fromFILE','comments']] = ['NOT SELECTED', 0, 'Non Refertabile from FILE']
            elif tipo == 'OMOLOGIA':
                df.loc[mask, ['Description','fromFILE','comments']] = ['OMOLOGIA', 0, 'Omologia from FILE']
            elif tipo == 'FUNCTIONAL IMPACT':
                cat = row['categoria'].lower()
                if 'patogenetica' in cat:
                    df.loc[mask, ['MAF<20','Description','fromFILE','comments']] = [1, 'SELECTED', 1, 'Patogenetica from FILE']
                elif 'benigna' in cat:
                    df.loc[mask, ['Description','fromFILE','comments']] = ['BENIGNA', 0, 'Benigna from FILE']


    df.loc[df['HGVS']=='chr12:121626873-121626874:-/A', ['Description', 'fromFILE', 'comments']] = ['NOT SELECTED', 0, 'Non Refertabile from FILE']
    df.loc[df['HGVS']=='chr6:42698426-42698426:G/C','Description', 'fromFILE', 'comments'] = ['NOT SELECTED', 0, 'Non Refertabile from FILE']
		
    return df


def tf_include_exceptions(sample_df: pd.DataFrame, exception_df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply exception overrides to a sample variant DataFrame.

    For any variant whose HGVS or rsID appears in the exceptions file, this function:
        - Resets `fromFILE` to 0
        - Marks `MAF<20` as 1
        - Sets `Description` to 'SELECTED'
        - Sets `MAF_selection` to 1
        - Adds a `comments` entry 'Selected from exception'

    Parameters:
    - sample_df: DataFrame with at least 'HGVS', 'ID', and the target columns.
    - exception_df: DataFrame with columns 'HGVS' and 'rsID'.

    Returns:
    - A new DataFrame with exception overrides applied.
    """
	
    df = sample_df.copy()
    # Initialize fromFILE for all variants
    df['fromFILE'] = 0

    exceptions = exception_df.copy()

    # Apply overrides for HGVS matches
    hgvs_mask = df['HGVS'].isin(exceptions['HGVS'])
    df.loc[hgvs_mask, ['MAF<20', 'Description', 'MAF_selection', 'comments']] = [1, 'SELECTED', 1, 'Selected from exception']

    # Apply overrides for rsID matches
    id_mask = df['ID'].isin(exceptions['rsID'])
    df.loc[id_mask, ['MAF<20', 'Description', 'MAF_selection', 'comments']] = [1, 'SELECTED', 1, 'Selected from exception']

    return df



def cut_MAF(SAMPLE,cols1,t=0.03):
	if len(SAMPLE) >=1:
		SAMPLE.loc[((SAMPLE['decisionmaf']<=t) & (SAMPLE['decisionINFO'] != 'ALL')),'MAF<20']=1
		SAMPLE.loc[((SAMPLE['decisionmaf']>t) & (SAMPLE['decisionINFO'] != 'ALL')  & (SAMPLE['MAF<20'] != 1)),'MAF<20']=0
		SAMPLE.loc[((SAMPLE['decisionmaf']>t) & (SAMPLE['decisionINFO'] != 'ALL')),'Description']='DISCARDED: MAF > 3%'
		# Include in the analysis variants with  unknown MAF
		SAMPLE['MAF<20'].fillna(value=1, inplace=True)
		SAMPLE = SAMPLE[SAMPLE['MAF<20']==1]
	else: SAMPLE = pd.DataFrame(columns=cols1)
	return SAMPLE



def MAF_t(SAMPLE,SOSPETTI,genes,ds,PREV_MAF,PREV,cols1):
	SAMPLE['MAF_selection']=0
	SAMPLE['MAF_selection_AD']=0
	SAMPLE['MAF_selection_AR']=0
	if len(SAMPLE) >=1:
		SAMPLE.loc[:,'controls_nhomalt'] = SAMPLE['gnomAD_exomes_controls_nhomalt'].astype(float)+SAMPLE['gnomAD_genomes_controls_nhomalt'].astype(float)
		mask1 = SAMPLE['controls_nhomalt'] == -1998.0
		mask2 = SAMPLE['controls_nhomalt'] == -999.0
		SAMPLE.loc[mask1,'controls_nhomalt'] = 0
		SAMPLE.loc[mask2,'controls_nhomalt'] = 0
		SAMPLE.reset_index(drop=True, inplace=True)
		for index, row in SAMPLE.iterrows():
				if row['GENE'] in genes:
					inheritance=SOSPETTI[SOSPETTI['GENE']==row['GENE']]['phenotypeInheritance'].unique()
					_allphenotype_ = SOSPETTI[SOSPETTI['GENE']==row['GENE']]['phenotype (phenotypeMimNumber)'].unique()
					try: prevalence=PREV[PREV['Sospetto diagnostico']==ds.upper()]['Prevalenza'].item()
					except: prevalence = '1:2000'
					temp=PREV_MAF[PREV_MAF['Frequenza']==prevalence]
					try:
						if (isinstance(inheritance[0], str)) and (inheritance[0].find('X-linked')==-1):
							if row['MAX_MAF']<=temp[temp['Ereditarietá ']=='AD']['MAF'].item():
								SAMPLE.loc[index,'MAF_selection_AD']=1
							if row['MAX_MAF']<=temp[temp['Ereditarietá ']=='AR']['MAF'].item():
								SAMPLE.loc[index,'MAF_selection_AR']=1
						elif (isinstance(inheritance[0], str)) and (inheritance[0].find('X-linked')!=-1):
							n1=float(prevalence.split(':')[0])
							n2=float(prevalence.split(':')[1])
							if row['MAX_MAF']<=(n1/n2):
								SAMPLE.loc[index,'MAF_selection_AD']=1
								SAMPLE.loc[index,'MAF_selection_AR']=1
								SAMPLE.loc[index,'comments']='X-linked'
						else:
							SAMPLE.loc[index,'comments']='hereditary model is missing!'
							if row['MAX_MAF']<=temp[temp['Ereditarietá ']=='AR']['MAF'].item():
								SAMPLE.loc[index,'MAF_selection_AR']=1
					except:
						SAMPLE.loc[index,'comments']='hereditary model is missing! Not Gene in Data!'
						if row['MAX_MAF']<=temp[temp['Ereditarietá ']=='AR']['MAF'].item():
							SAMPLE.loc[index,'MAF_selection_AR']=1
					try: inheritancefinal=SOSPETTI[SOSPETTI['GENE']==row['GENE']]['CGD'].unique()
					except: inheritancefinal=SOSPETTI[SOSPETTI['GENE']==row['GENE']]['phenotypeInheritance'].unique()
					SAMPLE=get_final_inheritance(inheritance,inheritancefinal, SAMPLE, index) # TODO:
					SAMPLE.loc[index,'allinheritance'] = str(sorted(inheritance))
					SAMPLE.loc[index,'allphenotype'] = str(sorted(_allphenotype_))
					for item in sorted(inheritance):
						if not item == 'None':
							if 'Autosomal dominant' in item:
								SAMPLE['MAF_selection']= np.where(SAMPLE['MAF_selection']==1,1,SAMPLE['MAF_selection_AD'])
								SAMPLE.loc[index,'inheritance']= 'AD'
								SAMPLE.loc[index,'prevalence'] = temp[temp['Ereditarietá ']=='AD']['MAF'].item()
							elif 'X-linked' in item:
								SAMPLE['MAF_selection']= np.where(SAMPLE['MAF_selection']==1,1,SAMPLE['MAF_selection_AR'])
								SAMPLE.loc[index,'inheritance']= 'XL'
								SAMPLE.loc[index,'prevalence'] = temp[temp['Ereditarietá ']=='AR']['MAF'].item()
							elif 'Autosomal recessive' in item:
								SAMPLE['MAF_selection']= np.where(SAMPLE['MAF_selection']==1,1,SAMPLE['MAF_selection_AR'])
								SAMPLE.loc[index,'inheritance']= 'AR'
								SAMPLE.loc[index,'prevalence'] = temp[temp['Ereditarietá ']=='AR']['MAF'].item()
								break
							else:
								SAMPLE.loc[index,'inheritance']= 'AD/AR' #np.nan
								SAMPLE.loc[index,'prevalence'] = 0.02200
					SAMPLE.loc[(SAMPLE['MAF_selection']==0) & (pd.isnull(SAMPLE['Description'])),'Description']='MAF above threshold'
				else:
					print(row['GENE'],'errore gene')

	else: SAMPLE = pd.DataFrame(columns=cols1)
	try:
		SAMPLE['inheritance'].fillna('AD/AR',inplace=True)
		SAMPLE['inheritance'].fillna(0.02200,inplace=True)
	except:
		SAMPLE['inheritance'] = 'AD/AR'
		SAMPLE['prevalence'] = 0.02200

	return SAMPLE


def tf_MAF_t(sample: pd.DataFrame, sospetti: pd.DataFrame, gene_panel: set[str], diagnosis: str, maf_table: pd.DataFrame, prevalence_df: pd.DataFrame, columns = list[str]) -> pd.DataFrame:
	
    """
    Apply prevalence-aware MAF thresholds and inheritance logic.

    Parameters
    ----------
    sample : DataFrame
        Variants for a single sample (already merged & pre-filtered).
        Must contain at least the columns:
        ['GENE', 'MAX_MAF', 'decisionmaf', 'decisionINFO',
         'gnomAD_exomes_controls_nhomalt', 'gnomAD_genomes_controls_nhomalt']
    suspects : DataFrame
        OMIM gene table with at least 'GENE', 'phenotypeInheritance', 'CGD'.
    gene_panel : set[str]
        Genes targeted by the assay (BED file).
    diagnosis : str
        Diagnostic suspicion (e.g. “RETINITIS PIGMENTOSA”).
    maf_table : DataFrame
        Lookup with cols ['Frequenza', 'Ereditarietá ', 'MAF'].
    prevalence_df : DataFrame
        Lookup with cols ['Sospetto diagnostico', 'Prevalenza'].
    columns : list[str]
        Column template used elsewhere in the pipeline.

    Returns
    -------
    DataFrame
        Same rows/cols as *sample* with the fields
        ['MAF_selection', 'MAF_selection_AD', 'MAF_selection_AR',
         'inheritance', 'final_inheritance', 'prevalence', 'comments']
        updated in-place.
    """
    if sample.empty:
        return pd.DataFrame(columns=columns)
	
    sample = sample.copy()

    for col in ("MAF_selection", "MAF_selection_AD", "MAF_selection_AR"):
        sample[col] = 0

    # consolidate control allele counts, treating magic “missing” sentinels
    sample["controls_nhomalt"] = (
        sample["gnomAD_exomes_controls_nhomalt"]
        .replace({-1998, -999}, 0)
        .astype(float)
        + sample["gnomAD_genomes_controls_nhomalt"]
        .replace({-1998, -999}, 0)
        .astype(float)
    )
	
    # determine population prevalence string (e.g. “1:5000”) for this diagnosis
    try:
        prevalence_str: str = (
            prevalence_df.loc[
                prevalence_df["Sospetto diagnostico"].str.upper() == diagnosis.upper(),
                "Prevalenza",
            ]
            .iloc[0]
        )
    except IndexError:
        prevalence_str = "1:2000"  # sensible default
		
    # build quick MAF-threshold lookup  ->  {"AD": 0.0004, "AR": 0.0020}
    maf_limits = (
        maf_table.loc[maf_table["Frequenza"] == prevalence_str]
        .set_index("Ereditarietá ")["MAF"]
        .to_dict()
    )
	
    # helper -------------------------------------------------------------
    def _inheritance_flags(gene: str) -> tuple[str, float]:
        """
        Translate gene-level inheritance strings to a tag (AD/AR/XL/AD/AR)
        and the matching threshold.
        """
        inherit_raw = suspects.loc[suspects["GENE"] == gene, "phenotypeInheritance"]
        model = next((h for h in inherit_raw if isinstance(h, str)), "AR")

        if "X-linked" in model:
            n1, n2 = map(float, prevalence_str.split(":"))
            return "XL", n1 / n2
        if "Autosomal dominant" in model:
            return "AD", maf_limits.get("AD", 0.0004)
        return "AR", maf_limits.get("AR", 0.0220)

    # ── row-wise evaluation (vectorised where it matters) ────────────────────────
    for idx, row in sample.iterrows():
        gene = row["GENE"]

        if gene not in gene_panel:
            sample.at[idx, "comments"] = "gene not in panel"
            continue

        model, threshold = _inheritance_flags(gene)
        passes_maf = row["MAX_MAF"] <= threshold

        # fill columns based on inheritance model
        if model in {"AD", "XL"}:
            sample.at[idx, "MAF_selection_AD"] = int(passes_maf)
        if model in {"AR", "XL"}:
            sample.at[idx, "MAF_selection_AR"] = int(passes_maf)

        sample.at[idx, "MAF_selection"] = max(
            sample.at[idx, "MAF_selection_AD"],
            sample.at[idx, "MAF_selection_AR"],
        )
        sample.at[idx, "inheritance"] = model
        sample.at[idx, "prevalence"] = threshold

        if not passes_maf and pd.isna(row.get("Description")):
            sample.at[idx, "Description"] = "MAF above threshold"

    # default / tidy-up
    sample["inheritance"].fillna("AD/AR", inplace=True)

    return sample



def make_predictions(SAMPLE,traduttore,traduttore2,PRED,cols1):
	SAMPLE['predictors_decision'] = 'NULL'
	SAMPLE['predictors_B'] = 'NULL'
	SAMPLE['predictors_D'] = 'NULL'

	if len(SAMPLE) >=1:
		P_TABLE=pd.DataFrame()
		for predictor in PRED:
			if predictor+'_pred' in SAMPLE.columns:
				P_TABLE[predictor+'_result']=0
				for index, row in SAMPLE.iterrows():
					item=str(row[predictor+'_pred']).split('&')[0]
					r1 = re.findall(r"^\w",str(row[predictor+'_pred']))
					r2 = re.findall(r"^.&\w",str(row[predictor+'_pred']))

					if r1: item = r1[0]
					elif r2: item = r2[0].split('&')[1]
					else: item = 'n'

					if not str(item): P_TABLE.loc[index,predictor+'_result']='NULL'
					elif str(item)=='nan': P_TABLE.loc[index,predictor+'_result']='NULL'
					elif str(item)=='n': P_TABLE.loc[index,predictor+'_result']='NULL'
					else:
						if predictor == 'MutationTaster':
							if traduttore2[item]=='D': P_TABLE.loc[index,predictor+'_result']='D'
							elif traduttore2[item]=='B': P_TABLE.loc[index,predictor+'_result']='B'
							else: P_TABLE.loc[index,predictor+'_result']='NULL'
						else:
							if traduttore[item]=='D': P_TABLE.loc[index,predictor+'_result']='D'
							elif traduttore[item]=='B': P_TABLE.loc[index,predictor+'_result']='B'
							else: P_TABLE.loc[index,predictor+'_result']='NULL'
			else:
				P_TABLE[predictor+'_result'] = SAMPLE[predictor+'_rankscore'].astype(float)
				SAMPLE[predictor+'_rankscore']=SAMPLE[predictor+'_rankscore'].astype(float)
				SAMPLE.loc[SAMPLE[predictor+'_rankscore']==float(-999),predictor+'_rankscore']=-999

				# if predictor == 'REVEL':
				# 	P_TABLE.loc[(SAMPLE[predictor+'_rankscore']>=0.7),predictor+'_result']='D'
				# 	P_TABLE.loc[(SAMPLE[predictor+'_rankscore']<0.7),predictor+'_result']='B'
				# 	P_TABLE.loc[(SAMPLE[predictor+'_rankscore']==-999),predictor+'_result']='NULL'
				# 	P_TABLE.loc[(SAMPLE[predictor+'_rankscore']==np.NaN),predictor+'_result']='NULL'
				# 	P_TABLE.loc[(SAMPLE[predictor+'_rankscore'].isnull()),predictor+'_result']='NULL'
				# elif predictor == 'CADD':
				# 	P_TABLE.loc[(SAMPLE[predictor+'_rankscore']>=0.6),predictor+'_result']='D'
				# 	P_TABLE.loc[(SAMPLE[predictor+'_rankscore']<0.6),predictor+'_result']='B'
				# 	P_TABLE.loc[(SAMPLE[predictor+'_rankscore']==-999),predictor+'_result']='NULL'
				# 	P_TABLE.loc[(SAMPLE[predictor+'_rankscore']==np.NaN),predictor+'_result']='NULL'
				# 	P_TABLE.loc[(SAMPLE[predictor+'_rankscore'].isnull()),predictor+'_result']='NULL'
				#else:
				P_TABLE.loc[(SAMPLE[predictor+'_rankscore']>=0.644),predictor+'_result']='D'
				P_TABLE.loc[(SAMPLE[predictor+'_rankscore']<=0.290),predictor+'_result']='B'
				P_TABLE.loc[(SAMPLE[predictor+'_rankscore']<0.644)&(SAMPLE[predictor+'_rankscore']>0.290),predictor+'_result']='NULL'
				P_TABLE.loc[(SAMPLE[predictor+'_rankscore']==-999),predictor+'_result']='NULL'
				P_TABLE.loc[(SAMPLE[predictor+'_rankscore']==np.NaN),predictor+'_result']='NULL'
				P_TABLE.loc[(SAMPLE[predictor+'_rankscore'].isnull()),predictor+'_result']='NULL'

		for index,row in P_TABLE.iterrows():
			benign = 0
			deleterius = 0
			null = 0

			counts = pd.DataFrame(row.value_counts()).T
			try: deleterius=counts['D'].values[0]
			except: deleterius=0
			try: benign=counts['B'].values[0]
			except: benign=0
			try: null=counts['NULL'].values[0]
			except: null=0
			total = deleterius+benign+null
			if int(total)==int(len(PRED)):
				SAMPLE.loc[index,'predictors_B']=(str(benign)+'/'+str(total))
				SAMPLE.loc[index,'predictors_D']=(str(deleterius)+'/'+str(total))
			else: print ('NON TORNANO I CONTI DEI PREDITTORI!!! Verificare!!!')

		SAMPLE.loc[(SAMPLE['predictors_D'].str.split('/').str.get(0).astype(int)>=8),'predictors_decision']='D'
		SAMPLE.loc[(SAMPLE['predictors_B'].str.split('/').str.get(0).astype(int)>=8),'predictors_decision']='B'
		SAMPLE['predictors_decision'].fillna('NULL',inplace=True)

	else: SAMPLE = pd.DataFrame(columns=cols1)
	return SAMPLE



class VariantSelection(Pipe):

    def __init__(self):
        pass

    def process(self, **kwargs):
        
        begin_time = datetime.datetime.now()

        self.sample = kwargs.get("sample")
        self.panel = kwargs.get("panel")
        self.dest = kwargs.get("dest")

        sample = self.sample

        if self.dest == 'r':
            dest = 'rovereto'
            path_django = '/home/magi/VIRTUAL/MAGIS/NGS_RESULT/annot'
            database = '/home/magi/VIRTUAL/MAGIS/magis.db'
        elif self.dest == 'b':
            dest = 'bolzano'
            path_django = '/home/magi/VIRTUAL/EUREGIO/NGS_RESULT/annot'
            database = '/home/magi/VIRTUAL/EUREGIO/euregio.db'
        elif self.dest == 'z':
            dest = 'ricerca'
            path_django = '/home/magi/VIRTUAL/RICERCA/NGS_RESULT/annot'
            database = '/home/magi/VIRTUAL/RICERCA/ricerca.db'

        """ Call Loaders """
        SOSPETTI = io_load_omim(io_select_omim(self.panel))
		
        PROBLEMATIC_REGIONS = io_load_problematic_region(iconfig.regioniproblematiche)

        EXCEPTIONS = io_load_exception(iconfig.exception)
		
        PREV = io_load_prevalenza(iconfig.prevalenza_path)
		
        PREV_MAF = io_load_prev_maf(iconfig.prevalzenza_freq_maf_path)

        SAMPLE_DS = io_load_sample_ds(os.path.join(folder_pheno,'phenotype'))


        sample_x = str(sample.name)
        print ('----->'+sample_x.split('_')[0]+'<------')

        bed = pd.read_csv(sample.bed, sep ='\t') 
        genes = (bed['GENE'].drop_duplicates()).values #.str.translate()

        pheno_predict = pd.read_csv(sample.pheno_predict, dtype=str, sep="\t")
        pheno_predict.fillna('-999',inplace=True)

        other_annot=pd.read_csv(sample.other_annot, dtype=str, sep='\t')

        ds = SAMPLE_DS[SAMPLE_DS['sample']==sample_x.split('_')[0]]['malattia'].values
        sample_diagnosis = ds[0]
		
        _stat_cov = os.path.join(folder_coverage,str(sample_x.split('_')[0]),str(sample_x.split('_')[0])+'_stat_cov.csv')
        stat_cov = pd.read_csv(_stat_cov, sep='\t')
        sesso=stat_cov.loc[0,'sesso']
		
        """ Merge pheno_predict and other_annot load sample files and create a unique dataframe. """
        SAMPLE = load_SAMPLE_data(pheno_predict, other_annot)

        SAMPLE['MAF<20'] = np.NaN
        SAMPLE['probsanger']=SAMPLE['probsanger'].astype(float)
        SAMPLE['STATO'] = np.where(SAMPLE['probsanger']<=0.5, 'In Lavorazione', 'NOSANGER')

        """ 1) ::: FILTERS ::: """
        # Invert variants with MAF > 0.8
        SAMPLE=invert_MAF(SAMPLE, iconfig.cols1)
		
        # Include variants from regioni problematiche
        SAMPLE=tf_exclude_problemregion(SAMPLE, PROBLEMATIC_REGIONS)
		
        # Include variants from exception file
        SAMPLE=tf_include_exceptions(SAMPLE, EXCEPTIONS)
		
        # Exclude variants with MAF > 0.3
        SAMPLE=cut_MAF(SAMPLE, iconfig.cols1)
		
        # Exclude variants above MAF treshold
        SAMPLE=tf_MAF_t(SAMPLE, SOSPETTI, genes, sample_diagnosis, PREV_MAF, PREV, iconfig.cols1)
		
        """ 2) ::: PREDICTORS ::: """
        # assign predictions to variants
        SAMPLE=make_predictions(SAMPLE, iconfig.traduttore, iconfig.traduttore2, iconfig.PRED, iconfig.cols1)

    

    def main(self, sample):
        pass

       