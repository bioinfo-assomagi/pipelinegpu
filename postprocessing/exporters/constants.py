dict_vusstatus = {
    'Cold':'FREDDA',
    'Middle': 'TIEPIDA',
    'Hot':'CALDA',
    'Unkown':'NON NOTO'
}

dict_vusstatus_it = {
    'FREDDA': 'Cold',
    'TIEPIDA': 'Middle',
    'CALDA':'Hot',
    'NON NOTO':'Unkown'
}

dict_BA1 = {
    'Stand-Alone': 'Stand Alone',
    'Stand-alone': 'Stand Alone',
    'Stand alone':'Stand Alone'
}

diction_cause = {
	'BA1_cause_final' : 'Allele frequency is >5% in Exome Sequencing Project, 1000 Genomes Project, or Exome Aggregation Consortium.',
	'BS1_cause_final': 'Allele frequency is greater than expected for disorder.',
	'BS2_cause_final': 'Observed in a healthy adult individual for a recessive (homozygous), dominant (heterozygous), or X-linked (hemizygous)\
 disorder, with full penetrance expected at an early age.',
	'BS3_cause_final': 'Well-established in vitro or in vivo functional studies show no damaging effect on protein function or splicing.',
	'BS4_cause_final': 'Lack of segregation in affected members of a family.',
	'BP1_cause_final': 'Missense variant in a gene for which primarily truncating variants are known to cause disease.',
	'BP2_cause_final': 'Observed in trans with a pathogenic variant for a fully penetrant dominant gene/disorder or observed in cis with a\
 pathogenic variant in any inheritance pattern.',
 	'BP3_cause_final': 'In-frame deletions/insertions in a repetitive region without a known function.',
	'BP4_cause_final' :'Multiple lines of computational evidence suggest no impact on gene or gene product (conservation, evolutionary,\
 splicing impact, etc).',
	'BP5_cause_final': 'Variant found in a case with an alternate molecular basis for disease.',
	'BP6_cause_final': 'Reputable source recently reports variant as benign, but the evidence is not available to the laboratory to\
 perform an independent evaluation.',
	'BP7_cause_final': 'A synonymous (silent) variant for which splicing prediction algorithms predict no impact to the splice consensus\
 sequence nor the creation of a new splice site AND the nucleotide is not highly conserved.',
	'PVS1_cause_final':'Null variant (nonsense, frameshift, canonical ±1 or 2 splice sites, initiation codon, single or multiexon deletion)\
 in a gene where LOF is a known mechanism of disease.',
 	'PS1_cause_final':'Same amino acid change as a previously established pathogenic variant regardless of nucleotide change..',
	'PS2_cause_final':'De novo (both maternity and paternity confirmed) in a patient with the disease and no family history.',
 	'PS3_cause_final':'Well-established in vitro or in vivo functional studies supportive of a damaging effect on the gene or gene product.',
	'PS4_cause_final':'The prevalence of the variant in affected individuals is significantly increased compared with the prevalence in controls.',
	'PM1_cause_final':'Located in a mutational hot spot and/or critical and well-established functional domain (e.g., active site of an enzyme)\
 without benign variation.',
	'PM2_cause_final':'Absent from controls (or at extremely low frequency if recessive) in Exome Sequencing Project, 1000 Genomes Project, or\
 Exome Aggregation Consortium.',
	'PM3_cause_final': 'For recessive disorders, detected in trans with a pathogenic variant.',
	'PM4_cause_final':'Protein length changes as a result of in-frame deletions/insertions in a non-repeat region or stop-loss variants.',
	'PM5_cause_final':'Novel missense change at an amino acid residue where a different missense change determined to be pathogenic has been seen\
 before.',
	'PM6_cause_final': 'Assumed de novo, but without confirmation of paternity and maternity.',
	'PP1_cause_final':'Cosegregation with disease in multiple affected family members in a gene definitively known to cause the disease.',
	'PP2_cause_final':'Missense variant in a gene that has a low rate of benign missense variation and in which missense variants are a common\
 mechanism of disease.',
	'PP3_cause_final':'Multiple lines of computational evidence support a deleterious effect on the gene or gene product (conservation, evolutionary,\
 splicing impact, etc.).',
	'PP4_cause_final':'Patient’s phenotype or family history is highly specific for a disease with a single genetic etiology.',
	'PP5_cause_final':'Reputable source recently reports variant as pathogenic, but the evidence is not available to the laboratory to perform an\
 independent evaluation.'
}

columns_criteri_cause_final = ['BA1_cause_final', 'BS1_cause_final', 'BS2_cause_final', 'BS3_cause_final', 'BS4_cause_final', 'BP1_cause_final',
									'BP2_cause_final', 'BP3_cause_final', 'BP4_cause_final', 'BP5_cause_final', 'BP6_cause_final', 'BP7_cause_final', 'PVS1_cause_final',
									'PS1_cause_final', 'PS2_cause_final', 'PS3_cause_final', 'PS4_cause_final', 'PM1_cause_final', 'PM2_cause_final', 'PM3_cause_final',
									'PM4_cause_final', 'PM5_cause_final', 'PM6_cause_final', 'PP1_cause_final', 'PP2_cause_final', 'PP3_cause_final', 'PP4_cause_final',
									'PP5_cause_final']

dict_choice_interpretation = {
	'VUS': 'Uncertain Significance',
	'Patogenetica': 'Pathogenic',
	'Probabilmente Patogenetica': 'Likely Pathogenic',
	'Probabilmente Benigna' : 'Likely Benign',
	'Benigna': 'Benign'
}
