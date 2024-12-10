#!/home/magi/miniconda3/envs/PY310/bin/python
#22/08/2015

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed MAY 18 2015

@author: Giuseppe Marceddu
"""

import sqlite3
import pandas as pd
from os import listdir, system

def get_disease(samples,dest):
	if dest == 'b':
		system('rsync -avz -e ssh bioinfo@192.168.1.120:/home/bioinfo/VIRTUAL/EUREGIO/euregio.db /home/magi/PROJECT/diagnosys/bin/DATABASE')
		data = '/home/magi/PROJECT/diagnosys/bin/DATABASE/euregio.db'
	elif dest == 'r':
		system('rsync -avz -e ssh bioinfo@192.168.1.120:/home/bioinfo/VIRTUAL/MAGIS/magis.db /home/magi/PROJECT/diagnosys/bin/DATABASE')
		data = '/home/magi/PROJECT/diagnosys/bin/DATABASE/magis.db'
	elif dest == 'z':
		system('rsync -avz -e ssh bioinfo@192.168.1.120:/home/bioinfo/VIRTUAL/RICERCA/ricerca.db /home/magi/PROJECT/diagnosys/bin/DATABASE')
		data = '/home/magi/PROJECT/diagnosys/bin/DATABASE/ricerca.db'
		#data = '/home/bioinfo/VIRTUAL/RICERCA/ricerca.db'

	db = sqlite3.connect(data)
	cursor = db.cursor()

	a_ = []
	for sample in samples:
		a = cursor.execute('''SELECT sample,panel_id,fenotipo_id FROM acept_sample WHERE sample=?''', (sample,)).fetchall()
		for x in a:
			#print x
			a_.append(x)

	b = cursor.execute('''SELECT id,malattia FROM acept_malattia''').fetchall()
	c = cursor.execute('''SELECT id,gene FROM acept_geni''').fetchall()
	d = cursor.execute('''SELECT malattia_id,geni_id FROM acept_malattia_gene_list''').fetchall()
	e = cursor.execute('''SELECT id, pannello FROM acept_pannelli ''').fetchall()

	df_a = pd.DataFrame(a_,columns=['sample','panel_id','malattia_id'])
	#print df_a
	df_b = pd.DataFrame(b,columns=['malattia_id','malattia'])
	#print df_b
	df_c = pd.DataFrame(d,columns=['malattia_id','geni_id'])
	#print df_c
	df_d = pd.DataFrame(c,columns=['geni_id','gene'])
	#print df_d
	df_e = pd.DataFrame(e,columns=['panel_id','panel'])
	#print df_e
	#print df_a,df_b,df_e
	df_ab = pd.merge(df_a,df_b,on='malattia_id',how='left')
	#print df_ab
	df_abc = pd.merge(df_ab,df_c,on='malattia_id',how='left')
	#print len(df_abc)
	df_abcd = pd.merge(df_abc,df_d,on='geni_id',how='left')
	df_abcde = pd.merge(df_abcd,df_e,on='panel_id',how='left')
	df_abcde.dropna(subset=['gene'],axis=0,how='all',inplace=True)
	df_abcde = df_abcde[['sample','panel','malattia','gene']]

	db.close()
	return df_abcde

#def add_gene():
#	dest = 'b'
#	if dest == 'r':
#		data = '/home/bioinfo/VIRTUAL/MAGI/magi.db'
#	elif dest == 'b':
#		data = '/home/bioinfo/VIRTUAL/EUREGIO/euregio.db'
#
#
##	gene_list = '/home/bioinfo/dataset/MAGIS/phenotype/uniq_genes'
#	gene_list = '/home/bioinfo/dataset/MAGIS/phenotype/lista_unique_geni_trusight'
#	gene_df = pd.read_csv(gene_list,sep='\t',header=True)
#
#	db = sqlite3.connect(data)
#	cursor = db.cursor()
#
#	print gene_df
#	for index, gene_ in gene_df.iterrows():
#		gene = gene_.values[0]
#		try:
#			a = cursor.execute('''INSERT INTO acept_geni(gene,timestamp_gene,updated_gene) VALUES(?,?,?)''',
#					(gene,'2015-08-24 08:47:57.612954','2015-08-24 08:47:57.613043'))
#		except sqlite3.IntegrityError:
#			pass
#
#		db.commit()

#if __name__=="__main__":
#	sample = ['OIVR20509.2020']
#	sample_malattia = get_disease(sample,'r')
	#print sample_malattia
	#gene = add_gene()
	#print gene
#	print sample_malattia
