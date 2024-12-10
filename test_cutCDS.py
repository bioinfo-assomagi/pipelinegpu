import pandas as pd
import mygene

import os

mg = mygene.MyGeneInfo()

phenotype = pd.read_csv("/home/magi/PROJECT/diagnosys/RESULT/validation_INFER_INFERTILITA/pheno/phenotype", sep="\t")
gene_list = list(phenotype["gene"])

a = mg.querymany(gene_list, scopes="symbol",fields=["symbol","exons"], species="human", as_dataframe=True,returnall=True)
print(a["out"])




def cutCDS(GENE):
        #SAMPLE = sample_obj.name # TODO: fix naming

        # lastFOLDER = FOLDER.split('/')[-2]
        #print lastFOLDER
        # lastFolder is equal to principal_directory name
        # if DEST == 'r':
        # 	dest = 'ROVERETO'
        # 	download_folder = '/home/magi/VIRTUAL/EUREGIO/DOWNLOADS/NGSINFO/BED_INFO/%s/%s/' % (dest,lastFOLDER)
        # 	if not os.path.exists(download_folder):
        # 		os.makedirs(download_folder)
        # elif DEST == 'b':
        # 	dest = 'BOLZANO'
        # 	download_folder = '/home/magi/VIRTUAL/EUREGIO/DOWNLOADS/NGSINFO/BED_INFO/%s/%s/' % (dest,lastFOLDER)
        # 	if not os.path.exists(download_folder):
        # 		os.makedirs(download_folder)
        # elif DEST == 's':
        # 	dest = 'SANFELICE'
        # 	download_folder = '/home/magi/VIRTUAL/SANFELICE/DOWNLOADS/NGSINFO/BED_INFO/%s/%s/' % (dest,lastFOLDER)
        # 	if not os.path.exists(download_folder):
        # 		os.makedirs(download_folder)
        # elif DEST == 'z':
        # 	dest = 'RICERCA'
        # 	download_folder = '/home/magi/VIRTUAL/RICERCA/DOWNLOADS/NGSINFO/BED_INFO/%s/' % (lastFOLDER)
        # 	if not os.path.exists(download_folder):
        # 		os.makedirs(download_folder)

        mg = ''
        mv = ''
        #mv = myvariant.MyVariantInfo()
        mg = mygene.MyGeneInfo()
        a = []

        CHROM = []
        hgvs38 = []
        hgvs37 = []
        START = []
        END = []
        START37 = []
        END37 = []
        geno = []
        exone = []
        tipo = []
        gene = []
        STRAND = []
        variazione = []
        annotazione = []
        CDSSTART=[]
        CDSEND=[]
        refseq = []
        INFO =[]
        ufficial_gene=[]
        snplist = []
        a = []
        TABLE = ''
        gene_DF = pd.DataFrame()
        PREsnp_DF3 = pd.DataFrame()
        TABLE_SNP = pd.DataFrame()

        a = mg.querymany(list(GENE), scopes="symbol",fields=["symbol","exons"], species="human", as_dataframe=True,returnall=True)
        try:
            gene_DF1 = a['out']
            gene_DF1['symbol'].fillna('ciao',inplace=True)
            gene_DF1 = gene_DF1[gene_DF1['symbol']!='ciao']
            gene_DF1['info'] =  'OfficialName'
        except: 
            gene_DF1 = pd.DataFrame()


        gene_DF = gene_DF1
        try:
            gene_DF.dropna(subset=['exons'],axis=0,how='all',inplace=True)
            gene_DF.reset_index(inplace=True)
        except: gene_DF=pd.DataFrame()

        gene_DF.dropna(subset=['exons'],axis=0,how='all',inplace=True)
        gene_DF.reset_index(inplace=True)

        if len(gene_DF)>=1:
                i=0
                z = dict(gene_DF['exons'])
                for DICT in z.keys():
                    refseq_rna = z[DICT]
                    for RNA in refseq_rna:
                        try:
                            chromosome = RNA['chr']
                            if len(chromosome.split('_')) > 1:
                                chromosome = RNA[1]['chr']
                                strand = RNA[1]['strand']
                                cdsstart = RNA[1]['cdsstart']
                                cdsend = RNA[1]['cdsend']
                                exons_number = len(RNA[1]['position'])
                                if int(strand) == -1:
                                    _exon_ = max(range(exons_number))+1
                                elif int(strand) == 1:
                                    _exon_ = min(range(exons_number))+1
                                for exons in RNA[1]['position']:
                                    _CHROM_ = 'chr'+str(chromosome)
                                    start = exons[0]
                                    end = exons[1]
                                    CHROM.append('chr'+str(chromosome))
                                    START.append(str(start))
                                    END.append(str(end))
                                    CDSSTART.append(str(cdsstart))
                                    CDSEND.append(str(cdsend))
                                    refseq.append(str(RNA['transcript']))
                                    geno.append('ND')
                                    exone.append(str(_exon_))
                                    tipo.append('FRAGMENT')
                                    gene.append(str(gene_DF['symbol'][i]))
                                    INFO.append(str(gene_DF['info'][i]))
                                    STRAND.append(str(strand))
                                    variazione.append(str(gene_DF['symbol'][i])+':'+str(RNA['transcript'])+':ex'+str(_exon_))
                                    annotazione.append('')
                                    if int(strand) == -1:
                                        _exon_ = _exon_-1
                                    elif int(strand) == 1:
                                        _exon_ = _exon_+1
                            else:
                                chromosome = RNA[0]['chr']
                                strand = RNA[0]['strand']
                                cdsstart = RNA[0]['cdsstart']
                                cdsend = RNA[0]['cdsend']
                                exons_number = len(RNA[0]['position'])

                                if int(strand) == -1:
                                    _exon_ = max(range(exons_number))+1
                                elif int(strand) == 1:
                                    _exon_ = min(range(exons_number))+1
                                for exons in RNA[0]['position']:
                                    _CHROM_ = 'chr'+str(chromosome)
                                    start = exons[0]
                                    end = exons[1]
                                    CHROM.append('chr'+str(chromosome))
                                    START.append(str(start))
                                    END.append(str(end))
                                    CDSSTART.append(str(cdsstart))
                                    CDSEND.append(str(cdsend))
                                    refseq.append(str(RNA['transcript']))
                                    geno.append('ND')
                                    exone.append(str(_exon_))
                                    tipo.append('FRAGMENT')
                                    gene.append(str(gene_DF['symbol'][i]))
                                    INFO.append(str(gene_DF['info'][i]))
                                    STRAND.append(str(strand))
                                    variazione.append(str(gene_DF['symbol'][i])+':'+str(RNA['transcript'])+':ex'+str(_exon_))
                                    annotazione.append('')
                                    if int(strand) == -1:
                                        _exon_ = _exon_-1
                                    elif int(strand) == 1:
                                        _exon_ = _exon_+1
                        except KeyError:
                            chromosome = RNA['chr']
                            strand = RNA['strand']
                            cdsstart = RNA['cdsstart']
                            cdsend = RNA['cdsend']
                            exons_number = len(RNA['position'])
                            if int(strand) == -1:
                                _exon_ = max(range(exons_number))+1
                            elif int(strand) == 1:
                                _exon_ = min(range(exons_number))+1
                                
                            for exons in RNA['position']:
                                _CHROM_ = 'chr'+str(chromosome)
                                start = exons[0]
                                end = exons[1]
                                CHROM.append('chr'+str(chromosome))
                                START.append(str(start))
                                END.append(str(end))
                                CDSSTART.append(str(cdsstart))
                                CDSEND.append(str(cdsend))
                                refseq.append(str(RNA['transcript']))
                                geno.append('ND')
                                exone.append(str(_exon_))
                                tipo.append('FRAGMENT')
                                gene.append(str(gene_DF['symbol'][i]))
                                INFO.append(str(gene_DF['info'][i]))
                                STRAND.append(str(strand))
                                variazione.append(str(gene_DF['symbol'][i])+':'+str(RNA['transcript'])+':ex'+str(_exon_))
                                annotazione.append('')
                                if int(strand) == -1:
                                    _exon_ = _exon_-1
                                elif int(strand) == 1:
                                    _exon_ = _exon_+1
                        else:
                            print ('GENE NON TROVATO!!!',str(gene_DF['symbol'][i]))
                            pass
                    i+=1
        else:
                    CHROM.append('chr'+str(0))
                    START.append(str(0))
                    END.append(str(0))
                    CDSSTART.append(str(0))
                    CDSEND.append(str(0))
                    refseq.append(str('Unknwon'))
                    geno.append('ND')
                    exone.append(str(0))
                    tipo.append('FRAGMENT')
                    gene.append(str('Unknown'))
                    INFO.append(str('Unknown'))
                    STRAND.append(str(0))
                    variazione.append(str('Unknown'))
                    annotazione.append('')

        TABLE = pd.DataFrame({'#CHROM':pd.Series(CHROM),'cdsstart':pd.Series(CDSSTART),
                'cdsend':pd.Series(CDSEND),'GENE':pd.Series(gene),'refseq':pd.Series(refseq),
                'START':pd.Series(START),'END':pd.Series(END),'INFO':pd.Series(INFO),
                'exone':pd.Series(exone),'geno':pd.Series(geno),'tipo':pd.Series(tipo),
                'strand':pd.Series(STRAND),'variazione':pd.Series(variazione),'annotazione':pd.Series(annotazione)})
        TABLE['NUMCHROM'] = TABLE['#CHROM'].apply(lambda x: len(x.split('_')))

        print(TABLE)


#         TABLE = TABLE[TABLE['NUMCHROM']==1]
        
#         TABLE.drop(['NUMCHROM'], axis=1,inplace=True)
#         TABLE2_A = pd.merge(TABLE,hgmd_df,on=['GENE','refseq'],how='left')
#         TABLE2_A.dropna(subset=['hgmd'],inplace=True)
#         TABLE2_B = pd.merge(TABLE,eccezioni_df,on=['GENE','refseq'],how='left')
#         TABLE2_B.dropna(subset=['hgmd'],inplace=True)
#         TABLE2 = pd.concat([TABLE2_A,TABLE2_B]).sort_values(by=['hgmd','GENE','exone'])
#         if len(TABLE2) == 0:
#             TABLE2 = TABLE
#             TABLE2['hgmd'] = 'hgmd'
#         TABLE2.reset_index(inplace=True)
#         TABLE2.drop(['index'],axis=1,inplace=True)
#         TABLE2['START'] = TABLE2['START'].astype(int)
#         TABLE2['END'] = TABLE2['END'].astype(int)
#         TABLE2['cdsstart'] = TABLE2['cdsstart'].astype(int)
#         TABLE2['cdsend'] = TABLE2['cdsend'].astype(int)
#         mask1 = TABLE2['START'] < TABLE2['cdsstart']
#         mask2 = TABLE2['END'] < TABLE2['cdsstart']
#         mask3 = TABLE2['START'] > TABLE2['cdsend']
#         mask4 = TABLE2['END'] > TABLE2['cdsend']
#         print (TABLE2)
# ###############################################################################
#         TABLE2.loc[mask1,'CODING_START'] = 'NoNcoding'
#         TABLE2.loc[mask2,'CODING_END'] = 'NoNcoding'
#         TABLE2.loc[mask3,'CODING_START'] = 'NoNcoding'
#         TABLE2.loc[mask4,'CODING_END'] = 'NoNcoding'
#         TABLE2['CODING_START'].fillna('coding',inplace=True)
#         TABLE2['CODING_END'].fillna('coding',inplace=True)

#         pattern = 'NR_'
#         maskNR = TABLE2['refseq'].str.contains(pattern)
#         TABLE2.loc[maskNR,'CODING_START'] = 'coding'
#         TABLE2.loc[maskNR,'CODING_END'] = 'coding'

#         mask1 = TABLE2['CODING_START']=='NoNcoding'
#         mask2 = TABLE2['CODING_END']=='NoNcoding'
#         TABLE2.loc[mask1,'START'] = TABLE2['cdsstart']
#         TABLE2.loc[mask2,'END'] = TABLE2['cdsend']
# ###############################################################################
#         TABLE2.sort_values(by=['hgmd','GENE','exone'],inplace=True)
#         TABLE_GROUP = TABLE2.groupby(['GENE'])

#         df = pd.DataFrame()
#         dfx = pd.DataFrame()
#         for index, group in TABLE_GROUP:
#             if len(group) <= 2:
#                 group.drop_duplicates(subset=['#CHROM','START','END'],keep='last',inplace=True)
#             elif len(group) > 2:
#                 if (((group['CODING_START']=='NoNcoding').all())&((group['CODING_END']=='NoNcoding').all())):
#                     if group['hgmd'].isin(['hgmd']).any():
#                         group = group[group['hgmd'].isin(['hgmd'])]
#                         group.drop_duplicates(subset=['#CHROM','START','END'],keep='last',inplace=True)
#                     else:
#                         group.drop_duplicates(subset=['#CHROM','START','END'],keep='first',inplace=True)
#                 else:
#                     group.drop_duplicates(subset=['#CHROM','START','END'],keep='first',inplace=True) # keeps exons from the hgmd transcript only
#             else: pass
#             df = pd.concat([df,group])
#             dfx = pd.concat([df,group])
#         #print df[df['GENE']=='GRIN2A']
# ###############################################################################
#         TABLE2 = df
#         TABLE2['START'] = TABLE2['START'].astype(int)
#         TABLE2['END'] = TABLE2['END'].astype(int)
#         TABLE2['strand'] = TABLE2['strand'].astype(int)
#         TABLE2['exone'] = TABLE2['exone'].astype(int)
#         TABLE2['length'] = (TABLE2['END'])-(TABLE2['START'])
#         TABLE2['length'] = TABLE2['length'].astype(int)

#         #TABLE2['START'] = TABLE2['START']-15
#         #TABLE2['END'] = TABLE2['END']+15
#         TABLE2['START'] = TABLE2['START']-5
#         TABLE2['END'] = TABLE2['END']+5
# ###############################################################################################################
#         group = TABLE2.groupby(['GENE'])
#         size = pd.DataFrame(group.size()).reset_index()
#         size['size'] = size[0]
#         gene_unique = size[size['size']==1]
#         for _gene_ in list(gene_unique['GENE']): # the exons that are the only exons for a gene will be marked as coding in both start and end
#             mask = TABLE2['GENE'] == _gene_
#             TABLE2.loc[mask,'CODING_START'] = 'coding'
#             TABLE2.loc[mask,'CODING_END'] = 'coding'

#         #pattern = 'NM_'
#         #TABLE2 = TABLE2[TABLE2['refseq'].str.contains(pattern)]
#         TABLE2.drop_duplicates(subset=['#CHROM','START','END'],inplace=True)
#         TABLE2B = TABLE2[~((TABLE2['CODING_START']=='NoNcoding')&(TABLE2['CODING_END']=='NoNcoding'))]
#         TABLE3B = TABLE2B #[['#CHROM','START','END','GENE','exone','length','strand','refseq','hgmd']]

#         _col_ = ['#CHROM','START','END','GENE','exone','refseq',
#             'hgmd','strand','variazione','length',
#             'INFO','annotazione','cdsend','cdsstart',
#             'geno','tipo','CODING_START','CODING_END']

#         if TABLE3B['GENE'].isin(['CBS']).any():
#             print ("TROVATO!!!"),
#             MORETABLE = TABLE3B[TABLE3B['refseq']=='NM_001178008'].sort_values(by=['GENE','START'],ascending=[True,True])
#             MORETABLE.drop_duplicates(subset=['GENE','exone'],keep='last',inplace=True)
#             _TABLE_ = TABLE3B[TABLE3B['GENE']!='CBS']
#             TABLE3 = pd.concat([_TABLE_,MORETABLE])[_col_].sort_values(by=['#CHROM','START'],ascending=[True,True])
#             TABLE3B = pd.concat([_TABLE_,MORETABLE])[_col_].sort_values(by=['#CHROM','START'],ascending=[True,True])
#         else:
#             TABLE3 = TABLE3B[_col_].sort_values(by=['#CHROM','START'],ascending=[True,True])

#         if TABLE3B['GENE'].isin(['CFTR']).any():
#             mask1 = ((TABLE3B['GENE']=='CFTR') & (TABLE3B['exone']==10))
#             TABLE3B.loc[mask1,'START'] = 117548600
#             TABLE3 =  TABLE3B[_col_].sort_values(by=['#CHROM','START'],ascending=[True,True])
#             TABLE3B = TABLE3B[_col_].sort_values(by=['#CHROM','START'],ascending=[True,True])
#         else:
#             TABLE3 = TABLE3B[_col_].sort_values(by=['#CHROM','START'],ascending=[True,True])

#         TABLE4_38MINUS = TABLE3[['#CHROM','START','END','GENE','exone','length','strand','refseq','hgmd']].sort_values(by=['#CHROM','START','END'])
#         if 'HRURF' in list(GENE): # add this region (exon) for HRURF if it is present in the original gene list
#             TABLE_PLUS = pd.DataFrame([['chr8','22130597','22130712','HRURF',1,106,-1,'NM_001394132','hgmd']],
#                     columns=['#CHROM','START','END','GENE','exone','length','strand','refseq','hgmd'])
#             TABLE4_38A = pd.concat([TABLE4_38MINUS,TABLE_PLUS])
#         else:TABLE4_38A=TABLE4_38MINUS
#         #print (list(GENE))
#         if 'MITOCONDRIO' in list(GENE): # add these regions if MITOCONDRION is present in the original gene list
#             #print ('WWWWWWWWWWWWWWWWWWWWWWWWWWWWWW')
#             TABLE_MITOCONDRIO = pd.DataFrame([['chrMT','1','16569','MITOCONDRIO',0,16569,1,'NC_012920','hgmd']],
#                 columns=['#CHROM','START','END','GENE','exone','length','strand','refseq','hgmd'])
#             TABLE4_38A = pd.concat([TABLE4_38MINUS,TABLE_MITOCONDRIO]).sort_values(by=['#CHROM','START','END'])
#         else:TABLE4_38A=TABLE4_38MINUS

#         #TABLE4_38A = TABLE3.sort_values(by=['#CHROM','START','END'])
#         TABLE4_38A.drop_duplicates(subset=['#CHROM','START'],keep='last',inplace=True)
#         TABLE4_38A.drop_duplicates(subset=['#CHROM','END'],keep='first',inplace=True)

#         if len(TABLE_SNP) > 0:
#             TABLE_SNP.drop_duplicates(subset=['#CHROM','END'],keep='first',inplace=True)
#             TABLE4_38B = pd.concat([TABLE4_38A,TABLE_SNP]).reset_index()
#         else:
#             TABLE4_38B = TABLE4_38A

#         TABLE3 = TABLE4_38B.sort_values(by=['GENE','START'])
#         TABLE4 = TABLE4_38B[['#CHROM','START','END','GENE','exone','length','strand','refseq','hgmd']].sort_values(by=['GENE','START'])
#         file_ = (os.path.join(FOLDER,'TOTAL_BED_'+str(SAMPLE)))
#         TABLE3.to_csv(file_,sep='\t',index=False)
#         file_ = (os.path.join(FOLDER,'bed_'+str(SAMPLE)))
#         sample_obj.BED = file_ # TODO: fix this
#         TABLE4.to_csv(file_,sep='\t',index=False)
#         TABLE4.to_csv(file_download,sep='\t',index=False)
# #######################################################################################################################################################
#         TABLE2X = dfx
#         TABLE2X['START'] = TABLE2X['START'].astype(int)
#         TABLE2X['END'] = TABLE2X['END'].astype(int)
#         TABLE2X['strand'] = TABLE2X['strand'].astype(int)
#         TABLE2X['exone'] = TABLE2X['exone'].astype(int)
#         TABLE2X['length'] = (TABLE2X['END'])-(TABLE2X['START'])
#         TABLE2X['length'] = TABLE2X['length'].astype(int)
#         TABLE2X['START'] = TABLE2X['START']-15
#         TABLE2X['END'] = TABLE2X['END']+15
# #######################################################################################################################################################
#         group = TABLE2X.groupby(['GENE'])
#         size = pd.DataFrame(group.size()).reset_index()
#         size['size'] = size[0]
#         gene_unique = size[size['size']==1]
#         for _gene_ in list(gene_unique['GENE']):
#             mask = TABLE2X['GENE'] == _gene_
#             TABLE2X.loc[mask,'CODING_START'] = 'coding'
#             TABLE2X.loc[mask,'CODING_END'] = 'coding'

#         TABLE2X.drop_duplicates(subset=['#CHROM','START','END'],inplace=True)
#         TABLE2BX = TABLE2X[~((TABLE2X['CODING_START']=='NoNcoding')&(TABLE2X['CODING_END']=='NoNcoding'))]
#         TABLE3BX = TABLE2BX #[['#CHROM','START','END','GENE','exone','length','strand','refseq','hgmd']]

#         _col_ = ['#CHROM','START','END','GENE','exone','refseq',
#             'hgmd','strand','variazione','length',
#             'INFO','annotazione','cdsend','cdsstart',
#             'geno','tipo','CODING_START','CODING_END']

#         if TABLE3BX['GENE'].isin(['CBS']).any():
#             print ("TROVATO!!!"),
#             MORETABLEX = TABLE3BX[TABLE3BX['refseq']=='NM_001178008'].sort_values(by=['GENE','START'],ascending=[True,True])
#             MORETABLEX.drop_duplicates(subset=['GENE','exone'],keep='last',inplace=True)
#             _TABLE_X = TABLE3BX[TABLE3BX['GENE']!='CBS']
#             TABLE3X = pd.concat([_TABLE_X,MORETABLEX])[_col_].sort_values(by=['#CHROM','START'],ascending=[True,True])
#             TABLE3BX = pd.concat([_TABLE_X,MORETABLEX])[_col_].sort_values(by=['#CHROM','START'],ascending=[True,True])
#         else:
#             TABLE3X = TABLE3BX[_col_].sort_values(by=['#CHROM','START'],ascending=[True,True])

#         if TABLE3BX['GENE'].isin(['CFTR']).any():
#             mask1 = ((TABLE3BX['GENE']=='CFTR') & (TABLE3BX['exone']==10))
#             TABLE3BX.loc[mask1,'START'] = 117548600
#             TABLE3X =  TABLE3BX[_col_].sort_values(by=['#CHROM','START'],ascending=[True,True])
#             TABLE3BX = TABLE3BX[_col_].sort_values(by=['#CHROM','START'],ascending=[True,True])
#         else:
#             TABLE3X = TABLE3BX[_col_].sort_values(by=['#CHROM','START'],ascending=[True,True])

#         if 'HRURF' in list(GENE):
#             TABLE_PLUS = pd.DataFrame([['chr8','22130597','22130712','HRURF',1,106,-1,'NM_001394132','hgmd']],
#                     columns=['#CHROM','START','END','GENE','exone','length','strand','refseq','hgmd'])
#             TABLE3X = pd.concat([TABLE3BX[TABLE3BX['GENE']!='HRURF'],TABLE_PLUS])
#             TABLE3BX = pd.concat([TABLE3BX[TABLE3BX['GENE']!='HRURF'],TABLE_PLUS])
#         else:
#             TABLE3X = TABLE3BX[_col_].sort_values(by=['#CHROM','START'],ascending=[True,True])

#         TABLE4_38AX = TABLE3X[_col_].sort_values(by=['#CHROM','START','END'])
#         TABLE4_38AX[_col_].drop_duplicates(subset=['#CHROM','START'],keep='last',inplace=True)
#         TABLE4_38AX[_col_].drop_duplicates(subset=['#CHROM','END'],keep='first',inplace=True)

#         if len(TABLE_SNP) > 0:
#             TABLE_SNP.drop_duplicates(subset=['#CHROM','END'],keep='first',inplace=True)
#             TABLE4_38BX = pd.concat([TABLE4_38AX,TABLE_SNP]).reset_index()
#         else:
#             TABLE4_38BX = TABLE4_38AX

#         TABLE3X = TABLE4_38BX.sort_values(by=['GENE','START'])
#         TABLE4X = TABLE4_38BX[['#CHROM','START','END','GENE','exone','length','strand','refseq','hgmd']].sort_values(by=['GENE','START'])
#         TABLE4X_XY = pd.DataFrame([['chrX','11296910','11296920','AMELX',5,3,1,'NM_001142','hgmd'],
#                     ['chrY','6869900','6869920','AMELY',5,9,-1,'NM_001143','hgmd'],
#                     ['chrY','2787170','2787270','SRY',1,94,-1,'NM_003140','hgmd']],
#                     columns=['#CHROM','START','END','GENE','exone','length','strand','refseq','hgmd'])

#         TABLE4X = pd.concat([TABLE4X,TABLE4X_XY])
#         file_ = (join(FOLDER,'TOTAL_BEDX_'+str(SAMPLE)))
#         TABLE3X.to_csv(file_,sep='\t',index=False)
#         file_ = (join(FOLDER,'bedX_'+str(SAMPLE)))
#         file_download = (join(config.DOWNLOAD_FOLDER,'bed_'+str(SAMPLE)))

#         TABLE4X.to_csv(file_,sep='\t',index=False)
#         TABLE4X.to_csv(file_download,sep='\t',index=False)
# ############################################################################################################
#         # TABLE4 is the BED file
#         vertical,verticalX = create_vertical(TABLE4,TABLE4X,SAMPLE,FOLDER)
#         file_2 = (join(FOLDER,'vertical_'+str(SAMPLE)))
#         vertical.to_csv(file_2,sep='\t',index=False)
#         file_2X = (join(FOLDER,'verticalX_'+str(SAMPLE)))
#         verticalX.to_csv(file_2X,sep='\t',index=False)
#         print ('FINITO!!! TUTTO OK!!!')

#         sample_obj.vertical = file_2
#         sample_obj.verticalX = file_2X
#         sample_obj.bedX = file_
#         sample_obj.bed = (join(FOLDER,'bed_'+str(SAMPLE)))
#         sample_obj.saveJSON()

#         return vertical,verticalX,TABLE4


cutCDS(gene_list)