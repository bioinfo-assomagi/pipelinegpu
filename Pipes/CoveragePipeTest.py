
from Pipes.Pipe import Pipe
import utils
import os
import pandas as pd
import csv
import glob
import logging

# silence all of urllib3’s DEBUG logs
logging.getLogger("urllib3").setLevel(logging.WARNING)

# if you only want to silence the connection‐pool chatter:
logging.getLogger("urllib3.connectionpool").setLevel(logging.WARNING)


import dir_tree

cols = ["#CHROM", "POS", "C%", "G%", "T%", "A%", "ins%", "del%", "sum", "DEPTH", "GENE", "exone", "length", "strand", "refseq", "hgmd", "sample"]

class CoveragePipeTest(Pipe):

    def __init__(self) -> None:
        super().__init__()

    def process(self, **kwargs):
        self.panel = kwargs.pop("panel")
        self.sample = kwargs.pop("sample")
        self.sample.temp_dir = os.path.join(dir_tree.principal_directory.temp.path, str(self.sample.name))
        dest = kwargs.pop("dest")
        self._logger = logging.getLogger(__name__)
        
        input_phenotype = os.path.join(dir_tree.principal_directory.pheno.path, "phenotype")

        if os.path.isfile(input_phenotype):
            self.phenotype = pd.read_csv(input_phenotype, sep="\t", header=0, dtype=str)
        else:
            raise Exception("No phenotype present. CoveragePipe terminating ...")

        genelist = list(self.phenotype["gene"][self.phenotype["sample"] == self.sample.name])
        self.vertical_df, self.verticalX_df, self.BED = self.cutCDS(genelist, dest) # cutCDS should be decoupled from CoverageAnalysis

        self.load_prereq()
        self.count_unbalance()
        count_disease, count_sex = self.merge_chunks()
        self.summary_disease(count_disease)
        self.summary_sex(count_sex)
        self.summary_macroarea()
        self.clean_temp_dir()

        #TODO: save results paths to sample.json

        kwargs.update({"panel": self.panel, "sample": self.sample, "dest": dest}) # only if donwstream pipes of the same thread will be present
        return kwargs


    def clean_temp_dir(self):
        pass

    def load_prereq(self):
        """
        The load_prereq function is responsible for specifically loading the vertical file - which again, is just an 
        exploded version of the BED file. The BED file loaded here, is the Panel specific BED file, as the disease
        specific BED file is already stored in the `Sample` object. At the same time, a file named buchiartificiali
        is loaded, that contains custom positions that are assumed to fall in regions of low coverage. 
        buchiartificiali in the current implementation, following the conventions of the old one is located in: /PROJECT/diagnosys/bin
        """
        #self.vertical_df, self.verticalX_df, self.vertical_macro_df = self.get_vertical_dataframes()
        _vertical_macro_path = utils.get_vertical_macro(self.panel)
        if _vertical_macro_path is not None:
            self.vertical_macro_df = pd.read_csv(_vertical_macro_path, sep="\t")
        else:
            self._logger.debug("No vertical_macro for panel {}".format(self.panel))
            self.vertical_macro_df = None
        self.buchiartificali = utils.get_buchiartificiali()             
        

    def count_unbalance(self):
        """
        Simple, call ``samtools mpileup`` and the _to_count file is produced, containing for each position in the
        covered regions of the genome, the piled bases, for each read which is aligned (covers) into that position.

        Then, due to the output being so large, the data is processed in chunks, and the bases our counted, after filtering the chunk with the BED
        files loaded from the load_prereq() function.

        TODO: make this functions receive inputs, instead of taking everything implicitly from the Sample object.
        Extract the data from the Sample object, and feed the data as arguments into these functions, to better 
        document them.
        """
        sample_name = str(self.sample.name)
        try:
            bam = self.sample.bam
        except Exception as e:
            raise("In coverage analysis, no bam file for sample {}".format(sample_name))
        
        # if self.sample.bam is None:
        #     raise Exception("No bam file for sample {}".format(sample_name))
        
        bam = self.sample.bam
        mpileup_out = os.path.join(dir_tree.principal_directory.temp.path, sample_name + "_to_count")
        folder_to_count = dir_tree.principal_directory.temp.to_count.path
        folder_to_macroarea = dir_tree.principal_directory.temp.to_macroarea.path
        #file_vcf = join(folder, "temp/", sample_name + "_samt.vcf")

        os.system(" ".join(["samtools", "mpileup", "-Q 0 -q 0 -d10000000 -L 100000 -A", bam, ">", mpileup_out]))

        mpileup_out_chunks = pd.read_csv(
                mpileup_out,
                sep="\t",
                header=None,
                quoting=csv.QUOTE_NONE,
                encoding="utf-8",
                low_memory=False,
                on_bad_lines="skip",
                names=["CHROM", "POS", "info", "DEPTH", "CALL", "quality"],
                chunksize=40 * 100024,
                
            )

              
        for index, chunk in enumerate(mpileup_out_chunks):
            chunk["sample"] = sample_name
            chunk.rename(columns={"CHROM": "#CHROM"}, inplace=True)


            # INNER JOIN WITH VERTICAL
            chunky_vertical = pd.merge(chunk, self.vertical_df, on=["#CHROM", "POS"], how='inner')
            chunky_verticalX = pd.merge(chunk, self.verticalX_df, on=["#CHROM", "POS"], how='inner')


            coverage_vertical = self.count_unbalance_per_filtered_chunk(chunky_vertical)
            coverage_verticalX = self.count_unbalance_per_filtered_chunk(chunky_verticalX)

            #coverage_vertical[["DEPTH", "sum", "selection"]] = coverage_vertical[["DEPTH", "sum", "selection"]].fillna(0, inplace=True)
            #coverage_verticalX[["DEPTH", "sum", "selection"]] = coverage_verticalX[["DEPTH", "sum", "selection"]].fillna(0, inplace=True)

            coverage_vertical.to_csv(os.path.join(folder_to_count, sample_name + "_" + str(index)) + "_disease", sep="\t", header=False, index=False)
            coverage_verticalX.to_csv(os.path.join(folder_to_count, sample_name + "_" + str(index)) + "_SEX", sep="\t", header=False, index=False)
            
            if self.vertical_macro_df is not None:
                chunky_vertical_macro = pd.merge(chunk, self.vertical_macro_df, on=["#CHROM", "POS"], how='inner')

                coverage_vertical_macro = self.count_unbalance_per_filtered_chunk(chunky_vertical_macro)
                #coverage_vertical_macro[["DEPTH", "sum", "selection"]] = coverage_vertical_macro[["DEPTH", "sum", "selection"]].fillna(0, inplace=True)
                coverage_vertical_macro.to_csv(os.path.join(folder_to_macroarea, sample_name + "_" + str(index) + "_macroarea"), sep="\t", index=False)


    def count_unbalance_per_filtered_chunk(self, filtered_chunk):
        """ Function which does the actual count of each nucleotide, over all reads aligned on a position covered by at least one read. This way counts are computed only on the
        filtered chunks. 

            Args:
                filtered_chunk (pandas.DataFrame): a chunk from the mpilesup output that is inner joined with the vertical

            Returns:
                pandas.DataFrame: a modified dataframe containing new columns representing the counts, and percentage counts, of each nucleotide.

        """
        #coverage = chunk[["CHROM", "POS", "DEPTH", "CALL"]]
    
        # col 5 of samtools mpileup output, each row represents the bases in the reads,
        # covering that position in the reference genome
        read_base_patterns = {".": r"[.]", 
                            ",": r"[,]",
                            "*": r"[*]",
                            "$": r"[$]",  
                            "C": r"[cC]", 
                            "G": r"[gG]", 
                            "T": r"[tT]", 
                            "A": r"[aA]", 
                            "ins": r"\+[0-9]+[ACGTNacgtn]+", 
                            "del": r"\-[0-9]+[ACGTNacgtn]+", 
                            "other": r"[*><$^]"
                            }
        
        data = []
        for key, value in read_base_patterns.items():
            data.append(filtered_chunk["CALL"].str.count(value).rename(key))
        
        coverage_counts = pd.concat(data, axis=1)

        coverage_counts['sum'] = coverage_counts[[".", ",", "C", "G", "T", "A", "ins", "del"]].sum(axis = 1)
        df = pd.concat([filtered_chunk[["#CHROM", "POS", "DEPTH", "CALL", "GENE", "exone", "length", "strand", "refseq", "hgmd"]], coverage_counts], axis = 1) # chunk contains the columns from mpileup (POS, DEPTH, C, G, etc.), as well as those from vertical (GENE, exone, length, refseq, etc.)
    
        df = df[["#CHROM", "POS", "DEPTH", "C", "G", "T", "A", "ins", "del", "sum", 'GENE',
                                'exone','length','strand','refseq','hgmd']]
        
        df[["DEPTH", "sum"]] = df[["DEPTH", "sum"]].fillna(0)
        df = df.astype({
            "sum":   "int64",
            "POS":   "int64",
            "DEPTH": "int64",
        })
        df["selection"] = 1
        

        for col in ["C", "G", "T", "A", "ins", "del"]:
            pct = (df[col].astype(float) / df["sum"]).round(3)
            df[f"{col}%"] = pct.map("{:.3f}".format)

        # df["C%"] = df["C"].astype(float) / df["sum"].astype(float)
        # df["G%"] = df["G"].astype(float) / df["sum"].astype(float)
        # df["T%"] = df["T"].astype(float) / df["sum"].astype(float)
        # df["A%"] = df["A"].astype(float) / df["sum"].astype(float)
        # df["ins%"] = df["ins"].astype(float) / df["sum"].astype(float)
        # df["del%"] = df["del"].astype(float) / df["sum"].astype(float)
        # df["C%"] = df["C%"].map("{:,.3f}".format)
        # df["G%"] = df["G%"].map("{:,.3f}".format)
        # df["T%"] = df["T%"].map("{:,.3f}".format)
        # df["A%"] = df["A%"].map("{:,.3f}".format)
        # df["ins%"] = df["ins%"].map("{:,.3f}".format)
        # df["del%"] = df["del%"].map("{:,.3f}".format)

        # df['sum'] = df['sum'].astype(int)
        # df['POS'] = df['POS'].astype(int)
        # df['DEPTH'] = df['DEPTH'].astype(int)
        # df['selection'] = df['selection'].astype(int)
        # df['exone'] = df['exone'].astype(int)
        # df['length'] = df['length'].astype(int)
        # df['strand'] = df['strand'].astype(int)

        df = df[['#CHROM','POS','C%','G%','T%','A%','ins%','del%','sum','DEPTH','GENE',
                                'exone','length','strand','refseq','hgmd',
                                'selection']].sort_values(by=['#CHROM','POS'],ascending=[True,True])

        return df


    def merge_chunks(self):
        #TODO: add full paths to the *_disease and *_SEX
        temp_folder = dir_tree.principal_directory.temp.path
        to_count_folder = dir_tree.principal_directory.temp.to_count.path

        final_disease_filename = os.path.join(temp_folder, "{}_final_disease".format(str(self.sample.name)))
        disease_filenames = os.path.join(to_count_folder, "{}*_disease".format(self.sample.name))

        final_sex_filename = os.path.join(temp_folder, "{}_final_sex".format(str(self.sample.name)))
        sex_filenames = os.path.join(to_count_folder, "{}*_SEX".format(self.sample.name))

        os.system(" ".join(["cat", disease_filenames, ">", final_disease_filename]))
        count_disease = pd.read_csv(final_disease_filename, sep="\t", header=None, names=cols)
        count_disease.fillna(0, inplace=True)

        os.system(" ".join(["cat", sex_filenames, ">", final_sex_filename]))
        count_sex = pd.read_csv(final_sex_filename, sep="\t", header=None, names=cols)
        count_sex.fillna(0, inplace=True)
        return count_disease, count_sex

    
    def summary_disease(self, COUNT):
        """ Return summary statistics about the coverage file filtered by the panel (the file named with the suffix of _final_disease, build by merging the respective _disease files.)
        
            Args:
                COUNT (pandas.DataFrame): the dataframe containing the coverage counts (produced by count_unbalance) which is filtered by the vertical
            
            Returns:
                pandas.DataFrame: a dataframe containing the buchi (positions with coverage less than 10) named with the suffix of _buchi
                pandas.DataFrame: a dataframe containing the marked (positions with coverage >= 10) named with the suffix of _marked
                pandas.DataFrame: a dataframe containing hte not marked (positions with coverage = 0) named with the suffix of _only_0
                pandas.DataFrame: a dataframe containing all the positions (slightly modified version of the input) named with the suffix _all #TODO: what are the `slight` modifications?
        """
        folder = dir_tree.principal_directory.coverage.path
        phenotype = self.phenotype
        BED = self.sample.get_bed_df()
        vertical = self.sample.get_vertical_df()
        buchiartificiali = utils.get_buchiartificiali()

        sample_name = self.sample.name
        
        folder_coverage = os.path.join(folder,sample_name)
        os.makedirs(folder_coverage, exist_ok = True)

        result = os.path.join(folder_coverage, sample_name + "_all")
        result_buchi = os.path.join(folder_coverage, sample_name + "_buchi")
        result_non_buchi = os.path.join(folder_coverage, sample_name + "_marked")
        result_not_marked = os.path.join(folder_coverage, sample_name + "_only_0")

        COUNT =  COUNT[['#CHROM','POS','C%','G%','T%','A%','ins%',
				'del%','sum','DEPTH']].sort_values(by=['#CHROM','POS'],ascending=[True,True])
        phenotype_ = phenotype[phenotype['sample'].astype(str) == str(sample_name)]
        a = phenotype_[['malattia','gene']].drop_duplicates()
        x1=pd.DataFrame(BED['GENE'].drop_duplicates())
        #x1 = pd.DataFrame(a['gene'])
        x2= pd.DataFrame({'GENE':pd.Series(['CTD-3074O7.11','TM4SF2'])})
        x = x1._append(x2)
        FILTER = vertical[vertical['GENE'].isin(x['GENE'])] # so merge two times with the vertical? NOTE: error on merging is thrown here. #CHROM seems to be different type between the two dataframes.
        #print("Colnames of Vertical = {}".format(FILTER.columns.values.tolist()))
        b = pd.merge(COUNT,FILTER,on=['#CHROM','POS'],how='right')
        for index,row in buchiartificiali.iterrows():
            x = row['#CHROM']
            y = row['START']
            z = row['END']
            mask1 =  ((b['#CHROM'] == x) & (b['POS'] >= y) & (b['POS'] <= z))
            b.loc[mask1,'DEPTH'] = 0
        b['DEPTH'].fillna(0,inplace=True)
        b['sum'].fillna(0,inplace=True)
        try: b['filt'].fillna(0,inplace=True)
        except: b['filt']=0
        b['sum'] = b['sum'].astype(int)
        b['POS'] = b['POS'].astype(int)
        b['DEPTH'] = b['DEPTH'].astype(int)
        b['filt'] = b['filt'].astype(int)
        b['exone'] = b['exone'].astype(int)
        b['length'] = b['length'].astype(int)
        b['strand'] = b['strand'].astype(int)

        COUNT = b
        COUNT['sample'] = sample_name
        COUNT.fillna(0,inplace=True)
        COUNT.drop(['filt'], axis=1,inplace=True)
        #print (len(COUNT))



        cov = len(COUNT[COUNT['DEPTH'] >= 10])
        buchi = len(COUNT[COUNT['DEPTH'] < 10])
        # cov = len(COUNT[COUNT['DEPTH'] >= 20])
        # buchi = len(COUNT[COUNT['DEPTH'] < 20])
        tot = cov+buchi
        COUNT['sample'] = str(sample_name)

        # buchi = COUNT[COUNT['DEPTH'] < 20]
        # marked = COUNT[COUNT['DEPTH'] >= 20]
        buchi = COUNT[COUNT['DEPTH'] < 10]
        marked = COUNT[COUNT['DEPTH'] >= 10]
        not_marked = COUNT[COUNT['DEPTH'] == 0]

        buchi.to_csv(result_buchi,sep='\t',index=False)
        marked.to_csv(result_non_buchi,sep='\t',index=False)
        not_marked.to_csv(result_not_marked,sep='\t',index=False)
        COUNT.to_csv(result,sep='\t',index=False)

        self._logger.debug('Len COV > 10: ' + str(cov))
        self._logger.debug('% Sequence COVERED: ' + '{:,.1f}'.format(float(cov)/float(tot)*100) + '%')
        
        try: self._logger.debug('% Sequence COVERED: ' + '{:,.1f}'.format(float(cov)/float(tot)*100) + '%')
        except ZeroDivisionError: self._logger.debug (0, 'debug', name="CoveragePipeTest") 

        self._logger.debug ('Len disease: ' + str(len(COUNT)))


    def summary_macroarea(self):
        """ Similar to summary_disease, but for macroarea files, that are coverage counts filtered by the vertical_macro. The final dataframes are written in the 'coverage' directory

            Args:
                macro_chunks (list): while not provided as arguments right now, a list of coverage counts files filtered by vertical macro.

            Returns:
                pandas.DataFrame: a dataframe containing the buchi (positions with coverage less than 10) named with the suffix of _buchi
                pandas.DataFrame: a dataframe containing the marked (positions with coverage >= 10) named with the suffix of _marked
                pandas.DataFrame: a dataframe containing hte not marked (positions with coverage = 0) named with the suffix of _only_0
                pandas.DataFrame: a dataframe containing all the positions (slightly modified version of the input) named with the suffix _all #TODO: what are the `slight` modifications?
                
        """

        folder_to_macroarea = dir_tree.principal_directory.temp.to_macroarea.path
        macro_chunks = glob.glob(os.path.join(folder_to_macroarea, "{}_*_macroarea".format(self.sample.name)))
        macro_count = pd.DataFrame(columns=["#CHROM", "POS", "C%", "G%", "T%", "A%", "ins%", "del%", "sum", "DEPTH", "GENE", "exone", "length", "strand", "refseq", "hgmd", "selection"])

        data_types = {
		"#CHROM": str,
		"POS": int,
		"C%": float,
		"G%": float,
		"T%": float,
		"A%": float,
		"ins%": float,
		"del%": float,
		"DEPTH": int,
		"GENE": str,
		"exone": int,
		"length": int,
		"strand": int,
		"refseq": str,
		"hgmd": str,
		"selection": int
	    }

        folder_coverage = os.path.join( dir_tree.principal_directory.coverage.path, self.sample.name)
        os.makedirs(folder_coverage, exist_ok = True)

        macro_data = []
        for macro_chunk in macro_chunks:
            macro_df = pd.read_csv(macro_chunk, sep="\t")
            macro_data.append(macro_df)
        
        

        df = pd.concat(macro_data, axis=0, ignore_index=True)
        macro_count = macro_count._append(df)
        macro_count = macro_count.sort_values(by=['#CHROM','POS'],ascending=[True,True])

        b = macro_count

        for index, row in self.buchiartificali.iterrows():
            x = row['#CHROM']
            y = row['START']
            z = row['END']

            mask1 =  ((b['#CHROM'] == x) & (b['POS'] >= y) & (b['POS'] <= z))
            b.loc[mask1,'DEPTH'] = 0

        b.fillna(0, inplace=True)

        b['DEPTH'].fillna(0,inplace=True)
        b['sum'].fillna(0,inplace=True)
        b['sum'] = b['sum'].astype(int)
        b['POS'] = b['POS'].astype(int)
        b['DEPTH'] = b['DEPTH'].astype(int)
        b['exone'] = b['exone'].astype(int)
        b['length'] = b['length'].astype(int)
        b['strand'] = b['strand'].astype(int)

        b['sample'] = str(self.sample.name)
        b.fillna(0, inplace=True)

        result = os.path.join(folder_coverage, self.sample.name + '_all_macro')
        result_buchi = os.path.join(folder_coverage, self.sample.name + '_buchi_macro')
        result_non_buchi = os.path.join(folder_coverage, self.sample.name + '_marked_macro')
        result_not_marked = os.path.join(folder_coverage, self.sample.name + '_only_0_macro')

        marked = b[b['DEPTH'] >= 10]
        buchi = b[b['DEPTH'] < 10]
        not_marked = b[b['DEPTH'] == 0]

        
        tot_length = len(marked) + len(buchi)

        buchi.to_csv(result_buchi,sep='\t',index=False)
        marked.to_csv(result_non_buchi,sep='\t',index=False)
        not_marked.to_csv(result_not_marked,sep='\t',index=False)
        b.to_csv(result,sep='\t',index=False)

        self._logger.debug ('Len BUCHI: ' + str(len(buchi)))
        try: 
            self._logger.debug ('% Sequence COVERED: ' + '{:,.1f}'.format(float(len(marked))/float(tot_length)*100) + '%')
        except ZeroDivisionError: self._logger.debug (0)

        self._logger.debug ('Len disease: ' + str(len(b)))


    def summary_sex(self, count_sex):
        """ In this case, summary_sex is equal to count_sex. """
        coverage_dir = dir_tree.principal_directory.coverage.path
        sample_dir = str(self.sample.name)
        
        # utils.print("Sending {}_finalsex to coverage dir ...".format(str(self.sample.name)), 'debug', name="CoveragePipeTest")
        count_sex.to_csv(os.path.join(coverage_dir, sample_dir, "{}_final_sex".format(str(self.sample.name))), sep="\t", index=False)
    

    def get_bam_filename(self):
        sample = self.sample
        


    def cutCDS(self, genelist, server_id):
        import cutCDS_jurgen 
        sample = self.sample
        folder_pheno = dir_tree.principal_directory.pheno.path
        vertical, verticalX, BED = cutCDS_jurgen.cutCDS(
            genelist, sample, folder_pheno, server_id
        )
        # Additional processing ... TODO

        sospetto = self.phenotype["malattia"][self.phenotype["sample"].astype(str) == str(self.sample.name)].unique()[0]
        if sospetto == "CONNETTIVOPATIE":
            region = pd.DataFrame({'#CHROM':pd.Series(['chr9','chr9','chr9','chr9','chr9','chr9','chr9','chr9','chr9','chr9','chr9','chr9','chr9','chr9','chr9']),
                                    'START':pd.Series([86584322,86585077,86585812,86585652,86586188,86586587,86586797,86587759,86588201,86588817,86589432,86590377,86591910,86592604,86593110]),
                                    'END':pd.Series([86584355,86585246,86585827,86585734,86586271,86586641,86587104,86587887,86588314,86588888,86589504,86590420,86591966,86592701,86593167]),
                                    'GENE':pd.Series(['int','int','int','int','int','int','int','int','int','int','int','int','int','int','int']),
                                    'length':pd.Series([34,170,16,83,84,55,308,129,114,72,73,44,57,98,58]),
                                    'exone':pd.Series([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]),'strand':pd.Series([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]),
                                    'refseq':pd.Series(['N','N','N','N','N','N','N','N','N','N','N','N','N','N','N']),
                                    'hgmd':pd.Series(['int','int','int','int','int','int','int','int','int','int','int','int','int','int','int'])})
        elif sospetto == "CROMATINOPATIE":
            region = pd.DataFrame({'#CHROM':pd.Series(['chr11','chr12','chr15','chr15','chr15','chr15','chr16','chr16','chr22','chr3','chr3','chr3','chr5','chr5','chr7','chr7','chr7','chr7','chr7','chr7','chr3','chr22','chr22','chr22','chr22','chr22','chr22','chrX']),
                                    'START':pd.Series([119077253,112910723,66995991,66995821,67001018,96875310,30134513,55515763,19770411,20202572,20202683,8591535,78280735,78280904,140624418,150693467,150700214,150700331,150700720,150700214,20216072,19748403,19753887,19754000,19754265,19754305,19444417,136649004]),
                                    'END':pd.Series([119077347,112910816,66996065,66995870,67001101,96875443,30134555,55515815,19770570,20202743,20202743,8591605,78280783,78280957,140624528,150693561,150700516,150700413,150700808,150700494,20216167,19748627,19754144,19754071,19754311,19754311,19444466,136649122]),
                                    'GENE':pd.Series(['CBL','PTPN11','int','SMAD6','int','NR2F2','MAPK3','int','int','int','int','int','int','ARSB','BRAF','NOS3','int','int','int','int','SGOL1','TBX1','TBX1','int','int','int','UFD1L','ZIC3']),
                                    'length':pd.Series([95,94,75,50,84,134,43,53,160,172,61,71,49,54,111,95,303,83,89,281,96,225,258,72,47,7,50,119]),
                                    'exone':pd.Series([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]),
                                    'strand':pd.Series([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]),
                                    'refseq':pd.Series(['N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N']),
                                    'hgmd':pd.Series(['int','int','int','int','int','int','int','int','int','int','int','int','int','int','int','int','int','int','int','int','int','int','int','int','int','int','int','int'])})



        return vertical, verticalX, BED
