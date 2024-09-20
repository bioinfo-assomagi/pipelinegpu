
from Pipes.Pipe import Pipe
import os
import config
import sys
import utils
import glob
import pandas as pd
import multiprocessing
import dir_tree
from Entities.Sample import Sample

from DBContext import DBContext

""" Pipe that sits in between two Pipelines. Processes the results of Pipeline1. """
class SampleListFam(Pipe):

    def __init__(self) -> None:
        super().__init__()

    def process(self, **kwargs):
       
        db_path = kwargs.pop('db_path', None)
        

        #sample_list_files = glob.glob("{}/sample_list_*.csv".format(dir_tree.principal_directory.path)) # TODO: sample_list files can also be passed in the queue
        #merged_samples_df = utils.merge_csv(sample_list_files)


        #merged_samples_df.to_csv(os.path.join(dir_tree.principal_directory.path, 'sample_list.csv'), sep='\t', index=False, encoding='utf-8')
        #sample_list = merged_samples_df['name'].tolist()

        sample_jsons = glob.glob(os.path.join(dir_tree.principal_directory.sample_data.path, "*.json"))
        samples = [Sample.fromJSON(json_file) for json_file in sample_jsons]
        sample_list = [sample.name for sample in samples]

        print("Building phenotype for samples: {}".format(sample_list))

        dbContext = DBContext(db_path)
        accetazione_df = dbContext.get_sample_familiari(sample_list)
        accetazione_df.to_csv("{}/sample_list_FAM.csv".format(dir_tree.principal_directory.path), sep='\t', index=False)

        phenotype_path = os.path.join(dir_tree.principal_directory.path, 'pheno', 'phenotype')
        pheno = DBContext(db_path).get_disease(sample_list) # could also use resynced_sample_list_names
        pheno.to_csv(phenotype_path, sep='\t', index=False, encoding='utf-8')

        #kwargs.update({'merged_samples_df': merged_samples_df, 'db_path': db_path})
        return kwargs


    



