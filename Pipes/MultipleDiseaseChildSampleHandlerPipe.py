import os
import sys
import utils
import glob
import pandas as pd
import dir_tree
import logging

from Entities.Sample import Sample


from DBContext import DBContext
from db_repo import DBRepo
import Pipeline
from Pipes.Pipe import Pipe
from pathlib import Path

""" In this pipe, the child samples id are taken. This comes with the new design of the business sytem, as multiple sample_ids (spawned due to the fact that a sample may be tested for multiple diseases),
are linked with a reference_id (parent sample, or the 'real' sample which will be tested for multiple diseases).
Each spawned sample will share the same .fastq, .bam, and .vcf with their parent sample, but will be analysed for different diseases, hence different genes, hence different BED files. 
In few words, a spawned child sample, is just the orginal sample (identified by the reference_id property) but corresponding to another disease. 

1) Get the samples (parents) you've processed so far and store them into a dict with (key, value) = (sample_name: str, ample: Sample(Object)))
2) Search the DB for the spawns of every sample in sample_dict
3) Save the retreived spawns as new JSONs in the sample_data directory

"""

class MultipleDiseaseSampleHandlerPipe(Pipe):

    def __init__(self) -> None:
        super().__init__()
 
    def process(self, **kwargs):

        print("PROGRESS_FLAG:{} - Running MultipleDiseaseSamplehandler ... ".format('60%'), flush=True)
        self._logger = logging.getLogger(__name__)
        self.panel = kwargs['panel']

        try:

            # 1) Get the samples (parents) you've processed so far and store them into a dict with (key, value) = (sample_name: str, ample: Sample(Object)))
            self.sample_dict = self.init_sample_dict()

            # 2) Search the DB for the spawns of every sample in sample_dict
            spawn_mapping = self.get_spawns_mapping(self.sample_dict.keys())
    
            # 3) Save the retreived spawns as new JSONs in the sample_data directory
            self.create_spawn_jsons(spawn_mapping)

            # 4) Create the sample_list.csv (complying to the old version of the pipeline)
            self.create_sample_list_file(self.sample_dict.values())

            # 5) Create the phenotype file containing the disease (SOSPETTO) of all samples (including the spawns) and the GENES belonging to each disease
            self.create_phenotype_file()

        except Exception as e:
            
            print(str(e), flush=True)
            sys.exit(1)
            

        return kwargs
    

    def init_sample_dict(self):
        d = {}
        sample_jsons = glob.glob(os.path.join(dir_tree.principal_directory.sample_data.path, "*.json"))

        samples = [Sample.fromJSON(json_file) for json_file in sample_jsons]
        for sample in samples:
            d[sample.name] = sample

        return d

    def get_spawns_mapping(self, samples: list) -> pd.DataFrame:
        """ 
        {"ref_id": "sample_id"} or {"parent_id": "spawn_id"}
        e.g.
        sample_mapping = {
            #"E366.2024" : "E99.2024",
            #"1.2017" : "1.2017"
        }
        """

        repo = DBRepo(DBContext("dummy_path"))
        mappings = []
        for parent_id in samples:
            sample_mapping_df = repo.get_spawns(parent_id)
            mappings.append(sample_mapping_df)

        df = pd.concat(mappings, ignore_index=True)

        return df

    """
    |sample_id  | origin_id |
    |-----------|-----------|
    | 9997.2024 | 9997.2024 |
    | 9997B.2024| 9997.2024 |
    """
    def create_spawn_jsons(self, spawn_mapping: pd.DataFrame) -> None:
        """ Additionally, it adds the spawwns to the class sample_dict. """
        # Create the json for the spawned childs. Up to this point, the spawns will be identical copies of the parents. They will share the same values for their attributes.
        for parent, spawn in spawn_mapping.itertuples(index=False, name=None):
            parent_id = str(parent)
            spawn_id = str(spawn)
            
            #spawn_id_object = copy.deepcopy(sample_dict[parent_id])
            spawn_object = Sample.fromDict(self.sample_dict[parent_id].__dict__)
            spawn_object.name = spawn_id
            spawn_object.saveJSON()
            self._logger.debug("Wrote spawn json: {}".format(spawn_object.name))
            self.sample_dict[spawn_id] = spawn_object    

    
    def create_sample_list_file(self, samples):
        """ Save updated sample_list (again, complying with the old flow of the pipeline). """

        sample_list_df = pd.DataFrame(columns=['name', 'forward', 'reverse'])
        
        for sample in samples:
            sr_row = pd.Series(data={'name': sample.name, 'forward': sample.forward, 'reverse': sample.reverse})
            sample_list_df = pd.concat([sample_list_df, sr_row.to_frame().T], ignore_index=True)

        sample_list_df.to_csv(os.path.join(dir_tree.principal_directory.path, 'sample_list.csv'), sep='\t', index=False, encoding='utf-8')


    def create_phenotype_file(self):

        """ Rereading the jsons here again; we should expect the new samples (spawns) to be there. 
        Reading from sample_data redone inside func: `create_phenotype_file` to keep it the function decoupled.
        """
        # TODO: put it into a seperate pipe

        sample_jsons = glob.glob(os.path.join(dir_tree.principal_directory.sample_data.path, "*.json"))
        samples = [Sample.fromJSON(json_file) for json_file in sample_jsons]
        sample_list = [sample.name for sample in samples]

        self._logger.debug("Building phenotype for samples: {}".format(sample_list))

        dbContext = DBContext("dummy_path")
        accetazione_df = dbContext.get_sample_familiari(sample_list)
        accetazione_df.to_csv("{}/sample_list_FAM.csv".format(dir_tree.principal_directory.path), sep='\t', index=False)

        phenotype_path = os.path.join(dir_tree.principal_directory.path, 'pheno', 'phenotype')
        pheno = DBContext("dummy_path").get_disease(sample_list) # could also use resynced_sample_list_names
       
        if pheno.iloc[:,3].isnull().any():
            self._logger.warning("Null GENES found in phenotype file ... ")
            self._logger.warning("Dropping the rows with NULL genes ... ")
            
            pheno.to_csv(Path(phenotype_path).parent / "phenotype_null_genes.csv", sep='\t', index=False, encoding='utf-8')
            self._logger.warning("Phenotype with nulls written at {}".format(Path(phenotype_path).parent / "phenotype_null_genes.csv"))
            
            pheno = pheno.dropna(subset=['gene'], axis=0, how='all')
    
        # Check if all samples are in phenotype, if not, stop, rerun
        # Group pheno by sample_id and than count and compare the count against sample_list
        try:
            pheno_grouped = pheno.groupby("sample")
        except Exception as e:
            self._logger.error("Error grouping phenotype by sample_id ... ")
            self._logger.error(str(e))
            sys.exit(1)

        if len(pheno_grouped) != len(sample_list):
            self._logger.error("Not all samples are in phenotype file. Probably some samples had only null genes. ")
            sys.exit(1)

        # check if any of the samples is in the wrong panel - TODO: for now stop the app; in the future maybe we skip this samples and remove their .json
        # for sample in sample_list:
        #     print(sample)
        #     if pheno[pheno["sample"] == sample]["panel"][0] != self.panel:
                
        #         self._logger.error("Sample {} is in the wrong panel ... ".format(sample))
        #         sys.exit(1)
        
        #print(pheno.loc[pheno["panel"] != self.panel, "sample"])
        if (pheno["panel"] != self.panel).any():
            s = pheno.loc[pheno["panel"] != self.panel, "sample"]
            val = s.iloc[0]  
            sample_str = str(val)
            self._logger.error("Sample {} in the wrong panel, will not abort, but be cautious ... ".format(sample_str))

            #sys.exit(1)

        pheno.to_csv(phenotype_path, sep='\t', index=False, encoding='utf-8')
        self._logger.debug("Phenotype file created at: {}".format(phenotype_path))
    



if __name__ == "__main__":

    dir_tree.build("/home/magi/PROJECT/diagnosys/RESULTS/TEST_3_Apr_2025_OCULARE_spawns/")
    MultipleDiseaseSampleHandlerPipe().process(dummy_arg = "dummy")