import json
import os
import pandas as pd

class Sample():

    name = None
    forward = None
    reverse = None
    bam = None
    vcf_haplotypecaller = None
    vcf_deepvariant = None
    vertcial_path = None

    def __init__(self, name, forward=None, reverse=None, bam=None, filepath=None) -> None:
        self.name = name
        self.forward = forward
        self.reverse = reverse
        self.bam = bam
        self.filepath = filepath

    
    @classmethod
    def fromJSON(cls, json_file):
        # NOTE: doesn't handle nested objects

        self = cls.__new__(cls)
        with open(json_file) as json_file:
            d = json.load(json_file)
        self.__dict__.update(d)
        return self
    
    @classmethod
    def fromDict(cls, d):
        self = cls.__new__(cls)
        self.__dict__.update(d)
        return self


    def toJSON(self):
        return json.dumps(
            self,
            default=lambda o: o.__dict__,
            sort_keys = True,
            indent = 4
        )
    
    def saveJSON(self):
        if self.filepath is None:
            raise("Cannot save, no filepath specified.")
        json_filename = os.path.join(self.filepath, "{}.json".format(self.name))
        with open(json_filename, 'w', encoding='utf-8') as f:
            json.dump(self.__dict__, f, ensure_ascii=False, indent=4)
    
    def set_filepath(self, filepath):
        self.filepath = filepath

    """TODO: implement dataframe getters"""

    def get_bed_df(self):
        try:
            return pd.read_csv(self.bed, sep="\t")
        except AttributeError:
            raise Exception("Bed attribute not found in the sample. Chec if the corresponding bed file exists.")

    def get_vertical_df(self):
        try:
            return pd.read_csv(self.vertical, sep="\t")
        except AttributeError:
            raise Exception("Vertical attribute not found in the sample. Check if the corresponding vertical file exists.")

    def get_vcf(self):
        pass


if __name__ == "__main__":
    sample = Sample("jurgen", forward="/home/magi/dummy/forward.fq", reverse="/home/magi/dummy/reverse/fq")
    #json = sample.toJSON()
    sample2 = Sample.fromJSON(sample.toJSON())
    sample2.name = "Gjergj"
    print(sample2.bam)