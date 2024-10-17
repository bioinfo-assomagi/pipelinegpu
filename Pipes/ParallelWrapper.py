import multiprocessing
import glob
import os
import dir_tree

from Pipeline import Pipeline
from Pipes.Pipe import Pipe
from Entities.Sample import Sample



""" Use this class when you want to make a Pipeline run in parallel. The Pipeline will be executed in parallel, but the pipes
in that Pipeline will be executed sequentially. You can also provide, of course, a Pipeline having only one Pipe. Since the ParallelWrapper
is a Pipe itself, you can build a "parent" Pipeline containing these parallel pipelines as pipes. """
class ParallelWrapper(Pipe):

    def __init__(self, pipeline) -> None:
        self.pipeline = pipeline


    def process(self, **kwargs):

        with multiprocessing.Pool(32) as pool:
            pool.map(self.worker, self.prepare_args(**kwargs))
            pool.close()
            pool.join()

        return kwargs

    def worker(self, kwargs):
        self.pipeline.start(**kwargs)

    def prepare_args(self, **kwargs):
        """ Build a Sample object from sample data located in the preset sample_data directory, and add them to the argument list that will be processed by the pipeline. """
        l = []
        # [sample1, sample2, ..., sampleN]
        sample_jsons = glob.glob(os.join(dir_tree.principal_directory.sample_data.path, "*.json"))
        for json_file in sample_jsons:
            sample = Sample.fromJSON(json_file)
            arg_d = kwargs.copy()
            arg_d.update({"sample": sample})
            l.append(arg_d)
        
        return l    

