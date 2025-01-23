import multiprocessing
import glob
import dir_tree
import os

from Pipes.Pipe import Pipe
from Pipeline import Pipeline
from Pipes.CoverageStatisticsPipe import CoverageStatisticsPipe
from Entities.Sample import Sample

class CoverageStatisticsWrapperPipe(Pipe):

    def process(self, **kwargs):
        # args_list = self.prepare_args(**kwargs)
        # for kwargs in args_list:
        #     Pipeline(VariantFilterPipe()).start(**kwargs)
        

        with multiprocessing.Pool(32) as pool:
            pool.map(self.worker, self.prepare_args(**kwargs))
            pool.close()
            pool.join()

        return kwargs

    def worker(self, kwargs): # TODO: let's add the subsequent pipes here? And change the name from VariantFilterWrapper to something else? And launch or the remianing steps
        # as a Pipeline itself, that will lead to multiple Pipelines being launched in parallel.
        Pipeline(CoverageStatisticsPipe()).start(**kwargs)

    def prepare_args(self, **kwargs):
        m = multiprocessing.Manager()
        lock = m.Lock()
        l = []
        # [sample1, sample2, ..., sampleN]
        sample_jsons = glob.glob(os.path.join(dir_tree.principal_directory.sample_data.path, "*.json"))
        for json_file in sample_jsons:
            sample = Sample.fromJSON(json_file)
            arg_d = kwargs.copy()
            arg_d.update({"sample": sample, "lock": lock})
            l.append(arg_d)
        
        return l