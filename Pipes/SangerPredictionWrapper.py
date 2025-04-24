import multiprocessing
import glob
import os

from Pipes.Pipe import Pipe
from Pipeline import Pipeline
from Pipes.SangerPredictionPipe import SangerPredictionPipe
from Entities.Sample import Sample

class SangerWrapper(Pipe):

    def process(self, **kwargs):
        print("PROGRESS_FLAG: Running SangerWrapper…")
        # usa il numero di thread passato da kwargs, altrimenti il massimo disponibile
        n_procs = int(kwargs.get('threads', multiprocessing.cpu_count()))
        args_list = self.prepare_args(**kwargs)
        with multiprocessing.Pool(n_procs) as pool:
            pool.map(self.worker, args_list)
        return kwargs

    def worker(self, arg):
        sample = arg['sample']
        print(f"Processing sample {sample.name}")
        Pipeline(SangerPredictionPipe()).start(**arg)

    def prepare_args(self, **kwargs):
        # prendi il principal_directory da kwargs
        proj_dir = kwargs['principal_directory']
        samples_dir = os.path.join(proj_dir, 'sample_data')
        json_files = glob.glob(os.path.join(samples_dir, "*.json"))

        args_list = []
        for jf in json_files:
            sample = Sample.fromJSON(jf)
            # copia TUTTI i kwargs, incluso principal_directory, threads, bwa, gatk, ecc.
            arg_d = kwargs.copy()
            arg_d['sample'] = sample
            args_list.append(arg_d)
        return args_list
