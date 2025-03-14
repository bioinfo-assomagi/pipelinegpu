

from Pipes.Pipe import Pipe
import os
import glob
import config
import utils
import multiprocessing
import queue
import subprocess
import dir_tree
from Entities.Sample import Sample

class Fq2BamPipe(Pipe):

    def __init__(self) -> None:
        super().__init__()

    def process(self, **kwargs):
        print("PROGRESS_FLAG:{} - Running BWA (alignment) ... ".format('20%'), flush=True)
        # principal_directory = kwargs.pop("principal_directory", None)
        principal_directory = dir_tree.principal_directory.path
        parabricks_input_directory = dir_tree.principal_directory.parabricks_input.path # TODO: careful with 'parabricks_input' here, store it somewhere globally, or pass it with kwargs
        # resynced_samples = kwargs.pop("resynced_samples", None)
        # resynced_samples_list = []

        # if type(resynced_samples) is list:
        #     resynced_samples_list = resynced_samples
        #     print("resynced_samples_list: {}".format(resynced_samples_list))
        # if type(resynced_samples) is multiprocessing.queues.Queue:
        #     print("queue type: {}".format(type(resynced_samples)))
        #     print("Queue size: {}".format(resynced_samples.qsize()))
        #     while resynced_samples.qsize() > 0:
        #         resynced_samples_list.append(resynced_samples.get())
        # if type(resynced_samples) is queue.Queue:
        #     while resynced_samples.qsize() > 0:
        #         resynced_samples_list.append(resynced_samples.get())
        

        # '-w', "{}".format(config.DOCKER_WORKDIR),

        # '--in-fq', "{}/parabricks_input.txt".format(config.DOCKER_WORKDIR),

        """ NOTE: putting the inputs in a text file to use with --in-fq-list gives an error, issues regarding read groups. Read groups are required in --in-fq-list, but not in --in-fq.
        So let's do it with --in-fq. In Parbricks doc, --in-fq option should be explicitly written for each sample."""
        input_fq = []
        docker_input_parabricks = os.path.join(config.DOCKER_WORKDIR, 'parabricks_input')
        # samples = utils.group_samples(glob.glob("{}/*".format(parabricks_input_directory)))

        # for sample, data  in samples.items():
        #     option = "--in-fq {} {}".format(data['forward'], data['reverse'])
        #     option = option.replace(principal_directory, config.DOCKER_WORKDIR)
        #     input_fq.append(option)

        # _input_option = ' '.join(input_fq)

        bam_samples = []
        samples_list = [Sample.fromJSON(sample_json) for sample_json in glob.glob(os.path.join(dir_tree.principal_directory.sample_data.path, "*.json"))]
        for sample in samples_list:
            sample_name = sample.name
            pairs_forward = sample.pairs_forward_filename
            pairs_reverse = sample.pairs_reverse_filename

            # NOTE: The respective pair forward and reverse fq should be found inside the parabricks_input directory

            _input_option = "--in-fq {} {}".format(os.path.join(docker_input_parabricks, pairs_forward), os.path.join(docker_input_parabricks, pairs_reverse))

            #print(_input_option)
            
            os.system(' '.join(['docker', 'run', '--gpus', 'all', '--rm', 
                                '--volume', "{}/:{}".format(config.REF, config.DOCKER_REFDIR), 
                                '--volume', "{}/:{}".format(principal_directory, config.DOCKER_WORKDIR), 
                                '--volume', "{}/:{}".format(os.path.join(principal_directory, "bam"), config.DOCKER_OUTPUTDIR), 
                                "{}".format(config.PARABRICKS_VERSION), 
                            'pbrun', 'fq2bam', 
                            '--ref', "{}/{}".format(config.DOCKER_REFDIR, config.REF_GENOME_NAME), 
                            # '--num-gpus', '1',
                            _input_option,
                            '--out-bam', "{}/{}_final.bam".format(config.DOCKER_OUTPUTDIR, sample_name)]))
            
            
            #TODO: consider returning this
            # bai_path = "{}/{}_final.bam.bai".format(os.path.join(principal_directory, "bam"), sample_name)
            # new_bai_path = "{}/{}_final.bai".format(os.path.join(principal_directory, "bam"), sample_name)
            # os.rename(bai_path, new_bai_path)

        #     sample['bam'] = "{}/{}_final.bam".format(os.path.join(principal_directory, "bam"), sample_name)
        #     sample['bai'] = new_bai_path

        #     #os.system(' '.join(['rm', "{}/{}_final_chrs.txt".format(os.path.join(principal_directory, "bam"), sample_name)]))

            sample.bam = "{}/{}_final.bam".format(os.path.join(principal_directory, "bam"), sample_name)
            sample.bai = "{}/{}_final.bam.bai".format(os.path.join(principal_directory, "bam"), sample_name)
            sample.saveJSON()

        # kwargs.update({"principal_directory": principal_directory, "samples": resynced_samples_list})
        return kwargs