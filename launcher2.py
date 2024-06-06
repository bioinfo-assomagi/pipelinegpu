import argparse
import time
from multiprocessing import Process, Queue
import sys

import config as config
from PipelineAssembler import PipelineAssembler
import utils
import math
import coverage_async
import tests

def parseInput():
    parser = argparse.ArgumentParser(prog='MAGI EUREGIO DIAGNOSYS', description='Pipe from FASTQ to BAM',
                                     epilog='Need to BWA 0.7.17 and SAMTOOLS v1.14 and GATK 4.2 (SET on global PATH or INSERT your path [--bwa] [--samtools]')

    parser.add_argument('-p', '--path', metavar='PATH', default=config.PROJECT_PATH,
                        help='[Default = Pipe run in same folder when launch command]')

    parser.add_argument('-d', '--dest', metavar='destination', choices=['b', 'r', 'z'],
                        required=True,
                        help='Choices destination: b = bolzano; r = rovereto; z = ricerca,  - required ')

    parser.add_argument('-pan', '--panel', metavar='PANEL', required=True,
                        help='Pannelli Consigliati: 1:CLM2; 2:OSSNEW; 3:OCULARE; 4:OCULARE2; 5:SIFSR; 6:ALLGENE;'
                             '7:trusightone; 8:GENETEST1; 9:GENETEST2; 10:CANCER; 11:MALFORMATIONI VASCOLARI;'
                             '12:ANOMALIE VASCOLARI; 13:INFERTILITA; 14:OBESITA; 15:GENODERMATOSI\nRequired')

    parser.add_argument('-fq', '--fastq', metavar='FASTQ',
                        help='FULL PATH to find FASTQ FILE - [default FROM NAS]')

    parser.add_argument('-proj', '--project', metavar='PROJECT NAME', help='Insert Family name to create Folder')

    parser.add_argument('-g', '--genome', metavar='choose Assembling', choices=['geno37', 'geno38'],
                        default='geno38',
                        help='Choices: geno37, geno38 - Run Analysis with genome GRCh38 or GRCh37')

    parser.add_argument('-q', '--quality', metavar='PHRED QUALITY FILTER', default=18,
                        help=' Choose a quality threshold [Default = 18]')

    parser.add_argument('-Q', '--quality_perc', metavar='PERC QUALITY READ', default=97,
                        help='Percentual of quality for read [Default = 97]')

    parser.add_argument('-N', '--threads', metavar='THREADS', default='12',
                        help='Number of threads to run [20]')

    parser.add_argument('-m', '--memory', metavar='MERORY', default='15G',
                        help='max memory per thread; suffix K/M/G [defaul = 15G]')

    parser.add_argument('-b', '--bwa', metavar='BWA', default='bwa',
                        help='Insert BWA path')

    parser.add_argument('-s', '--samtools', metavar='SAMTOOLS', default='samtools',
                        help='Insert SAMTOOLS path')

    parser.add_argument('-t', '--samstat', metavar='SAMSTAT', default='samstat',
                        help='Insert SAMSTAT path')

    parser.add_argument('-gk', '--gatk', metavar='GATK', default=config.GATK, help='Insert gatk path')

    parser.add_argument('-v', '--bcftools', metavar='BCFTOOLS', default='bcftools',
                        help='Insert BCFTOOLS path')

    parser.add_argument('-ovr', '--over', action='store_true',
                        default=False,
                        help='Choices: ALLERT!!! Set on "False" overwrite old data [Default = False]')

    parser.add_argument('-vep', '--VEP', metavar='variant effect predictor', default=config.VEPS,
                        help='Insert VEP PATH')

    parser.add_argument('-del', '--delete', metavar='Delete temp file after Backup', choices=['False', 'True'],
                        required=True,
                        help='Default False = Not Delete, put True = Delete')

    parser.add_argument('-freq', '--runfrequency', metavar='Choose if Run Frequency Calculation or Not',
                        choices=['False', 'True'],
                        default='True',
                        help='Default False = Not Delete, put True = Delete')

    parser.add_argument('-core', '--runcore', metavar='Choose if Run Frequency Calculation or Not',
                        choices=['False', 'True'],
                        default='True',
                        help='Default False = Not Delete, put True = Delete')

    parser.add_argument('-importapp', '--importAPP', metavar='Choose if Run ImportAPP or Not',
                        choices=['False', 'True'],
                        default='True',
                        help='Default False = Not Delete, put True = Delete')

    parser.add_argument('-backup', '--backup', metavar='Choose if make BackUp or Not', choices=['False', 'True'],
                        default='True',
                        help='Default False = Not Backup, put True = Backup')

    return parser.parse_args()

def run_in_parallel(num_threads=3, **kwargs):
    
    queue = Queue()

    def get_fastq_files_thread(sample_dict, thread_id, thread_num_samples):
        fastq_files = []
        sample_names_list = list(sample_dict.keys())
        # get the samples that will be processed by thread with id=thread_id
        thread_samples = sample_names_list[thread_id * thread_num_samples : min((thread_id + 1) * thread_num_samples, len(sample_names_list))]
        # min can be removed, it is just to make sure that we are not going out of index

        for sample_name in thread_samples:
            fastq_files.append(sample_dict[sample_name]['forward'])
            fastq_files.append(sample_dict[sample_name]['reverse'])

        return fastq_files

    fastq_files = kwargs.pop("fastq_files")
    sample_dict = utils.group_samples(fastq_files)
    print(sample_dict)

    thread_pool = []

    thread_num_samples = math.ceil(len(sample_dict.keys()) / num_threads)
    if (thread_num_samples < 1):
        print("Max number of threads allowed is {}".format(len(sample_dict.keys())))
        sys.exit()

    for thread_id in range(num_threads):
        fastq_files_thread = get_fastq_files_thread(sample_dict, thread_id, thread_num_samples)
        kwargs.update({"fastq_files": fastq_files_thread, "thread_id": thread_id, "queue": queue})

        p = Process(target=pipelineAssembler.factory('process').start, kwargs=kwargs)
        #p = Process(target=dummy_function, kwargs=_kwargs)
        thread_pool.append(p)

    for p in thread_pool:
        p.start()

    for p in thread_pool:
        p.join()

    return queue

    """ 
    The p.join() method ensures that the main program waits for each process in the thread_pool to finish before moving on to the next line of code.

    The threads will still be executed in parallel due to the earlier use of p.start(), which initiates the concurrent execution of each thread. 
    
    The p.join() calls simply block the main program's execution until each respective thread completes.

    the main program will continue execution after the longest-running thread in the thread_pool has finished. 
    
    The use of p.join() doesn't prevent the threads from being executed in parallel; it just ensures that the main program waits for their completion before proceeding. """


def tests():
    args = parseInput()
    # from Pipes.InputPipe import Setup
    # setup = Setup()
    # kwargs = setup.process(**vars(args))
    testPipelineOut = PipelineAssembler().factory('test').start(**vars(args))

if __name__ == "__main__":
    # args = parseInput()
    
    print("Currently unavailable. Coming back soon ... ")

    tests()
    
    # pipelineAssembler = PipelineAssembler()
    # resource_pipeline_out = pipelineAssembler.factory('resource').start(**vars(args))
    # fastq_files = resource_pipeline_out['fastq_files']

    # kwargs = resource_pipeline_out
    # kwargs.update({"fastq_files": fastq_files})


    # start = time.time()
     
    # # TODO: redesign the concurrency, at least for the coverage pipeline. No need to assemble a seperate pipeline, put perform the parallelization inside the CoveragePipe

    # # TODO: run in parallel maybe should return a set of queues, with the different outputs (written files, directories, etc.) that can be used by the next steps of the pipeline
    # # For example, InterstagePipe gets the sample by using glob.glob on the pirnicpal_directory, but it would be better if InterstagePipe reads from the queue
    # # so the actual files are passed directly from firstPipe, to interstagePipe, without having to hardcode the name on the files to search with glob
    # if len(fastq_files) == 0:
    #     print("No fastq files were read from the server!")
    #     sys.exit()
   
    # samples_queue = run_in_parallel(num_threads=int(len(fastq_files) / 2), **kwargs)
    # #run_in_parallel(fastq_files=["test/E378.2023_R1_001.fastq", "test/E378.2023_R2_001.fastq", "test/E379.2023_R1_001.fastq", "test/E379.2023_R2_001.fastq", "test/E380.2023_R1_001.fastq", "test/E380.2023_R2_001.fastq"])
    
    
    # print("\033[92mEntering NVIDIA Parabricks ...")
    # """ fq2bam will run non-concurrently. fq2bam will start after fastx for every sample is completed. """
    # kwargs.update({"resynced_samples": samples_queue})
    # fq2bam_out = pipelineAssembler.factory('fq2bam').start(**kwargs)
    # print("\033[0m")
    
    # kwargs = fq2bam_out
    # coverage_queue = coverage_async.run(**kwargs)
    # while not coverage_queue.empty():
    #     if coverage_queue.empty():
    #         break
    #     print(coverage_queue.get())


