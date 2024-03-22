import argparse
import time
from multiprocessing import Process, Queue
import sys

import config
from PipelineAssembler import PipelineAssembler
import utils
import math

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

    #num_files_thread = int(len(fastq_files) / num_threads)
    # for i in range(num_threads):
    #     _kwargs = kwargs
    #     current_files = fastq_files[i * num_files_thread: (i + 1) * num_files_thread]
    #     _kwargs.update({"fastq_files": current_files, "thread_id": i})

    #     p = Process(target=pipelineAssembler.factory('process').start, kwargs=kwargs)
    #     thread_pool.append(p)

    """ 
    The p.join() method ensures that the main program waits for each process in the thread_pool to finish before moving on to the next line of code.

    The threads will still be executed in parallel due to the earlier use of p.start(), which initiates the concurrent execution of each thread. 
    
    The p.join() calls simply block the main program's execution until each respective thread completes.

    the main program will continue execution after the longest-running thread in the thread_pool has finished. 
    
    The use of p.join() doesn't prevent the threads from being executed in parallel; it just ensures that the main program waits for their completion before proceeding. """


def dummy_function(**kwargs):
    print("Running thread-{} on files {}".format(kwargs["thread_id"], kwargs["fastq_files"]))
    

# TODO: haplotype caller parabricks. merge the two vcfs. merge haplotype caller vcf with deepvariant vcf. then filter according to the bed
def run_coverage_in_parallel(**kwargs):
    queue = Queue()
     
    def get_samples(thread_id, samples, num_threads):
        num_samples_per_thread = math.ceil(len(samples) / num_threads)
        samples_for_thread_id = samples[thread_id * num_samples_per_thread : min((thread_id + 1) * num_samples_per_thread, len(samples))]
        
        return samples_for_thread_id

    samples = kwargs.pop("samples")
    assert(len(samples) > 0)
    num_threads = len(samples)
    thread_pool = []
    
    for thread_id in range(num_threads):
        kwargs.update({"samples": get_samples(thread_id, samples, num_threads), "thread_id": thread_id, "queue": queue})
        p = Process(target=pipelineAssembler.factory('coverage').start, kwargs=kwargs)
        thread_pool.append(p)

    for p in thread_pool:
        p.start()

    for p in thread_pool:
        p.join()
    
    return queue


""" NOTE: principal_directory will be created by the first pipeline, in the setup pipe. This is the place were the files for the current project will be written. 
    However, keeping track of files do not necessarily requires principal_directory, since in each step of the pipeline, each output is saved in the sample dictionary. 
    Therefore, a coupling between principal_directory and the files we are working with is introduced, implying that the pipeline should be executed as a whole. 
    Consider the case where we execute it up to the alignmentpipe. We will have the bam files located inside the bam folder of the principal_directory. Then,
    if we want to execute the second half of the pipeline (Variant Calling), and put them in a different workplace (i.e. different principal_directory), we will not be able to do so,
    given the current code. This is because the current code relies on principal_directory as the place where are files can be found. i.e. in our sample E719.2023 we have saved the alignment files
    inisde 23March/bam/ directory, and they are present in sample["bam"] as 23March/bam/E719.2023_final.bam. If we execute the VariantCalling pipe with sample E719.2023 as an input, and principal
    _directory 24March, instead of 23March, as an input, we will have an error. Therefore, each step of the pipeline should have one common principal_directory. In our case, we would have to run 
    the VariantCall Pipe using 23March as the principal_directory.
"""
if __name__ == "__main__":
    args = parseInput()

    #test_stuff()
    pipelineAssembler = PipelineAssembler()
    resource_pipeline_out = pipelineAssembler.factory('resource').start(**vars(args))
    fastq_files = resource_pipeline_out['fastq_files']
    principal_directory = resource_pipeline_out['principal_directory']

    kwargs = resource_pipeline_out
    kwargs.update({"fastq_files": fastq_files})
  
    start = time.time()
    
    # TODO: run in parallel maybe should return a set of queues, with the different outputs (written files, directories, etc.) that can be used by the next steps of the pipeline
    # For example, InterstagePipe gets the sample by using glob.glob on the pirnicpal_directory, but it would be better if InterstagePipe reads from the queue
    # so the actual files are passed directly from firstPipe, to interstagePipe, without having to hardcode the name on the files to search with glob
    if len(fastq_files) == 0:
        print("No fastq files were read from the server!")
        sys.exit()
   
    samples_queue = run_in_parallel(num_threads=int(len(fastq_files) / 2), **kwargs)
    #run_in_parallel(fastq_files=["test/E378.2023_R1_001.fastq", "test/E378.2023_R2_001.fastq", "test/E379.2023_R1_001.fastq", "test/E379.2023_R2_001.fastq", "test/E380.2023_R1_001.fastq", "test/E380.2023_R2_001.fastq"])
    
    
    print("\033[92mEntering NVIDIA Parabricks ...")
    """ fq2bam will run non-concurrently. fq2bam will start after fastx for each sample is finished"""
    kwargs.update({"resynced_samples": samples_queue})
    fq2bam_out = pipelineAssembler.factory('fq2bam').start(**kwargs)
    print("\033[0m")

    kwargs = fq2bam_out
    queue_of_samples = run_coverage_in_parallel(**kwargs)
    #coverage_out = pipelineAssembler.factory('coverage').start(**kwargs)

    
    
    """ Refactor this. """
    #interstage_pipeline_out = pipelineAssembler.factory('interstage').start(**kwargs)
    end = time.time()
    print("Runtime - {}".format(str(end-start)))

    
    #print("fq2bam_out = {}".format(fq2bam_out))

    """ Get the output of fq2bam"""
    
    
# TODO: sample_list.csv will only contain the last list of samples processed (i.e. the output of the longest thread). We can fix that.