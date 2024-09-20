import argparse
import time
from multiprocessing import Process, Queue
import sys

import config as config
from PipelineAssembler import PipelineAssembler
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
    
    # parser.add_argument('-type', '--pipelinetype', metabar='What type of pipeline do you want to run?')

    return parser.parse_args()

    """ 
    The p.join() method ensures that the main program waits for each process in the thread_pool to finish before moving on to the next line of code.

    The threads will still be executed in parallel due to the earlier use of p.start(), which initiates the concurrent execution of each thread. 
    
    The p.join() calls simply block the main program's execution until each respective thread completes.

    the main program will continue execution after the longest-running thread in the thread_pool has finished. 
    
    The use of p.join() doesn't prevent the threads from being executed in parallel; it just ensures that the main program waits for their completion before proceeding. """


def tests():
    from Pipeline import Pipeline
    from Pipes.DiagnosysCorePipe import CoverageWrapperPipe
    from Pipes.VariantFilterWrapperPipe import VariantFilterWrapperPipe
    from Pipes.VariantCallPipe import VariantCallPipe
    from Pipes.InputPipe import Setup
    from Pipes.DiagnosysCorePipe import ResyncDBPipe
    from Pipes.InterStage import SampleListFam

    args = parseInput()
    # from Pipes.InputPipe import Setup
    # setup = Setup()
    # kwargs = setup.process(**vars(args))

    pipeline_type = config.PIPELINETYPE

    if pipeline_type == 'setup':
        print("Im in setup...")
        Pipeline(Setup()).start(**vars(args))
    elif pipeline_type == 'bamstart':
        PipelineAssembler().factory('bamstart').start(**vars(args))
    elif pipeline_type == 'test':
        PipelineAssembler().factory('test').start(**vars(args))
    elif pipeline_type == 'coverage':
        print("Im in coverage ...")
        Pipeline(Setup()).assemblePipe(ResyncDBPipe()).assemblePipe(SampleListFam()).assemblePipe(CoverageWrapperPipe()).start(**vars(args))
        
    #Pipeline(Setup()).assemblePipe(VariantCallPipe()).start(**vars(args))
    #Pipeline(Setup()).assemblePipe(VariantFilterWrapperPipe()).start(**vars(args))

    #processPipelineOut = PipelineAssembler().factory('process').start(**vars(args))
    

if __name__ == "__main__":
    # args = parseInput()
    
    print("Currently unavailable. Coming back soon ... ")

    tests()
    
    


