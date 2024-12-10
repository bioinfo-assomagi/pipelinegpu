import multiprocessing
from PipelineAssembler import PipelineAssembler
from AsyncRunner import AsyncRunner

def worker(kwargs):
    pipelineAssembler = PipelineAssembler()
    # pipelineAssembler.factory('coverage').start(**kwargs)
    pipelineAssembler.independent_pipeline('coverage').start(**kwargs)

def worker2(sample, kwargs):
    print('processing: ', sample['name'])
    kwargs.update({"sample": sample})
    pipelineAssembler = PipelineAssembler()
    pipelineAssembler.factory('coverage').start(**kwargs)

def run(**kwargs):
    q = AsyncRunner(**kwargs).run(worker)
    return q
    

if __name__ == "__main__":
    fq2bam_out_test = {'quality': 18, 'quality_perc': 97, 'threads': '12', 'memory': '15G', 'bwa': 'bwa', 
     'samtools': 'samtools', 'samstat': 'samstat', 'gatk': '/home/magi/tools/gatk-4.2.0.0/gatk', 
     'bcftools': 'bcftools', 'VEP': '/home/magi/tools/ensembl-vep-release-110/vep', 
     'delete': 'True', 'runfrequency': 'True', 'runcore': 'True', 
     'importAPP': 'True', 'backup': 'True', 'proj_name': '27_Mar_2024', 
     'panel': 'CANCER', 'path': '/home/magi/PROJECT/diagnosys/', 'genome': 'geno38', 'over': True, 
     'fastq_folder': '/home/magi/PROJECT/diagnosys/RESULTS_jurgen/27_Mar_2024_CANCER/fastq/', 
     'fastq': '../RESULT/14_Feb_2024_CANCER/fastq/', 'dest': 'b', 
     'fastq_files': ['/home/magi/PROJECT/diagnosys/RESULTS_jurgen/27_Mar_2024_CANCER/fastq/E719.2023_R2_001.fastq.gz', '/home/magi/PROJECT/diagnosys/RESULTS_jurgen/27_Mar_2024_CANCER/fastq/E719.2023_R1_001.fastq.gz'],
     'db_path': '/home/magi/PROJECT/diagnosys/bin/DATABASE/euregio.db', 'principal_directory': '/home/magi/PROJECT/diagnosys/RESULTS_jurgen/27_Mar_2024_CANCER', 
     'samples': [{'name': 'E719.2023', 'pairs_forward': 'E719.2023_R1_001_new.fastq_pairs_R1.fastq', 'pairs_reverse': 'E719.2023_R2_001_new.fastq_pairs_R2.fastq', 'bam': '/home/magi/PROJECT/diagnosys/RESULTS_jurgen/27_Mar_2024_CANCER/bam/E719.2023_final.bam', 'bai': '/home/magi/PROJECT/diagnosys/RESULTS_jurgen/27_Mar_2024_CANCER/bam/E719.2023_final.bai'}]}
    
    kwargs = fq2bam_out_test
    run(**kwargs)

