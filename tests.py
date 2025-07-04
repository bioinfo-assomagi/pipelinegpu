import pandas as pd
import os
import config

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
     'sample': {'name': 'E719.2023', 'pairs_forward': 'E719.2023_R1_001_new.fastq_pairs_R1.fastq', 'pairs_reverse': 'E719.2023_R2_001_new.fastq_pairs_R2.fastq', 'bam': '/home/magi/PROJECT/diagnosys/RESULTS_jurgen/27_Mar_2024_CANCER/bam/E719.2023_final.bam', 'bai': '/home/magi/PROJECT/diagnosys/RESULTS_jurgen/27_Mar_2024_CANCER/bam/E719.2023_final.bai'}}


def dummy():
    print(genome)

def ParabricksHaplotypeCaller(samples, principal_directory):

    docker_input_parabricks = os.path.join(config.DOCKER_WORKDIR, 'bam')

    for sample in samples:
        sample_name = str(sample['name'])
        bam_filename = sample['bam'].split('/')[-1]

        os.system(' '.join(['docker', 'run', '--gpus', 'all', '--rm', 
                            '--volume', "{}/:{}".format(config.REF, config.DOCKER_REFDIR), 
                            '--volume', "{}/:{}".format(principal_directory, config.DOCKER_WORKDIR), 
                            '--volume', "{}/:{}".format(os.path.join(principal_directory, "temp"), config.DOCKER_OUTPUTDIR), 
                            "{}".format(config.PARABRICKS_VERSION), 
                        'pbrun', 'haplotypecaller', 
                        '--ref', "{}/Homo_sapiens_assembly38.fasta".format(config.DOCKER_REFDIR), 
                        # '--num-gpus', '1',
                        "--in-bam", os.path.join(docker_input_parabricks, bam_filename),
                        '--out-variants', "{}/{}_pb_gatk.vcf".format(config.DOCKER_OUTPUTDIR, sample_name)]))


def DeepVariant(samples, principal_directory):

    docker_input_parabricks = os.path.join(config.DOCKER_WORKDIR, 'bam')

    for sample in samples:
        sample_name = str(sample['name'])
        bam_filename = sample['bam'].split('/')[-1]

        os.system(' '.join(['docker', 'run', '--gpus', 'all', '--rm', 
                            '--volume', "{}/:{}".format(config.REF, config.DOCKER_REFDIR), 
                            '--volume', "{}/:{}".format(principal_directory, config.DOCKER_WORKDIR), 
                            '--volume', "{}/:{}".format(os.path.join(principal_directory, "temp"), config.DOCKER_OUTPUTDIR), 
                            "{}".format(config.PARABRICKS_VERSION_DEEPVARIANT), 
                        'pbrun', 'deepvariant', 
                        '--ref', "{}/Homo_sapiens_assembly38.fasta".format(config.DOCKER_REFDIR), 
                        "--in-bam", os.path.join(docker_input_parabricks, bam_filename),
                        '--out-variants', "{}/{}_pb_deepvariant.vcf".format(config.DOCKER_OUTPUTDIR, sample_name)]))


def test_variant_call():
    
    from Pipes.VariantCallPipe import VariantCallPipe
    from Pipes.CoveragePipe import CoveragePipe
    
    pipe = VariantCallPipe()
    pipe.process(principal_directory = principal_directory, samples = testdata.samples)


def test_annotation():
    import utils

    seq1 = utils.Sequence("AAAAACTGACCC", "REF")
    seq2 = utils.Sequence("AAAAACCC", "ALT")

    print(seq1.remove_common_prefix(seq2))
    #print(seq1.find_indel(seq2))


def test_spawns():
    import dir_tree
    from Pipes.MultipleDiseaseChildSampleHandlerPipe import MultipleDiseaseSampleHandlerPipe
    dir_tree.build("/home/magi/PROJECT/diagnosys/RESULT/TEST_3_Apr_2025_OCULARE_spawns/", create=False)
    MultipleDiseaseSampleHandlerPipe().process(dummy_arg = "dummy")

if __name__ == "__main__":
    # TODO: haplotype caller parabricks. merge the two vcfs. merge haplotype caller vcf with deepvariant vcf. then filter according to the bed
    test_spawns()


    # genome = 'hey'
    # from test_import import *
    # import testdata
    # dummy()
    # #mg_results['out'].to_csv('mygene_results_test.csv', sep='\t')
    # #gene_DF1 = mg_results['out']
    
    

    # principal_directory = "/home/magi/PROJECT/diagnosys/RESULT/26_Mar_2024_CANCER/"
    # genome_type = 'geno38'
    # panel = 'OCULARE'
    # path = '/home/magi/PROJECT/diagnosys/'
    
    # arr = [1]


    # def filter_for_vcf():
    #     vertical = pd.read_csv("/home/magi/PROJECT/diagnosys/RESULTS_jurgen/07_Mar_2024_OCULARE/pheno/vertical_281.2024", sep="\t")
    #     vertical['filt'] = 1
    #     vcf = pd.read_csv("/home/magi/PROJECT/diagnosys/RESULTS_jurgen/07_Mar_2024_OCULARE/temp/281.2024_gatk.vcf", sep='\t', comment='#',
	# 			names=['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT', "281.2024"])
    #     #print(vcf)
    #     FILTER = pd.merge(vcf, vertical, on=['#CHROM','POS'], how='left')

    #     FILTER_ONLY_CDS = FILTER[FILTER['filt']==1]
    #     FILTER_INTRONIC = FILTER
    #     FILTER_ONLY_CDS.drop('filt', axis=1,inplace=True)
    #     FILTER_INTRONIC.drop('filt', axis=1,inplace=True)
    #     FILTER_ONLY_CDS['control'] = 0

    #     FILTER_INNER = pd.merge(vcf, vertical, on=['#CHROM','POS'], how='inner')
    #     print(FILTER_ONLY_CDS)
    

    # import csv
    # def read_mpileup(sample_name, mpileup_out):
        
    #     a = pd.read_csv(
    #         mpileup_out,
    #         sep="\t",
    #         header=None,
    #         quoting=csv.QUOTE_NONE,
    #         encoding="utf-8",
    #         low_memory=False,
    #         on_bad_lines="skip",
    #         names=["CHROM", "POS", "info", "DEPTH", "CALL", "quality"],
    #         chunksize=40 * 100024,
    #     )

    #     i = 0
    #     for chunk in a:
    #         i += 1
    #         chunk["sample"] = sample_name
    #         print(type(chunk))
            
    #         print(i)
    #filter_for_vcf()

    #ParabricksHaplotypeCaller(samples, principal_directory)
    
    #DeepVariant(samples, principal_directory)       

    #test_variant_call()

  
    #test_annotation()

    #read_mpileup("E530.2019", "/home/magi/PROJECT/diagnosys/RESULTS_jurgen/19_Apr_2024_LYMPHOBESITY/temp/E530.2019_to_count")
    



    # CHECK THE BACKUP.py to find the old vcf.

    """
    FOLDER = principal_directory/pheno/
    lastFOLDER = pirncipal_directory
    
    """