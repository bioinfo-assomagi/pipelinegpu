#!/home/magi/miniconda3/envs/PY310/bin/python
#17/06/2023
# -*- coding: utf-8 -*-
"""
Created on Wed MAY 18 2023
@author: Giuseppe Marceddu
"""

from multiprocessing import Process, Lock
import argparse
import csv
import datetime
import glob
import os
import datetime
import subprocess
import sys
import pandas as pd
import numpy as np
from os import listdir, system
from os.path import isfile, join
import time

pd.options.mode.chained_assignment = None
import first_diagnosys_core_V3_jurgen as step

path = os.getcwd()
geno37 = join('/home/magi/', 'dataset/GENOME/37/all_chr37.fa')
geno38 = join('/home/magi/', 'dataset/GENOME/38/all_chr38.fa')
##########################################################################################
dbsnp144_37 = join('/home/magi/', 'dataset/dbsnp144/37/common_all_20150605_2.vcf')
dbsnp144_38 = join('/home/magi/', 'dataset/dbsnp144/38/common_all_20150603_2.vcf')
indel_37 = join('/home/magi/', 'dataset/dbsnp144/37/common_all_20150605_indel.vcf')
indel_38 = join('/home/magi/', 'dataset/dbsnp144/38/common_all_20150603_indel.vcf')
clinvar = join('/home/magi/', 'dataset/dbsnp144/38/clinvar_20140929_2.vcf')
clinvar_indel = join('/home/magi/', 'dataset/dbsnp144/38/clinvar_20140929_indel.vcf')
#############################################################################################
BUCHIARTIFICIALI = join('/home/magi/', 'PROJECT/diagnosys/bin/BUCHIARTIFICIALI.txt')


##############################################################################################
def InputPar():
    ####Introducing arguments
    parser = argparse.ArgumentParser(prog='MAGI EUREGIO DIAGNOSYS', description='Pipe from BAM to VCF ',
                                     epilog='Need to SAMTOOLS v1.9 and in global path or put your path [--samtools] [--bcftools]')

    parser.add_argument('-p', '--path', metavar='PATH', default=path,
                        help='[Default = Pipe run in same folder when launch command]')

    parser.add_argument('-proj', '--project', metavar='PROJECT NAME',
                        help='Insert Family name to create Folder')

    parser.add_argument('-g', '--genome', metavar='choose Assembling', choices=['geno37', 'geno38'],
                        default='geno38',
                        help='Choices: geno37, geno38 Run Analysis with genome GRCh38 or GRCh37')

    parser.add_argument('-q', '--quality', metavar='QUALITYFILTER', default=18,
                        help=' Choose a quality threshold [Default = 18]')

    parser.add_argument('-N', '--threads', metavar='THREADS', default='32',
                        help='Number of threads to run [8]')

    parser.add_argument('-s', '--samtools', metavar='SAMTOOLS', default='samtools',
                        help='Insert SAMTOOLS path')

    parser.add_argument('-v', '--bcftools', metavar='BCFTOOLS', default='bcftools',
                        help='Insert BCFTOOLS path')

    parser.add_argument('-t', '--samstat', metavar='SAMSTAT', default='samstat',
                        help='Insert SAMSTAT path')

    parser.add_argument('-gk', '--gatk', metavar='GATK', default='gatk', help='Insert gatk path')

    parser.add_argument('-ovr', '--over', metavar='New/OLD project', choices=['True', 'False'],
                        default='False',
                        help='Choices: Attention option False overwrite old data [Default = False] ')

    parser.add_argument('-d', '--dest', metavar='destination', choices=['b', 'r', 's', 'z', 'p'],
                        required=True,
                        help='Choices destination: b = bolzano; r = rovereto; s = sanfelice; z = ricerca; p = privato,  - required ')

    parser.add_argument('-pan', '--panel', metavar='PANEL', required=True,
                        help='Pannelli Consigliati: 1:CLM2; 2:OSSNEW; 3:OCULARE; 4:OCULARE2; 5:SIFSR; 6:ALLGENE;'
                             '7:trusightone; 8:GENETEST1; 9:GENETEST2; 10:CANCER; 11:MALFORMATIONI VASCOLARI;'
                             '12:ANOMALIE VASCOLARI; 13:INFERTILITA; 14:OBESITA; 15:GENODERMATOSI\nRequired')

    #	parser.add_argument('-o','--name', metavar='Output ', required=True,help='Choose an output filename (Required)')
    return parser.parse_args()


class Writer:

    def __init__(self, stdout, filename):
        self.stdout = stdout
        self.logfile = open(filename, 'a')

    def write(self, text):
        self.stdout.write(text)
        self.logfile.write(text)

    def close(self):
        self.stdout.close()
        self.logfile.close()


def PrintLog(command, folder):
    ###Print all commands in a log file
    path = join(folder, 'Commands2.log')
    cl = open(path, 'a')
    cl.write(command)
    cl.write('\n')
    cl.close()


##########################################################################
def vcf_from_samtools(param, genome, vertical, verticalX, folder, bam, phenotype, vertical_macro):
    print('Trovato:', bam)
    a = bam.split("/")
    if '-' in a[-1]:
        x = a[-1].split("-")
        x = (x[0] + '.' + x[1])
    else:
        x = a[-1].split(".")
        x = (x[0] + '.' + x[1])

    x1 = x.split("_")
    sample = x1[0]
    print('-------->' + sample + '<------------')
    file_vcf = join(folder, 'temp/', sample + '_samt.vcf')

    print(' '.join(['samtools', 'mpileup', '-C50', '-q', str(param.quality),
                    '-Q', str(param.quality), '-d 100000', '-L 100000', '-t DP', '-Buf',
                    genome, bam, '|', 'bcftools', 'call', '-O', 'v', '-mv', '>', file_vcf]))

    # TODO: substitue with parabricks
    system(' '.join(['samtools', 'mpileup', '-C50', '-q', str(param.quality),
                     '-Q', str(param.quality), '-d 100000', '-L 100000', '-t DP', '-Buf',
                     genome, bam, '|', 'bcftools', 'call', '-O', 'v', '-mv', '>', file_vcf]))

    print('load vcf!!!')
    vcf = pd.read_csv(file_vcf, sep='\t', comment='#',
                      names=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', sample])
    ########################################################################
    vcf_ = filter_for_vcf(verticalX, folder, vcf, sample, buchiartificiali)  # TODO: deepvariant
    ########################################################################
    print('------------------------')
    return str(sample)


########################################################################
########################################################################
def filter_for_vcf(vertical, folder, VCF, sample, buchiartificiali):
    vertical['filt'] = 1
    name = join(folder, 'vcf/', sample)

    FILTER = pd.merge(VCF, vertical, on=['#CHROM', 'POS'], how='left')

    FILTER_ONLY_CDS = FILTER[FILTER['filt'] == 1]
    FILTER_INTRONIC = FILTER
    FILTER_ONLY_CDS.drop('filt', axis=1, inplace=True)
    FILTER_INTRONIC.drop('filt', axis=1, inplace=True)

    FILTER_ONLY_CDS['control'] = 0
    for index, row in buchiartificiali.iterrows():
        x = row['#CHROM']
        y = row['START']
        z = row['END']
        mask1 = ((FILTER_ONLY_CDS['#CHROM'] == x) & (FILTER_ONLY_CDS['POS'] >= y) & (FILTER_ONLY_CDS['POS'] <= z))
        FILTER_ONLY_CDS.loc[mask1, 'control'] = 1
    FILTER_ONLY_CDS2 = FILTER_ONLY_CDS[FILTER_ONLY_CDS['control'] == 0]
    FILTER_ONLY_CDS2.drop('control', axis=1, inplace=True)

    FILTER_INTRONIC.to_csv(name + '_samt_all.vcf', sep='\t', index=False)
    FILTER_ONLY_CDS2.to_csv(name + '_samt_only_CDS.vcf', sep='\t', index=False)
    return FILTER_ONLY_CDS


#####################################################################################
def GATK_unfied(param, folder, bam, dbsnp, vertical, buchiartificiali):
    print('Trovato:', bam)
    a = bam.split("/")
    if "-" in a[-1]:
        x = a[-1].split("-")
        x = (x[0] + '.' + x[1])
    else:
        x = a[-1].split(".")
        x = (x[0] + '.' + x[1])

    x1 = x.split("_")
    sample = str(x1[0])
    print('-------->' + str(sample) + '<------------')

    file_vcf_gatk = join(folder, 'temp/', sample + '_unfied.vcf')
    # TODO: haplotype caller parabricks. merge the two vcfs. merge haplotype caller vcf with deepvariant vcf. then filter according to the bed
    system(' '.join([param.gatk, '--java-options -Xmx32G HaplotypeCaller', '-I', bam, '-O', file_vcf_gatk,
                     '-R', genome, '--native-pair-hmm-threads', str(param.threads), '--minimum-mapping-quality',
                     str(param.quality)]))
    VCF = pd.read_csv(file_vcf_gatk, sep='\t', comment='#',
                      names=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', str(sample)])

    vertical['filt'] = 1
    name = join(folder, 'vcf/', sample)

    print('-------FILTERING---GATK----')
    FILTER = pd.merge(VCF, vertical, on=['#CHROM', 'POS'], how='left')

    FILTER_ONLY_CDS = FILTER[FILTER['filt'] == 1]
    FILTER_INTRONIC = FILTER
    FILTER_ONLY_CDS.drop('filt', axis=1, inplace=True)
    FILTER_INTRONIC.drop('filt', axis=1, inplace=True)
    FILTER_ONLY_CDS['control'] = 0
    for index, row in buchiartificiali.iterrows():
        x = row['#CHROM']
        y = row['START']
        z = row['END']
        mask1 = ((FILTER_ONLY_CDS['#CHROM'] == x) & (FILTER_ONLY_CDS['POS'] >= y) & (FILTER_ONLY_CDS['POS'] <= z))
        FILTER_ONLY_CDS.loc[mask1, 'control'] = 1
    FILTER_ONLY_CDS2 = FILTER_ONLY_CDS[FILTER_ONLY_CDS['control'] == 0]
    FILTER_ONLY_CDS2.drop('control', axis=1, inplace=True)

    FILTER_INTRONIC.to_csv(name + '_unfied_all.vcf', sep='\t', index=False)
    FILTER_ONLY_CDS2.to_csv(name + '_unfied_only_CDS.vcf', sep='\t', index=False)
    return FILTER_ONLY_CDS


########################################################################
#########################################################################
def plot(param, folder, sample):
    print('in plot')
    x = sample.split("_")
    sample = str(x[0])
    sample_sam_ = join(folder, 'temp/', sample + '.sam')
    sample_bam_ = join(folder, 'bam/', sample + '_final.bam')
    vcf_sample_ = join(folder, 'temp/', sample + '_samt.vcf')
    vcf_sample_unfied = join(folder, 'temp/', sample + '_unfied.vcf')
    fastqc_folder = join(folder, 'fastQC')
    if param.dest == 'r':
        dest = 'rovereto'
        django_folder = '/home/magi/VIRTUAL/MAGIS/MAGIS/templates/ngs/stat/%s/' % (sample)
        if not os.path.exists(django_folder):
            os.makedirs(django_folder)
    elif param.dest == 'b':
        dest = 'bolzano'
        django_folder = '/home/magi/VIRTUAL/EUREGIO/EUREGIO/templates/ngs/stat/%s/' % (sample)
        if not os.path.exists(django_folder):
            os.makedirs(django_folder)
    elif param.dest == 's':
        dest = 'sanfelice'
        django_folder = '/home/magi/VIRTUAL/SANFELICE/SANFELICE/templates/ngs/stat/%s/' % (sample)
        if not os.path.exists(django_folder):
            os.makedirs(django_folder)
    elif param.dest == 'z':
        dest = 'ricerca'
        django_folder = '/home/magi/VIRTUAL/RICERCA/RICERCA/templates/ngs/stat/%s/' % (sample)
        if not os.path.exists(django_folder):
            os.makedirs(django_folder)

    if os.path.exists(join(fastqc_folder, sample + '_R1_001_fastqc.html')):
        system(' '.join(['cp', join(fastqc_folder, sample + '*_R1_001_fastqc.html'), django_folder]))
    if os.path.exists(join(fastqc_folder, sample + '_R2_001_fastqc.html')):
        system(' '.join(['cp', join(fastqc_folder, sample + '_R2_001_fastqc.html'), django_folder]))
    else:
        print('FASTQC not found!!!')
    return


##############################################################################
##############################################################################
if __name__ == "__main__":
    lock = Lock()
    args = InputPar()

    if args.over == 'False':
        folder_name = step.principal_folder(args, step.create_folder(), over='False')
    if args.over == 'True':
        folder_name = step.principal_folder(args, step.create_folder(), over='True')
        step.path_creation(args, folder_name)
        print('NOW YOU HAVE TO COPY BAM FILES IN FOLDER (BAM)\nAND RUN COMMAND AGAIN\nWHITOUT (--OVER)')
        sys.exit

    folder_pheno = join(folder_name, 'pheno')
    input_phenotype = join(folder_name, 'pheno/phenotype')
    print('here')

    print_args = vars(args)  #take args in a dict and print in a log
    for k, v in print_args.items():
        strategy_out = '='.join([str(k), str(v)])
        path = join(folder_name, 'log')
        PrintLog(strategy_out, path)
    #####################################################################################
    ####################################################################################
    if args.genome == 'geno38':
        genome = geno38
        dbsnp = dbsnp144_38
    elif args.genome == 'geno37':
        genome = geno37
        dbsnp = dbsnp144_37
    path_ = join(folder_name, 'bam/*.bam')
    path_bai = join(folder_name, 'bam/*.bai')
    print(path_)
    files_bai = glob.glob(path_bai)
    for bai2 in files_bai:
        a = bai2.split("/")
        if '-' in a[-1]:
            x = a[-1].split("-")
            x = (x[0] + '.' + x[1])
        else:
            x = a[-1].split(".")
            x = (x[0] + '.' + x[1])
        x1 = x.split("_")
        name = x1[0]
        bai = join(folder_name, 'bam/', name + '_final.bai')
        os.rename(bai2, bai)

    files_bam = glob.glob(path_)
    for bam2 in files_bam:
        a = bam2.split("/")
        if '-' in a[-1]:
            x = a[-1].split("-")
            x = (x[0] + '.' + x[1])
        else:
            x = a[-1].split(".")
            x = (x[0] + '.' + x[1])

        x1 = x.split("_")
        name = x1[0]
        bam = join(folder_name, 'bam/', name + '_final.bam')
        os.rename(bam2, bam)
        print(name)
        log_file2 = join(folder_name, 'log', name + '_vcf.log')
        writer2 = Writer(sys.stdout, log_file2)
        sys.stdout = writer2
        phenotype = pd.read_csv(input_phenotype, sep='\t', header=0)  #,encoding='utf-8')
        phenotype['sample'] = phenotype['sample'].astype('str')
        phenotype['sample'].replace(r'\.202$', '.2020', inplace=True, regex=True)
        phenotype['sample'] = phenotype['sample'].astype('str')
        GENELIST = []
        GENELIST = list(phenotype['gene'][phenotype['sample'].astype(str) == str(name)])
        print('NAME:', name)
        ####################################################################################
        sospetto = phenotype['malattia'][phenotype['sample'].astype(str) == str(name)].unique()[0]
        print(sospetto)

        bed_ = join(folder_pheno, 'bed_' + str(name))
        bed = pd.read_csv(bed_, sep='\t')
        bedx_ = join(folder_pheno, 'bedX_' + str(name))
        bedx = pd.read_csv(bedx_, sep='\t')
        vertical_ = join(folder_pheno, 'vertical_' + str(name))
        vertical = pd.read_csv(vertical_, sep='\t')
        verticalx_ = join(folder_pheno, 'verticalX_' + str(name))
        verticalX = pd.read_csv(verticalx_, sep='\t')

        if not vertical.empty:
            print('ifnot')
            buchiartificiali = pd.read_csv(BUCHIARTIFICIALI, sep='\t', header=0)
            if 'OCULAR' in args.panel:
                _vertical_macro = join(args.path, 'bin', 'VERTICAL', 'vertical_ONA20509.2020')
            elif 'VASCULAR' in args.panel:
                _vertical_macro = join(args.path, 'bin', 'VERTICAL', 'vertical_VASCNA20509.2020')
            elif 'NEUROLOGY' in args.panel:
                _vertical_macro = join(args.path, 'bin', 'VERTICAL', 'vertical_NEURNA20509.2020')
            elif 'INTEGRACARDIOSTANCHEZZA' in args.panel:
                _vertical_macro = join(args.path, 'bin', 'VERTICAL', 'vertical_INTEGRACARDIOSTANCHEZZA.2023')
            elif 'MIXED' in args.panel:
                _vertical_macro = join(args.path, 'bin', 'VERTICAL', 'vertical_MIX120509.2020')
            elif 'LYMPHOBESITY' in args.panel:
                _vertical_macro = join(args.path, 'bin', 'VERTICAL', 'vertical_LIMPHNA20509.2020')
            elif 'INFERTILIT' in args.panel:
                _vertical_macro = join(args.path, 'bin', 'VERTICAL', 'vertical_INA20509.2020')
            elif 'CHERATOCONO' in args.panel:
                _vertical_macro = join(args.path, 'bin', 'VERTICAL', 'vertical_CHERATOCONO.2020')
            elif 'PCDH19' in args.panel:
                _vertical_macro = join(args.path, 'bin', 'VERTICAL', 'vertical_PCDH19')
            elif 'GENEOBNEW' in args.panel:
                _vertical_macro = join(args.path, 'bin', 'VERTICAL', 'vertical_GENEOBNEW')
            elif 'CUSTOM' in args.panel:
                _vertical_macro = join(folder_pheno, 'vertical_' + name)
            elif 'CANCER' in args.panel:
                print('CANCERRRRRRRRRRRRRRRRRRRRRRRR', join(args.path, 'bin', 'VERTICAL', 'vertical_CANC20509.2022'))
                _vertical_macro = join(args.path, 'bin', 'VERTICAL', 'vertical_CANC20509.2022')
            else:
                _vertical_macro = ''  #'' #print('VERTICAL not found')

            try:
                vertical_macro = pd.read_csv(_vertical_macro, sep='\t')
            except:
                print('ERRORE!!!')  #verti
            sample = vcf_from_samtools(args, genome, vertical, verticalX, folder_name, bam, phenotype, vertical_macro)
            name = join(folder_name, 'temp/', sample + '_final_disease')
            GATK_unfied(args, folder_name, bam, dbsnp, verticalX, buchiartificiali)
            plot(args, folder_name, sample)
        else:
            pass
