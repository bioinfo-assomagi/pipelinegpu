#!/home/magi/miniconda3/envs/PY310/bin/python
#27/06/2015

# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 2023

@author: Giuseppe Marceddu
"""

from multiprocessing import Process, Lock
import multiprocessing
import argparse
import re
import csv
import datetime
import glob
import os
import subprocess
import sys
import pandas as pd
#import sys, csv, re
from os import listdir, system
from os.path import isfile, join
import time
import connect_database

############ PERCORSI FILE (DA MODIFICARE)###############################################
path = os.getcwd()
path_data= path+'/FASTQ/'

geno37 = join('/home/magi/','dataset/GENOME/37/all_chr37.fa')
geno38 = join('/home/magi/','dataset/GENOME/38/all_chr38.fa')

dbsnp144_37 = join('/home/magi/','dataset/dbsnp144/37/common_all_20150605_2.vcf')
dbsnp144_38 = join('/home/magi/','dataset/dbsnp144/38/common_all_20150603_2.vcf')

#indel_37 = join('/home/magi/','dataset/old_dbsnp/Mills_and_1000G_gold_standard.indels.b37.sites.vcf')
indel_37 = join('/home/magi/','dataset/dbsnp144/37/common_all_20150605_indel.vcf')
indel_38 = join('/home/magi/','dataset/dbsnp144/38/common_all_20150603_indel.vcf')

clinvar = join('/home/magi/','dataset/dbsnp144/38/clinvar_20140929_2.vcf')
clinvar_indel = join('/home/magi/','dataset/dbsnp144/38/clinvar_20140929_indel.vcf')

#GATK = join('/home/magi/','tools/GenomeAnalysisTK36.jar')
# GATK = join('/home/magi/','tools/gatk-4.1.4.0/gatk')
GATK = join('/home/magi/','tools/gatk-4.2.0.0/gatk')

def InputPar():
	####Introducing arguments
	parser = argparse.ArgumentParser(prog='MAGI EUREGIO DIAGNOSYS',description='Pipe from FASTQ to BAM',
			epilog='Need to BWA 0.7.12-r1039 and SAMTOOLS v1.3 (SET on global PATH or INSERT your path [--bwa] [--samtools]')

	parser.add_argument('-p','--path', metavar='PATH', default=path,
			help='[Default = Pipe run in same folder when launch command]')

	parser.add_argument('-fq','--fastq', metavar='FASTQ',
			help='FULL PATH to find FASTQ FILE - [default FROM NAS]')

	parser.add_argument('-d','--dest', metavar='destination', choices=['b','r','s','z','p'],
	                    required=True,
	                    help='Choices destination: b = bolzano; r = rovereto; s = sanfelice; z = ricerca; p = privato,  - required ')

	parser.add_argument('-pan','--panel',metavar='PANEL',required=True,
				help='Pannelli Consigliati: 1:CLM2; 2:OSSNEW; 3:OCULARE; 4:OCULARE2; 5:SIFSR; 6:ALLGENE;'
					'7:trusightone; 8:GENETEST1; 9:GENETEST2; 10:CANCER; 11:MALFORMATIONI VASCOLARI;'
					'12:ANOMALIE VASCOLARI; 13:INFERTILITA; 14:OBESITA; 15:GENODERMATOSI\nRequired')
	parser.add_argument('-proj','--project',metavar='PROJECT NAME', help='Insert Family name to create Folder')
	parser.add_argument('-g','--genome', metavar='choose Assembling', choices=['geno37','geno38'],
			default='geno38',
			help='Choices: geno37, geno38 - Run Analysis with genome GRCh38 or GRCh37')
	parser.add_argument('-q','--quality',  metavar='PHRED QUALITY FILTER',default=20,
			help=' Choose a quality threshold [Default = 20]')
	parser.add_argument('-Q','--quality_perc',  metavar='PERC QUALITY READ',default=97,
			help='Percentual of quality for read [Default = 97]')

	parser.add_argument('-N','--threads', metavar='THREADS', default='8',
			help='Number of threads to run [20]')

	parser.add_argument('-m','--memory', metavar='MERORY', default='10G',
			help='max memory per thread; suffix K/M/G [defaul = 12G]')

	parser.add_argument('-b','--bwa', metavar='BWA', default='bwa',
			help='Insert BWA path')

	parser.add_argument('-s','--samtools', metavar='SAMTOOLS', default='samtools',
			help='Insert SAMTOOLS path')

	parser.add_argument('-gk','--gatk', metavar='GATK', default=GATK, help='Insert gatk path')

	parser.add_argument('-ovr','--over', metavar='New/OLD project', choices=['True','False'],
			default='True',
			help='Choices: ALLERT!!! Set on "False" overwrite old data [Default = True]')

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

#import subprocess
#proc = subprocess.Popen(["cat", "/etc/services"], stdout=subprocess.PIPE, shell=True)
#(out, err) = proc.communicate()
#print "program output:", out
def PrintLog(command,folder):
	###Print all commands in a log file
	path = join(folder,'Commands.log')
	cl = open(path,'a')
	cl.write(command)
	cl.write('\n')
	cl.close()

def create_folder():
	today = datetime.date.today()
	name = "{:%d_%b_%Y}".format(today)
	return name

def principal_folder(param, name, over=None):

	if over == 'True':

		if param.panel == 'trusightone':
			panel = param.panel
		else:
			panel = param.panel.upper()

		if param.project: name_folder = join(param.path,'RESULT',param.project+'_'+panel)
		else: name_folder = join(param.path,'RESULT',name+'_'+panel)

		try:
			os.makedirs(name_folder)
		except OSError:
			sys.exit("Folder exist just - Choose another project name, please!!! You can run with <-proj> option")
		finally:
			print ('-----------------------------------------')
		return name_folder

	elif over == 'False':

		if param.panel == 'trusightone':
			panel = param.panel
		else:
			panel = param.panel.upper()

		if param.project: name_folder = join(param.path,'RESULT',param.project+'_'+panel)
		else: name_folder = join(param.path,'RESULT',name+'_'+panel)

		return name_folder


def path_creation(param,folder):
	print ('built tree folder...')
	spec_fastq = join(folder,'fastq/')
	spec_fastqfilter = join(folder,'fastqfiltered/')
	spec_temp = join(folder,'temp/')
	spec_control = join(folder,'control/')
	spec_vcf = join(folder,'vcf/')
	spec_bam = join(folder,'bam/')
	spec_baminanalysis = join(folder,'bam/inanalysis')
	spec_plot = join(folder,'plot/')
	spec_annot = join(folder,'annotation/')
	spec_pheno = join(folder,'pheno/')
	spec_cov = join(folder,'coverage/')
	folder_report = join(folder,'report')
#	spec_cov_fig = join(folder,'coverage/figure')
	spec_qc = join(folder,'fastQC/')
	spec_final = join(folder,'final/')
	spec_indel = join(folder,'indel/')
	#spec_json = join(folder,'json/')
	#spec_json_stat = join(folder,'json/','cov_stat')
	#spec_json_annot = join(folder,'json/','annot')
	spec_log = join(folder,'log')
	if not os.path.exists(spec_qc):
		os.makedirs(spec_qc)
	if not os.path.exists(spec_temp):
		os.makedirs(spec_temp)
	if not os.path.exists(spec_vcf):
		os.makedirs(spec_vcf)
	if not os.path.exists(spec_indel):
		os.makedirs(spec_indel)
	if not os.path.exists(spec_bam):
		os.makedirs(spec_bam)
	if not os.path.exists(spec_baminanalysis):
		os.makedirs(spec_baminanalysis)
	if not os.path.exists(spec_plot):
		os.makedirs(spec_plot)
	if not os.path.exists(spec_fastq):
		os.makedirs(spec_fastq)
	if not os.path.exists(spec_fastqfilter):
		os.makedirs(spec_fastqfilter)
	if not os.path.exists(spec_annot):
		os.makedirs(spec_annot)
	if not os.path.exists(spec_pheno):
		os.makedirs(spec_pheno)
	if not os.path.exists(spec_cov):
		os.makedirs(spec_cov)
	if not os.path.exists(spec_final):
		os.makedirs(spec_final)
	if not os.path.exists(spec_control):
		os.makedirs(spec_control)
	#if not os.path.exists(spec_json):
	#	os.makedirs(spec_json)
	#if not os.path.exists(spec_json_stat):
	#	os.makedirs(spec_json_stat)
	#if not os.path.exists(spec_json_annot):
	#	os.makedirs(spec_json_annot)
	if not os.path.exists(spec_log):
		os.makedirs(spec_log)
#	if not os.path.exists(spec_cov_fig):
#		os.makedirs(spec_cov_fig)
	if not os.path.exists(folder_report):
		os.makedirs(folder_report)
	return

def copy_fastq(param,folder):
	print ('copy files...')
	fastq_folder = join(folder,'fastq/')
	pheno_folder = join(folder,'pheno/')
	fastqc_folder = join(folder,'fastQC/')

	if param.fastq:
		files_fq = glob.glob(join(param.fastq,'*'))
		for files in files_fq:
			if '.fastq.gz' in files:
				system(' '.join(['cp',files,fastq_folder]))
			elif '.fq.gz' in files:
				system(' '.join(['cp',files,fastq_folder]))
			elif '.fastaq.gz' in files:
				system(' '.join(['cp',files,fastq_folder]))
			elif '.fastq' in files:
				system(' '.join(['cp',files,fastq_folder]))
			elif '.fastaq' in files:
				system(' '.join(['cp',files,fastq_folder]))
			elif '.fq' in files:
				system(' '.join(['cp',files,fastq_folder]))
			else:
				if 'phenotype' in files: 'Copying phenotype file...ever is better check inside this file....'
				else: print (files,'is not a VALID file and will not copy')
	else:
		# files_fq = system(' '.join(['scp server@192.168.1.201:/media/4e955bfb-88f0-4cc7-a824-27ee0b4bf6e2/NGS/tmp/analysis/germinal/*',fastq_folder]))
		if param.dest == 'r':
			#files_fq = system(' '.join(['scp root@192.168.2.188:/sharedfolders/NGS/tmp/analysis/germinal/*',fastq_folder])) #added
			files_fq = system(' '.join(['scp root@192.168.1.51:/home/NGS/tmp/analysis/germinal/*',fastq_folder])) #added
		else:
			files_fq = system(' '.join(['scp root@192.168.1.51:/home/NGS/tmp/analysis/germinal/*',fastq_folder])) #added
		#files_fq = system(' '.join(['scp server@192.168.1.201:/media/4e955bfb-88f0-4cc7-a824-27ee0b4bf6e2/NGS/tmp/analysis/germinal2/*',fastq_folder]))
		#files_fq = system(' '.join(['scp server@192.168.1.201:/media/4e955bfb-88f0-4cc7-a824-27ee0b4bf6e2/NGS/tmp/analysis/somatic1/*',fastq_folder]))
	destin_file = glob.glob(fastq_folder+'*')

	for files2 in destin_file:
		a = files2.split("/")
		if '-' in a[-1]:
			x = a[-1].split("-")
			x = (x[0]+'.'+x[1])
		else:
			x = a[-1].split(".")
			x = (x[0]+'.'+x[1])
		x1 = x.split("_")
		name = x1[0]
		files3 = files2.replace('-','.')
		files = re.sub(r"(_S\d+_L\d+|_L\d+)",'', files3)
		print (files)
		os.rename(files2,files)
		if 'new' not in files:
			if '.gz' in files:
				print ('unzip files...')
				system(' '.join(['gunzip',files]))

	list_file = glob.glob(fastq_folder+'*')
	for files in list_file:
		system(' '.join(['fastqc',files,'-t',param.threads,files,'-o',fastqc_folder]))
		if 'new' not in files:
			print ('####################################')
			print (files)
			print ('####################################')

			system(' '.join(['fastq_quality_filter -v -Q33 -q',str(param.quality),
					 '-p',str(param.quality_perc),'-i',files,'|','fastx_trimmer -Q33 -t 5 -m 20 -v',
					 '-o',join(files[:-6]+'_new.fastq')]))

			PrintLog(' '.join(['fastq_quality_filter -v -Q33 -q',str(param.quality),
					 '-p',str(param.quality_perc),'-i',files,'|','fastx_trimmer -Q33 -t 5 -m 20 -v',
					 '-o',join(files[:-6]+'_new.fastq')]),join(folder,'log'))

			# if '_R1_' in files:
			#  	print ('FORWARDDDDDDDD!!!')
			#  	system(' '.join(['fastx_trimmer -Q33 -t 5 -m 20 -v','-i',files,'|',
			# 			'fastq_quality_filter -v -Q33 -q',str(param.quality),'-p',str(97),
			#  			'-o',join(files[:-6]+'_new.fastq')]))
			#
			#  	PrintLog(' '.join(['fastx_trimmer -Q33 -t 5 -m 20 -v','-i',files,'|',
			# 			'fastq_quality_filter -v -Q33 -q',str(param.quality),'-p',str(97),
			#  			'-o',join(files[:-6]+'_new.fastq')]),join(folder,'log'))
			#
			# elif '_R2_' in files:
			#  	print ('REVERSEEEEEEE!!!')
			#  	system(' '.join(['fastx_trimmer -Q33 -t 5 -m 20 -v','-i',files,'|',
			# 			'fastq_quality_filter -v -Q33 -q',str(param.quality),'-p',str(97),
			#  			'-o',join(files[:-6]+'_new.fastq')]))
			#
			#  	PrintLog(' '.join(['fastx_trimmer -Q33 -t 5 -m 20 -v','-i',files,'|',
			# 			'fastq_quality_filter -v -Q33 -q',str(param.quality),'-p',str(97),
			#  			'-o',join(files[:-6]+'_new.fastq')]),join(folder,'log'))

	list_file = glob.glob(fastq_folder+'*_new.fastq')
	return list_file

def set_samples(param,samples,folder):
	forward = []
	reverse = []
	name = []
	for sample in samples:
		a = sample.split("/")
		x = a[-1]
		x1 = x.split("_")
		name.append(str(x1[0]))
		strand = x1[1]
		strand2 = x1[-2]
		strand3 = x1[-3]
		if strand == 'R1':
			forward.append(sample)
		elif strand2 == 'R1':
			forward.append(sample)
		elif strand3 == 'R1':
			forward.append(sample)
		elif strand == 'R2':
			reverse.append(sample)
		elif strand2 == 'R2':
			reverse.append(sample)
		elif strand3 == 'R2':
			reverse.append(sample)

	name_ID = list(set(name))
	print (name_ID)
	pheno_folder = join(folder,'pheno','phenotype')
	dest = param.dest
	print (dest)
	pheno = connect_database.get_disease(name_ID,dest)
	print (pheno)
	pheno.to_csv(pheno_folder,sep='\t',index=False,encoding='utf-8')

	try:
		date_frw = pd.DataFrame({'forward':pd.Series(forward)})
		date_rev = pd.DataFrame({'reverse':pd.Series(reverse)})
		date_frw['name'] = date_frw['forward'].str.split('/').str.get(-1).str.split('_').str.get(0)
		date_rev['name'] = date_rev['reverse'].str.split('/').str.get(-1).str.split('_').str.get(0)
		date_df = pd.merge(date_frw,date_rev,on=['name'],how='outer')
		date_df = date_df[['name','forward','reverse']]
		date_df.sort_values(by=['name'],inplace=True)
		date_df.to_csv(join(folder,'sample_list.csv'),sep='\t',index=False,encoding='utf-8')
		return date_df
	except AttributeError:
		print ('ATTENTION!!!', args.fastq)
		return sys.exit

### Allineament #####param.name
def PreAlignment(lock,param,folder,name,forward,reverse):
#	lock.acquire()
	name_process = multiprocessing.current_process().name
	print (name)
	try:
		print (name_process, 'Starting PRE_ALIGNMENT...', name)
		print ('##################################################')
		if param.genome == 'geno38':
			system(' '.join(['/home/magi/miniconda3/envs/PY270/bin/python', '/home/magi/PROJECT/diagnosys/bin/resync_fastq.py', forward,reverse]))
		if param.genome == 'geno37':
			system(' '.join(['/home/magi/miniconda3/envs/PY270/bin/python', '/home/magi/PROJECT/diagnosys/bin/resync_fastq.py',forward,reverse]))
		system(' '.join(['rm',forward]))
		system(' '.join(['rm',reverse]))

	except ExceptionType as Argument:
		print ("You need to run.....!!!")
		print ("Exception: %s" % str(Argument))
		sys.exit(1)

	finally:
#		lock.release()
#	print name_process, 'Finished', name
		print ('##################################################')
	return


def Alignment(lock,param,folder,name,forward,reverse):
#	lock.acquire()
	name_process = multiprocessing.current_process().name
	print (name_process, 'Starting BWA..', name)
	try:
		print ('#################################################')
		sample_sam_ = join(folder,'temp/',name+'.sam')
		header = '"@RG\\tID:group1\\tSM:%s\\tLB:%s\\tPL:Illumina"' % (name,name)
		if param.genome == 'geno38':
			system(' '.join([param.bwa, 'mem','-S','-w 150','-R',header,
				'-t', param.threads, geno38,forward,reverse,
				'>', sample_sam_]))

			PrintLog(' '.join([param.bwa, 'mem','-S','-w 150','-R',header,
				'-t', param.threads, geno38,forward,reverse,
				'>', sample_sam_]),join(folder,'log'))

		if param.genome == 'geno37':
			system(' '.join([param.bwa, 'mem','-R',header,
				'-t', param.threads, geno37,forward,reverse,
				'>', sample_sam_]))

			PrintLog(' '.join([param.bwa, 'mem','-R',header,
				'-t', param.threads, geno37,forward,reverse,
				'>', sample_sam_]),join(folder,'log'))

	except ExceptionType as Argument:
		print ("You need to run.....!!!")
		print ("Exception: %s" % str(Argument))
		sys.exit(1)

	finally: pass
#		lock.release()
#	print name_process, 'Finished', name
	return

### FROM SAM TO BAM AND FIXMATE - REALIGNMENT #####
def FromSAMtoBAM(lock,param,folder,name):
#	lock.acquire()
	try:
		print ('##################################################')
		name_process = multiprocessing.current_process().name
		print (name_process, 'Starting')

		sample_sam_ = join(folder,'temp/',name+'.sam')
		sample_bam_ = join(folder,'temp/',name+'.bam')

		if param.genome == 'geno37':
			system(' '.join(['samtools','faidx',geno37]))
		elif param.genome == 'geno38':
			system(' '.join(['samtools','faidx',geno38]))
		system(' '.join(['samtools','view','-bS',sample_sam_,'>',sample_bam_]))
	finally: pass
	return
######### sort BAM #############
def sortBAM(lock,param,folder,name):
#	lock.acquire()
	try:
		print ('##################################################')
		sample_bam_ = join(folder,'temp/',name+'.bam')
		sample_bam_fixmate_ = join(folder,'temp/',name+'_fixmate.bam')
		sample_sorted_ = join(folder,'temp/',name+'_sorted.bam')

		print ('Fixmate...')
		system(' '.join(['samtools','fixmate','-r',sample_bam_,sample_bam_fixmate_]))
		#system(' '.join([param.samtools,'index',sample_bam_fixmate_]))
		print ('-----------------------------------------------------')

#		system(' '.join(['/home/magi/tools/samtools-1.6/sambamba_v0.6.7', 'sort','-m',param.memory,'-t',param.threads,
#					sample_bam_fixmate_,'-o', sample_sorted_+'.bam']))
		print ('Sort...')
		system(' '.join(['samtools', 'sort','-m','4G','-@','4',
					'-T SAMPLE' ,'-o', sample_sorted_,sample_bam_fixmate_]))

		system(' '.join(['samtools', 'index', sample_sorted_]))
	finally: pass
#		lock.release()
	return

#### remove duplicates ####
def remove_duplicate(lock,param,folder,name):
#	lock.acquire()
	try:
		print ('##################################################')
		print ('remove duplicates...')
		sample_sorted_ = join(folder,'temp/',name+'_sorted.bam')
		sample_bam_calmd_ = join(folder,'temp/',name+'_calmd.bam')
		sample_bam_rmdup = join(folder,'temp/',name+'_rmdup.bam')

		######  Tolto perche troppo lento###RIVALUTARE IN SEGUITO
#		if param.genome == 'geno38':
#			print 'Calmd...'
#			system(' '.join([param.samtools,'calmd','-Arb',sample_sorted_, geno38,'>',sample_bam_calmd_]))
#			system(' '.join([param.samtools, 'index', sample_bam_calmd_]))
#
#		if param.genome == 'geno37':
#			print 'Calmd...'
#			system(' '.join([param.samtools,'calmd','-Arb',sample_sorted_, geno37,'>',sample_bam_calmd_]))
#			system(' '.join([param.samtools, 'index', sample_bam_calmd_]))
#
#		system(' '.join(['sambamba_v0.6.4','markdup','-r','-t',param.threads,sample_bam_calmd_, sample_bam_rmdup]))
		print ('rmdup...')
		system(' '.join(['/home/magi/tools/samtools-1.9/sambamba_v0.6.7','markdup','-r','-t',param.threads,sample_sorted_, sample_bam_rmdup]))

		PrintLog(' '.join(['/home/magi/tools/samtools-1.9/sambamba_v0.6.7','markdup','-r','-t',
					param.threads,sample_sorted_, sample_bam_rmdup]),join(folder,'log'))

		system(' '.join(['samtools', 'index', sample_bam_rmdup]))
	finally: pass
#		lock.release()
	return


def gatk(lock,param,folder,name):
#	lock.acquire()
	try:
		print ('##################################################')
		lane = join(folder,'temp/'+name+'_lane.intervals')
		lane2 = join(folder,'temp/'+name+'_lane2.intervals')

		sample_bam_rmdup = join(folder,'temp/'+name+'_rmdup.bam')
		sample_bam_int = join(folder,'temp/'+name+'_int.bam')
		sample_bam_int2 = join(folder,'temp/'+name+'_int2.bam')

		sample_bam_final = join(folder,'bam/'+name+'_final.bam')

		if param.genome == 'geno38':
			print ('GATK......FIRST step.....!!!!')

			#system(' '.join(['java -Xmx10g -jar', param.gatk, '-T RealignerTargetCreator','-nt',param.threads,
			#		 '-R', geno38,'-known',clinvar_indel,'-I', sample_bam_rmdup ,'-o',lane]))

			#PrintLog(' '.join(['java -Xmx10g -jar', param.gatk, '-T RealignerTargetCreator','-nt',param.threads,
			#		 '-R', geno38,'-known',clinvar_indel,'-I', sample_bam_rmdup ,'-o',lane]),join(folder,'log'))

			#system(' '.join(['java -Xmx10g -jar', param.gatk, '-T IndelRealigner',
			#		 '-R', geno38,'-I', sample_bam_rmdup,'-targetIntervals', lane,'-o', sample_bam_int]))

			#PrintLog(' '.join(['java -Xmx10g -jar', param.gatk, '-T IndelRealigner',
			#		 '-R', geno38,'-I', sample_bam_rmdup,'-targetIntervals',
			#		 lane,'-o', sample_bam_int]),join(folder,'log'))


			#system(' '.join([param.gatk,'--java-options -Xmx10G LeftAlignIndels','-I',sample_bam_rmdup,'-R',geno38,'--OUTPUT',sample_bam_int]))
			system(' '.join([param.gatk,'--java-options -Xmx10G LeftAlignIndels','-I',sample_bam_rmdup,'-R',geno38,'-O',sample_bam_int]))

			print ('GATK.......SECOND step.....!!!!')

			system(' '.join([param.gatk,'--java-options -Xmx10G BaseRecalibrator','-R',geno38,
					'--known-sites','/home/magi/dataset/GATKRESOURCE/1000G_omni2.5.hg38.vcf.gz',
					'--known-sites','/home/magi/dataset/GATKRESOURCE/Homo_sapiens_assembly38.variantEvalGoldStandard.vcf.gz',
					'--known-sites','/home/magi/dataset/GATKRESOURCE/Homo_sapiens_assembly38.dbsnp.vcf.gz',
					'--known-sites','/home/magi/dataset/GATKRESOURCE/Homo_sapiens_assembly38.known_indels.vcf.gz',
					'--known-sites','/home/magi/dataset/GATKRESOURCE/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz',
					'-I',sample_bam_int,
					'-O', lane2]))

			system(' '.join([param.gatk,'--java-options -Xmx10G ApplyBQSR','-bqsr',lane2,
					 '-I',sample_bam_int,'-O',sample_bam_final ]))

			#system(' '.join(['gatk --java-options -Xmx10G GatherBamFiles',
					 #'-I',sample_bam_int2,'-O',sample_bam_final ]))


		if param.genome == 'geno37':
			print ('GATK.......first step.....!!!!')

			system(' '.join(['java -Xmx10g -jar', param.gatk, '-T RealignerTargetCreator','-nt',param.threads,
					'-R', geno37,'-I', sample_bam_rmdup ,'-o',lane]))

			PrintLog(' '.join(['java -Xmx10g -jar', param.gatk, '-T RealignerTargetCreator','-nt',param.threads,
					'-R', geno37,'-I', sample_bam_rmdup ,'-o',lane]),join(folder,'log'))

			system(' '.join(['java -Xmx10g -jar', param.gatk, '-T IndelRealigner',
					'-R', geno37,'-I', sample_bam_rmdup,'-targetIntervals', lane,'-o', sample_bam_int]))

			PrintLog(' '.join(['java -Xmx10g -jar', param.gatk, '-T IndelRealigner',
					'-R', geno37,'-I', sample_bam_rmdup,'-targetIntervals',
					 lane,'-o', sample_bam_int]),join(folder,'log'))

			print ('GATK.......second step.....!!!!')

			system(' '.join(['java -Xmx10g -jar',param.gatk,'-T BaseRecalibrator',
					'-R',geno37,'-knownSites', dbsnp144_37,'-I',sample_bam_int,'-o', lane2]))

			PrintLog(' '.join(['java -Xmx10g -jar',param.gatk,'-T BaseRecalibrator',
					'-R',geno37,'-knownSites', dbsnp144_37,'-I',sample_bam_int,'-o', lane2]) ,join(folder,'log'))

			system(' '.join(['java -Xmx10g -jar',param.gatk,'-T PrintReads',
					'-R',geno37,'-I',sample_bam_int,'-BQSR',lane2, '-o',sample_bam_final ]))

			PrintLog(' '.join(['java -Xmx10g -jar',param.gatk,'-T PrintReads',
					'-R',geno37,'-I',sample_bam_int,'-BQSR',lane2, '-o',sample_bam_final ]),join(folder,'log'))

#		system(' '.join([param.samtools, 'index', sample_bam_final]))
	finally: pass
#		lock.release()
	return

############################################################################
############################################################################
def remove_temp_file(param,folder):
	return os.remove(join(folder,'temp/*'))
############################################################################
if __name__=="__main__":
	lock = Lock()
	args=InputPar()

	if args.over == 'True':
		folder = create_folder()
		folder_name = principal_folder(args,folder,over=args.over)
		path_creation(args,folder_name)
	elif args.over == 'False':
		folder = create_folder()
		folder_name = principal_folder(args,folder,over=args.over)

	###folder_name = '/home/magi/PROJECT/diagnosys/RESULT/24_Jun_2016_CLM2'
	###sample_list = glob.glob('/home/magi/PROJECT/diagnosys/RESULT/24_Jun_2016_CLM2/fastq/'+'*_new.fastq')
	#print_args=vars(args)#take args in a dict and print in a log
	#for k,v in print_args.iteritems():
	#	strategy_out='='.join([str(k),str(v)])
	#	path = join(folder_name,'log')
	#	PrintLog(strategy_out,path)
	###print sample_list

	if args.over == 'True':
		sample_list = copy_fastq(args,folder_name)
		#print sample_list
		data = set_samples(args,sample_list,folder_name)
	elif args.over == 'False':
		sample_list = glob.glob('/home/magi/PROJECT/diagnosys/RESULT/19_Jul_2023_LYMPHOBESITY/fastq/'+'*')
		print (sample_list)
		data = set_samples(args,sample_list,folder_name)
		#data = pd.read_csv(join(folder_name,'sample_list.csv'),sep='\t',header=0,encoding='utf-8')
###################################################################################################
###################################################################################################
	jobs0 = []
	jobs1 = []
	jobs2 = []
	jobs3 = []
	jobs4 = []
	jobs5 = []
	try:
		for index, sample in data.iterrows():
			try:
				name = str(sample['name'])
				log_file=join(folder_name,'log',name+'_core.log')
				writer = Writer(sys.stdout,log_file)
				sys.stdout = writer
				forward = sample['forward']
				reverse = sample['reverse']
				if args.over == 'True':
					p0 = PreAlignment(lock,args,folder_name,name,forward,reverse)
					forward = sample['forward'] + '_pairs_R1.fastq'
					reverse = sample['reverse'] + '_pairs_R2.fastq'
					p1 = Alignment(lock,args,folder_name,name,forward,reverse)
					p2 = FromSAMtoBAM(lock,args,folder_name,name)
					p3 = sortBAM(lock,args,folder_name,name)
					p4 = remove_duplicate(lock,args,folder_name,name)
					p5 = gatk(lock,args,folder_name,name)
				elif args.over == 'False':
					#p2 = FromSAMtoBAM(lock,args,folder_name,name)
					p3 = sortBAM(lock,args,folder_name,name)
					p4 = remove_duplicate(lock,args,folder_name,name)
					p5 = gatk(lock,args,folder_name,name)
			except: pass
		spec_fastq = join(folder_name,'fastq/')
		spec_fastqfilter = join(folder_name,'fastqfiltered/')
		system(' '.join(['mv',spec_fastq+'*pairs*',spec_fastqfilter]))
		system(' '.join(['rm',spec_fastq+'*single*']))
##############################################################################################
#############################################################################################
	except AttributeError:
		print ('Check your FASTQ FOLDER and try again!!!')
		sys.exit
