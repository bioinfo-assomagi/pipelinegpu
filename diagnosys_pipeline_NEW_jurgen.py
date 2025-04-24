#!/home/magi/miniconda3/envs/PY310/bin/python
#27/06/2015

# -*- coding: utf-8 -*-
"""
Created on Wed MAY 16 2023
@author: Giuseppe Marceddu
"""
from multiprocessing import Process, Lock
import multiprocessing
import argparse
from os import listdir, system, getcwd,getenv
from os.path import isfile, join, exists
import datetime
import subprocess
import time
from os import getenv
import time,os,sys

path = getcwd()
#GATK = join('/home/magi/','tools/GenomeAnalysisTK36.jar')
GATK = join('/home/magi/','tools/gatk-4.2.0.0/gatk')
#VEPS = join('/home/magi/','tools/ensembl-vep-release-104/vep')
VEPS = join('/home/magi/','tools/ensembl-vep-release-110/vep')

#if exists(join(getcwd(),'nohup.out')):
def InputPar():
    ####Introducing arguments
	parser = argparse.ArgumentParser(prog='MAGI EUREGIO DIAGNOSYS',description='Pipe from FASTQ to BAM',
			epilog='Need to BWA 0.7.17 and SAMTOOLS v1.14 and GATK 4.2 (SET on global PATH or INSERT your path [--bwa] [--samtools]')

	parser.add_argument('-p','--path', metavar='PATH', default=path,
	                    help='[Default = Pipe run in same folder when launch command]')

	parser.add_argument('-d','--dest', metavar='destination', choices=['b','r','z'],
	                    required=True,
	                    help='Choices destination: b = bolzano; r = rovereto; z = ricerca,  - required ')

	#parser.add_argument('-pan','--panel',metavar='PANEL',type=int, choices=range(1,10),
	#			required=True,
	#			help='Choose Panel from 1 to 7: 1:CLM2; 2:OSSNEW; 3:OCULARE; 4:OCULARE2; 5:SIFSR; 6:ALLGENE; '
	#				'7:TruSight One; 8:GENETEST1; 9:GENETEST2; \nRequired')
	parser.add_argument('-pan','--panel',metavar='PANEL',required=True,
				help='Pannelli Consigliati: 1:CLM2; 2:OSSNEW; 3:OCULARE; 4:OCULARE2; 5:SIFSR; 6:ALLGENE;'
					'7:trusightone; 8:GENETEST1; 9:GENETEST2; 10:CANCER; 11:MALFORMATIONI VASCOLARI;'
					'12:ANOMALIE VASCOLARI; 13:INFERTILITA; 14:OBESITA; 15:GENODERMATOSI\nRequired')

	parser.add_argument('-fq','--fastq', metavar='FASTQ',
	                    help='FULL PATH to find FASTQ FILE - [default FROM NAS]')

	parser.add_argument('-proj','--project',metavar='PROJECT NAME', help='Insert Family name to create Folder')

	parser.add_argument('-g','--genome', metavar='choose Assembling', choices=['geno37','geno38'],
	                    default='geno38',
	                    help='Choices: geno37, geno38 - Run Analysis with genome GRCh38 or GRCh37')

	parser.add_argument('-q','--quality',  metavar='PHRED QUALITY FILTER',default=18,
	                    help=' Choose a quality threshold [Default = 18]')

	parser.add_argument('-Q','--quality_perc',  metavar='PERC QUALITY READ',default=97,
	                    help='Percentual of quality for read [Default = 97]')

	parser.add_argument('-N','--threads', metavar='THREADS', default='12',
	                    help='Number of threads to run [20]')

	parser.add_argument('-m','--memory', metavar='MERORY', default='15G',
	                    help='max memory per thread; suffix K/M/G [defaul = 15G]')

	parser.add_argument('-b','--bwa', metavar='BWA', default='bwa',
	                    help='Insert BWA path')

	parser.add_argument('-s','--samtools', metavar='SAMTOOLS', default='samtools',
	                    help='Insert SAMTOOLS path')

	parser.add_argument('-t','--samstat', metavar='SAMSTAT', default='samstat',
	                    help='Insert SAMSTAT path')

	parser.add_argument('-gk','--gatk', metavar='GATK', default=GATK, help='Insert gatk path')

	parser.add_argument('-v','--bcftools', metavar='BCFTOOLS', default='bcftools',
	                    help='Insert BCFTOOLS path')

	parser.add_argument('-ovr','--over', metavar='New/OLD project', choices=['True','False'],
	                    default='True',
	                    help='Choices: ALLERT!!! Set on "False" overwrite old data [Default = True]')

	parser.add_argument('-vep','--VEP', metavar='variant effect predictor', default=VEPS,
                            help='Insert VEP PATH')

	parser.add_argument('-del','--delete', metavar='Delete temp file after Backup',choices=['False','True'],  required=True,
                            help='Default False = Not Delete, put True = Delete')

	parser.add_argument('-freq','--runfrequency', metavar='Choose if Run Frequency Calculation or Not',choices=['False','True'],
	                    default='True',
                            help='Default False = Not Delete, put True = Delete')

	parser.add_argument('-core','--runcore', metavar='Choose if Run Frequency Calculation or Not',choices=['False','True'],
	                    default='True',
                            help='Default False = Not Delete, put True = Delete')

	parser.add_argument('-importapp','--importAPP', metavar='Choose if Run ImportAPP or Not',choices=['False','True'],
	                    default='True',
                            help='Default False = Not Delete, put True = Delete')

	parser.add_argument('-backup','--backup', metavar='Choose if make BackUp or Not',choices=['False','True'],
	                    default='True',
                            help='Default False = Not Backup, put True = Backup')

	return parser.parse_args()

def create_folder():
	today = datetime.date.today()
	name = "{:%d_%b_%Y}".format(today)
	return name

if __name__=="__main__":
	args = InputPar()
	proj_name = create_folder()
	if args.panel == 'trusightone':
		panel = 'trusightone' #args.panel
	else:
		panel = args.panel.upper()

	if args.dest == 'r':
		dest = 'rovereto'
	elif args.dest == 'b':
		dest = 'bolzano'
	elif args.dest == 'z':
		dest = 'ricerca'

	if (args.dest == 'b') | (args.dest == 'z') | (args.dest == 'r'):
		if args.genome == 'geno37': genome = 'geno37'
		elif args.genome == 'geno38': genome = 'geno38'

		if args.over == 'True': over = 'True'
		elif args.over == 'False': over = 'False'
		if args.project:
			proj_name = args.project
			name_folder= join(path,'RESULT',proj_name+'_'+panel)
			print ('\nIL NOME DEL PROGETTO E\': '+proj_name+'_'+panel+' ---> PER: '+dest+'\n')
		else:
			name_folder= join(path,'RESULT',proj_name+'_'+panel)
			print ('\nIL NOME DEL PROGETTO E\': '+proj_name+'_'+panel+' ---> PER: '+dest+'\n')

		#name_folder= join(path,'RESULT',proj_name+'_'+panel)
		print (name_folder+'!!!')

		if ((os.path.exists(name_folder)) & (args.over == 'True')):
			#pass
			sys.exit("Folder exist just - Choose another project name, please!!! You can run with <-proj> option")
		elif ((os.path.exists(name_folder)) & (args.over == 'False')):
			args.runcore = 'False'
		else: pass
#########################################################################################################################
#########################################################################################################################
		if args.runcore == 'False': print ('CORE NOT RUN!!!')
		if args.runcore == 'True':
			if args.fastq:
				print (' '.join(['/home/magi/PROJECT/diagnosys/bin_jurgen/first_diagnosys_core_V3_jurgen.py','-p',args.path,'-proj',proj_name,'-pan',str(panel),
					'-d',args.dest,'-g',genome,'-ovr',over,'-fq',args.fastq,'-q',str(args.quality),'-Q',str(args.quality_perc),
					'-N',str(args.threads),'-m',str(args.memory),'-b',args.bwa,'-s',args.samtools,'-gk',args.gatk]))
				time.sleep(3)
				system(' '.join(['/home/magi/PROJECT/diagnosys/bin_jurgen/first_diagnosys_core_V3_jurgen.py','-p',args.path,'-proj',proj_name,'-pan',str(panel),
					'-d',args.dest,'-g',genome,'-ovr',over,'-fq',args.fastq,'-q',str(args.quality),'-Q',str(args.quality_perc),
					'-N',str(args.threads),'-m',str(args.memory),'-b',args.bwa,'-s',args.samtools,'-gk',args.gatk]))
			else:
				print (' '.join(['/home/magi/PROJECT/diagnosys/bin_jurgen/first_diagnosys_core_V3_jurgen.py','-p',args.path,'-proj',proj_name,'-pan',str(panel),
					'-d',args.dest,'-g',genome,'-ovr',over,'-q',str(args.quality),'-Q',str(args.quality_perc),
					'-N',str(args.threads),'-m',str(args.memory),'-b',args.bwa,'-s',args.samtools,'-gk',args.gatk]))
				time.sleep(3)
				system(' '.join(['/home/magi/PROJECT/diagnosys/bin_jurgen/first_diagnosys_core_V3_jurgen.py','-p',args.path,'-proj',proj_name,'-pan',str(panel),
					'-d',args.dest,'-g',genome,'-ovr',over,'-q',str(args.quality),'-Q',str(args.quality_perc),
					'-N',str(args.threads),'-m',str(args.memory),'-b',args.bwa,'-s',args.samtools,'-gk',args.gatk ]))
#############################################################################################################################
#############################################################################################################################
		print ('\n')
		over = 'False'
		print (' '.join(['/home/magi/PROJECT/diagnosys/bin_jurgen/second_A_diagnosys_COV_NEW_jurgen.py','-p',args.path,'-proj',proj_name,'-pan',str(panel),
			'-d',args.dest,'-g',genome,'-ovr',over,'-q',str(args.quality),
			'-N',str(args.threads),'-v',args.bcftools,'-t',args.samstat,'-s',args.samtools,'-gk',args.gatk]))
		time.sleep(3)
		system(' '.join(['/home/magi/PROJECT/diagnosys/bin_jurgen/second_A_diagnosys_COV_NEW_jurgen.py','-p',args.path,'-proj',proj_name,'-pan',str(panel),
			'-d ',args.dest,'-g',genome,'-ovr',over,'-q',str(args.quality),
			'-N',str(args.threads),'-v',args.bcftools,'-t',args.samstat,'-s',args.samtools,'-gk',args.gatk]))
########################################################################################################################
		time.sleep(3)
########################################################################################################################
		print ('\n')
		over = 'False'
		print (' '.join(['/home/magi/PROJECT/diagnosys/bin_jurgen/second_B_diagnosys_VCF_NEW_jurgen.py','-p',args.path,'-proj',proj_name,'-pan',str(panel),
			'-d',args.dest,'-g',genome,'-ovr',over,'-q',str(args.quality),
			'-N',str(args.threads),'-v',args.bcftools,'-t',args.samstat,'-s',args.samtools,'-gk',args.gatk]))
		time.sleep(3)
		system(' '.join(['/home/magi/PROJECT/diagnosys/bin_jurgen/second_B_diagnosys_VCF_NEW_jurgen.py','-p',args.path,'-proj',proj_name,'-pan',str(panel),
			'-d ',args.dest,'-g',genome,'-ovr',over,'-q',str(args.quality),
			'-N',str(args.threads),'-v',args.bcftools,'-t',args.samstat,'-s',args.samtools,'-gk',args.gatk]))
########################################################################################################################
		time.sleep(3)
########################################################################################################################
		print ('\n')
		print (' '.join(['/home/magi/PROJECT/diagnosys/bin/third_a_diagnosys_annotation_VEP111_jurgen.py','-p',args.path,'-proj',
			proj_name,'-pan',str(panel),'-d',args.dest,'-g',genome,'-ovr',over,'-vep',args.VEP]))
		time.sleep(3)
		system(' '.join(['/home/magi/PROJECT/diagnosys/bin/third_a_diagnosys_annotation_VEP111_jurgen.py','-p',args.path,'-proj',
			proj_name,'-pan',str(panel),'-d ',args.dest,'-g',genome,'-vep',args.VEP]))
#########################################################################################################################
		#time.sleep(3)
###########################################################################################################################
		#print( '\n')
		#print (' '.join(['/home/magi/PROJECT/diagnosys/bin/third_b_diagnosys_coverage.py','-p',args.path,'-proj',
		#			proj_name,'-pan',str(panel),'-d',args.dest,'-g',genome]))
		#time.sleep(3)
		#system(' '.join(['/home/magi/PROJECT/diagnosys/bin/third_b_diagnosys_coverage.py','-p',args.path,'-proj',
		#			proj_name,'-pan',str(panel),'-d ',args.dest,'-g',genome]))
##########################################################################################################################
		#time.sleep(3)
##########################################################################################################################
		#print ('\n')
		#print (' '.join(['/home/magi/PROJECT/diagnosys/bin/kinship_test.py','-p',args.path,'-proj',proj_name,'-pan',str(panel),
		#			'-d',args.dest,'-g',genome]))
		#time.sleep(3)
		#system(' '.join(['/home/magi/PROJECT/diagnosys/bin/kinship_test.py','-p',args.path,'-proj',proj_name,'-pan',str(panel),
		#			'-d ',args.dest,'-g',genome]))
##########################################################################################################################
		#time.sleep(3)
##########################################################################################################################
		#if args.runfrequency == 'False': print ('FREQUENCY NOT RUN!!!')
		#elif args.runfrequency == 'True':
		#	print (' '.join(['/home/magi/PROJECT/diagnosys/bin/frequency_analysis_multifileversion.py','-p',args.path,'-proj',
		#		proj_name,'-pan',str(panel),'-d',args.dest,'-g',genome]))
		#	time.sleep(3)
		#	system(' '.join(['/home/magi/PROJECT/diagnosys/bin/frequency_analysis_multifileversion.py','-p',args.path,'-proj',
		#		proj_name,'-pan',str(panel),'-d',args.dest,'-g',genome]))
##########################################################################################################################
		#time.sleep(3)
##########################################################################################################################
		#print ('\n')
		#print (' '.join(['/home/magi/PROJECT/diagnosys/bin/analisi_TRIO_NEW.py','-p',args.path,'-proj',proj_name,
		#	'-pan',str(panel),'-d',args.dest]))
		#time.sleep(3)
		#system(' '.join(['/home/magi/PROJECT/diagnosys/bin/analisi_TRIO_NEW.py','-p',args.path,'-proj',proj_name,
		#	'-pan',str(panel),'-d ',args.dest]))
##########################################################################################################################
		#time.sleep(3)
##########################################################################################################################
		#print ('\n')
		#print (' '.join(['/home/magi/PROJECT/diagnosys/bin/third_d_diagnosys_nosanger_NEW.py','-p',args.path,'-proj',
		#		proj_name,'-pan',str(panel),'-d',args.dest,'-g',genome]))
		#time.sleep(3)
		#system(' '.join(['/home/magi/PROJECT/diagnosys/bin/third_d_diagnosys_nosanger_NEW.py','-p',args.path,'-proj',
		#	proj_name,'-pan',str(panel),'-d ',args.dest,'-g',genome]))
##########################################################################################################################
		#time.sleep(3)
###########################################################################################################################
		#print ('\n')
		#print (' '.join(['/home/magi/PROJECT/diagnosys/bin/copy_on_igv.py','-p',args.path,'-proj',proj_name,'-pan',str(panel),
		# 	'-d',args.dest,'-g',genome]))
		#time.sleep(3)
		#system(' '.join(['/home/magi/PROJECT/diagnosys/bin/copy_on_igv.py','-p',args.path,'-proj',proj_name,'-pan',str(panel),
		# 	'-d',args.dest,'-g',genome]))
##########################################################################################################################
		#time.sleep(3)
##########################################################################################################################
		#print ('\n')
		#print (' '.join(['/home/magi/PROJECT/diagnosys/bin/third_c_diagnosys_indel_NEW.py','-p',args.path,'-proj',
		# 	proj_name,'-pan',str(panel),'-d',args.dest,'-g',genome]))
		#time.sleep(3)
		#system(' '.join(['/home/magi/PROJECT/diagnosys/bin/third_c_diagnosys_indel_NEW.py','-p',args.path,'-proj',
		# 	proj_name,'-pan',str(panel),'-d ',args.dest,'-g',genome]))
##########################################################################################################################
		#time.sleep(3)
##########################################################################################################################
		#print ('\n')
		#print (' '.join(['/home/magi/PROJECT/diagnosys/bin/third_c_diagnosys_indel_INTERPRETATION.py','-p',args.path,'-proj',
		# 		proj_name,'-pan',str(panel),'-d',args.dest,'-g',genome]))
		#time.sleep(3)
		#system(' '.join(['/home/magi/PROJECT/diagnosys/bin/third_c_diagnosys_indel_INTERPRETATION.py','-p',args.path,'-proj',
		# 		proj_name,'-pan',str(panel),'-d ',args.dest,'-g',genome]))
##########################################################################################################################
		#time.sleep(3)
#########################################################################################################################
###############################OLDDDDDDDDDDDDDD###########################################################################
		#print ('\n')
		#print (' '.join(['/home/magi/PROJECT/diagnosys/bin/fourth_diagnosys_selezioneautomatica_NEW.py','-p',args.path,'-proj',
		#	proj_name,'-pan',str(panel),'-d',args.dest]))
		#time.sleep(3)
		#system(' '.join(['/home/magi/PROJECT/diagnosys/bin/fourth_diagnosys_selezioneautomatica_NEW.py','-p',args.path,'-proj',
		#			proj_name,'-pan',str(panel),'-d ',args.dest]))
##########################################################################################################################
		#time.sleep(3)
##########################################################################################################################
		#print ('\n')
		#print (' '.join(['/home/magi/PROJECT/diagnosys/bin/fourth_diagnosys_selezioneautomatica_FAMILY.py','-p',args.path,'-proj',proj_name,'-pan',str(panel),
		#		'-d',args.dest]))
		#time.sleep(3)
		#system(' '.join(['/home/magi/PROJECT/diagnosys/bin/fourth_diagnosys_selezioneautomatica_FAMILY.py','-p',args.path,'-proj',proj_name,'-pan',str(panel),
		#		'-d ',args.dest]))
##########################################################################################################################
		#time.sleep(3)
##########################################################################################################################
		#print ('\n')
		#print (' '.join(['/home/magi/PROJECT/diagnosys/bin/fourth_b_diagnosys_variantinterpretation.py','-p',args.path,'-proj',proj_name,'-pan',str(panel),
		#		'-d',args.dest]))
		#time.sleep(3)
		#system(' '.join(['/home/magi/PROJECT/diagnosys/bin/fourth_b_diagnosys_variantinterpretation.py','-p',args.path,'-proj',proj_name,'-pan',str(panel),
		#		'-d ',args.dest]))
##########################################################################################################################
		#time.sleep(3)
##########################################################################################################################
		#print ('\n')
		#print (' '.join(['/home/magi/PROJECT/diagnosys/bin/fifth_diagnosys_automatic_reporting_IT.py','-p',args.path,'-proj',proj_name,'-pan',str(panel),
		#		'-d',args.dest]))
		#time.sleep(3)
		#system(' '.join(['/home/magi/PROJECT/diagnosys/bin/fifth_diagnosys_automatic_reporting_IT.py','-p',args.path,'-proj',proj_name,'-pan',str(panel),
		#		'-d ',args.dest]))
##########################################################################################################################
		#time.sleep(3)
##########################################################################################################################
		#print ('\n')
		#if args.importAPP == 'False': print ('IMPORT APP NOT RUN!!!')
		#if args.importAPP == 'True':
		#	if args.dest == 'b':
		#		print(' '.join(['sh','/home/magi/PROJECT/diagnosys/bin/rsyncdata_BOLZANO.sh']))
		#		system(' '.join(['sh','/home/magi/PROJECT/diagnosys/bin/rsyncdata_BOLZANO.sh']))
		#	elif args.dest == 'z':
		#		print(' '.join(['sh','/home/magi/PROJECT/diagnosys/bin/rsyncdata_RICERCA.sh']))
		#		system(' '.join(['sh','/home/magi/PROJECT/diagnosys/bin/rsyncdata_RICERCA.sh']))
		#	elif args.dest == 'r':
		#		print(' '.join(['sh','/home/magi/PROJECT/diagnosys/bin/rsyncdata_ROVERETO.sh']))
		#		system(' '.join(['sh','/home/magi/PROJECT/diagnosys/bin/rsyncdata_ROVERETO.sh']))
#########################################################################################################################
		#print ('\n')
		#if args.backup == 'False': print ('BACKUP NOT RUN!!!')
		#if args.backup == 'True':
		#	print (' '.join(['/home/magi/PROJECT/diagnosys/bin/backup_diagnosys.py','-p',args.path,'-proj',proj_name,'-pan',str(panel),
		#				'-d',args.dest,'-g',genome,'-ovr',over,'-del',str(args.delete)]))
		#	time.sleep(3)
		#	system(' '.join(['/home/magi/PROJECT/diagnosys/bin/backup_diagnosys.py','-p',args.path,'-proj',proj_name,'-pan',str(panel),
		#				'-d',args.dest,'-g',genome,'-del',str(args.delete)]))
##########################################################################################################################
