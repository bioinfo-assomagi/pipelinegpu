from Pipes.Pipe import Pipe
from Pipes.ParallelPipe import ParallelPipe
import os

import config
import sys
import utils
from os.path import join
import glob
import re
import shutil
import tools
from DBContext import DBContext

"""
3 Pipelines contained here - PrincipalFolderPipe, ReadFastQFilesPipe, ProcessFastQFilesPipe

"""

"""
['/home/bioinfo/PROJECT/diagnosys/bin/first_diagnosys_core_V2.py','-p',args.path,'-proj',proj_name,'-pan',str(panel),
					'-d',args.dest,'-g',genome,'-ovr',over,'-fq',args.fastq,'-q',str(args.quality),'-Q',str(args.quality_perc),
					'-N',str(args.threads),'-m',str(args.memory),'-b',args.bwa,'-s',args.samtools,'-gk',args.gatk])"""


""" Creates the principal directory of the sequencing project. """
class PrincipalFolderPipe(Pipe):

    def __init__(self) -> None:
        super().__init__()

    def process(self, **kwargs):
        print("Inside PrincipalFolderPipe ... ")
    
        over = kwargs.pop("over", None)  # if True, overwrites old data
        principal_directory = kwargs.pop('name_folder', None)
        dest = kwargs.pop('dest', None)

        if os.path.exists(principal_directory):
            if not over:
                sys.exit("Folder exist just - Choose another project name, please!!! You can run with <-proj> option")
            else:
                shutil.rmtree(principal_directory)

        os.makedirs(principal_directory)
        #config.update("DEFAULTS", "principal_directory", principal_directory)

        self.create_paths(principal_directory)

        kwargs.update({'over': over, 'principal_directory': principal_directory, 'dest': dest})
        print("Principal Directory: {}".format(principal_directory))

        return kwargs
        

    """ Builds the project tree inside the principal directory. TODO: write this to config """
    def create_paths(self, folder_name):
        print('Building tree folder...')  # TODO: check if there is a general way to build directory trees
        # programmatically
        directories = ["fastq", "fastqfiltered", "temp", "control", "vcf", "bam", "bam/inalaysis", "plot", "annotation",
                       "pheno", "coverage", "report", "fastQC", "final", "indel"]

        parabricks_input = join(folder_name, 'parabricks_input/')
        parabricks_output = join(folder_name, 'parabricks_output/')
        os.makedirs(parabricks_input)
        os.makedirs(parabricks_output)
        #config.update("PARABRICKS", "parabricks_input", parabricks_input)
        #config.update("PARABRICKS", "parabricks_output", parabricks_output)


        spec_fastq = join(folder_name, 'fastq/')
        spec_fastqfilter = join(folder_name, 'fastqfiltered/')
        spec_temp = join(folder_name, 'temp/')
        spec_control = join(folder_name, 'control/')
        spec_vcf = join(folder_name, 'vcf/')
        spec_bam = join(folder_name, 'bam/')
        spec_baminanalysis = join(folder_name, 'bam/inanalysis')
        spec_plot = join(folder_name, 'plot/')
        spec_annot = join(folder_name, 'annotation/')
        spec_pheno = join(folder_name, 'pheno/')
        spec_cov = join(folder_name, 'coverage/')
        folder_report = join(folder_name, 'report/')
        # spec_cov_fig = join(folder,'coverage/figure')
        spec_qc = join(folder_name, 'fastQC/')
        spec_final = join(folder_name, 'final/')
        spec_indel = join(folder_name, 'indel/')
        # spec_json = join(folder,'json/')
        # spec_json_stat = join(folder,'json/','cov_stat')
        # spec_json_annot = join(folder,'json/','annot')
        spec_log = join(folder_name, 'log')
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
        # if not os.path.exists(spec_json):
        #	os.makedirs(spec_json)
        # if not os.path.exists(spec_json_stat):
        #	os.makedirs(spec_json_stat)
        # if not os.path.exists(spec_json_annot):
        #	os.makedirs(spec_json_annot)
        if not os.path.exists(spec_log):
            os.makedirs(spec_log)
        #	if not os.path.exists(spec_cov_fig):
        #		os.makedirs(spec_cov_fig)
        if not os.path.exists(folder_report):
            os.makedirs(folder_report)

    

class ResyncDBPipe(Pipe):

    def __init__(self) -> None:
        super().__init__()
    
    def process(self, **kwargs):
        
        server_id = kwargs.pop("dest", None)

        print("Resyncing DB ...")
        db_server = utils.get_db_server(server_id)
        local_db_path = utils.get_db_path(server_id)
        #os.system("rsync -avz -e ssh {} {}".format(db_server, config.DB_PATH))

        kwargs.update({"dest": server_id, "db_path": local_db_path})
        return kwargs

""" Transfer the FastQ files from the external server to the project directory. """
class ReadFastQFilesPipe(Pipe):

    def __init__(self) -> None:
        super().__init__()

    def process(self, **kwargs):

        principal_directory = kwargs.pop('principal_directory')
        fastq = kwargs.pop('fastq', None)
        dest = kwargs.pop('dest', None)

        fastq_folder = self.copy_fastq_files(fastq, dest, principal_directory)
        fastq_files = glob.glob(fastq_folder + '*')

        kwargs.update({"fastq_folder": fastq_folder, "principal_directory": principal_directory, "fastq": fastq, "dest": dest, "fastq_files": fastq_files})
        return kwargs

    def copy_fastq_files(self, fastq_source, server_id, destination):
        
        print('Copying fastq files ... ')
        fastq_folder = join(destination, 'fastq/')

        if fastq_source:
            files_fq = glob.glob(join(fastq_source, '*'))
            for file in files_fq:
                if '.fastq.gz' in file or '.fq.gz' in file or '.fastaq.gz' in file or '.fastq' in file or '.fastaq' in file or '.fq' in file:
                    os.system(' '.join(['cp', file, fastq_folder]))
                elif 'phenotype' in file:
                    'Copying phenotype file...ever is better check inside this file....'
                else:
                    print(file, 'is not a VALID file and will not copy')
        else:
            if server_id == 'r':
                #files_fq = os.system(' '.join(['scp root@192.168.2.188:/sharedfolders/NGS/tmp/analysis/germinal/*', fastq_folder])) #added
                #files_fq = os.system(' '.join(['s *', fastq_folder]))  # added
                files_fq = os.system(
                    ' '.join(['scp root@192.168.1.51:/home/NGS/tmp/analysis/germinal/*', fastq_folder]))
            else:
                files_fq = os.system(
                    ' '.join(['scp root@192.168.1.51:/home/NGS/tmp/analysis/germinal/*', fastq_folder]))  # added
            # files_fq = system(' '.join(['scp server@192.168.1.201:/media/4e955bfb-88f0-4cc7-a824-27ee0b4bf6e2/NGS
            # /tmp/analysis/germinal2/*',fastq_folder])) files_fq = system(' '.join(['scp
            # server@192.168.1.201:/media/4e955bfb-88f0-4cc7-a824-27ee0b4bf6e2/NGS/tmp/analysis/somatic1/*',
            # fastq_folder]))

        return fastq_folder


class UnzipFastQFilesPipe(Pipe):

    def __init__(self) -> None:
        super().__init__()

    def process(self, **kwargs):
        self.thread_id = kwargs.pop('thread_id')
        fastq_folder = kwargs.pop("fastq_folder")
        # principal_directory = kwargs.pop("name_folder")
        fastq_files = kwargs.pop("fastq_files")
        utils.thread_print(self.thread_id, "Inside UnzipFastQFilesPipe pipe. List of fastq_files = {}".format(fastq_files))
        unizzped_fastq_files = self.unzip_fastq(fastq_files)
        utils.thread_print(self.thread_id, "Inside UnzipFastQFilesPipe pipe. List of unizzped_fastq_files = {}".format(unizzped_fastq_files))
        kwargs.update({"fastq_folder": fastq_folder, "fastq_files": unizzped_fastq_files, "thread_id": self.thread_id})
        return kwargs

    def unzip_fastq(self, fastq_files):
        # TODO: run this only if there are zipped files
        # TODO: however the preprocessing filename part should be kept, thus decouple these two functionalities

        # fastq_files = glob.glob(fastq_folder+'*')

        # Rename the files
        utils.thread_print(self.thread_id, "Renaming FastQ files ... ")

        unzipped_files = []
        for file in fastq_files:
            utils.thread_print(self.thread_id, "Old filename: {}".format(file))
            a = file.split("/")
            if '-' in a[-1]:
                x = a[-1].split("-")
                x = (x[0] + '.' + x[1])
            else:
                x = a[-1].split(".")
                x = (x[0] + '.' + x[1])
            x1 = x.split("_")
            name = x1[0]
            files3 = file.replace('-', '.')
            filename_new = re.sub(r"(_S\d+_L\d+|_L\d+)", '', files3)
            utils.thread_print(self.thread_id, "New filename: {}".format(filename_new))
            os.rename(file, filename_new)

            if 'new' not in filename_new: #NOTE: 'new' is not related to filename_new. 'new' is appended after fastq_quality_control and fastx_trimmer. filename_new is just a processed filename according to the rules above
                if '.gz' in filename_new:
                    utils.thread_print(self.thread_id, 'unzip files...')
                    os.system(' '.join(['gunzip', filename_new]))
                    unzipped_files.append(
                        filename_new[:-3])  # since the unzipped file will not contain .gz extension in its filename NOTE: this name is also automatically generated by 'gunzip'?
                else:
                    unzipped_files.append(filename_new)

        return unzipped_files


class ProcessFastQFilesPipe2(Pipe):
    def __init__(self) -> None:
        super().__init__()

    def process(self, **kwargs):
        fastq_folder = kwargs.pop("fastq_folder")
        principal_directory = kwargs.pop("principal_directory")
        fastq_files = kwargs.pop("fastq_files")

        threads = kwargs.pop("threads")
        quality = kwargs.pop("quality")
        quality_perc = kwargs.pop("quality_perc")
        self.thread_id = kwargs.pop("thread_id")

        fastqc_folder = join(principal_directory, 'fastQC/')

        porecessed_fastq_files = self.process_fastq(fastq_files, fastqc_folder, threads, quality, quality_perc)

        kwargs.update(
            {"principal_directory": principal_directory, "fastq_folder": fastq_folder, "fastqc_folder": fastqc_folder,
             "fastq_files": porecessed_fastq_files, "threads": threads, "quality": quality,
             "quality_perc": quality_perc, "thread_id": self.thread_id})  # TODO: remove args if not used anymore
        
        utils.thread_print(self.thread_id, "ProcessFastQFilesPipe finished excecuting FastQC and FastXTrimmer!")
        return kwargs

    def process_fastq(self, fastq_files, fastqc_folder, threads, quality, quality_perc):
        utils.thread_print(self.thread_id, "Inside process_fastq function. List of fastq_files = {}".format(fastq_files))
        # list_file = glob.glob(fastq_folder+'*')
        processed_files = []
        for fastq_file in fastq_files:
            utils.thread_print(self.thread_id, "Processing fastq_file ... {}".format(fastq_file))
            if not fastq_file.endswith('new'):
                tools.FastQC(fastq_file, t=threads, o=fastqc_folder)
                quality_filter_res = tools.FastQualityFilter(v=True, Q33=True, q=str(quality), p=str(quality_perc),
                                                             i=fastq_file, capture_output=True)
                #print("quality_filter_res = {}".format(quality_filter_res))
                tools.FastXTrimmer(Q33=True, t=5, m=20, v=True, o=join(fastq_file[:-6] + '_new.fastq'),
                                   input=quality_filter_res[1])
                processed_file = "{}_new.fastq".format(fastq_file[:-6])
                processed_files.append(processed_file)

        utils.thread_print(self.thread_id, "Processed fastq_files: {}".format(processed_files))
        return processed_files


class SetSamplesPipe(ParallelPipe):

    def __init__(self) -> None:
        super().__init__()

    def process(self, **kwargs):

        principal_directory = kwargs.pop("principal_directory")
        fastq_files = kwargs.pop("fastq_files")
        dest = kwargs.pop("dest")
        fastq = kwargs.pop("fastq")
        self.thread_id = kwargs.pop("thread_id")

        df_samples = self.set_samples(principal_directory, fastq_files, dest)

        kwargs.update(
            {"principal_directory": principal_directory, "fastq_files": fastq_files, "dest": dest, "fastq": fastq,
             "samples_dataframe": df_samples, "thread_id": self.thread_id})
        
        self.thread_print("SetSamples Pipe finished execution!")
        return kwargs

    def set_samples(self, principal_directory, fastq_files, dest):
        import pandas as pd

        sample_dict = {}
        for fastq_file in fastq_files:
            fastq_name = fastq_file.split('/')[-1]  # get the name of the fastq_file without the absolute path
            sample_name = fastq_name.split('_')[0]  # get only the sample name e.g. E380.2023

            if sample_name not in sample_dict:
                sample_dict[sample_name] = {'name': sample_name, 'forward': None, 'reverse': None}

            if 'R1' in fastq_name:
                sample_dict[sample_name]['forward'] = fastq_file
            elif 'R2' in fastq_name:
                sample_dict[sample_name]['reverse'] = fastq_file

        #utils.thread_print(self.thread_id, "sample_dict = {}".format(sample_dict))

        df_samples = pd.DataFrame(list(sample_dict.values()))
        df_samples.sort_values(by=['name'], inplace=True)
        df_samples.to_csv(join(principal_directory, "sample_list_{}.csv".format(self.thread_id)), sep='\t', index=False, encoding='utf-8')

        self.thread_print("Samples {} written to: {}".format(sample_dict.keys(), join(principal_directory, "sample_list_{}.csv".format(self.thread_id))))

        return df_samples


class PreAlignmentPipe(Pipe):

    def __init__(self) -> None:
        super().__init__()

    def process(self, **kwargs):
        genome = kwargs.pop('genome')
        df_samples = kwargs.pop('samples_dataframe')
        principal_directory = kwargs.pop('principal_directory')
        self.thread_id = kwargs.pop("thread_id", None)
        queue = kwargs.pop("queue", None)

        resynced_samples = []
        for index, sample in df_samples.iterrows():
            resynced_sample = self.prealignment(principal_directory, genome, sample)
            resynced_samples.append(resynced_sample)

            if queue is not None:
                queue.put(resynced_sample)
                
        print("QUEUE SIZE: {}".format(queue.qsize()))
        kwargs.update({"queue": queue, "samples_dataframe": df_samples, "principal_directory": principal_directory, "resynced_samples": resynced_samples, "genome": genome, "thread_id": self.thread_id})
        return kwargs

    def prealignment(self, principal_directory, genome, sample):
        name = str(sample['name'])
        forward = sample['forward']
        reverse = sample['reverse']

        #utils.thread_print(self.thread_id, sample)

        try:
            # if genome == 'geno38':
            #     os.system(' '.join(['python', '/home/magi/PROJECT/diagnosys/bin_jurgen/resync_fastq.py', forward, reverse]))
            # if genome == 'geno37':
            #     os.system(' '.join(['python', '/home/magi/PROJECT/diagnosys/bin_jurgen/resync_fastq.py', forward, reverse]))

            if genome == 'geno38':
                os.system(' '.join(['/home/magi/miniconda3/envs/PY270/bin/python', '/home/magi/PROJECT/diagnosys/bin/resync_fastq.py', forward, reverse]))
            if genome == 'geno37':
                os.system(' '.join(['/home/magi/miniconda3/envs/PY270/bin/python', '/home/magi/PROJECT/diagnosys/bin/resync_fastq.py', forward, reverse]))


            """ Probably will be moved to another pipe. The following code just manipulates filenames and file organization to remove temp files, as well as to prepare it for Parabricks. """
            #fastq_directory = os.path.join(principal_directory, 'fastq/')
            parabricks_input_directory = os.path.join(principal_directory, 'parabricks_input/')

            """ Get the outputs of resync """
            pairs_forward = forward.split(".")[0] + '.' + forward.split(".")[1] + '.fastq_pairs_R1.fastq'
            pairs_reverse = reverse.split(".")[0] + '.' + reverse.split(".")[1] + '.fastq_pairs_R2.fastq'

            """ Get only the filenames, then build the proper directories that will be accessed by nvidia parabricks and docker. """
            pairs_forward_filename = pairs_forward.split('/')[-1]
            pairs_reverse_filename = pairs_reverse.split('/')[-1]

            singles_forward = forward.split(".")[0] + '.' + forward.split(".")[1] + '.fastq_singles.fastq'
            
            os.system(' '.join(['mv', pairs_forward, parabricks_input_directory]))
            os.system(' '.join(['mv', pairs_reverse, parabricks_input_directory]))

            os.system(' '.join(['rm', forward]))
            os.system(' '.join(['rm', reverse]))
            os.system(' '.join(['rm', singles_forward]))
            
            """NOTE: POTENTIAL FOR HIGH COUPLING BETWEEN DOCKER AND LOCAL DIRECTORY NAMES. THINK OF A BETTER SOLUTION! """
            docker_input_parabricks = os.path.join(config.DOCKER_WORKDIR, 'parabricks_input')
            #read_group = '"@RG\\tID:group1\\tSM:%s\\tLB:%s\\tPL:Illumina"' % (name,name)
            
            read_group = '"@RG\\tID:%s\\tLB:%s\\tPL:Illumina\\tSM:%s\\tPU:%s"' % (name, name, name, name)
            #read_group = '"@RG\\tID:%s\\tSM:%s\\tLB:%s\\tPL:Illumina"' % (name, name, name)
            utils.log_pair_end_FASTQ(os.path.join(principal_directory, 'parabricks_input.txt'), 
                                     os.path.join(docker_input_parabricks, pairs_forward_filename), 
                                     os.path.join(docker_input_parabricks, pairs_reverse_filename),
                                     read_group)

            utils.thread_print(self.thread_id, "Resync completed succesfully!")

            return {"name": name, "pairs_forward": pairs_forward_filename, "pairs_reverse": pairs_reverse_filename}

        except Exception as e:
            utils.thread_print(self.thread_id, "Resync Failed.....!!!")
            utils.thread_print(self.thread_id, "Exception: %s" % str(e))
            sys.exit(1)
