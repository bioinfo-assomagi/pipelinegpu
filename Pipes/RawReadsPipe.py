from Pipes.Pipe import Pipe


"""
TODO:
first_diagnosys_core_V2.py - 
get the copy_fastq function, germinalprot; files_fq = system(' '.join(['scp root@192.168.1.51:/home/NGS/tmp/analysis/germinalprot/*',fastq_folder])) in germinalprot there are 3-4 fastq, each fastq has two files, R1 and R2
formation of filenames
fastq_quality_filter, fastx_trimmer

connect_database to get the list of genes
resync_fastq 

then BWA
"""

class RawReadsPipe(Pipe):
    def __init__(self):
        print("RawReadsPipe init!")
    def process(self, **kwargs):
        name_folder = kwargs.pop("name_folder", None)
        if name_folder:
            print("RawReadPipe is reading the {} directory ...".format(name_folder))

        return {"name_folder" : name_folder, "dummy_var" : 'dummy_var1'}
