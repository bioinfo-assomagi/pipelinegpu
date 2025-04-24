class AlignmentPipe(Pipe):

    def __init__(self) -> None:
        super().__init__()

    def process(self, **kwargs):
        genome = kwargs.pop('genome')
        df_samples = kwargs.pop('samples_dataframe')
        principal_directory = kwargs.pop('principal_directory')
        threads = kwargs.pop('threads')

        for index, sample in df_samples.iterrows():
            self.Alignment(genome, threads, principal_directory, sample)

        kwargs.update({'genome': genome, 'samples_dataframe': df_samples, 'principal_directory': principal_directory})
        return kwargs

    def Alignment(self, genome, threads, principal_directory, sample):
        name = str(sample['name'])
        forward = sample['forward']
        reverse = sample['reverse']

        print("Running BWA for sample {}".format(name))

        try:
            sample_sam_ = join(principal_directory, 'temp/', name + '.sam')
            header = '"@RG\\tID:group1\\tSM:%s\\tLB:%s\\tPL:Illumina"' % (name, name)

            if genome == 'geno38':
                # TODO: use tools.bwa, nvidia Parabricks will be used

                os.system(' '.join([os.path.join(config.BWA, 'bwa'), 'mem', '-S', '-w 150', '-R', header,
                                    '-t', threads, genome, forward, reverse,
                                    '>', sample_sam_]))

            if genome == 'geno37':
                os.system(' '.join([os.path.join(config.BWA, 'bwa'), 'mem', '-R', header,
                                    '-t', threads, genome, forward, reverse,
                                    '>', sample_sam_]))



        except Exception as e:
            print("You need to run.....!!!")
            print("Exception: %s" % str(e))
            sys.exit(1)



class SAMToBAMPipe(Pipe):

    def __init__(self) -> None:
        super().__init__()

    def process(self, **kwargs):
        genome = kwargs.pop('genome')
        df_samples = kwargs.pop('samples_dataframe')
        principal_directory = kwargs.pop('principal_directory')

        for index, sample in df_samples.iterrows():
            self.samtobam(genome, principal_directory, sample)

        return kwargs.update({'genome': genome, 'samples_dataframe': df_samples, 'principal_directory': principal_directory})

    def samtobam(self, genome, principal_directory, sample):
        name = str(sample['name'])
        sample_sam_ = join(principal_directory, 'temp/', name + '.sam')
        sample_bam_ = join(principal_directory, 'temp/', name + '.bam')

        os.system(' '.join([os.path.join(config.SAMTOOLS, 'samtools'), 'faidx', genome]))
        os.system(' '.join([os.path.join(config.SAMTOOLS, 'samtools'), 'view', '-bS', sample_sam_, '>', sample_bam_]))







def test_stuff():
    # DBContext().setup_database('r')
    # db = DBContext().database
    from DBContext import DBContext
    #samples = ['home/S1_R1_new.fastq', 'home/S1_R2_new.fastq', 'home/S2_R1_new.fastq', 'home/S2_R2_new.fastq']
    db1 = DBContext('')
    db1.setup_database('r')
    results = db1.get_disease(['105.2016', '436.2015'])
    print(results)

def test_setsample():
    import pandas as pd
    samples = ['home/S1_R1_new.fastq', 'home/S1_R2_new.fastq', 'home/S2_R1_new.fastq', 'home/S2_R2_new.fastq']

    sample_dict = {}

    for sample in samples:
        fastq_name = sample.split('/')[-1]  # get the name of the fastq_file without the absolute path
        sample_name = fastq_name.split('_')[0]  # get only the sample name e.g. E380.2023

        if sample_name not in sample_dict:
            sample_dict[sample_name] = {'name': sample_name, 'forward': 'na', 'reverse': 'na'}

        if 'R1' in fastq_name:
            sample_dict[sample_name]['forward'] = sample
        elif 'R2' in fastq_name:
            sample_dict[sample_name]['reverse'] = sample

    df_samples = pd.DataFrame(list(sample_dict.values()))

    for index, sample in df_samples.iterrows():
        print(sample['name'])
        print(sample['forward'])
        print(sample['reverse'])