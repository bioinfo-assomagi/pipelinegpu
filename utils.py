import logging
import config # TODO: utils should not know about config, fix it
import os
import json


def get_buchiartificiali():
    import pandas as pd
    return pd.read_csv(config.BUCHIARTIFICIALI, sep="\t", header=0)

def remove_buchi(df):
    buchiartificiali = get_buchiartificiali() # TODO: provide it as input?
    for index, row in buchiartificiali.iterrows():
        chrom_num = row['#CHROM']
        start = row['START']
        end = row['END']
        mask1 = ((df['#CHROM'] == chrom_num) & (df['POS'] >= start) & (df['POS'] <= end))
        df.loc[mask1,'DEPTH'] = 0
    
    return df

def get_vertical_macro(panel):
    if "OCULAR" in panel:
        _vertical_macro = config.OCULAR
    elif "VASCULAR" in panel:
        _vertical_macro = config.VASCULAR
    elif "NEUROLOGY" in panel:
        _vertical_macro = config.NEUROLOGY
    elif "INTEGRACARDIOSTANCHEZZA" in panel:
        _vertical_macro = config.INTEGRACARDIOSTANCHEZZA
    elif "MIXED" in panel:
        _vertical_macro = config.MIXED
    elif "LYMPHOBESITY" in panel:
        _vertical_macro = config.LYMPHOBESITY
    elif "INFERTILIT" in panel:
        _vertical_macro = config.INFERTILIT
    elif "CHERATOCONO" in panel:
        _vertical_macro = config.CHERATOCONO
    elif "PCDH19" in panel:
        _vertical_macro = config.PCDH19
    elif "GENEOBNEW" in panel:
        _vertical_macro = config.GENEOBNEW
    elif "CUSTOM" in panel:
        _vertical_macro = config.CUSTOM
    elif "CANCER" in panel:
        _vertical_macro = config.CANCER
    else:
        _vertical_macro = None

    return _vertical_macro 


def sample_json(sample_json_path):
    with open(sample_json_path) as f:
        d = json.load(f)
    return d
    

def get_db_path(server_id):
    if server_id == "b":
        return config.DB_EUREGIO
    elif server_id == "r":
        return config.DB_MAGI
    elif server_id == "z":
        return config.DB_RICERCA


def get_db_server(server_id):
    if server_id == "b":
        return config.SERVER_EUREGIO
    elif server_id == "r":
        return config.SERVER_MAGIS
    elif server_id == "z":
        return config.SERVER_RICERCA


def createProjectName():
    import datetime

    today = datetime.date.today()
    return "{:%d_%b_%Y}".format(today)


def log_pair_end_FASTQ(filename, r1, r2, read_group):
    logpath = filename
    logger = logging.getLogger("log")
    logger.setLevel(logging.INFO)
    ch = logging.FileHandler(logpath)
    ch.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(ch)
    logger.info("{} {} {}".format(r1, r2, read_group))


def thread_print(thread_id, msg):
    # print(multiprocessing.current_process().name)
    if thread_id is not None:
        print("{}THREAD-{}{}: {}".format("\033[91m", thread_id, "\033[0m", msg))
    else:
        print(msg)


def group_samples(fastq_files):

    sample_dict = {}

    for fastq_file in fastq_files:
        sample_name = fastq_file.split("/")[-1].split("_")[0]

        if sample_name not in sample_dict:
            sample_dict[sample_name] = {
                "name": str(sample_name),
                "forward": None,
                "reverse": None,
            }

        if "R1" in fastq_file:
            sample_dict[sample_name]["forward"] = fastq_file
        elif "R2" in fastq_file:
            sample_dict[sample_name]["reverse"] = fastq_file

    return sample_dict


def merge_csv(csv_list):
    import pandas as pd

    df_concat = pd.concat(
        [pd.read_csv(f, dtype=str, sep="\t") for f in csv_list], ignore_index=True
    )

    return df_concat


def unwrap_queue(queue):
    l = []
    while queue.qsize() > 0:
        l.append(queue.get())

    return l


if __name__ == "__main__":
    # tests
    print(get_vertical_macro("CANCER"))


class Sequence:

    def __init__(self, sequence: str, name=None) -> None:
        self.sequence = sequence
        self.name = name
        self.length = len(sequence)

    def remove_common_prefix(self, seq2):
        commonprefix = os.path.commonprefix([self.sequence, seq2.sequence])

        noprefix_seq1 = self.sequence[len(commonprefix) :]
        noprefix_seq2 = seq2.sequence[len(commonprefix) :]

        commonprefix_reverse = os.path.commonprefix(
            [noprefix_seq1[::-1], noprefix_seq2[::-1]]
        )
        nosuffix_seq1 = noprefix_seq1[: len(noprefix_seq1) - len(commonprefix_reverse)]
        nosuffix_seq2 = noprefix_seq2[: len(noprefix_seq2) - len(commonprefix_reverse)]

        return {self.name: nosuffix_seq1, seq2.name: nosuffix_seq2}

    def find_indel(self, seq2):
        start_index = None
        end_index = None

        for i in range(len(self.sequence)):

            if start_index == None:
                if self.sequence[i] != seq2.sequence[i]:
                    start_index = i
                if i + 1 >= seq2.length:
                    start_index = i + 1

        if end_index == None:
            if self.sequence[self.length - i - 1] != seq2.sequence[seq2.length - i - 1]:
                end_index = i

        if start_index is not None and end_index is not None:
            return {
                self.name: self.sequence[start_index : self.length - end_index],
                seq2.name: seq2.sequence[start_index : seq2.length - end_index],
            }




