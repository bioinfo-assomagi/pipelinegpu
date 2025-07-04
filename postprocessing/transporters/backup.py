import os
from os import system
from pathlib import Path

ROVERETO_RAW        = "admin@192.168.1.202:/share/CACHEDEV1_DATA/homes/storage/NGSDATA/NGSDATA_GPU/RAW/ROVERETO"
ROVERETO_RESULTS    = "admin@192.168.1.202:/share/CACHEDEV1_DATA/homes/storage/NGSDATA/NGSDATA_GPU/result/ROVERETO"

BOLZANO_RAW         = "admin@192.168.1.202:/share/CACHEDEV1_DATA/homes/storage/NGSDATA/NGSDATA_GPU/RAW/BOLZANO"
BOLZANO_RESULTS     = "admin@192.168.1.202:/share/CACHEDEV1_DATA/homes/storage/NGSDATA/NGSDATA_GPU/result/BOLZANO"

RICERCA_RAW         = "admin@192.168.1.202:/share/CACHEDEV1_DATA/homes/storage/NGSDATA/NGSDATA_GPU/RAW/RICERCA"
RICERCA_RESULTS     = "admin@192.168.1.202:/share/CACHEDEV1_DATA/homes/storage/NGSDATA/NGSDATA_GPU/result/RICERCA"

def run_backup_results(destination: str, res_dir: str):
    """Backup result folders (vcf, pheno, indel, final, control, annotation)
    contained in *res_dir* to the NAS result share corresponding to
    *destination*."""
    if destination == "r":
        destination_path = ROVERETO_RESULTS
    elif destination == "b":
        destination_path = BOLZANO_RESULTS
    elif destination == "z":
        destination_path = RICERCA_RESULTS
    else:
        raise ValueError(f"Unknown destination {destination}")

    res_dir_path = Path(res_dir).resolve()
    parent_name = res_dir_path.name  # e.g. TEST_SPAWN_NOFAM_OCULARE
    remote_parent_dir = f"{destination_path}/{parent_name}/"

    for folder in ["vcf", "pheno", "indel", "final", "control", "annotation"]:
        src_folder = res_dir_path / folder
        # Copy each folder into the dedicated parent directory on the NAS, keeping
        # the folder name (vcf, pheno, ...).
        command = f'rsync -avzhe ssh "{src_folder}" "{remote_parent_dir}"'
        print(command)
        system(command)

def run_backup_fastq(destination: str, fq_dir: str):
    """Backup FASTQ (and fastq.gz) files found in *fq_dir* to the RAW NAS share
    corresponding to *destination*."""
    if destination == "r":
        destination_path = ROVERETO_RAW
    elif destination == "b":
        destination_path = BOLZANO_RAW
    elif destination == "z":
        destination_path = RICERCA_RAW
    else:
        raise ValueError(f"Unknown destination {destination}")

    # Copy all FASTQ files (and any other contents) from the source directory into
    # the remote RAW destination, keeping the same filenames.
    # A trailing slash after the source directory copies only its content.
    p = Path(fq_dir)
    command = f'rsync -avzhe ssh "{p}" "{destination_path}/"'
    print(command)
    system(command)


def run_backup(dest, result_dir, fq_dir):
    run_backup_results(dest, result_dir)
    run_backup_fastq(dest, fq_dir)
    


if __name__ == "__main__":
    run_backup("r", "/share/CACHEDEV1_DATA/homes/storage/NGSDATA/NGSDATA_GPU/RAW/ROVERETO", "/share/CACHEDEV1_DATA/homes/storage/NGSDATA/NGSDATA_GPU/RAW/ROVERETO")