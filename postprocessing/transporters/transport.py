#!/usr/bin/env python3
"""
postprocessing/transporters/transport.py

Move pipeline outputs from RESULT_DIR to the staging directory
used by downstream parsing-and-DB-import scripts.
"""

import os
import sys
import shutil
import logging
from pathlib import Path

from ...LogFormatter import ColorFormatter

# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------
# You could also pull these in from a global settings module
logger = logging.getLogger(__name__)


# logging.basicConfig(
#     level=getattr(logging, LOG_LEVEL),
#     format="%(asctime)s.%(msecs)03d %(levelname)-8s [%(name)s:%(lineno)d] %(message)s",
#     datefmt="%Y-%m-%d %H:%M:%S"
# )

# -----------------------------------------------------------------------------
# MAIN FUNCTIONS
# -----------------------------------------------------------------------------
def discover_result_files(result_dir: Path) -> list[Path]:
    """
    Scan RESULT_DIR and return a list of files to transport.
    Could filter by extension, age, naming pattern, etc.
    """
    # TODO: refine filtering logic

    # get files of interest in annot
    final_files = [f for f in (result_dir/'final').iterdir() if f.is_file()]

    # get files of interest in coverage

    sample_coverage_dirs = [d for d in (result_dir/'coverage').iterdir() if d.is_dir()]
    coverage_files = []
    for cov_dir in sample_coverage_dirs:
        coverage_files.extend([f for f in cov_dir.iterdir() if f.is_file()])


    # get files of interest in indel
    sample_indel_dirs = [d for d in (result_dir/'indel').iterdir() if d.is_dir()]
    indel_files = []
    for indel_dir in sample_indel_dirs:
        indel_files.extend([f for f in indel_dir.iterdir() if f.is_file()])


    # get files of interest from report
    report_files = [f for f in (result_dir/'report').iterdir() if f.is_file()]

    result_dict = {'final': final_files, 'coverage': coverage_files, 'indel': indel_files, 'report': report_files}
    return result_dict


def filter_files(orig_filepath_list: list, suffixes: tuple) -> list:
    """
    Given a list of filepaths, return a list of filepaths filtered by suffix.
    """
    filtered_list = []
    for filepath in orig_filepath_list:

        if any(filepath.name.endswith(suffix) for suffix in suffixes):
            filtered_list.append(filepath)

    return filtered_list

def determine_destination(filepaths_dict: dict, staging_dir: Path) -> dict:
    """
    Given a source file, decide where in STAGING_DIR it should go.
    E.g. group by file suffix or prefix.
    """

    filepaths_final = filepaths_dict.get('final', [])
    filepaths_coverage = filepaths_dict.get('coverage', [])
    filepaths_indel = filepaths_dict.get('indel', [])
    filepaths_report = filepaths_dict.get('report', [])
    
    print()
    logger.info("Found %d final files", len(filepaths_final))
    logger.info("Found %d coverage files", len(filepaths_coverage))
    logger.info("Found %d indel files", len(filepaths_indel))
    logger.info("Found %d report files", len(filepaths_report))


    # keep only files with suffix pheno_predict.csv, variant_selection.csv, pheno_annot.csv, family_annot.csv, other_annot.csv
    filtered_final = filter_files(filepaths_final, ("pheno_predict.csv", "variant_selection.csv", "pheno_annot.csv", "family_annot.csv", "other_annot.csv"))

    # keep only files with suffix stat_cov.csv, gene_cov.csv
    filtered_coverage = filter_files(filepaths_coverage, ('stat_cov.csv', 'gene_cov.csv'))
    
    # keep only files with suffix final_indel.csv
    filtered_indel = filter_files(filepaths_indel, ('final_indel.csv',))

    # take everything - only one should be there, others are inside sample dict which don't interest us for now
    filtered_report = filter_files(filepaths_report, ('.csv',))

    # keep only files with suffix interpretation_varsome.csv
    filtered_checklist = filter_files(filepaths_final, ('interpretation_varsome.csv',))

    print()
    logger.info("Maintained %d final files", len(filtered_final))
    logger.info("Maintained %d coverage files", len(filtered_coverage))
    logger.info("Maintained %d indel files", len(filtered_indel))
    logger.info("Maintained %d report files", len(filtered_report))
    logger.info("Maintained %d checklist files", len(filtered_checklist))


    destination_dict = {Path('annot'): filtered_final, Path('coverage'): filtered_coverage, Path('INDEL'): filtered_indel, Path('REFERTAZIONE'): filtered_report, Path('CHECKLIST'): filtered_checklist}

    
    return destination_dict

def transport_file(dst_dict: dict, destination_dir: Path, overwrite: bool = True) -> None:
    """
    Actually copy (or move) the file from src to dst.
    Ensures target directory exists.
    """
    print()
    for dst, files in dst_dict.items():
        (destination_dir/(dst.name)).mkdir(parents=True, exist_ok=True)
        for file in files:
            if overwrite and (destination_dir/(dst.name)/file.name).exists():
                (destination_dir/(dst.name)/file.name).unlink()

            shutil.copy2(file, destination_dir/(dst.name))
            
            logger.info("Copied %s to directory %s ... ", file.name, destination_dir/(dst.name))

def run_transport(RESULT_DIR: Path, STAGING_DIR: Path):
    """
    Main entrypoint: discover files, transport them, and report.
    """
    result_dict = discover_result_files(RESULT_DIR)
    if not result_dict:
        logger.warning("No files found in %s", RESULT_DIR)
        return

    try:
        dst_dict = determine_destination(result_dict, STAGING_DIR)
        transport_file(dst_dict, STAGING_DIR)
    except Exception as e:
        logger.error("Failed to transport %s: %s", result_dict, e)


def tests():
    pass
    #determine_destination(discover_result_files(RESULT_DIR), STAGING_DIR)
    #run_transport(RESULT_DIR, STAGING_DIR)

# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    run_transport(RESULT_DIR, STAGING_DIR)
