import os
from os.path import join as pathjoin
import logging
from glob import glob
import re
import shutil as sh
from .core import transfer_run_folder

REGEX_COMPLETED_FOLDER = re.compile(r"\d{6}_[A-Z0-9]+_\d{3}_[A-Z0-9]{9}\/hcv-nf.+$")


def get_completed_runs(analysis_dir):
    run_folders = glob(pathjoin(analysis_dir, "*", "hcv-nf-*"))

    completed_runs = []

    for run_folder in run_folders:
        # print(run_folder)

        # simplify to the relevant base path
        run_folder = "/".join(run_folder.split("/")[-2:])

        if not REGEX_COMPLETED_FOLDER.search(run_folder):
            continue

        analysis_json_path = pathjoin(
            analysis_dir, run_folder, "analysis_complete.json"
        )

        if os.path.isfile(analysis_json_path):
            completed_runs.append(run_folder)
        else:
            logging.info(
                f"Skipping {run_folder}. analysis_complete.json does not exist."
            )

    return set(completed_runs)


def get_transferred_runs(transfer_path):
    report_paths = glob(pathjoin(transfer_path, "hcv-nf*", "*", "*_summary_report.csv"))
    transferred_runs = [x.split("/")[-2] + "/" + x.split("/")[-3] for x in report_paths]
    transferred_runs = [x for x in transferred_runs if REGEX_COMPLETED_FOLDER.search(x)]
    return set(transferred_runs)


def read_exclude_list(path):
    with open(path, "r") as infile:
        entries = {x.strip() for x in infile.readlines()}
    return entries


def run_auto_detect(analysis_path, transfer_path):
    logging.debug(f"Detecting runs")

    completed_runs = get_completed_runs(analysis_path)
    transferred_runs = get_transferred_runs(transfer_path)

    logging.debug(f"Reading exclude list")

    if not os.path.isfile("exclude_list.txt"):
        open("exclude_list.txt", "a").close()

    exclude_list = read_exclude_list("exclude_list.txt")

    logging.info("Starting scan...")

    runs_to_transfer = completed_runs - transferred_runs - exclude_list

    logging.debug(f"Runs to transfer: {runs_to_transfer}")

    for run_name, pipeline_name in [x.split("/") for x in runs_to_transfer]:
        logging.info(
            f"Detected run {run_name} for pipeline {pipeline_name}. Starting transfer."
        )
        transfer_run_folder(
            pathjoin(analysis_path, run_name, pipeline_name),
            pathjoin(transfer_path, pipeline_name, run_name),
        )
