import datetime
import glob
import json
import logging
import os
import re
import shutil
import subprocess
import uuid

from typing import Iterator, Optional

import auto_cpo.fastq as fastq

def find_fastq_dirs(config, check_symlinks_complete=True):
    miseq_run_id_regex = "\d{6}_M\d{5}_\d+_\d{9}-[A-Z0-9]{5}"
    nextseq_run_id_regex = "\d{6}_VH\d{5}_\d+_[A-Z0-9]{9}"
    fastq_by_run_dir = config['fastq_by_run_dir']
    subdirs = os.scandir(fastq_by_run_dir)
    analysis_base_outdir = config['analysis_output_dir']
    for subdir in subdirs:
        run_id = subdir.name
        analysis_outdir = os.path.abspath(os.path.join(analysis_base_outdir, run_id))
        matches_miseq_regex = re.match(miseq_run_id_regex, run_id)
        matches_nextseq_regex = re.match(nextseq_run_id_regex, run_id)
        if check_symlinks_complete:
            ready_to_analyze = os.path.exists(os.path.join(subdir.path, "symlinks_complete.json"))
        else:
            ready_to_analyze = True
        analysis_not_already_initiated = not os.path.exists(analysis_outdir)
        conditions_checked = {
            "is_directory": subdir.is_dir(),
            "matches_illumina_run_id_format": ((matches_miseq_regex is not None) or (matches_nextseq_regex is not None)),
            "ready_to_analyze": ready_to_analyze,
            "analysis_not_already_initiated": analysis_not_already_initiated,
        }
        conditions_met = list(conditions_checked.values())
        pipeline_parameters = {}

        if all(conditions_met):
            logging.info(json.dumps({"event_type": "fastq_directory_found", "sequencing_run_id": run_id, "fastq_directory_path": os.path.abspath(subdir.path)}))
            pipeline_parameters['fastq_input'] = os.path.abspath(subdir.path)
            yield pipeline_parameters
        else:
            logging.debug(json.dumps({"event_type": "directory_skipped", "fastq_directory_path": os.path.abspath(subdir.path), "conditions_checked": conditions_checked}))
            yield None
    

def scan(config: dict[str, object]) -> Iterator[Optional[dict[str, object]]]:
    """
    Scanning involves looking for all existing runs and storing them to the database,
    then looking for all existing symlinks and storing them to the database.
    At the end of a scan, we should be able to determine which (if any) symlinks need to be created.

    :param config: Application config.
    :type config: dict[str, object]
    :return: A run directory to analyze, or None
    :rtype: Iterator[Optional[dict[str, object]]]
    """
    logging.info(json.dumps({"event_type": "scan_start"}))
    for symlinks_dir in find_fastq_dirs(config):    
        yield symlinks_dir


def analyze_run(config, run):
    """
    Initiate an analysis on one directory of fastq files.
    """
    base_analysis_outdir = config['analysis_output_dir']
    base_analysis_work_dir = config['analysis_work_dir']
    no_value_flags_by_pipeline_name = {
        "BCCDC-PHL/taxon-abundance": [],
        "BCCDC-PHL/routine-assembly": ["unicycler", "hybrid", "bakta", "prokka"],
        "BCCDC-PHL/mlst-nf": [],
        "BCCDC-PHL/plasmid-screen": [],
    }
    if 'notification_email_addresses' in config:
        notification_email_addresses = config['notification_email_addresses']
    else:
        notification_email_addresses = []
    for pipeline in config['pipelines']:
        pipeline_parameters = pipeline['pipeline_parameters']
        pipeline_short_name = pipeline['pipeline_name'].split('/')[1]
        pipeline_minor_version = ''.join(pipeline['pipeline_version'].rsplit('.', 1)[0])
        analysis_output_dir_name = '-'.join([pipeline_short_name, pipeline_minor_version, 'output'])
        outdir = os.path.abspath(os.path.join(config['analysis_output_dir'], os.path.basename(run['fastq_input']), analysis_output_dir_name))
        pipeline_parameters['outdir'] = outdir
        if pipeline['pipeline_name'] == 'BCCDC-PHL/taxon-abundance':
            run_fastq_files = glob.glob(os.path.join(run['fastq_input'], '*.f*q.gz'))
            first_fastq = None
            if len(run_fastq_files) > 0:
                first_fastq = run_fastq_files[0]
            else:
                first_fastq = None
                logging.error(json.dumps({"event_type": "find_fastq_files_failed", "fastq_directory_path": os.path.abspath(run['fastq_input'])}))
                return None

            read_length = fastq.estimate_read_length(fastq.get_first_n_reads(first_fastq, 100))
            pipeline_parameters['read_length'] = str(read_length)
        
        analysis_timestamp = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
        analysis_run_id = os.path.basename(run['fastq_input'])
        analysis_work_dir = os.path.abspath(os.path.join(base_analysis_work_dir, 'work-' + analysis_run_id + '-' + analysis_timestamp))
        analysis_report_path = os.path.abspath(os.path.join(outdir, analysis_run_id + '_report.html'))
        analysis_trace_path = os.path.abspath(os.path.join(outdir, analysis_run_id + '_trace.txt'))
        analysis_timeline_path = os.path.abspath(os.path.join(outdir, analysis_run_id + '_timeline.html'))
        pipeline_command = [
            'nextflow',
            'run',
            pipeline['pipeline_name'],
            '-r', pipeline['pipeline_version'],
            '-profile', 'conda',
            '--cache', os.path.join(os.path.expanduser('~'), '.conda/envs'),
            '-work-dir', analysis_work_dir,
            '-with-report', analysis_report_path,
            '-with-trace', analysis_trace_path,
            '-with-timeline', analysis_timeline_path,
        ]
        if 'send_notification_emails' in config and config['send_notification_emails']:
            pipeline_command += ['-with-notification', ','.join(notification_email_addresses)]
        for flag, config_value in pipeline_parameters.items():
            if config_value is None and flag not in no_value_flags_by_pipeline_name[pipeline['pipeline_name']]:
                value = run[flag]
                pipeline_command += ['--' + flag, value]
            elif config_value is None and flag in no_value_flags_by_pipeline_name[pipeline['pipeline_name']]:
                pipeline_command += ['--' + flag]
            else:
                value = config_value
                pipeline_command += ['--' + flag, value]
            
        logging.info(json.dumps({"event_type": "analysis_started", "sequencing_run_id": analysis_run_id, "pipeline_command": " ".join(pipeline_command)}))
        try:
            subprocess.run(pipeline_command, capture_output=True, check=True)
            logging.info(json.dumps({"event_type": "analysis_completed", "sequencing_run_id": analysis_run_id, "pipeline_command": " ".join(pipeline_command)}))
            shutil.rmtree(analysis_work_dir, ignore_errors=True)
            logging.info(json.dumps({"event_type": "analysis_work_dir_deleted", "sequencing_run_id": analysis_run_id, "analysis_work_dir_path": analysis_work_dir}))
        except subprocess.CalledProcessError as e:
            logging.error(json.dumps({"event_type": "analysis_failed", "sequencing_run_id": analysis_run_id, "pipeline_command": " ".join(pipeline_command)}))
        except OSError as e:
            logging.error(json.dumps({"event_type": "delete_analysis_work_dir_failed", "sequencing_run_id": analysis_run_id, "analysis_work_dir_path": analysis_work_dir}))
