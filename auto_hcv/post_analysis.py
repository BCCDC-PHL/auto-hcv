import csv
import datetime
import glob
import json
import logging
import os
from os.path import join as pathjoin
import shutil


# Default formats of sample files to be transferred
DEFAULT_SAMPLE_FILES = [
	("depth_plot","{}_depth_plots.png"),
	("tree_core","RAxML_bestTree.{}_core"),
	("tree_ns5b","RAxML_bestTree.{}_ns5b"),
	("genotype_calls","{}_genotype_calls_nt.csv"),
	("consensus_seqs","{}_consensus_seqs.fa"),
]

# Function to transfer a file from source to destination
def transfer_file(src_path, dest_path):
	try: 
		#logging.info(f'Copying {src_path} to {dest_path}')
		shutil.copy2(src_path, dest_path)
	except FileNotFoundError as e :
		# Log a warning if the file does not exist
		logging.warning(json.dumps({"event_type": "transfer_file_does_not_exist", "file": src_path}))
		

# Function to transfer HCV results
def transfer_hcv_results(config, pipeline, run, fstring_list=DEFAULT_SAMPLE_FILES):
	# Get the pipeline details from the configuration
	pipeline_short_name = pipeline['pipeline_name'].split('/')[1]
	pipeline_minor_version = ''.join(pipeline['pipeline_version'].rsplit('.', 1)[0])
	pipeline_path_name = '-'.join([pipeline_short_name, pipeline_minor_version, 'output'])

	# Set the source and destination paths for file transfer
	src_path = pathjoin(config['analysis_output_dir'], run['run_id'], pipeline_path_name)
	dest_path = pathjoin(config['analysis_report_dir'], pipeline_path_name, run['run_id'], datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))

	# Check if the run already exists in the summary folder
	if os.path.isfile(pathjoin(dest_path, "transfer_complete.json")):
		logging.warning(json.dumps({"event_type": "transfer_hcv_folder_exists_complete","sequencing_run_id": run['run_id'],"pipeline_name": pipeline['pipeline_name']}))
		return
	elif os.path.isdir(dest_path):
		logging.error(json.dumps({"event_type": "transfer_hcv_folder_exists_incomplete","sequencing_run_id": run['run_id'],"pipeline_name": pipeline['pipeline_name']}))
		return

	# Record the start time of file transfer
	transfer_complete = {'timestamp_transfer_start': datetime.datetime.now().isoformat()}
	logging.info(json.dumps({"event_type": "transfer_hcv_results_start","sequencing_run_id": run['run_id'],"pipeline_name": pipeline['pipeline_name']}))
	
	try: 
		# Create the destination directory
		os.makedirs(dest_path)

		# Transfer the summary report for the run
		logging.debug(json.dumps({"event_type": "transfer_hcv_summary_report","sequencing_run_id": run['run_id'],"pipeline_name": pipeline['pipeline_name']}))

		transfer_file(pathjoin(src_path, run['run_id'] + "_run_summary_report.csv"), pathjoin(dest_path, run['run_id'] + "_run_summary_report.csv"))
		
		# Transfer sample folders for the run
		logging.debug(json.dumps({"event_type": "transfer_hcv_sample_folders","sequencing_run_id": run['run_id'],"pipeline_name": pipeline['pipeline_name']}))
		sample_dirs = [x for x in os.listdir(src_path) if os.path.isdir(pathjoin(src_path, x))]

		# Iterate through each sample folder 
		for sample_name in sample_dirs:
			os.makedirs(pathjoin(dest_path, sample_name))

			# Copy each individual file
			for str_name, fstring in fstring_list:
				if str_name == 'depth_plot':
					name = sample_name.replace('-','o')
				else:
					name = sample_name
			
				filename = fstring.format(name)

				transfer_file(pathjoin(src_path, sample_name, filename), pathjoin(dest_path, sample_name, filename))

		# Record the completion time of file transfer
		transfer_complete['timestamp_transfer_complete'] = datetime.datetime.now().isoformat()

		# Write the transfer completion details to a JSON file
		with open(pathjoin(dest_path, 'transfer_complete.json'), 'w') as f:
			json.dump(transfer_complete, f, indent=2)

		logging.info(json.dumps({"event_type": "transfer_hcv_results_complete","sequencing_run_id": run['run_id'],"pipeline_name": pipeline['pipeline_name']}))
	
	except FileExistsError as e:
		logging.error(json.dumps({"event_type": "transfer_hcv_directory_exists_error","sequencing_run_id": run['run_id'],"pipeline_name": pipeline['pipeline_name'], 'error': e}))
		

def post_analysis_hcv_nf(config, pipeline, run):
	logging.debug(json.dumps({
				"event_type": "post_analysis_hcv_nf_start",
				"sequencing_run_id": run['run_id'],
				"pipeline_name": pipeline['pipeline_name']
			}))

	transfer_hcv_results(config, pipeline, run)

def find_latest_glob(path):
	list_dirs = glob.glob(path)
	if len(list_dirs) > 0:
		return list_dirs[-1]
	else:
		return None
	
def post_analysis(config, pipeline, run):
	"""
	Perform post-analysis tasks for a pipeline.

	:param config: The config dictionary
	:type config: dict
	:param pipeline: The pipeline dictionary
	:type pipeline: dict
	:param run: The run dictionary
	:type run: dict
	:return: None
	"""

	pipeline_name = pipeline['pipeline_name']
	pipeline_short_name = pipeline_name.split('/')[1]
	pipeline_version = pipeline['pipeline_version']
	sequencing_run_id = run['run_id']
	base_analysis_work_dir = config['analysis_work_dir']

	# The work_dir includes a timestamp, so we need to glob to find the most recent one
	work_dir_glob = os.path.join(base_analysis_work_dir, 'work-' + sequencing_run_id + '_' + pipeline_short_name + '_' + '*')
	work_dir = find_latest_glob(work_dir_glob)

	# Remove the working directory tree
	if work_dir:
		try:
			shutil.rmtree(work_dir, ignore_errors=True)
			logging.info(json.dumps({
				"event_type": "analysis_work_dir_deleted",
				"sequencing_run_id": sequencing_run_id,
				"analysis_work_dir_path": work_dir
			}))
		except OSError as e:
			logging.error(json.dumps({
				"event_type": "delete_analysis_work_dir_failed",
				"sequencing_run_id": sequencing_run_id,
				"analysis_work_dir_path": work_dir
			}))
	else:
		logging.warning(json.dumps({
			"event_type": "analysis_work_dir_not_found",
			"sequencing_run_id": sequencing_run_id,
			"analysis_work_dir_glob": work_dir_glob
		}))

	# a dictionary mapping pipeline names to their functions
	post_analysis_fn_map = {
		'BCCDC-PHL/hcv-nf': post_analysis_hcv_nf
	}

	if pipeline_name in post_analysis_fn_map :
		logging.info(json.dumps({
			"event_type": f"post_analysis_{pipeline_name}",
			"sequencing_run_id": sequencing_run_id,
			"pipeline_name": pipeline_name
		}))
		return post_analysis_fn_map[pipeline_name](config, pipeline, run)

	else:
		logging.warning(json.dumps({
			"event_type": "post_analysis_not_implemented",
			"sequencing_run_id": sequencing_run_id,
			"pipeline_name": pipeline_name
		}))
		return None
