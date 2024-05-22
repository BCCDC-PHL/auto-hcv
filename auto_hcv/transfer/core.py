import os, sys 
from os.path import join as pathjoin
import logging
import json 
import shutil as sh
from datetime import datetime

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
		sh.copy2(src_path, dest_path)
	except FileNotFoundError as e :
        # Log a warning if the file does not exist
		logging.warning(f"Failed to copy {src_path}. File does not exist. {e}")
	except Exception as e :
		logging.error(f"ERROR: An unexpected error occurred: {e}")

# Function to transfer HCV results
def transfer_hcv_results(config, run, fstring_list=DEFAULT_SAMPLE_FILES):

    # Get the pipeline details from the configuration
	pipeline = config['pipelines'][0]
	pipeline_name = "-".join([pipeline['pipeline_name'].split("/")[-1], pipeline['pipeline_version'], 'output'])

    # Set the source and destination paths for file transfer
	src_path = pathjoin(config['analysis_output_dir'], run['run_id'], pipeline_name)
	dest_path = pathjoin(config['analysis_summary_dir'], pipeline_name, run['run_id'])

    # Check if the run already exists in the summary folder
	if os.path.isdir(dest_path) or os.path.isfile(pathjoin(dest_path, "transfer_complete.json")):
		logging.info(f"Run {run['run_id']} already exists in summary folder. Skipping.")
		return

    # Record the start time of file transfer
	transfer_complete = {'timestamp_transfer_start': datetime.now().isoformat()}
	logging.info(f"Beginning file transfer for run {run['run_id']} at {datetime.now().isoformat()}")
	
	# Create the destination directory

	os.makedirs(pathjoin(dest_path))

    # Transfer the summary report for the run
	logging.info(f"Transferring summary report for run {run['run_id']}")
	transfer_file(pathjoin(src_path, run['run_id'] + "_run_summary_report.csv"), pathjoin(dest_path, run['run_id'] + "_run_summary_report.csv"))
	
	# Transfer sample folders for the run
	logging.info(f'Transferring sample folders for run {run['run_id']}')
	sample_dirs = [x for x in os.listdir(src_path) if os.path.isdir(pathjoin(src_path, x))]

	for sample_name in sample_dirs:
		os.makedirs(pathjoin(dest_path, sample_name))

		for str_name, fstr in fstring_list:
			if str_name == 'depth_plot':
				name = sample_name.replace('-','o')
			else:
				name = sample_name
		
			filename = fstr.format(name)

			transfer_file(pathjoin(src_path, sample_name, filename), pathjoin(dest_path, sample_name, filename))

	# Record the completion time of file transfer

	transfer_complete['timestamp_transfer_complete'] = datetime.now().isoformat()

	# Write the transfer completion details to a JSON file

	with open(os.path.join(dest_path, 'transfer_complete.json'), 'w') as f:
		json.dump(transfer_complete, f, indent=2)

# Function to transfer an entire run folder (not implemented yet)
def transfer_run_folder(src_path, dest_path):
	raise NotImplementedError()


