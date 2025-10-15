import pandas as pd
from pathlib import Path
import base64
import os
import json
import yaml
from datetime import datetime
import logging
from os.path import join as pathjoin
import glob

def build_report_html(config,pipeline,run):
    # Get the current local date and time
    current_datetime = datetime.now()
    pipeline_short_name = pipeline['pipeline_name'].split('/')[1]
    pipeline_minor_version = ''.join(pipeline['pipeline_version'].rsplit('.', 1)[0])
    pipeline_path_name = '-'.join([pipeline_short_name, pipeline_minor_version, 'output'])

    sequencing_run_id = os.path.join(run['run_id'],pipeline_path_name)
    analysis_run_output_dir = os.path.join(config['analysis_output_dir'], sequencing_run_id)

    sample_dirs = [x for x in os.listdir(analysis_run_output_dir) if os.path.isdir(pathjoin(analysis_run_output_dir, x))]

    # Iterate through each sample folder 
    for sample_name in sample_dirs:
        #print("Current date and time:", current_datetime)
        # Define paths to your files
        demix_tsv = Path(os.path.join(analysis_run_output_dir,sample_name,"demix", sample_name+"_demixing_results.tsv"))
        consensus_tsv = Path(os.path.join(analysis_run_output_dir,sample_name,sample_name + "_consensus_seqs_report.tsv"))
        print(sample_name)
        pattern = os.path.join(analysis_run_output_dir,sample_name,sample_name + "_[0-9]*_provenance.yml")
        
        files = glob.glob(pattern)

        if not files:
            print("no yml files")
            provenance_yml = os.path.join(analysis_run_output_dir,sample_name,sample_name + "_[0-9]*_provenance.yml")
        else:
            provenance_yml = files[0]

        core_plot_png = Path(os.path.join(analysis_run_output_dir,sample_name,sample_name + "_core_db_depth_plots.png" ))
        ns5b_plot_png = Path(os.path.join(analysis_run_output_dir,sample_name,sample_name + "_ns5b_db_depth_plots.png" ))
        depth_plot_name = sample_name.replace('-','o')
        print(os.path.join(analysis_run_output_dir,sample_name,depth_plot_name + "_depth_plots.png" ))
        depth_plots = Path(os.path.join(analysis_run_output_dir,sample_name,depth_plot_name + "_depth_plots.png" ))
        core_tree = Path(os.path.join(analysis_run_output_dir,sample_name,sample_name + "_core_tree.png" ))
        ns5b_tree = Path(os.path.join(analysis_run_output_dir,sample_name,sample_name + "_ns5b_tree.png" ))
        core_subtype_tree = Path(os.path.join(analysis_run_output_dir,sample_name,sample_name + "_core_subtype_tree.png" ))
        ns5b_subtype_tree = Path(os.path.join(analysis_run_output_dir,sample_name,sample_name + "_ns5b_subtype_tree.png" ))
        blastn_result = Path(os.path.join(analysis_run_output_dir,sample_name,sample_name + "_blast_results_prefilter.csv" ))
        genotype_csv = Path(os.path.join(analysis_run_output_dir,sample_name,sample_name + "_genotype_calls_nt.csv" ))
        genome_result_csv = Path(os.path.join(analysis_run_output_dir,sample_name,sample_name + "_parsed_genome_results.csv"))
        #igv_report_path = os.path.join(analysis_run_output_dir,sample_name, sample_name + "_igv_report.html")
        #print(igv_report_path)
        # Helper to safely create image tags
        if not ( consensus_tsv.exists() or core_plot_png.exists() or blastn_result.exists() or depth_plots.exists() or genotype_csv.exists() or demix_tsv.exists() ):
            continue
         
        def img_tag_if_exists(img_path, height=None, alt="image"):
            if img_path.exists():
                data_uri = base64.b64encode(img_path.read_bytes()).decode('utf-8')
                height_attr = f' height="{height}px"' if height else ""
                return f'<img src="data:image/png;base64,{data_uri}"{height_attr} alt="{alt}">'
            return f'<span style="color:#888;">Missing: {img_path.name}</span>'

        # Helper to safely create html tables
        def table_if_exists(table_path, **kwargs):
            if table_path.exists():
                try:
                    df = pd.read_csv(table_path, **kwargs)
                    if table_path == genotype_csv:  
                        #df = pd.read_csv(genotype_csv, index_col=0)
                        df = df.drop(['subject_strand','e_value'],axis=1)
                        df['amplicon'] = df.apply(lambda row: row['query_seq_id'].split('|')[1], axis=1)
                        df = df.sort_values(['amplicon', 'bitscore'], ascending=[True, False])
                        df = df.groupby('amplicon').head(10).reset_index(drop=True)
                    
                    if table_path == blastn_result:
                        #df = pd.read_csv(blastn_result, index_col=0)
                        df = df.sort_values(['amplicon', 'bitscore'], ascending=[True, False])
                        df = df.groupby('amplicon').head(10).reset_index(drop=True)

                    return df.to_html(index=False, classes='data-table', border=0)
                except Exception as e:
                    return f"<div style='color:red'>Error reading {table_path.name}: {e}</div>"
            return f"<div style='color:#888'>Table file missing: {table_path.name}</div>"

        # Helper to safely read YAML
        if os.path.exists(provenance_yml):
            with open(provenance_yml, 'r') as f:
                provenance_data = yaml.safe_load(f)

            html_content = "<ul>\n"

            for item in provenance_data:
                for key, value in item.items():
                    html_content += f"  <li>{key}: {value}</li>\n"

            html_content += "</ul>\n"
        else:
            html_content = "<div style='color:#888'>Provenance file missing.</div>"

        # Usage in template
        consensus_html = table_if_exists(consensus_tsv, sep='\t')
        demix_html = table_if_exists(demix_tsv, sep='\t')
        blastn_html = table_if_exists(blastn_result, index_col=0)
        genotype_html = table_if_exists(genotype_csv, index_col=0)
        genome_results_html = table_if_exists(genome_result_csv, index_col = 0)
        img_tag_mapref_core = img_tag_if_exists(core_plot_png, alt="Reads mapped to core db")
        img_tag_mapref_ns5b = img_tag_if_exists(ns5b_plot_png, alt="Reads mapped to ns5b db")
        img_tag_depth = img_tag_if_exists(depth_plots, alt="Core/NS5B depth")
        img_tag_tree_core = img_tag_if_exists(core_tree, alt="Core genotype tree")
        img_tag_tree_ns5b = img_tag_if_exists(ns5b_tree, alt="NS5B genotype tree")
        img_tag_tree_core_subtype = img_tag_if_exists(core_subtype_tree, alt="Core subtype tree")
        img_tag_tree_ns5b_subtype = img_tag_if_exists(ns5b_subtype_tree, alt="NS5B subtype tree")

        #with open(provenance_yml,'r') as file:
        #    provenance_data = yaml.safe_load(file)

        #print(provenance_data)
        # Read tabular data into dataframes
        #demix_df = pd.read_csv(demix_tsv, sep='\t', nrows=20)         # reading first 20 rows as preview
        #consensus_df = pd.read_csv(consensus_tsv, sep='\t', nrows=20)

        #genotype_df = pd.read_csv(genotype_csv, index_col = 0)
        #genotype_df = genotype_df.drop(['subject_strand','e_value'],axis=1)
        #genotype_df['amplicon'] = genotype_df.apply(lambda row: row['query_seq_id'].split('|')[1], axis=1)
        #genotype_sorted = genotype_df.sort_values(['amplicon', 'bitscore'], ascending=[True, False])
        #genotype_sorted_top10 = genotype_sorted.groupby('amplicon').head(10).reset_index(drop=True)

        #blastn_df = pd.read_csv(blastn_result,index_col=0)
        #blastn_df = blastn_df.sort_values(['amplicon', 'bitscore'], ascending=[True, False])
        #blastn_df_top10 = blastn_df.groupby('amplicon').head(10).reset_index(drop=True)


        #html_content = "<ul>\n"
        #for item in provenance_data:
        #    for key,value in item.items():
        #        html_content += f"  <li>{key}: {value}<//li>\n"

        #html_content += "</ul>\n"

        # Start the HTML report string
        html = f"""
        <!DOCTYPE html>
        <html lang="en">
        <head>
        <meta charset="UTF-8" />
        <title>R3510056725 HCV Typing and Analysis Results</title>
        <meta name="viewport" content="width=device-width, initial-scale=1">

        <style>
            html, body {{
            font-family: 'Segoe UI', Arial, sans-serif;
            background: #f7fafc;
            color: #222;
            margin: 0;
            padding: 0;
            min-height: 100vh;
            }}
            body {{
            margin: 0;
            padding: 0 0 50px 0;
            }}
            .container {{
            max-width: 1100px;
            margin: 30px auto 60px auto;
            padding: 30px 30px 20px 30px;
            background: #fff;
            border-radius: 16px;
            box-shadow: 0 2px 20px rgba(30,40,80,0.13), 0 1.5px 6px rgba(0,0,0,0.04);
            }}
            h1 {{
            font-size: 2.2em;
            margin-bottom: 15px;
            }}
            h2 {{
            font-size: 1.35em;
            color: #29598c;
            border-bottom: 1px solid #e1e7ef;
            padding-bottom: 3px;
            margin-top: 40px;
            }}
            section {{
            margin-bottom: 36px;
            padding-bottom: 14px;
            border-bottom: 1px solid #eaeaea;
            }}
            .data-table {{
            width: 100%;
            border-collapse: separate;
            border-spacing: 0;
            background: #fafbfc;
            margin: 10px 0 20px 0;
            font-size: 0.98em;
            overflow-x: auto;
            display: block;
            max-width: 100%;
            box-shadow: 0 0.5px 1.2px rgba(80,130,170,0.05);
            }}
            .data-table th, .data-table td {{
            border: 1px solid #dde4ea;
            padding: 7px 8px;
            text-align: left;
            font-size: 0.97em;
            }}
            .data-table th {{
            background: #e4ecf6;
            font-weight: 500;
            color: #263e5a;
            }}
            .data-table tbody tr:nth-child(even) {{
            background: #f5faff;
            }}
            img {{
            width: auto;
            max-width: 97%;
            margin: 0 0 18px 0;
            border-radius: 8px;
            border: 1.2px solid #d1dbe6;
            box-shadow: 0 1px 8px rgba(80,120,140,0.10);
            display: block;
            }}
            .miniimg {{
            max-height: 120px;
            margin-right: 18px;
            display: inline-block;
            border-radius: 4px;
            box-shadow: none;
            border: 1px solid #eee;
            }}
            iframe {{
            border-radius: 10px;
            border: 1.2px solid #d5dbe7;
            margin-top: 10px;
            background: #f7fafc;
            }}
            ul {{
            padding-left: 18px;
            margin: 8px 0;
            }}
            li {{
            margin-bottom: 4px;
            font-size: 1em;
            }}
            .provenance-list {{
            background: #f4f8fb;
            border-radius: 8px;
            padding: 10px 18px 10px 18px;
            font-size: 0.98em;
            color: #2d2c2e;
            border: 1px solid #e2e6ec;
            margin-top: 7px;
            }}
            a {{
            color: #2077c7;
            text-decoration: none;
            border-bottom: 1px dashed #2077c7;
            }}
            a:hover {{
            color: #0a3157;
            border-bottom: 1px solid #0a3157;
            }}
            @media (max-width: 700px) {{
            .container {{ padding: 8px; border-radius: 0; box-shadow: none; }}
            h1 {{ font-size: 1.4em; }}
            section {{ margin-bottom: 20px; }}
            .data-table th, .data-table td {{ padding: 5px 5px; font-size: 0.91em; }}
            }}
        </style>
        </head>
        <body>
        <div class="container">
            <h1>{sample_name} HCV Typing &amp; Analysis Results</h1>
            <p>Report generated on: {current_datetime} </p>
            <section>
            <h2>Consensus Sequences Report</h2>
            <p>Note the segment coverage was calculated based on the full length hcv genome that has a length of an approximate 9600 nucleotides.
            But we are only sequencing core and ns5b amplicons, 3.53 (core) and 3.35 (ns5b) are full coverage. 
            </p>
            {consensus_html}
            </section>
            <section>
            <h2>Alignment Statistics</h2>
            {genome_results_html}
            </section>
            <section>
            <h2>Freyja Mixture Analysis</h2>
            {demix_html}
            </section>
            
            <section>
            <h2>Blast Results (HCV Reference Database)</h2>
            <div>
                <span style="font-size:1em;">References downloaded from: 
                <a href="https://ictv.global/sg_wiki/flaviviridae/hepacivirus/hcv_files" target="_blank">
                    ictv.global HCV files
                </a>
                </span>
            </div>
            {blastn_html}
            </section>
            
            <section>
            <h2>Blast Results (Core_nt databases, top 10 per amplicon)</h2>
            {genotype_html}
            </section>
            
            <section>
            <h2>Depth Plots</h2>
            <div style="margin-bottom:12px;">
            {img_tag_depth}
            </div>

            </section>
            
            <section>
            <h2>Phylogenetic Analysis</h2>
            <div style="display: flex; flex-wrap: wrap; gap: 24px 16px; margin-bottom: 13px;">
                <div>
                <div style="font-weight:600;color:#364e73;font-size:1.04em;margin-bottom:4px;">Core genotype clustering</div>
                {img_tag_tree_core}
                </div>
                <div>
                <div style="font-weight:600;color:#364e73;font-size:1.04em;margin-bottom:4px;">Core subtype clustering</div>
                {img_tag_tree_core_subtype}
                </div>
                <div>
                <div style="font-weight:600;color:#364e73;font-size:1.04em;margin-bottom:4px;">NS5B genotype clustering</div>
                {img_tag_tree_ns5b}
                </div>
                <div>
                <div style="font-weight:600;color:#364e73;font-size:1.04em;margin-bottom:4px;">NS5B subtype clustering</div>
                {img_tag_tree_ns5b_subtype}
                </div>
            </div>
            </section>
            
            <section>
            <h2>Reads Mapped to the Reference DB</h2>
            <p>Here we take the raw reads and map directly onto the 237 HCV references. The coverage of the reference with 
            the most reads mapped is plotted. Note the coverage plots here do not represent the real coverage as the current database
            used do not capture all the variability of the actual HCV samples. This is why the assembly mode is better.
            The results here give confirmation to the reported results above and in case of missing assembly, suggest what the genotype/subtypes
            might be. We expect these to have comparable coverage as the HCV database grows. 
                
            </p>
            <div style="display: flex; flex-wrap: wrap; gap: 28px 16px; align-items: center;">
                {img_tag_mapref_core}
                {img_tag_mapref_ns5b}
            </div>
            </section>
            
            <section>
            <h2>Provenance</h2>
            <div class="provenance-list">
                {html_content}
            </div>
            </section>
        </div>
        </body>
        </html>
        """

        # Write the HTML report to disk
        report_file =  os.path.join(analysis_run_output_dir,sample_name,sample_name + '_report.html')
        with open(report_file, 'w') as f:
            f.write(html)

        print(f"HTML report generated: {report_file}")