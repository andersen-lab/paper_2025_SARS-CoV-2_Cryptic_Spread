import sys
sys.path.append('../bjorn/src/')
import shutil
import argparse
import pandas as pd
import re
import bjorn_support as bs



if __name__ == "__main__":
    # grab user args
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",
                            type=str,
                            required=True,
                            help="Folder containing output analysis results from A-labs SARS-CoV-2 informatics pipeline")
    parser.add_argument("-o", "--out-fp",
                            type=str,
                            required=True,
                            help="Output filepath storing QC information for each sample")
    args = parser.parse_args()
    analysis_fp = args.input
    out_fp = args.out_fp
    depth_cols = ['SAMPLE_ID', 'COVERAGE', 'AVG_DEPTH', 'PATH_depth']
    mapped_cols = ['SAMPLE_ID', 'mapped', 'unmapped', 'PATH_mapped']
    # get depth information from samtools output
    depth_fps = bs.get_filepaths(analysis_fp, data_type='trimmed_bams/illumina/reports',
                                data_fmt='tsv',
                                generalised=False,
                                return_type='list')
    depth_df = pd.concat([pd.read_csv(fp, sep='\t') for fp in depth_fps])
    depth_df.fillna('NA', inplace=True)
    depth_df.rename(columns={'SAMPLE': 'PATH_depth'}, inplace=True)
    depth_df['SAMPLE_ID'] = depth_df['PATH_depth'].apply(lambda x: x.split('/')[-1].split('_')[0])
    depth_df = depth_df[depth_cols]
    # get mapped vs unmapped read count information
    mapped_cols = ['SAMPLE_ID', 'mapped', 'unmapped', 'PATH_mapped']
    mapped_fps = bs.get_filepaths(analysis_fp, data_type='merged_aligned_bams/illumina/reports',
                                data_fmt='tsv',
                                generalised=False,
                                return_type='list')
    mapped_df = pd.concat([pd.read_csv(fp, sep='\t') for fp in mapped_fps])
    mapped_df.fillna('NA', inplace=True)
    mapped_df.rename(columns={'SAMPLE': 'PATH_mapped'}, inplace=True)
    mapped_df['SAMPLE_ID'] = mapped_df['PATH_mapped'].apply(lambda x: x.split('/')[-1].split('_')[0])
    # fuse depth and mapped information to create temporary QC dataframe
    qc_df_tmp1 = pd.merge(depth_df, mapped_df, on='SAMPLE_ID')
    # get trim information from iVar's log files
    trim_fps = bs.get_filepaths(analysis_fp,data_type='logs/trimmed',
                             data_fmt='log',
                             generalised=False,
                             return_type='list')
    trim_df = pd.DataFrame(columns=['SAMPLE_ID', 'trimmed_pct', 'quality_pct', 
                                    'trimmed_count', 'quality_count', 'PATH_trim'])
    for fp in trim_fps:
        with open(fp, 'r') as fh:
            sample_data = {}
            sample_data['PATH_trim'] = fp
            sample_data['SAMPLE_ID'] = fp.split('/')[-1].split('_')[0]
            data = fh.readlines()
            try:
                trim_line = [l for l in data if 'Trimmed primers' in l][0]
                sample_data['trimmed_pct'] = re.findall('(\d+(?:\.\d+)?)', trim_line)[0]
                sample_data['trimmed_count'] = re.findall('\((\d+)\)', trim_line)[0]
                quality_line = [l for l in data if 'quality trimmed' in l][0]
                sample_data['quality_pct'] = re.findall('(\d+(?:\.\d+)?)', quality_line)[0]
                sample_data['quality_count'] = re.findall('\((\d+)\)', quality_line)[0]
            except:
                print(f"Not able to collect trim information from log file: {fp}")
            trim_df = trim_df.append(pd.Series(sample_data), ignore_index=True)
    # Fuse trim information to create final QC dataframe
    qc_df = pd.merge(qc_df_tmp1, trim_df, on='SAMPLE_ID')
    # save QC dataframe to file
    qc_df.to_csv(out_fp, index=False)