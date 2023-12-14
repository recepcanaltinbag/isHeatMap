#Recep Can Altınbağ, 2023

import os
import glob
from Bio import SeqIO
import csv
import subprocess
import shutil
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import ticker
from scipy.cluster.hierarchy import linkage, dendrogram


#isHeatMap programm will use the isescan in the future, this function is not usefeul for now
def isescan_driver(path_of_isescan, genome_file, out_folder, n_thread=2):
    #isescan.py --seqfile saponiphilaDSM9751.fna --output resultSap --nthread 12
    args = f"{path_of_isescan} --seqfile {genome_file} --output {out_folder} --nthread {n_thread} "
    my_process = subprocess.run(args, shell=True, executable='/bin/bash', text=True, check=True,  stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return 0


#Blast related functions, 
def delete_blast_db(db_dir):
    try:
        shutil.rmtree(db_dir)
    except OSError as e:
        print("ERROR: %s - %s." % (e.filename, e.strerror))


def make_blast_db(path_of_makeblastdb, db_input, db_output):
    args = f"{path_of_makeblastdb} -in  {db_input} -dbtype nucl -out {db_output} "
    my_process = subprocess.run(args, shell=True, executable='/bin/bash', text=True, check=True,  stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def run_blast(path_of_blastn, query, database, output):
    extras = '"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen"'
    args = f'{path_of_blastn} -db {database} -query {query} -out {output} -outfmt {extras}'
    my_process = subprocess.run(args, shell=True, executable='/bin/bash', text=True, check=True,  stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def blast_driver(path_of_makeblastdb, path_of_blastn, out_blast_folder, db_path, blast_query, q_vs_s_score_list):
    cycle_file_name = os.path.basename(db_path)
    db_output = os.path.join(out_blast_folder,"blast_temp", cycle_file_name)
    db_dir = out_blast_folder + "/blast_temp"
    db_input = db_path
    result_path = os.path.join(out_blast_folder,"blast_files", os.path.basename(blast_query) + ".out")

    if not os.path.exists(db_dir):
        os.mkdir(db_dir)
    if not os.path.exists(out_blast_folder + "/blast_files"):
        os.mkdir(out_blast_folder + "/blast_files")
    if not os.path.exists(blast_query):
        print('WARNING: No available Blast Query file')
        return False

    try:
        if os.path.exists(result_path):
            print(f'WARNING: File existed{result_path}, skipping the blast run!')
        else:
            make_blast_db(path_of_makeblastdb, db_input, db_output)
            run_blast(path_of_blastn, blast_query, db_output, result_path)
    except Exception as e:
        print('ERROR: Blast Error')
        raise UserWarning('Blast Error')

    q_vs_s_score_list.extend(parsing_blast_file(result_path))
    delete_blast_db(db_dir)
    return q_vs_s_score_list


#According to blast file, distance matrix identities of insertion sequences are determined. There can be partial alignments, 
#So, the score can be different from a->b and b->a in the scores.  
#the function returns the list: (query is id (string), subject is id (string), score (float), partial info (Boolean))
def parsing_blast_file(blast_result_file, name_of_query=''):
    list_of_cds = []

    try:
        blast_result = pd.read_table(blast_result_file, header=None)
    except pd.errors.EmptyDataError:
        print('WARNING: Empty blast file')
        return list_of_cds

    default_outfmt6_cols = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen'.strip().split(' ')
    blast_result.columns = default_outfmt6_cols
    df_filtered = blast_result[(blast_result['pident'] >= 50.0) & (blast_result['evalue'] < 0.1)]
    q_vs_s_score_list = []

    if len(df_filtered) > 0:
        if name_of_query == '':
            name_of_query = (df_filtered.iloc[0]['sseqid'])
        
        df_filtered.sort_values(by=['qstart'], ascending=False)
        groups = df_filtered.groupby('sseqid')
        #print(groups['sseqid'])

        for qseqid, frame in groups:
            #print(frame.iloc[0]['sstart'])
            #print(frame.iloc[0]['qseqid'], frame.iloc[0]['sseqid'], frame.iloc[0]['pident'], frame.iloc[0]['qstart'], frame.iloc[0]['qend'])

            the_best_list = []
            temp_b = 0
            the_identity_score = 0
            the_cov_score = 0
            partial = False
            if len(frame) > 1:
                partial = True
            
            for trv_in_group in range(0,len(frame)):
                the_identity_score =  frame.iloc[trv_in_group]['pident']/100 + the_identity_score
                the_a = frame.iloc[trv_in_group]['qstart']
                the_b = frame.iloc[trv_in_group]['qend']
                the_cov_score = the_b - the_a

                if temp_b != 0 and temp_b > the_a:
                    the_cov_score = the_cov_score - abs(temp_b - the_a)       
                temp_b = the_b

            the_final_score = round((the_identity_score/len(frame)) * (the_cov_score/frame.iloc[0]['qlen']),2)
            if the_final_score < 0:
                the_final_score = the_final_score * -1.0
            #print(the_final_score, 'for', frame.iloc[0]['qseqid'], ' and ', frame.iloc[0]['sseqid'])
            q_vs_s_score_list.append((frame.iloc[0]['qseqid'], frame.iloc[0]['sseqid'], the_final_score, partial))

    return q_vs_s_score_list


#list of records other than the given the_record_id
def other_fasta_records(records, the_record_id):
    out_records = []
    for record in records:
        if record.id != the_record_id:
            out_records.append(record)
    return out_records


#finds the scores according to the id1 and id2, returns the score
def find_id_score_from_list_with_ids(the_list, the_id1, the_id2):
    for element in the_list:
        id1 = element[0]
        id2 = element[1]
        id_score = element[2]

        if the_id1 == id1 and the_id2 == id2:
            if element[3] == True:
                return id_score
            return id_score

    for element in the_list:
        id1 = element[0]
        id2 = element[1]
        id_score = element[2]

        if the_id1 == id2 and the_id2 == id1:
            if element[3] == True:
                return id_score
            return id_score

    return 0.0

def font_scale_ratio(number_of_elements):
    font_scale = number_of_elements * (-0.0012) + 0.55
    return font_scale
