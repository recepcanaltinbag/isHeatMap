
from heat_functions import *
import sys
import argparse
import datetime
import os
import textwrap

version = '1.0.0'

the_red_code = "\033[91m"
the_normal = "\033[0m"

the_logo =  f'''
         {the_red_code}_{the_normal}       _    _            _     __  __              
        {the_red_code}(_){the_normal}     | |  | |          | |   |  \/  |            
        {the_red_code} _ ___{the_normal}  | |__| | ___  __ _| |_  | \  / | __ _ _ __  
        {the_red_code}| / __| {the_normal}|  __  |/ _ \/ _` | __| | |\/| |/ _` | '_ \ 
        {the_red_code}| \__ \ {the_normal}| |  | |  __/ (_| | |_  | |  | | (_| | |_) |
        {the_red_code}|_|___/ {the_normal}|_|  |_|\___|\__,_|\__| |_|  |_|\__,_| .__/ 
                                                     | |    
                                                     |_|    '''

if __name__ == "__main__":

    description = the_logo + ''' \nisHeatMap is created to find dominant insertion sequences in the genome and the relation between them. It uses the output of ISEScan program and creates a heatmap.
    '''

    parser = argparse.ArgumentParser(prog='isHeatMap', description= textwrap.dedent(description), formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--version', action='version', version='%(prog)s' + ' ' + version)
    parser.add_argument('--in_folder', required=True, default='', help='Folder that has input files, you can test the program with entering test_data folder')
    parser.add_argument('--threads', required=False, type=int, default=1, help='Number of CPU cores to use')
    parser.add_argument('--out_folder', required=True, default='output', help='Output Folder to save the final figure')

    args = parser.parse_args()

    input_folder = args.in_folder
    number_of_threads = args.threads
    output_folder = args.out_folder

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    fna_file_extension = "*is.fna"
    csv_file_extension = "*.csv"

    blast_folder = "blast_files"
    blast_input = os.path.join(blast_folder, "blast_in")
    blast_output = os.path.join(blast_folder, "blast_out")

    if not os.path.exists(blast_folder):
        os.makedirs(blast_folder)
    if not os.path.exists(blast_input):
        os.makedirs(blast_input)
    if not os.path.exists(blast_output):
        os.makedirs(blast_output)
    
    fna_files = glob.glob(os.path.join(input_folder, fna_file_extension)) 
    csv_files = glob.glob(os.path.join(input_folder, csv_file_extension))

    print(the_logo)
    print('\n\nisHeatMap starts at ', datetime.datetime.now().ctime())

    if len(fna_files) != len(csv_files):
        print('ERROR: The number of is fna files and csv files are not match, quitting from the program. Please control the input files in the input folder!')
    
    else:
        if len(fna_files) == 0:
            print('There is no input document, please control the folder or the folder name')
        for i in range(0, len(fna_files)):
            fna_file = fna_files[i]
            csv_file = glob.glob((fna_file).split('.')[0] + '*.csv')[0]

            print('Processing ', fna_file)

            records = list(SeqIO.parse(fna_file, "fasta"))
            record_dict = SeqIO.index(fna_file, "fasta") #indexing the fasta file
            q_vs_s_score_list_final = []
            list_ises = []
            is_dict = {}

            with open(csv_file, newline='') as csvfile:
                reader = csv.DictReader(csvfile)
                for row in reader:
                    seqIDinFastaFile = "_".join((row['seqID'],row['isBegin'], row['isEnd'], row['strand']))
                    #print(seqIDinFastaFile)
                    #print(record_dict[seqIDinFastaFile])
                    list_ises.append(seqIDinFastaFile)

                    row_seq_file = os.path.join(blast_input, seqIDinFastaFile)
                    if not os.path.exists(row_seq_file):
                        os.makedirs(row_seq_file)
                    
                    SeqIO.write(other_fasta_records(records, seqIDinFastaFile), os.path.join(row_seq_file,"subject.fasta"), "fasta")
                    SeqIO.write(record_dict[seqIDinFastaFile], os.path.join(row_seq_file,"query.fasta"), "fasta")
                    
                    q_vs_s_score_list_final.extend(blast_driver("./isheatmap/executables/makeblastdb", "./isheatmap/executables/blastn", row_seq_file, os.path.join(row_seq_file,"subject.fasta"),  os.path.join(row_seq_file,"query.fasta"), []))
                    is_dict[seqIDinFastaFile] = row['family']


            N = len(list_ises)
            distance_matrix = [ [0] * N for _ in range(N)]

            for ind1 in range(0, len(list_ises)):
                for ind2 in range(0, len(list_ises)):
                    if ind1 == ind2:
                        distance_matrix[ind1][ind2] = 1.0
                    else:
                        distance_matrix[ind1][ind2] = find_id_score_from_list_with_ids(q_vs_s_score_list_final, list_ises[ind1], list_ises[ind2])
            
            distance_df = pd.DataFrame(distance_matrix, index=list_ises, columns=list_ises)
            row_total_scores = distance_df.mean(axis=1)
            sorted_distance_df = distance_df.loc[row_total_scores.sort_values(ascending=True).index, row_total_scores.sort_values(ascending=True).index]
            sorted_distance_df = sorted_distance_df.rename(index=is_dict, columns=is_dict)

            # Create a mask to hide the lower triangle
            mask = np.triu(np.ones_like(sorted_distance_df, dtype=bool))

            # Create a custom color map transitioning from red to green
            cmap = sns.diverging_palette(220, 120, as_cmap=True)
            cmap2 = sns.color_palette("coolwarm", as_cmap=True)
            cmap3 = sns.diverging_palette(240, 10, s=75,l=60, as_cmap=True)
            sns.set(font_scale=font_scale_ratio(len(list_ises)))
            clustered = sns.clustermap(sorted_distance_df, metric="correlation", method="average", cmap=cmap3,  cbar_pos=(.1, .2, .03, .4),  figsize=(12, 12), col_cluster=True)
            col_order = clustered.dendrogram_col.reordered_ind


            df_reordered = sorted_distance_df.iloc[col_order, col_order]
            clustermap = sns.clustermap(df_reordered, cmap=cmap3, metric="correlation", cbar_pos=(.1, .2, .03, .4),  figsize=(12, 12), row_cluster=False, col_cluster=False)


            clustermap.fig.suptitle("ISes of " + str(os.path.basename(fna_file).split('.')[0]), fontsize=20)
            clustermap.ax_heatmap.tick_params(length=1, width=1)


            plt.savefig(output_folder + "/" + str(os.path.basename(fna_file).split('.')[0]) + ".pdf", dpi=360)
            print( '\n' + str(os.path.basename(fna_file).split('.')[0]) + ".pdf was saved in " + output_folder)
        
        print('END of the query, please check output folder')
        print('\n\nisHeatMap ends at ', datetime.datetime.now().ctime())

        print('Please cite if the program is helpful :)')