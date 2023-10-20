# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 13:33:49 2023

@author: rkerrid
"""
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import os

def create_report(df, path):
    # Count rows for each dataframe
    counts_df = df['Run'].value_counts().reset_index()
    counts_df.columns = ['Run', 'Counts']  # Naming the columns
    
    # Set up the PDF
    create_reports_directory(path)
    output_dir = path + '/reports'
    pdf_path = os.path.join(output_dir, 'protein_groups_report.pdf')

    with PdfPages(pdf_path) as pdf:
        # Title and introduction
        plt.figure(figsize=(11, 8))
        plt.axis('off')
        plt.text(0.5, 0.98, "Protein Ratios QC Report", ha='center', va='top', fontsize=15, fontweight='bold')
        
        # Creating the barplot
        plt.figure(figsize=(10, 7))
        barplot = sns.barplot(x='Run', y='Counts', data=counts_df, orient='v')
        plt.title('Counts per Run')
        plt.xlabel('Run')
        plt.ylabel('Counts')
        
        # Rotating x-axis labels
        barplot.set_xticklabels(barplot.get_xticklabels(), rotation=45, horizontalalignment='right')
        
        # Adding counts on the bars
        for index, row in counts_df.iterrows():
            barplot.text(row.name, row.Counts, row.Counts, color='black', ha="center")
        
        pdf.savefig()  # Saves the current figure into the PDF
        plt.close()

#create preprocessing directory for new files 
def create_reports_directory(path):
    # Combine the paths
    new_folder_path = os.path.join(path, 'reports')
    
    # Create the new folder
    if not os.path.exists(new_folder_path):
        os.makedirs(new_folder_path)
        print(f"Folder reports created successfully at {new_folder_path}")
    else:
        print(f"Folder reports already exists at {new_folder_path}")