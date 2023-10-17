# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 17:37:24 2023

@author: rkerrid
"""
import seaborn as sns
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def silac_precursor_qc(df, path):
    counts_df = df['Run'].value_counts()
    
    # Set up the PDF
    create_reports_directory(path)
    output_dir = path + '/reports'
    pdf_path = os.path.join(output_dir, 'silac_precursors_report.pdf')
    print(os.path.exists(pdf_path))
    df_grouped = df.groupby('Run')
    
    with PdfPages(pdf_path) as pdf:
        # Title and introduction
        plt.figure(figsize=(8, 11))
        plt.axis('off')
        plt.text(0.5, 0.98, "Silac Precursors QC Report", ha='center', va='top', fontsize=15, fontweight='bold')
        plt.text(0.5, 0.85, "test", ha='center', va='center', wrap=True)
                
        # Display table of counts
        table_data = [["Run", "DF Count"]]
        for run in counts_df.keys():
            table_data.append([run, counts_df.get(run, 0)])
        
        plt.table(cellText=table_data, cellLoc = 'center', loc='center', colWidths=[0.2,0.2])
                
        pdf.savefig()  # Saves the current figure into the PDF
        plt.close()

        # Get the list of unique runs
        runs = df_grouped.groups.keys()

        # For each run, plot the histograms and save to PDF
        for run in runs:
            plot_histograms_for_run(run, df_grouped, ['H intensity', 'M intensity', 'L intensity'])
            
            # )
            pdf.savefig()  # Saves the current figure into the PDF
            plt.close()
        # for run in runs:
        #     plot_scatter_for_run(run, df_grouped)
        #     pdf.savefig()  # Saves the current figure into the PDF
        #     plt.close()
            
            
def plot_histograms_for_run(run, df_grouped, labels):
    """Plot histograms for a given run from multiple dataframes."""
    colors = sns.color_palette("husl", 3)  # Get a color palette

    # Extract the DataFrame for the specific run
    df_run = df_grouped.get_group(run)

    for label, color in zip(labels, colors):
        # Calculate log2 values
        data = np.log2(df_run[label].dropna())  # Add 1 to avoid log2(0)
        
        # Filter out -inf and inf values
        data = data[np.isfinite(data)]
        plt.hist(data, alpha=0.5, label=label, bins=300, color=color)

    plt.title(f'Histogram for {run}')
    plt.xlabel("log2 Intensity")
    plt.ylabel('Frequency')
    plt.legend()

# def plot_scatter_for_run(run, df_grouped):
#     """Plot scatter for a given run from the dataframe."""
#     # Extract the DataFrame for the specific run
#     df_run = df_grouped.get_group(run)

#     # Calculate the required values
#     y_values = np.log10(df_run['H intensity'] + df_run['L intensity'])
    
#     # Avoid division by zero by adding a small constant
#     ratio = (df_run['H intensity'] + 1e-10) / (df_run['L intensity'] + 1e-10)
#     x_values = np.log2(ratio)

#     # Plotting
#     plt.scatter(x_values, y_values, alpha=0.5)

#     plt.title(f'Scatter Plot for {run}')
#     plt.xlabel("log2(H intensity / L intensity)")
#     plt.ylabel("log10(H intensity + L intensity)")

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