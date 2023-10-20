# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 15:08:12 2023

@author: rkerrid

        placeholder untill ready for protein intensities QC
"""

import seaborn as sns
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import warnings


def create_report(df, path):
    
    # Set up the PDF
    create_reports_directory(path)
    output_dir = path + '/reports'
    pdf_path = os.path.join(output_dir, 'filtering_report.pdf')

    with PdfPages(pdf_path) as pdf:
        # Title and introduction
        plt.figure(figsize=(8, 11))
        plt.axis('off')
        plt.text(0.5, 0.98, "Filtering QC Report", ha='center', va='top', fontsize=15, fontweight='bold')
       
        df_grouped = df.groupby("Run")
        # Get the list of unique runs
        runs = df_grouped.groups.keys()

        # For each run, plot the histograms and save to PDF
        # Suppress runtime warnings
        warnings.filterwarnings("ignore", category=RuntimeWarning) 
        for run in runs:
            plot_histograms_for_run(run, [df_grouped], ['df', 'contams', 'filtered_out'])
            # )
            pdf.savefig()  # Saves the current figure into the PDF
            plt.close()


def plot_histograms_for_run(run, dfs, labels, column='Precursor.Quantity'):
    """Plot histograms for a given run from multiple dataframes."""
    colors = sns.color_palette("husl", len(dfs))  # Get a color palette

    for df, label, color in zip(dfs, labels, colors):
        # Calculate log2 values
        data = np.log2(df.get_group(run)[column].dropna())
        # Filter out -inf and inf values
        data = data[np.isfinite(data)]
        plt.hist(data, alpha=0.5, label=label, bins=300, color=color)  # Removed edgecolor and added color parameter

    plt.title(f'Histogram for {run}')
    plt.xlabel(f"log2({column})")
    plt.ylabel('Frequency')
    plt.legend()

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
