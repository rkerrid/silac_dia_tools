<<<<<<< HEAD
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 13:33:49 2023

@author: rkerrid
"""
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import os
import json 


def create_intensities_report(df, path):

    # Construct the description string
    # Load filtering parameters from JSON
    CONFIG_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'configs')
    json_path = os.path.join(CONFIG_DIR, 'filtering_parameters.json')
    with open(json_path, 'r') as f:
        params = json.load(f)
    params_str = "\n".join([f"{key} {item['op']} {item['value']}" for key, item in params['apply_filters'].items()])
    description = f"Parameters used:\n{params_str}"
    
    # Set up the PDF
    create_reports_directory(path)
    output_dir = path + '/reports'
=======
from icecream import ic
import pandas as pd
import seaborn as sns
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import warnings
from silac_dia_tools.pipeline.utils import manage_directories


# def create_correlation_heatmap(df, pdf):
#     # Assume first column is 'Proteins' and the rest are samples
#     protein_col = df.columns.values.tolist()[0]
#     sample_cols = df.columns.values.tolist()[1:]

#     # Calculate the correlation matrix
#     corr = df[sample_cols].corr()

#     # Create a mask for the upper triangle
#     mask = np.triu(np.ones_like(corr, dtype=bool))

#     # Set up the matplotlib figure
#     f, ax = plt.subplots(figsize=(11, 9))

#     # Generate a custom diverging colormap
#     cmap = sns.diverging_palette(230, 20, as_cmap=True)

#     # Draw the heatmap with the mask and correct aspect ratio
#     sns.heatmap(corr, mask=mask, cmap=cmap, vmax=.3, center=0,
#                 square=True, linewidths=.5, cbar_kws={"shrink": .5})

#     plt.title('Sample Correlation Heatmap')

#     # Save the heatmap to the PDF
#     pdf.savefig(f)
#     plt.close(f)

def create_correlation_heatmap(file_path, pdf):

    # Load the dataset
    df = pd.read_csv(file_path, index_col=0)  # Assuming the first column is 'Protein Group'

    # Calculate the correlation matrix
    corr_matrix = df.corr(method='pearson')

    # Plot the heatmap
    f, ax = plt.subplots(figsize=(11, 9))
    sns.heatmap(corr_matrix, annot=True, cmap='coolwarm')
  
    # Add 1 to last_slash_index to start after the "/", and csv_index is the end point
    title = file_path[file_path.rfind('/') + 1 : file_path.find('.csv')]
     
    plt.title(f"Heatmap of Sample Correlations in {title}")
    # Save the heatmap to the PDF
    pdf.savefig(f)
    plt.close(f)


def create_report(path, params):
    # Construct the description string
    params_str = "\n".join([f"{key} {item['op']} {item['value']}" for key, item in params['apply_filters'].items()])
    description = f"Parameters used:\n{params_str}"

    # Set up the PDF
    manage_directories.create_directory(path, 'reports')
    output_dir = os.path.join(path, 'reports')
>>>>>>> silac_dia_tools_oop
    pdf_path = os.path.join(output_dir, 'protein_intensities_report.pdf')

    with PdfPages(pdf_path) as pdf:
        # Title and introduction
        plt.figure(figsize=(11, 8))
        plt.axis('off')
<<<<<<< HEAD
        plt.text(0.5, 0.98, "Protein Intensities QC Report", ha='center', va='top', fontsize=15, fontweight='bold')
        plt.text(0.5, 0.85, description, ha='center', va='center', wrap=True)
        
        pdf.savefig()  # Saves the current figure into the PDF
        plt.close()
        
        # # Creating the barplot
        # plt.figure(figsize=(10, 7))
        # barplot = sns.barplot(x='Run', y='Counts', data=counts_df, orient='v')
        # plt.title('Counts per Run')
        # plt.xlabel('Run')
        # plt.ylabel('Counts')
        
        # # Rotating x-axis labels
        # barplot.set_xticklabels(barplot.get_xticklabels(), rotation=45, horizontalalignment='right')
        
        # # Adding counts on the bars
        # for index, row in counts_df.iterrows():
        #     barplot.text(row.name, row.Counts, row.Counts, color='black', ha="center")
        
        # pdf.savefig()  # Saves the current figure into the PDF
        # plt.close()

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
=======
        plt.text(0.5, 0.98, "Intensities QC Report", ha='center', va='top', fontsize=15, fontweight='bold')
        plt.text(0.5, 0.85, description, ha='center', va='center', wrap=True)

        pdf.savefig()  # Saves the current figure into the PDF
        plt.close()

        path = f'{path}protein intensities/'
        file_list = [f'{path}light_href.csv', f'{path}light_unnorm.csv', f'{path}nsp_href.csv', f'{path}nsp_unnorm.csv']        
        for file in file_list:
            
            # Add correlation heatmap page
            create_correlation_heatmap(file, pdf)





# import pandas
# import seaborn as sns
# import os
# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.backends.backend_pdf import PdfPages
# import warnings
# from silac_dia_tools.pipeline.utils import manage_directories


# def create_correlation_heatmap(df, pdf):
#     # Assume first column is 'Proteins' and the rest are samples
#     protein_col = df.columns.values.tolist()[0]
#     sample_cols = df.columns.values.tolist()[1:]

#     # Calculate the correlation matrix
#     corr = df[sample_cols].corr()

#     # Create a mask for the upper triangle
#     mask = np.triu(np.ones_like(corr, dtype=bool))

#     # Set up the matplotlib figure
#     f, ax = plt.subplots(figsize=(11, 9))

#     # Generate a custom diverging colormap
#     cmap = sns.diverging_palette(230, 20, as_cmap=True)

#     # Draw the heatmap with the mask and correct aspect ratio
#     sns.heatmap(corr, mask=mask, cmap=cmap, vmax=.3, center=0,
#                 square=True, linewidths=.5, cbar_kws={"shrink": .5})

#     plt.title('Sample Correlation Heatmap')

#     # Save the heatmap to the PDF
#     pdf.savefig(f)
#     plt.close(f)



# def create_report(file_list, path, params):
#     # Construct the description string
#     params_str = "\n".join([f"{key} {item['op']} {item['value']}" for key, item in params['apply_filters'].items()])
#     description = f"Parameters used:\n{params_str}"

#     # Set up the PDF
#     # create_reports_directory(path)
#     manage_directories.create_directory(path,'reports')
#     output_dir = path + '/reports'
#     pdf_path = os.path.join(output_dir, 'protein_intensities_report.pdf')

#     with PdfPages(pdf_path) as pdf:
#         # Title and introduction
#         plt.figure(figsize=(8, 11))
#         plt.axis('off')
#         plt.text(0.5, 0.98, "Intensities QC Report", ha='center', va='top', fontsize=15, fontweight='bold')
#         plt.text(0.5, 0.85, description, ha='center', va='center', wrap=True)

#         pdf.savefig()  # Saves the current figure into the PDF
#         plt.close()
#         for file in file_list:
#             print(file)
#         # Add correlation heatmap page
#         # create_correlation_heatmap(df, pdf)



>>>>>>> silac_dia_tools_oop
