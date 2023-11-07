# # -*- coding: utf-8 -*-
# """
# Created on Thu Oct 19 15:08:12 2023

# @author: rkerrid

#         placeholder untill ready for protein intensities QC
# """

# import seaborn as sns
# import os
# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.backends.backend_pdf import PdfPages
# import warnings


# def create_report(df, path, params):
#     # Construct the description string
#     params_str = "\n".join([f"{key} {item['op']} {item['value']}" for key, item in params['apply_filters'].items()])
#     description = f"Parameters used:\n{params_str}"
    
#     # Set up the PDF
#     create_reports_directory(path)
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
        
# #create preprocessing directory for new files 
# def create_reports_directory(path):
#     # Combine the paths
#     new_folder_path = os.path.join(path, 'reports')
    
#     # Create the new folder
#     if not os.path.exists(new_folder_path):
#         os.makedirs(new_folder_path)
#         print(f"Folder reports created successfully at {new_folder_path}")
#     else:
#         print(f"Folder reports already exists at {new_folder_path}")


import seaborn as sns
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import warnings

def create_correlation_heatmap(df, pdf):
    # Assume first column is 'Proteins' and the rest are samples
    protein_col = df.columns[0]
    sample_cols = df.columns[1:]

    # Calculate the correlation matrix
    corr = df[sample_cols].corr()

    # Create a mask for the upper triangle
    mask = np.triu(np.ones_like(corr, dtype=bool))

    # Set up the matplotlib figure
    f, ax = plt.subplots(figsize=(11, 9))

    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(230, 20, as_cmap=True)

    # Draw the heatmap with the mask and correct aspect ratio
    sns.heatmap(corr, mask=mask, cmap=cmap, vmax=.3, center=0,
                square=True, linewidths=.5, cbar_kws={"shrink": .5})

    plt.title('Sample Correlation Heatmap')

    # Save the heatmap to the PDF
    pdf.savefig(f)
    plt.close(f)

def create_reports_directory(path):
    # Combine the paths
    new_folder_path = os.path.join(path, 'reports')

    # Create the new folder
    if not os.path.exists(new_folder_path):
        os.makedirs(new_folder_path)
        print(f"Folder 'reports' created successfully at {new_folder_path}")
    else:
        print(f"Folder 'reports' already exists at {new_folder_path}")

def create_report(df, path, params):
    # Construct the description string
    params_str = "\n".join([f"{key} {item['op']} {item['value']}" for key, item in params['apply_filters'].items()])
    description = f"Parameters used:\n{params_str}"

    # Set up the PDF
    create_reports_directory(path)
    output_dir = path + '/reports'
    pdf_path = os.path.join(output_dir, 'protein_intensities_report.pdf')

    with PdfPages(pdf_path) as pdf:
        # Title and introduction
        plt.figure(figsize=(8, 11))
        plt.axis('off')
        plt.text(0.5, 0.98, "Intensities QC Report", ha='center', va='top', fontsize=15, fontweight='bold')
        plt.text(0.5, 0.85, description, ha='center', va='center', wrap=True)

        pdf.savefig()  # Saves the current figure into the PDF
        plt.close()

        # Add correlation heatmap page
        create_correlation_heatmap(df, pdf)

# Now you can call create_report with the dataframe, path, and parameters
# Example usage:
# create_report(your_dataframe, 'path_to_directory', your_params_dictionary)
