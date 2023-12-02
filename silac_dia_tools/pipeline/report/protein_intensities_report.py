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
    # create_reports_directory(path)
    manage_directories.create_directory(path,'reports')
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
        path = f'{path}protein intensities/'
        file_list = [f'{path}light_href.csv', f'{path}light_unnorm.csv', f'{path}nsp_href.csv', f'{path}nsp_unnorm.csv']        
        for file in file_list:
            
            # Add correlation heatmap page
            create_correlation_heatmap(file, pdf)


