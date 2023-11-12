
import pandas
import seaborn as sns
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import warnings
from silac_dia_tools.pipeline.utils import manage_directories


def create_correlation_heatmap(df, pdf):
    # Assume first column is 'Proteins' and the rest are samples
    protein_col = df.columns.values.tolist()[0]
    sample_cols = df.columns.values.tolist()[1:]

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



def create_report(file_list, path, params):
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
        for file in file_list:
            print(file)
        # Add correlation heatmap page
        # create_correlation_heatmap(df, pdf)


