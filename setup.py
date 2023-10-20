from setuptools import setup, find_packages

setup(
    name='silac_dia_tools',
    version='0.1',
    packages=find_packages(where="silac_dia_tools"),
    include_package_data=True,
    package_dir={"": "silac_dia_tools"},
    
    install_requires=[
        'pandas>=1.0.0',  # specify the minimum version you require
        'matplotlib>=3.0.0', 
        'numpy>=1.18.0', 
        'spyder',  # this will install the latest version of spyder
        'icecream',
        'seaborn',
        'fpdf',# adding numpy as an example, you can add more packages as needed
        # add more packages here as needed
    ],
    dependency_links=[
        'git+https://github.com/MannLabs/alphapept.git',
        'git+https://github.com/MannLabs/directlfq.git',
    ],
)
