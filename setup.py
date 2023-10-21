from setuptools import setup, find_packages

setup(
    name='silac_dia_tools',
    version='0.1',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'pandas>=1.0.0',  # specify the minimum version you require
        'matplotlib>=3.0.0',
        'numpy>=1.18.0',
        'spyder',
        'icecream',
        'seaborn',
        'fpdf',
        'alphapept @ git+https://github.com/MannLabs/alphapept.git@master#egg=alphapept',  # specify the branch, tag, or commit
        'directlfq @ git+https://github.com/MannLabs/directlfq.git@master#egg=directlfq',  # specify the branch, tag, or commit
    ],
)
