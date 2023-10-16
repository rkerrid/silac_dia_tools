from setuptools import setup, find_packages

setup(
    name='silac_dia_tools',
    version='0.1',
    packages=find_packages(where="src"),
    package_dir={"": "src"},
)