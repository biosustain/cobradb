# -*- coding: utf-8 -*-

from os.path import abspath, dirname, isfile
from sys import path
from setuptools import setup, find_packages

# To temporarily modify sys.path
SETUP_DIR = abspath(dirname(__file__))

requirement_path = f"{SETUP_DIR}/requirements.txt"
install_requires = []
if isfile(requirement_path):
    with open(requirement_path) as f:
        install_requires = f.read().splitlines()

setup(
    name="cobradb",
    version="0.3.0",
    description="""COBRAdb loads genome-scale metabolic models and genome
                   annotations into a relational database.""",
    url="https://github.com/SBRG/cobradb",
    author="Zachary King",
    author_email="zaking@ucsd.edu",
    license="Apache-2.0",
    classifiers=[
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
    ],
    keywords="systems biology, genome-scale model",
    packages=find_packages(),
    package_data={
        "cobradb": [
            "data/escher_templates/*",
        ]
    },
    install_requires=install_requires,
)
