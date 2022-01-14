import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="atlas_anndata",
    version="0.0.1",
    author="Jonathan Manning",
    author_email="jmanning@ebi.ac.uk",
    description="Functions for preparing annData files for experiment inclusion in Single-cell Expression Atlas"
    license="Apache Software License",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ebi-gene-expression-group/atlas-anndata",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    include_package_data=True,
    install_requires=[
        "scanpy",
        "yaml",
        "click",
        "jsonschema"
    ],
    python_requires=">=3.6",
)
