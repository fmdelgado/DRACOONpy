from setuptools import setup, find_packages
with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='dracoon',
    version='0.1.0',
    url='https://github.com/fmdelgado/DRACOONpy',
    author='Fernando M. Delgado Chaves, PhD',
    author_email='fernando.miguel.delgado-chaves@uni-hamburg.de',
    description='DRaCOoN (Differential Regulation and CO-expression Networks) is a data-driven tool optimized for effectively retrieving differential relationships between genes across two distinct conditions. ',
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=['app'],  # automatically find all packages
    include_package_data=True,
    classifiers=[
        'Programming Language :: Python :: 3.10',
        'Intended Audience :: Science/Research',
        # "License :: OSI Approved :: MIT License",
        # "Operating System :: OS Independent",
    ],                                      # Information to filter the project on PyPi website
    python_requires='>=3.10',
    install_requires=[
        "matplotlib",
        "pyvis",
        "pandas",
        "networkx",
        "seaborn",
        "plotly",
        "statsmodels>=0.14.0",
        "numpy<=1.24.3",
        "fitter",
        "scipy",
        "tqdm",
        "streamlit-chat",
        "streamlit-pills",
        "numba"
    ]
)
