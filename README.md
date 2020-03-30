# Copoly

[![Build Status](https://travis-ci.org/smerz1989/copoly-analysis.svg?branch=master)](https://travis-ci.org/smerz1989/copoly-analysis)

## Description

A Python module for analyzing and deploying LAMMPS simulations of polymerization reaction.  This includes methods for encapsulating simulation timesteps as individual objects, `SimulationSnapshots`, which allows quick analysis of the graph structure of the polymers using the efficient iGraph module.  Deployment of LAMMPS simulations can be done over SSH connections to remote HPC servers using SLURM queuing systems.

## Getting Started

### Prerequisites

Copoly depends on the following Python packages which are installed when installing using pip:

* **matplotlib** [(homepage)](https://matplotlib.org/)
* **seaborn** [(homepage)](https://seaborn.pydata.org/)
* **tqdm** [(homepage)](https://tqdm.github.io/)
* **networkx** [(homepage)](https://networkx.github.io/)
* **numpy** [(homepage)](https://numpy.org/)
* **paramiko** [(homepage)](https://www.paramiko.org/)
* **pandas** [(homepage)](https://pandas.pydata.org/)
* **python_igraph** [(homepage)](https://igraph.org/python/)
* **moltemplate** [(homepage)](https://www.moltemplate.org/)

Additionally, simulations are run using [LAMMPS](https://lammps.sandia.gov/) as such in order to deploy simulations using Copoly LAMMPS must be installed on the target environment.  This can be locally or on a remote cluster depending on how one where one wishes to deploy the simulations.

### Quickstart

Copoly can be quickly installed on systems containing a working installation of Python 3 using pip running the following command from the copoly-analysis directory:

`pip install .`

## Running Tests

Unit tests are kept in `copoly/tests` and can be run individually or as a group with the command `pytest` executed within the folder `copoly/tests`
