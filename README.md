[![Build Status](https://travis-ci.org/smerz1989/copoly-analysis.svg?branch=master)](https://travis-ci.org/smerz1989/copoly-analysis)

# Copoly

## Description

A Python module for analyzing and deploying LAMMPS simulations of polymerization reaction.  This includes methods for encapsulating simulation timesteps as individual objects, `SimulationSnapshots`, which allows quick analysis of the graph structure of the polymers using the efficient iGraph module.  Deployment of LAMMPS simulations can be done over SSH connections to remote HPC servers using SLURM queuing systems.

## Installation

### Quickstart

Copoly can be quickly installed on systems containing a working installation of Python 3 using pip running the following command from the copoly-analysis directory:

`pip install .`
