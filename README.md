# labodock_binder
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/zenodo/10.5281/zenodo.10208365/?labpath=wedock.ipynb) 

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10208365.svg)](https://doi.org/10.5281/zenodo.10208365)

[![version](https://img.shields.io/badge/version-0.0.2-orange)](https://github.com/quantaosun/labodock_binder/releases)

[![license](https://img.shields.io/badge/license-MIT-yellow)](https://github.com/quantaosun/labodock_binder/blob/main/LICENSE)

# Alternatively binder launched directly from GitHub
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/quantaosun/labodock_binder/HEAD?labpath=wedock.ipynb)

## What happens once you click the launch binder badge ?

This Binder section is copied from https://discourse.jupyter.org/t/how-to-reduce-mybinder-org-repository-startup-time/4956 

Each time a user clicks a Binder link, these things happen:

A slot (called a “pod”) is reserved on one of the cloud machines. This takes 1-2 seconds.
Binder looks to see if a Docker image exists for that repository
If it doesn’t, Binder must first build the image for that repo using repo2docker (this takes time)
Binder looks for a built image on the machine the user will use
If it isn’t on the machine, Binder must first pull the image onto that machine (this takes time)
Binder launches the user’s session. This includes:
a small amount of time to start the “init pods to limit network access”,
a few seconds for the Jupyter process to start,
a few seconds for BinderHub to notice,
and finally, your browser needs to follow the redirect.

## Introduction
labodock_binder is an enhanced iteration of Labodock (https://github.com/RyanZR/labodock), with a primary focus on improved accessibility and user-friendliness. Notable modifications from Labodock version 2.0.0 include:

**Seamless Access**: Unlike the original Google Colab login requirement, labodock_binder allows free access from Chrome, Firefox, and Safari, without the need for signup or login. This is made possible through the implementation of Binder technology (https://mybinder.org/).

**Pre-built Images**: The labodock_binder image is pre-built. When a user clicks the link, the image is pulled from the cloud to the local browser, instead of installing all packages from scratch. This significantly accelerates the readiness for docking.

**Simplified Execution**: labodock_binder streamlines the process with one-click execution, eliminating the need for repetitive variable modifications.

These enhancements collectively contribute to a more user-friendly and accessible experience. labodock_binder emerges as a versatile and convenient tool for users worldwide, requiring only internet availability—no Google account necessary.

<img width="775" alt="image" src="https://github.com/quantaosun/labodock_binder/assets/75652473/d0333255-36f2-46da-a105-5290bd91a063">


## Usage
**PDB_ID** is a four-letter code from the PDB bank website, for example, 3HTB

**Native_lig** below means the native ligand of the co-crystal structure of the PDB bank, it can be found on the PDB website, it should always be a three-letter code, for example, JZ4

**smiles** is the small molecule you want to dock, you can obtain this from ChemDraw, from **Edit** > **Copy as smiles**

(In wedock2.ipynb, **smiles** variable was replaced by a CSV file containing multiple smiles)

The first three cells are the only thing you should modify. From the top left, click **Run** > **Run All Cells**, and the docking will start automatically. Your docked 3D interaction diagram will be at the bottom of the notebook.

**wedock.ipynb does not need any manual intervention, but for wedock2.ipynb, if you wish to overlay two docked small molecules, you do need to specify several variables manually.**

## Acknowledgments

- Credits to the [Labodock](https://github.com/RyanZR/labodock) project for developing an excellent and open-source suite for docking on Google Colab.
- Credits to the [Binder](https://mybinder.org/) project for developing the fantastic interactive Jupiter notebook used on web browsers.
