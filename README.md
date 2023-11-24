# labodock_binder
## Introduction
labodock_binder is an enhanced iteration of Labodock (https://github.com/RyanZR/labodock), with a primary focus on improved accessibility and user-friendliness. Notable modifications from Labodock version 2.0.0 include:

**Seamless Access**: Unlike the original Google Colab login requirement, Wedock allows free access from Chrome, Firefox, and Safari, without the need for signup or login. This is made possible through the implementation of Binder technology (https://mybinder.org/).

**Pre-built Images**: The Wedock image is pre-built. When a user clicks the link, the image is pulled from the cloud to the local browser, instead of installing all packages from scratch. This significantly accelerates the readiness for docking.

**Simplified Execution**: Wedock streamlines the process with one-click execution, eliminating the need for repetitive variable modifications.

These enhancements collectively contribute to a more user-friendly and accessible experience. Wedock emerges as a versatile and convenient tool for users worldwide, requiring only internet availability—no Google account necessary.

## Usage
**PDB_ID** is a four-letter code from the PDB bank website, for example, 3HTB

**Native_lig** below means the native ligand of the co-crystal structure of the PDB bank, it can be found on the PDB website, it should always be a three-letter code, for example, JZ4

**smiles** is the small molecule you want to dock, you can obtain this from ChemDraw, from Edit > Copy as smiles

The first three cells are the only thing you should modify. From the top left, click Run > Run All Cells, and the docking will start automatically. Your docked 3D interaction diagram will be at the bottom of the notebook.

**webdock.ipynb do not need any manual intervention, but for webdock2.ipynb, if you wish to overlay two docked small molecules, you do need to specify several variables manually.**
