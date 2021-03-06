# ptairMS

[![Travis build status](https://travis-ci.org/camilleroquencourt/ptairMS.svg?branch=master)](https://travis-ci.org/camilleroquencourt/ptairMS)


# Introduction
The _**ptairMS**_ package provides a workflow to process PTR-TOF-MS raw data in the open Hierarchical Data Format 5 ([HDF5](https://www.hdfgroup.org/); .h5 extension), and generate the peak table as an `ExpressionSet` object for subsequent data analysis with the many methods and packages available in [R](https://www.r-project.org/). Applications include the analysis of exhaled air, headspace or ambient air. The package offers several features to check the raw data and tune the few processing parameters. It also enables to include new samples in a study without re-processing all the previous data, providing a convenient management for cohort studies. 

# Volatolomics
Characterization of volatile organic compounds (VOCs) emitted by living organisms is of major interest in medicine, food sciences, and ecology. As an example, thousands of VOCs have been identified in the exhaled breath, resulting from normal metabolism or pathological processes (de lacy costello 2014). The main advantage of breath analysis in medicine is that the sampling is non-invasive (devillier, 2017). Methods based on mass spectrometry (MS) are the reference technologies for VOC analysis because of their sensitivity and large dynamic range.

# Installation 

`devtools::install_github("camilleroquencourt/ptairMS")`

# Acknowledgements
This research was developed within the Data Sciences for Deep Phenotyping and Precision Medicine team at CEA (C. Roquencourt, P. Zheng, P. Roger-Mele and E. Thévenot) in collaboration with the Exhalomics platform from the Foch Hospital (S. Grassin-Delyle, P. Devillier, H. Salvator, E. Naline and L.-J. Couderc) within the SoftwAiR project supported by the Agence Nationale de la Recherche (ANR-18-CE45-0017 grant).
