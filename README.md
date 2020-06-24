# Readme for GLM-BCI Repository
## Purpose and Publication
This README provides a brief documentation for the code and data provided in this repository that is linked with the following [publication](https://www.frontiersin.org/articles/10.3389/fnhum.2020.00030):
> A. von Lühmann, A. Ortega-Martinez, D.A. Boas and M.A. Yücel, 2020, *Using the General Linear Model to Improve Performance in fNIRS Single Trial Analysis and Classification: A Perspective*, Frontiers in Human Neuroscience, vol. 14, no. 30, doi: https://doi.org/10.3389/fnhum.2020.00030  

Please cite this work when using any of the presented methods or code.

## Directory Structure
**.../aux functions**

- auxilliary functions 

**.../Homer functions**

- processing functions from the homer2 toolbox that our scripts are dependent on. The full toolbox can be downloaded [here](https://www.nitrc.org/projects/homer2)

**.../lit rev data**

- the data and function for generating the pie-charts in the paper

**.../main**

- contains the analysis scripts and functions to generate the main results and corresponding figures from the manuscript. 

**.../sim HRF**

- contains synthetic HRF templates and function to augment fNIRS resting data to generate the ground truth data used for the evaluation.




## How To Reproduce Figures
The following list points to the analysis scripts that generate the results and figures of the publication:
- Figure 1:  .../lit rev data/create_piecharts.m
- Figure 2:  .../aux functions/generate_Fig2_signals.m 
- Figure 3:  .../lit rev data/create_piecharts.m
- Figure 7a: .../main/plot_Features_SSvsNo.m
- Figure 7b: .../aux functions/generate_Fig7_right_panel.m
- Figure 8:  .../lit rev data/create_piecharts.m
- Figure 9:  .../main/plot_Features_SSvsNo.m
- Figure 10:  .../lit rev data/create_piecharts.m
- Figure 11: .../main/plot_Features_SSvsNo.m
- Figure 12: .../main/classify_SSvsNo.m   - please note that this code depends on the BBCI toolbox that you can find [here](https://github.com/bbci/bbci_public)

