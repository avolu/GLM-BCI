# Readme for GLM-BCI Repository
## Purpose and Publication
This README provides a brief documentation for the code and data provided in this repository that is linked with the following [publication](https://www.frontiersin.org/articles/10.3389/fnhum.2020.00030):
> A. von Lühmann, A. Ortega-Martinez, D.A. Boas and M.A. Yücel, 2020, *Using the General Linear Model to Improve Performance in fNIRS Single Trial Analysis and Classification: A Perspective*, Frontiers in Human Neuroscience, vol. X, no. X, doi: https://doi.org/10.1016/j.neuroimage.2019.116472 

Please cite this work when using any of the presented methods or code.

## Directory Structure
**.../block design**
- scripts to estimate hrf for the visual stimulation data using tCCA GLM

**.../data**

- provides resulting data from CrossValidation pipeline described in section 2.3.3. and Figure 3. Needed to plot results figures. Data is in .mat format.

**.../external**

- external scripts and functions used for the generation of our figures that are not authored by us. Please see the README file in this folder and the LICENSE files in the subfolders

**.../Homer functions**

- processing functions from the homer2 toolbox that our scripts are dependent on

**.../power analysis**

- contains a script to generate figure 8 and corresponding results

**.../sim HRF**

- contains synthetic HRF templates and function to augment fNIRS resting data to generate the ground truth data used for the evaluation.

**.../tCCA**

- contains the analysis scripts and functions to generate the main results and corresponding figures from the manuscript. The main function that performs tCCA GLM analysis of ground truth augmented resting state data is *run_CCA.m*. 


## How To Reproduce Figures
The following list points to the analysis scripts that generate the results and figures of the publication:
- Figure 2A: .../tCCA/fig2A_hrf_comp.m
- Figure 2B: .../tCCA/run_CCA.m
- Figure 6: .../tCCA/eval_and_plot_CVresults.m
- Figure 7: .../tCCA/eval_and_plot_invAUXresults.m
- Figure 8: .../power analysis/run_power.m  
- Figure 9: .../tCCA/eval_and_plot_CVresults.m
- Figure 10: .../tCCA/eval_and_plot_CVresults.m
- Figure 11: .../block design/results_eval_block.m
- Figure C1 and C2: .../tCCA/eval_and_plot_CVresults.m

