# Readme for GLM-BCI Repository
## Purpose and Publication
This README provides a brief documentation for the code and data provided in this repository that is linked with the following [publication](https://www.frontiersin.org/articles/10.3389/fnhum.2020.00030):
> A. von Lühmann, A. Ortega-Martinez, D.A. Boas and M.A. Yücel, 2020, *Using the General Linear Model to Improve Performance in fNIRS Single Trial Analysis and Classification: A Perspective*, Frontiers in Human Neuroscience, vol. X, no. X, doi: XXX 

Please cite this work when using any of the presented methods or code.

## Directory Structure
**.../aux functions**
- scripts to estimate hrf for the visual stimulation data using tCCA GLM

**.../Homer functions**

- processing functions from the homer2 toolbox that our scripts are dependent on

**.../lit rev data**

- provides resulting data from CrossValidation pipeline described in section 2.3.3. and Figure 3. Needed to plot results figures. Data is in .mat format.

**.../main**

- contains the analysis scripts and functions to generate the main results and corresponding figures from the manuscript. The main function that performs tCCA GLM analysis of ground truth augmented resting state data is *run_CCA.m*. 

**.../sim HRF**

- contains synthetic HRF templates and function to augment fNIRS resting data to generate the ground truth data used for the evaluation.




## How To Reproduce Figures
The following list points to the analysis scripts that generate the results and figures of the publication:
- Figure 7: .../tCCA/eval_and_plot_invAUXresults.m
- Figure 9: .../tCCA/eval_and_plot_CVresults.m
- Figure 11: .../block design/results_eval_block.m
- Figure 12: .../tCCA/eval_and_plot_CVresults.m

