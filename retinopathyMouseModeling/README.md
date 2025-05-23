This repository contains all scripts and data needed to reproduce the in silico metabolic modeling analyses in Heinken et al., "Systems-level analysis of lipid metabolism in oxygen-induced retinopathy".

Please run the script ""RunAll" to execute all analyses and simulations that require MATLAB.
The simulations were run in MATLAB version 2020b, which is recommended.

Used model:
- iMM1865, retrieved from https://www.nature.com/articles/s41598-020-63235-w. 

Software requirements:
- COBRA Toolbox (https://github.com/opencobra/cobratoolbox)
- rFASTCORMICS (https://github.com/sysbiolux/rFASTCORMICS)
- IBM CPLEX solver
- Parallel Computing Toolbox, Bioinformatics Toolbox, Statistics and Machine Learning Toolbox in MATLAB
- donut function (https://uk.mathworks.com/matlabcentral/fileexchange/56833-donut)

To reproduce Figure 6a, afterwards run the R script "plotCorrelations". 
R version 4.3.0 with RStudio had been used in the paper.

- Almut Heinken, John M. Asara, Gopalan Gnanaguru, Charandeep Singh, September 2024

