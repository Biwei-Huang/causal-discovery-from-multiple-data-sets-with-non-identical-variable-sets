# causal-discovery-from-multiple-data-sets-with-non-identical-variable-sets

Copyright (c) 
 2019-2020 Biwei Huang

This package contains code to the paper:
"Biwei Huang, Kun Zhang, Mingming Gong, Clark Glymour. Causal Discovery from Multiple Data Sets with Non-Identical Variable Sets. AAAI’20"
The code is written in Matlab R2019a.


### IMPORTANT FUNCTIONS
function [A, mu, sigma, w, loglAll] = ifaEM_my(X, IDX0, numGauss, parsEM, filename)
EM algorithm to learn the mixing matrix over the intergrated set of variables
Inputs:
  *  X: data points
  *  IDX0: each column denotes the indices of variables in each data set; 
  *  numGauss: number of Gaussian for each component
  *  parsEM: parameters for EM
  *  filename: the file name to save the results

Outputs:
  *  A: mixing matrix
  *  mu: mixture of gaussian means
  *  sigma: mixture of gaussian stds for all components
  *  w: mixture of gaussian weights
  *  loglAll: log-likelihood



### EXAMPLE
example1.m gives an example of using this package.

### Notes
To avoid getting stuck on local optimum, one may randomly set different initial parameters, repeat the experiments, and finally choose the results with maximum log-likelihood. One may also set the initial parametres by leveraging the domain knowledge.


### CITATION
If you use this code, please cite the following paper:
"Biwei Huang, Kun Zhang, Mingming Gong, Clark Glymour. Causal Discovery from Multiple Data Sets with Non-Identical Variable Sets. AAAI’20"

If you have problems or questions, do not hesitate to send an email to biweih@andrew.cmu.edu
