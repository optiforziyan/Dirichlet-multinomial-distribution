# Inferring single- and multi-species distributional aggregation using quadrat sampling

## Simulation

###  **Folders**
#### functions
1. Dirichlet-Multinomial distribution-S1.R
* rdirichlet   - random variate generation
* rdirichlet1  - random variate generation (avoid NAs)
* dDir         - probability density (not mass!) function of Dirichlet distribution
* dMDir        - probability mass (or density) function of the multinomial-Dirichlet distribution
* dNBD         - probability mass function for independent negative binomial model 
* dNMD         - probability mass function for negative multinomial model
1.7 rMDir        - simulation of multinomial-Dirichlet distribution
1.8 rMDir1       - simulation of multinomial-Dirichlet distribution by avoiding NAs
1.9 rand         - random matrices generation
1.10 likelihood  - calculation of the negative log-likelihood function for the Dirichlet-Multinomial model
1.11 likelihood0 - calculation of the negative log-likelihood function for the null model: multinomial model
1.12 fit         - fitting of the SDM model
1.13 fitNBD      - fitting of the independent NBD model
1.14 fitNMD      - fitting of the NMD model
1.15 fitMD       - fitting of the ordinary multinomial model
2. Spatial distribution Index.R
(1) CECI()       - Clark and Evans Competition Index
(2) DC()         - Deviation Coefficient or Diffusion Coefficient 
(3) species.distribution() - Simulation of species distribution using Poisson cluster process
(4) r2.test      - calculate r2, RMSE, NRMSD based on actual values and predicted values


###  Code Availability
#### Software
1. R version 4.0.2 (2020-06-22)
#### Required R packages
1. MASS version 7.3-51.6
2. dirmult version 0.1.3-5
3. MCMCpack version 1.6-3
4. coda version 0.19-3
5. ggplot2 version 3.36
6. spatstat version 1.64-1
7. colorspace version 1.4-1
8. ggpmisc version 0.3.5
9. dgof version 1.2
10. vioplot version 0.3.7
11. doParallel version 1.0.15
## Citation

The American Naturalist (In Review)


## Author(s)

## References: 

