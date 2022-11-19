# Inferring single- and multi-species distributional aggregation using quadrat sampling

## Simulation

Here, we provide a self-explained R code (AMNAT-SDM-Simulations.R) for repoduce all the results presented in our manuscript. 

###  **Folders**
#### *functions* - Custom R Script
1. Dirichlet-Multinomial distribution-S1.R
* rdirichlet()   - random variate generation
* rdirichlet1()  - random variate generation (avoid NAs)
* dDir()         - probability density (not mass!) function of Dirichlet distribution
* dMDir()        - probability mass (or density) function of the multinomial-Dirichlet distribution
* dNBD()         - probability mass function for independent negative binomial model 
* dNMD()         - probability mass function for negative multinomial model
* rMDir()        - simulation of multinomial-Dirichlet distribution
* rMDir1()       - simulation of multinomial-Dirichlet distribution by avoiding NAs
* rand()         - random matrices generation
* likelihood()   - calculation of the negative log-likelihood function for the Dirichlet-Multinomial model
* likelihood0()  - calculation of the negative log-likelihood function for the null model: multinomial model
* fit()          - fitting of the SDM model
* fitNBD()       - fitting of the independent NBD model
* fitNMD()       - fitting of the NMD model
* fitMD()        - fitting of the ordinary multinomial model
2. Spatial distribution Index.R
* CECI()         - Clark and Evans Competition Index
* DC()           - Deviation Coefficient or Diffusion Coefficient 
* species.distribution() - Simulation of species distribution using Poisson cluster process
* r2.test()      - calculate r2, RMSE, NRMSD based on actual values and predicted values
#### *input*
* bci.spptable.rdata           - for more information, please visit: https://datadryad.org/stash/dataset/doi:10.15146/5xcp-0d46
* bci.treen.rdata (n = 1 to 8) - for more information, please visit: https://datadryad.org/stash/dataset/doi:10.15146/5xcp-0d46
#### *output*
* randomized_single_species_distribution.rds - the results of numeric simulation – Randomized spatial-point distribution data
* randomized_10_species_distribution.rds     - the results of Numeric simulation – Randomized spatial-point distribution data
* other files are the results of application SDM to BCI dataset.

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
Clark, P. J., and Evans, F. C. (1954). Distance to nearest neighbor as a measure of spatial relationships in populations. Ecology, 35, 445–453.

Blackman, G. E. (1942). Statistical and ecological studies in the distribution of species in plant communities: I. dispersion as a factor in the study of changes in plant populations. Annals of Botany, 6, 351–370. 
