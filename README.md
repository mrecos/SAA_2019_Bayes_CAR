# SAA 2019 - Bayesian Hierarchical CAR Model 
### April 11th, 2019 @ Albuquerque, MN


### Session: 
#### Novel Statistical Techniques in Archaeology II (QUANTARCH II) (Sponsored by SAA QUANTARCH Interest Group)


### Title: 
#### __*Estimating the Effect of Endogenous Spatial Dependency with a Hierarchical Bayesian CAR Model on Archaeological Site Location Data*__


### Authors: 
#### Matthew Harris & Mary Lennon


### Abstract:
This research presents a method to test the endogenous spatial correlation effect when modeling the landscape sensitivity for archaeological sites. The effects of endogenous spatial correlation are inferred using a Hierarchical Bayesian model with a Conditional Auto-Regressive (CAR) component to better understand the importance of modeling spatial cultural process. In current practice, effects of endogenous spatial autocorrelation are rarely explicitly incorporated into quantitative archaeological predictive models. This is due in part to the difficulties of measuring how cultural process relate across space and time, as well as accepting the assumption that geographically near sites are implicitly more related than distant sites. Typically these difficulties are side-stepped by including aspects of cultural processes as features and ignoring endogenous spatial correlation by assuming sites are spatially independent phenomena. While there are benefits to this approach, aside from convenience, the validity of either of these assumptions has not previously been tested. The approach developed here leads to better understanding the penalty for assuming spatial independence and the development of methods to model spatial cultural process.

*****

### TO DO:


#### Research

:black_square_button: - read sources below


#### Data Prep.

:black_square_button: - select study areas

:black_square_button: - get site data

:black_square_button: - create fishnet grid (regular size/shape or physio-shed based?)

:black_square_button: - calculate target (as count, area, and density)


##### Code

:black_square_button:


#### Writing

:black_square_button:


####  Poster/Graphics

:black_square_button:


### Resources:


- [Spatial Models in Stan: Intrinsic Auto-Regressive Models for Areal Data](https://mc-stan.org/users/documentation/case-studies/icar_stan.html)
  * Stan case study in spatial modeling with CAR and ICAR models, descriptions, and code

- [Exact sparse CAR models in Stan](https://mc-stan.org/users/documentation/case-studies/mbjoseph-CARStan.html)
  * Stan code and example for implimenting a sparse (and non-sparse) prior over Phi

- [Fast CAR: Two weird tricks for fast conditional autoregressive models in Stan](https://andrewgelman.com/2016/09/02/two-weird-tricks-for-fast-conditional-autoregressive-models-in-stan/)
  * Gelman blog post on sparse CAR model. Comments section has lots of relevant discussion.
  
- [Using rstan and spdep for Spatial Modelling](https://rpubs.com/chrisbrunsdon/carstan)
  * rpubs post on spatial models, CAR, Stan, and Bayes in general. Lots of intro info and code.
  
- [package brms](https://cran.r-project.org/web/packages/brms/brms.pdf)
  * brms package version of CAR with both sparse and otherwise models using cor_car() for spatial weight matrix






- [Bayesian analysis of a CAR model](https://www4.stat.ncsu.edu/~reich/SpatialStats/code/CAR.html)
  * Bayesian CAR model with homegrown MCMC function, examples, and description of the model.


- [Bayesian analysis of conditional autoregressive models](https://www.ism.ac.jp/editsec/aism/pdf/10463_2010_Article_298.pdf)
  * Detailed workup on CAR model statistics and inference. Discusses choice of priors.
  

  
- [CARBayes: An R Package for Bayesian Spatial Modeling with Conditional Autoregressive Priors](https://www.jstatsoft.org/article/view/v055i13/v55i13.pdf)
  * Article describing R package for CAR, the CAR model, and examples 
  



- [Neighborhood Dependence in Bayesian Spatial Models](https://pdfs.semanticscholar.org/60db/f7abf83011690dffd8ae62b805c475c04694.pdf)
  * In depth deiscussion of correlation between prior and posterior covariance matricies based on adjacency graph.
  
  
- [A close look at the spatial structure implied by the CAR and SAR models](https://www4.stat.ncsu.edu/~reich/CUSP/wall.car.sar.pdf)
  * A paper describing, testing, and somewhat critiquing the spatial structre matrix W in the SAR and CAR model
  * Good description of models and stats
  
  
- [NOTEBOOK FOR SPATIAL DATA ANALYSIS ](https://www.seas.upenn.edu/~ese502/NOTEBOOK/Part_III/6_Spatial_Regression_Models_for_Areal_Data_Analysis.pdf)
  * Class notes from Tony Smith at UPenn for sptail regression and touching in CAR
  
