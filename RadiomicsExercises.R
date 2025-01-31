###############################################################################
###############################################################################
## 
## Radiomics prediction Exercises
## On the basis of FMradio 1.1
## Course: Statistics for Omics: Radiomics
## Place:  Amsterdam, 28/06/2019
##
## Code:   Carel F.W. Peeters
##         Department of Epidemiology & Biostatistics
##         VU University medical center Amsterdam
##         Amsterdam, the Netherlands
##         cf.peeters@vumc.nl
## Date:   22/06/2019, Amsterdam, VUmc
##
###############################################################################
###############################################################################



#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Preliminaries**
#' **------------------------------------------------------------------------**




## Set working directory
setwd("/Users/raindy/Documents/GitHub/RadiomicsSession_Exercrise")



## Install necessary packages
pkgs <- rownames(installed.packages())
if (!"DandEFA" %in% pkgs) install.packages("DandEFA")
if (!"FMradio" %in% pkgs) install.packages("FMradio")
if (!"pec" %in% pkgs) install.packages("pec")
if (!"rags2ridges" %in% pkgs) install.packages("rags2ridges")
if (!"randomForestSRC" %in% pkgs) install.packages("randomForestSRC")
if (!"survival" %in% pkgs) install.packages("survival")

if (!"Biobase" %in% pkgs || !"RBGL" %in% pkgs) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  if (!"Biobase" %in% pkgs) {
    BiocManager::install("Biobase")
  }
  if (!"RBGL" %in% pkgs) {
    BiocManager::install("RBGL")
  }
}

## Load all packages
library(FMradio)
library(Biobase)
library(rags2ridges)
library(DandEFA)
library(randomForestSRC)
library(pec)
library(survival)
library(RBGL)

## Convenience functions
source("Convenience.R")




#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 1: Get acquainted with the data** [DONE]
#' **------------------------------------------------------------------------**

## Data packaged as expressionset
## Will invoke basic functions for looking at data

## Load data, get to know objects
load("AnonymousRadio.Rdata")

## Look at radiomic measurements
exprs(AnonymousRadio)[10:30,1:5]

## Look at Clinical (sample) information
head(pData(AnonymousRadio))


#' **---------------Answers to Exercise 1-------------------**
# This is a matrix containing measurements of 𝑝= 432 radiomic features (rows) on 𝑛=174 samples (columns),  
# and with 0 (no death observed) or 1 (death observed) as well as the follow-up time for each patient.



###############################################################################
###############################################################################
## 
## Section 1: Factor Analytic Projection
##
###############################################################################
###############################################################################

#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 2: Assess redundancy in the correlation matrix** [DONE]
#' **------------------------------------------------------------------------**

## Scale data and get raw correlation matrix
DATAscaled <- scale(t(exprs(AnonymousRadio)))
R          <- cor(DATAscaled)

## Redundancy visualization
radioHeat(R, diag = FALSE,     # correlation matrix
          threshold = TRUE,    # logical
          threshvalue = .95,   # logical
          labelsize = .01)     # value for thresholding
                               # Red indicates a positive correlation above the threshold (>0.95), 
                               # blue indicates a negative correlation below the threshold (<-0.95).


#' **---------------Answers to Exercise 2-------------------**
# Our aim is to project the high-dimensional (and highly collinear) radiomic feature-space onto a lower-dimensional (and non-collinear) latent meta-feature space.
# These latent meta-features may then be used in downstream classification and prediction.

# Since the features have different scales and as the variability of the features may differ substantially,
# we first scale the (transposed) data. This scaling is typically in high-dimensional prediction problems, but not a necessity for finding the correlation matrix.
# Then, we visualize feature-redundancy through a thresholded correlation heatmap.

# The resulting visualization gives a thresholded correlation heatmap of the radiomic features.
# To visualize only correlations equaling or exceeding an absolute marginal correlation threshold of 0.95,  
# we set those absolute correlations < 0.95 to zero.

# Colored regions (red & blue) then visualize blocks of information redundancy or strong collinearity, 
# which means the information contained in some features is almost completely represented by other features.



# Through the thresholded correlation heatmap, we noticed that there is quite some redundant information in the data. 
# Hence, one could consider redundancy-filtering (in Exercris 3) as a preprocessing step. 
# The goal would be to remove truly redundant features.

#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 3: Redundancy-filter the correlation matrix** [DONE]
#' **------------------------------------------------------------------------**

## Redundancy filtering
## And subsetting data
## 124 features remain
RFs         <- RF(R, t = .95)             # R: correlation matrix; t: thresholding value
DATAscaledS <- subSet(DATAscaled, RFs)    # DataScaledS: data matrix or expressionset; RFs: redundancy-filtered correlation matrix






#' **---------------Answers to Exercise 3-------------------**
# The `RF` function removes the minimal number of redundant features, 
# therefore no pair of features remains with an absolute marginal correlation that equals or exceeds t (0.95).
# The `subSet` function is a convenience function that subsets the data to those features that were retained after redundancy-filtering.

# Show the number remaining features in the matrix
dim(RFs)[2] ## 124 features remain

# Actually, these features will still be highly collinear in the sense that they share many high-correlations. 
# However, they are not redundant from the perspective of the specified argument t (0.95). 
# The filtered correlation matrix will be the main input for the next step (Exercise 4).



#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 4: Find an optimal regularized correlation matrix** [DONE]
#' **------------------------------------------------------------------------**

## Optimal penalty
set.seed(303)
OPT <- regcor(DATAscaledS, fold = 5)

## Look at optimal penalty-value
## Obtain regularized correlation matrix
## Conditioning can again be assessed with, e.g., CNplot from rag2ridges
OPT$optPen       # Optimal penalty parameter
Re = OPT$optCor  # Correlation estimate under optimal penalty parameter


#' **---------------Answers to Exercise 4-------------------**
## Q: Why should we find an optimal regularized correlation matrix after removing the redundant features?

## After removing the truly redundant features, we apply regularization (i.e., penalization) ensures that the resulting estimate will be numerically stable.
## This numerical stability is needed for the factor analytic projection.
## Also, the resulting regularized correlation matrix will be basis for the factor analysis procedure.


# The `regcor` function determines the optimal value of the penalty parameter, 
# by applying the Brent algorithm (1971) to the K-fold cross-validated negative log-likelihood score (using a "regularized ridge estimator" for the correlation matrix). 
# The search for the optimal value is automatic. The penalty value at which the K-fold cross-validated negative log-likelihood score is minimized is deemed optimal. 

# The output is an object of class `list` (in the upper right corner of Environment), 
# wit `$optPen` containing the optimal penalty value (0.05595778), 
# and `$optCor` containing the correlation matrix under the optimal value of the penalty parameter.



#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 5: Perform factor analysis on the regularized correlation matrix** [DONE]
#' **------------------------------------------------------------------------**

## Assess dimensionality factor solution / Calculate (upper-) bounds to the dimensionality of the latent vector
## 13 considered upper bound (First.lower-bound)
## Variance explained would suggest 8 factors (Second.lower-bound)
dimGB(Re)
dimVAR(Re, 15, graph = TRUE)

## Assessing solutions around 8
## 9th factor seems weak
## Will keep solution at 8
## ML factor analysis with Varimax rotation
fito <- mlFA(R = Re, m = 8)                                # mlFA: main function for performing a factor analysis ; R = input (the regularized correlation matrix) ; m = the number of latent meta-features (numeric indicating)
print(fito$Loadings, digits = 2, cutoff = .3, sort = TRUE)





## Visualizing solution using Dandelion plot
dandpal <- rev(rainbow(100, start = 0.4, end = 0.6))
dandelion(fito$Loadings, bound = .3, mcex = c(1,1), palet = dandpal)

## Export to pdf for inspection
pdf(file = "Dandelion.pdf", width = 11, height = 11)
dandelion(fito$Loadings, bound = .3, mcex = c(1,1), palet = dandpal)
dev.off()



#' **---------------Answers to Exercise 5-------------------**
#' Factor analysis is a technique that assumes that observed variables can be grouped based on their correlation into a lower-dimensional linear combination of latent features. 
#' Hence, a main goal of this technique is to project the information contained in a higher-dimensional observation space onto a lower-dimensional latent meta-feature space. 
#' This projection can be done orthogonally, such that the latent features are uncorrelated with each other.
#' Here, we perform a form of regularized factor analysis, for which one of the main inputs will be our previously obtained regularized correlation matrix.



#' Usually, the dimension of `m` (the number of latent meta-features) will be unknown a priori.
#' Therefore, we precede the factor analysis with the `dimGB` function, which calculates (upper-) bounds to the dimensionality of the latent vector.

#' The choice of the number of factors can be further assessed with the `dimVAR` function, 
#' which returns the total cumulative explained variance against the dimension of the factor solution possibly with visualization, 
#' informing if the result of the first bound should be accepted or be treated as an upper-bound.

#' In high-dimensional situations, Peeters et al. (2019) suggested that the first reported bound (13) provides a reliable upper-bound to the latent dimension.


#' The `mlFA` function performs a maximum likelihood factor analysis with an orthogonal simple structure rotation, which returns an object of class `list`. 

#' `fito$Loadings` contains a matrix whose elements represent the loadings of the observed features on the latent meta-features, 
#' which may be used to interpret the structure of the factor solution and the possible meaning of the latent meta-features. 

#' `fito$Uniqueness` contains a diagonal matrix whose elements represent the unique variances (variance not explain by the model).



#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 6: Obtain factor scores** [DONE]
#' **------------------------------------------------------------------------**

## Factor scores
Lambda <- fito$Loadings
Psi    <- fito$Uniqueness
Scores <- facScore(DATAscaledS, Lambda, Psi)       # DATAscaledS = (subsetted) data matrix; Lambda = loading matrix, Psi = uniquenesses matrix



#' **---------------Answers to Exercise 6-------------------**
#' Once a factor model is fitted one may desire an estimate of the score each object/individual would obtain on each of the latent factors.
#' Such scores are referred to as factor scores. These scores can then serve as the low-dimensional and non-collinear input into (standard) prediction procedures.

#' The function returns an object of class `data.frame`. Observations are represented in the rows. Each column represent a latent factor (8 factors in total).




#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 7: Assess factor scores** [DONE]
#' **------------------------------------------------------------------------**

## Determinacy factor scores
## Highly determinate
DF <- facSMC(Re, Lambda); DF      # Re = (regularized) correlation matrix, Lambda = loadings matrix




#' **---------------Answers to Exercise 7-------------------**
#' This step was performed to assess the quality of the factor extraction.
#' This can be probed through the `squared multiple correlations` between the observed features and the common latent factors. 
#' It indicates how well a common factor can be predicted by the observed features or, from another perspective, to what extent the factor scores are determinate.

#' 
#' `facSMC` is a function with which one may evaluate the factor scores,
#' returning a numeric vector indicating, for each factor, the squared multiple correlation between the observed features and the common latent factor.

#' The closer to unity (1), the better one is able to `uniquely` determine the factor scores. 
#' And the more useful these scores are in prediction and classification.

#'   Factor1   Factor2   Factor3   Factor4   Factor5   Factor6   Factor7   Factor8 
#' 0.9904875 0.9841129 0.9863895 0.9616917 0.9709467 0.9465867 0.9456656 0.9467385 



###############################################################################
###############################################################################
## 
## Section 2: Prediction
##
###############################################################################
###############################################################################


#' Our general aim would then be to predict survival probabilities based on the prognostic features available at baseline. 
#' These prognostic features are our factor scores. 
#' As the factor scores represent low-dimensional and non-collinear data, we may use standard modeling techniques (Cox proportional hazards model).


#' Here we have a right-censored survival (time to event) outcome. 
#' The response concerns the time from a starting point 𝑡= 0 to an event of interest or to censoring. 
#' Tied to this response is a status indicator that is 1 when the event of interest has occurred during the study-time, 
#' and 0 indicates when the event of interest did not occur during the study-time or when the subject was lost to follow-up (right-censoring).


#' We will also take interest in comparing the performance of the simple model with the projected data with a popular ML technique often used on radiomics data: random forests. 
#' The random survival forest is a learning method in which the predictor is an ensemble formed by combining many decision trees over the original observed radiomic features.


#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 8: Concatenate original and projected data** [DONE]
#' **------------------------------------------------------------------------**

## Combine original (scaled) data with projected meta-features
DAT <- cbind(DATAscaled, Scores)

## Include the survival information
Status  <- as.numeric(AnonymousRadio$Death_yesno) - 1
time    <- AnonymousRadio$Death_followuptime_months
DAT     <- cbind(time, Status, DAT)



#' **---------------Answers to Exercise 8-------------------**
#' Here we use `cbind` to conveniently column-wise collect all data of interest. 
#' The data of interest concern the observed radiomic features, the factor scores, as well as the time-to-event data and its associated event-indicator.



#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 9: Set up model comparisons** [DONE]
#' **------------------------------------------------------------------------**

## Formulating the model formula's
FitRSF     <- as.formula(paste("Surv(time, Status)~", 
                               paste(colnames(DAT)[3:434], collapse="+")))
FitMetaCox <- as.formula(paste("Surv(time, Status) ~", 
                               paste(colnames(DAT[,c(435:442)]), collapse="+")))

models <- list("MetaCox" = coxph(FitMetaCox, data = DAT, x = TRUE, y = TRUE),
               "RforestCox" = rfsrc(FitRSF, data = DAT))



#' **---------------Answers to Exercise 9-------------------**
#' Evaluation metrics used in the comparison

#' 1. Prediction error: Time-dependent Brier score; Can be seen as mean square error of prediction.

#' 2. R^2: Measure of (overall) explained residual variation; Can be calculated from the Brier score.

#' 3. Apparent and calibrated performance: 
#' Apparent error: prediction error obtained when the data used from training are re-used for validation.
#' Calibrated error: cross-validated prediction error.




#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 10: Compare models w.r.t. prediction error** [DONE]
#' **------------------------------------------------------------------------**

## Assessing prediction error
## Median follow-up time = 25.7
## (Averaged) repeated 5-fold cross-validation
set.seed(446464)
PredError <- pec(object = models,                     # list of prediction models
                 formula = Surv(time, Status) ~ 1,    # censoring model
                 data = DAT,                          # data frame
                 exact = TRUE,                        # logical determining time-points per curve estimation 
                 maxtime = median(time),              # max time for estimating per curve
                 cens.model = "marginal",             # method estimation inv.prob. censoring weights
                 splitMethod = "cv5",                 # calibration method
                 B = 50,                              # number of repeats/cycles
                 verbose = TRUE)

## Summary results:
## Apparent and cross-validated prediction error and R2
crps(PredError)          # produce a table with apparent and cross-validated errors
Or2(crps(PredError))


#' **---------------Answers to Exercise 10-------------------**
#' For model evaluation we focus on prediction error, which is a measure of inaccuracy and smaller measures are desired.
#' We use the time-dependent Brier score (Brier, 1950), which may be seen as a mean square error of prediction. 
#' It can also be used to calculate a measure of (overall) explained residual variation (𝑅^2).
#' 
#' Since we have only one dataset available, so we have to harness ourselves from overoptimistic prediction errors.
#' Hence, we will use a cross-validated prediction error, with 5-fold cross-validation to this end. 

#' To mitigate the dependency on the random-choice of folds, we repeat the 5-fold cross-validation many times (B = 50) and average the corresponding Brier scores over these repeats. 
#' This can be conveniently done with the pec function from the pec package.


#' Calling the `crps` function will produce a table with apparent and cross-validated errors. 
#' The apparent error is the prediction error obtained when the data used from training are re-used for validation, which will result in optimistic prediction errors.
#' And the cross-validation will harness against this overoptimism.

#' The Cox model with the latent features and the random forest with the original observed radiomic features the table will also give apparent and cross-validated errors for a reference model. 
#' The reference model is the Kaplan-Meier prediction rule, which acts as a null model by ignoring all feature information. One desires a model to do better than the reference model.



#' Looking at the results, we see that both models do better than the reference model. 
#' The simple model with the factor-analytic meta-features has a better cross-validated prediction error than the random forest model. 

#' Moreover, it has a much smaller gap between the apparent and cross-validated errors. 
#' These results are mirrored in the (apparent and cross-validated) residual explained varation which can be obtained with the Or2 function from the source file.



#'                                BI                  |           R2 
#'                       Apparent    Cross-validated  | Apparent     Cross-validated 
#'Reference model          .128          .130              –               -        
#'FMradio                  .098          `.108`           .236             `.170` 
#'Random survival forest   .061          .114            .526             .119


#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Exercise 11: Visualize the results** [DONE]
#' **------------------------------------------------------------------------**

## Visualize apparent prediction error
plot(PredError, what = "AppErr",
     xlab = "Time (months)",
     ylab = "Apparent prediction error",
     legend.cex = .9,
     legend.lty = 1,
     legend.lwd = 2,
     legend.legend = c("Reference model", 
                       "FMradio",
                       "Random survival forest"),
     add.refline = TRUE,
     lwd = 1.5,
     legend.y.intersp = 1.7)

## Visualize cross-validated prediction error
plot(PredError, what = "crossvalErr",
     xlab = "Time (months)",
     ylab = "Averaged cross-validated prediction error",
     legend.cex = .9,
     legend.lty = 1,
     legend.lwd = 2,
     legend.legend = c("Reference model", 
                       "FMradio",
                       "Random survival forest"),
     add.refline = TRUE,
     lwd = 1.5,
     legend.y.intersp = 1.7)

## Visualize apparent residual explained variation 
R2Table <- R2(PredError, times = seq(0,median(time),.01), reference = 1)
plotR2Table(R2Table, "AE", Xlab = "Time (months)")

## Visualize cross-validated residual explained variation 
plotR2Table(R2Table, "CV", Xlab = "Time (months)")


## Exporting
pdf("Internal.pdf", width = 15, height = 15)
par(mfrow=c(2,2))
## Visualize apparent prediction error
plot(PredError, what = "AppErr",
     xlab = "Time (months)",
     ylab = "Apparent prediction error",
     legend.cex = .9,
     legend.lty = 1,
     legend.lwd = 2,
     legend.legend = c("Reference model", 
                       "FMradio",
                       "Random survival forest"),
     add.refline = TRUE,
     lwd = 1.5,
     legend.y.intersp = 1.7)

## Visualize cross-validated prediction error
plot(PredError, what = "crossvalErr",
     xlab = "Time (months)",
     ylab = "Averaged cross-validated prediction error",
     legend.cex = .9,
     legend.lty = 1,
     legend.lwd = 2,
     legend.legend = c("Reference model", 
                       "FMradio",
                       "Random survival forest"),
     add.refline = TRUE,
     lwd = 1.5,
     legend.y.intersp = 1.7)

## Visualize apparent residual explained variation 
R2Table <- R2(PredError, times = seq(0,median(time),.01), reference = 1)
plotR2Table(R2Table, "AE", Xlab = "Time (months)")

## Visualize cross-validated residual explained variation 
plotR2Table(R2Table, "CV", Xlab = "Time (months)")
dev.off()


###############################################################################
###############################################################################
## 
## Section 3: Hidden Gems
##
###############################################################################
###############################################################################

## Really?
FMradio:::.Airwolf()

## You betcha!
FMradio:::.Airwolf2()


