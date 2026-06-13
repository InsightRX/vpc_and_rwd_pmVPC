# Script to create the figure with real-world dataset from PAGE poster 
# "Visual Predictive Checks for Real-World Data using Propensity-Score Matching",
# by Keizer, Bergstrand, Hughes.
#
# This script shows how to set up the "corrections" of the simulated data
# using propensity-score matching against the observed data and fitted subjects.
#
# To install vpc library from CRAN, or from source:
#
# `install.packages("vpc")`
#   or
# `remotes::install_github("ronkeizer/vpc")`
#
# Also install the packages for propensity-score matching, if you don't have yet:
#
# `install.packages(c("optmatch", "MatchIt"))`
#
# The datasets (also the observed) in this example were simulated to mimic
# real world data with treatment adaptation, either using dose-, interval-,
# or sampling-adaptation. The source code to run those simulations is not
# replicated here, but available here: https://github.com/InsightRX/vpc_and_rwd. 
# More details on these dataset are in the original CPT paper: 
# https://github.com/InsightRX/vpc_and_rwd
# 
# The example in this script walks through the example in the middle row
# of the poster ("Case 1"), showing biased VPCs caused by interval-adaptation in 
# the data.

library(vpc)
library(ggplot2)
library(dplyr)
library(patchwork)

## setwd("<path to the folder where you put this repo>")

## Case 1: First row in the poster, VPC in case of interval-adjustments
obs_1 <- read.csv(file = "data/obs_nm.csv")
## Run NONMEM to simulate data for regular VPC and pcVPC
## This will generate a `simtab1` file
system("cd nm && /opt/NONMEM/nm751/run/nmfe75 sim1.mod sim1.lst")
## Or, just unzip nm_output.zip, this will also create the simtab file.

## Then we can read in the simtab file
sim_1 <- vpc::read_table_nm("nm/simtab1") # simulated dataset for regular VPC

## Regular VPC: indicates lots of bias (but unwarranted, the model is unbiased)
vpc(
  sim = sim_1, 
  obs = obs_1, 
  show = list(obs_dv = TRUE)
) +
  scale_y_log10()

## Prediction-corrected VPC: doesn't correct for the bias from treatment adjustment,
## since it is interval-adjustment (and not dose-adjustments)
vpc(
  sim = sim_1, 
  obs = obs_1,
  pred_corr = TRUE,
  show = list(obs_dv = TRUE)
) +
  scale_y_log10()

## "Propensity-score matched" VPC: corrects for the bias properly.
##
## Now we need to resimulate the data in a step-wise manner
## 1. Fit the original model (MAXEVAL=0) and get the ETAs. It is better to
##    use ETAs rather than EBEs, because the bias is induced by factors not
##    explained by the structural model. EBEs will include also covariate
##    effects, and hence will likely dilute the signal during propensity-matching.
##    In our case the model had no covariates, so it doesn't matter.
## 2. Draw subjects from omega for the simulation
## 3. Perform matching based on propensity score (Draws from step 2 
##    vs the ETAs from step 1). Then pick the dose/sampling designs
##    from the observed subjects and construct simulation dataset for VPC./
## 4. Simulate PK (or PK/PD) based on the dataset constructed in step 3. 
##
## Step 1: Get the ETAs from a MAP fit in NONMEM (MAXEVAL=0). 
system("cd nm && /opt/NONMEM/nm751/run/nmfe75 run1.mod run1.lst")
## If you don't want to run NONMEM, get it from the zip file.
etas <- vpc::read_table_nm("nm/etatab1")

## Step 2-3: Draw n_sim x n individual parameter sets from the dataset
##    n_sim = number of simulations for the VPC dataset, e.g. m = 250
##    n = number of subjects in the dataset, e.g. n = 200 (in our case)
##    We'll do this in a loop. For efficiency purposes, we'll do step 3
##    (the actual propensity-score matching) in the same loop, otherwise
##    we'd have to store all the sampled etas in memory.
##    This step will take 20-60 seconds.
par_draws <- list()
n_etas <- 2 # number of ETAs
n_sim  <- 200
ids <- unique(etas$ID)
omega_matrix <- c(
  0.25000000,
  0.00050625, 0.09000000
)
parameters = list(CL = 8.66, V = 100)
for (j in 1:n_sim) {
  eta_draws <- PKPDsim::mvrnorm2(
    n = length(ids), # draw the same number of subjects as in the observed data
    mu = rep(0, n_etas),
    Sigma = PKPDsim::triangle_to_full(omega_matrix)
  ) |>
    as.data.frame()
  ## The names need to match up with the data.frame loaded in step 1:
  eta_names <- c(paste0("ETA", 1:n_etas))
  names(eta_draws) <- eta_names
  ## Then combine the simulated draws and observed ETAs into a single data.frame
  ## that we'll use in the next step for matching. The `comb_eta` data.frame
  ## will exactly have n rows with simulated subjects (type 0) and n rows with 
  ## "observed" (ETA) subjects.
  comb_eta <- bind_rows(
    eta_draws |> 
      dplyr::select(!!eta_names) |> 
      dplyr::mutate(type = 0), # Type 0 is simulated parameter sets
    etas |> 
      dplyr::select(!!eta_names) |> 
      dplyr::mutate(type = 1) # type 1 is observed datasets
  )

  ## Step 3: Match simulated with observed:
  formula_str <- paste("type ~", paste(eta_names, collapse = " + "))
  matched <- MatchIt::matchit(
    as.formula(formula_str),
    data = comb_eta,
    replace = FALSE, # without replacement, we want all subjects to match!
    method = "optimal", # we found these settings to work well, but please experiment!
    distance = "mahalanobis"
  )
  ## Turn the MatchIt object into a data.frame
  match_matrix <- matched$match.matrix |>
    data.frame() |>
    stats::setNames("obs") |>
    dplyr::mutate(obs = as.numeric(obs), sim = 1:n()) |>
    dplyr::mutate(id = 1:dplyr::n())
  ## Now we'll rearrange the order of the drawn parameter sets so that they match
  ## with their observed counterparts. In that way we can just use the observed input 
  ## dataset for simulation (this is slightly easier than rearranging the dataset).
  eta_draws <- eta_draws[match_matrix$obs,] |>
    dplyr::select_at(eta_names) |>
    dplyr::mutate(ID = ids, iteration = j)
  ## Add any fixed parameters (not needed if you just want to output the ETAs
  ## and handle simulation in NONMEM)
  par_draws[[j]] <- eta_draws  
}

## combine into single data.frame:
sim_params <- par_draws |> 
  bind_rows()

## Now we have a data.frame with the drawn individual parameter estimates 
## that we can simulate from. When doing the simulation in NONMEM, make sure 
## to set `$SIM ... NSIM=1` and not to e.g. NSIM=200, because we don't want
## NONMEM to repeat the simulation 200 times, it just has to simulate the entire 
## large input dataset once (which already has the repetitions).
obsdat <- read.csv(file = "data/obs_nm.csv")
simdat <- lapply(1:200, function(x) obsdat |> mutate(iteration = x)) |>
  bind_rows() |>
  merge(sim_params) |>
  rename(ITER = iteration) |>
  mutate(OID = ID) |>
  mutate(ID = ITER*1000 + OID) |> # ensure unique ID for NONMEM
  select(ID, TIME, AMT, EVID, DV, CMT, ETA1, ETA2, ITER, OID) |>
  arrange(ITER, ID, TIME)
write.csv(simdat, file = "./data/sim_nm.csv", quote=F, row.names=F)

## Run simulation in NONMEM
system("cd nm && /opt/NONMEM/nm751/run/nmfe75 sim1_pm.mod sim1_pm.lst")

## Now make a "regular" VPC with the simulated data, and you will
## get the pmVPC, which is (mostly) unbiased.
sim_1_pm <- read.table("nm/simtab1_pm") |>
  setNames(c("ID", "TIME", "EVID", "DV", "ITER", "OID")) |>
  dplyr::filter(EVID == 0) |>
  dplyr::mutate(ID = OID) |>
  dplyr::arrange(ITER, ID, TIME)
vpc(
  sim = sim_1_pm, 
  obs = obs_1,
  show = list(obs_dv = TRUE)
) +
  scale_y_log10()
