# Script to create the figure with real-world dataset from PAGE poster
# "Visual Predictive Checks for Real-World Data using Propensity-Score Matching",
# by Keizer, Bergstrand, Hughes.
#
# This is a dependency-light variant of `example_poster_nonmem.R`. It is
# identical except that the propensity-score matching no longer relies on the
# `MatchIt` / `optmatch` packages -- the matching is implemented with a few
# small base-R helpers (see the "Matching helpers" section below). The result
# of the "optimal" method is identical to the MatchIt-based original.
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
# No `MatchIt` or `optmatch` install is needed -- the matching helpers below are
# pure base R.
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

# -----------------------------------------------------------------------------
# Matching helpers (replacement for MatchIt / optmatch)
#
# The propensity-score matching we need is a 1:1 matching, without replacement,
# between the n drawn parameter sets and the n observed (fitted) subjects. The
# underlying maths is simple enough to write out directly. Each helper returns,
# for every observed subject (in input order), the index of the matched draw --
# exactly the vector we need to reorder the draws.
#
# `method` options:
#   "nearest" - greedy nearest-neighbour: subjects are processed in order and
#               each takes its closest still-available draw (Mahalanobis).
#   "optimal" - minimises the TOTAL Mahalanobis distance over all pairs (the
#               linear sum assignment problem), solved with the Hungarian
#               algorithm. This is what the poster used.
#   "rank"    - cheapest: rank subjects and draws by their Mahalanobis distance
#               from the centre and pair them by rank (1st with 1st, ...).
# -----------------------------------------------------------------------------

## Mahalanobis distance matrix between two sets of points: rows of `a` (observed)
## vs rows of `b` (draws). Uses the pooled covariance of the combined points,
## mirroring MatchIt's `distance = "mahalanobis"`. Returns an nrow(a) x nrow(b)
## matrix of distances.
mahalanobis_dist <- function(a, b, Sigma = NULL) {
  a <- as.matrix(a); b <- as.matrix(b)
  if (is.null(Sigma)) Sigma <- stats::cov(rbind(a, b))
  Sinv <- solve(Sigma)
  ## (ai - bj)' Sinv (ai - bj) = ai'Sinv ai + bj'Sinv bj - 2 ai'Sinv bj
  aa <- rowSums((a %*% Sinv) * a)
  bb <- rowSums((b %*% Sinv) * b)
  ab <- a %*% Sinv %*% t(b)
  d2 <- outer(aa, bb, "+") - 2 * ab
  d2[d2 < 0] <- 0          # guard against tiny negative round-off
  sqrt(d2)
}

## Greedy nearest-neighbour matching without replacement. `D` has observed
## subjects in rows, draws in columns.
match_nearest <- function(D) {
  n <- nrow(D)
  avail <- rep(TRUE, ncol(D))
  out <- integer(n)
  for (i in 1:n) {
    j <- which(avail)[which.min(D[i, avail])]
    out[i] <- j
    avail[j] <- FALSE
  }
  out
}

## Rank matching: order observed subjects and draws by their Mahalanobis
## distance from the centre, then pair them by rank.
match_rank <- function(obs, sim, Sigma = NULL) {
  obs <- as.matrix(obs); sim <- as.matrix(sim)
  if (is.null(Sigma)) Sigma <- stats::cov(rbind(obs, sim))
  Sinv <- solve(Sigma)
  score <- function(x) sqrt(rowSums((x %*% Sinv) * x))
  out <- integer(nrow(obs))
  out[order(score(obs))] <- order(score(sim))
  out
}

## Optimal assignment via the Hungarian algorithm: returns, for each row of the
## (square) cost matrix `cost`, the matched column, minimising the total cost.
hungarian <- function(cost) {
  n <- nrow(cost); m <- ncol(cost)           # requires n <= m
  INF <- .Machine$double.xmax
  u <- numeric(n + 1); v <- numeric(m + 1)   # dual potentials (index 0 = virtual)
  p <- integer(m + 1); way <- integer(m + 1) # p[j] = row matched to column j
  for (i in 1:n) {
    p[1] <- i
    j0 <- 0L
    minv <- rep(INF, m + 1)
    used <- rep(FALSE, m + 1)
    repeat {
      used[j0 + 1] <- TRUE
      i0 <- p[j0 + 1]
      cur  <- cost[i0, ] - u[i0 + 1] - v[2:(m + 1)]
      free <- !used[2:(m + 1)]
      better <- free & (cur < minv[2:(m + 1)])
      if (any(better)) {
        minv[c(FALSE, better)] <- cur[better]
        way[c(FALSE, better)]  <- j0
      }
      cand <- ifelse(free, minv[2:(m + 1)], INF)
      j1 <- which.min(cand)
      delta <- cand[j1]
      u[p[used] + 1] <- u[p[used] + 1] + delta
      v[used] <- v[used] - delta
      minv[!used] <- minv[!used] - delta
      j0 <- j1
      if (p[j0 + 1] == 0) break
    }
    repeat {                                  # augment along the found path
      j1 <- way[j0 + 1]
      p[j0 + 1] <- p[j1 + 1]
      j0 <- j1
      if (j0 == 0) break
    }
  }
  out <- integer(n)
  for (j in 1:m) {
    r <- p[j + 1]
    if (r >= 1 && r <= n) out[r] <- j
  }
  out
}

## Unified entry point. `obs_eta` / `sim_eta` are matrices/data.frames of ETAs
## (one row per subject). Returns, for each observed subject, the matched draw
## index.
match_draws <- function(obs_eta, sim_eta, method = "optimal", Sigma = NULL) {
  if (method == "rank") return(match_rank(obs_eta, sim_eta, Sigma))
  D <- mahalanobis_dist(obs_eta, sim_eta, Sigma)   # rows = obs, cols = draws
  switch(
    method,
    nearest = match_nearest(D),
    optimal = hungarian(D),
    stop("unknown matching method: ", method)
  )
}

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
match_method <- "optimal" # "optimal", "nearest", or "rank" (see helpers above)
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

  ## Step 3: Match simulated draws with observed (fitted) subjects. `match_draws`
  ## returns, for each observed subject, the index of the matched draw -- a 1:1
  ## matching without replacement, so every subject gets exactly one draw. We
  ## use Mahalanobis distance for the "optimal" and "nearest" methods.
  obs_to_sim <- match_draws(
    obs_eta = etas[, eta_names],
    sim_eta = eta_draws[, eta_names],
    method  = match_method
  )
  ## Now we'll rearrange the order of the drawn parameter sets so that they match
  ## with their observed counterparts. In that way we can just use the observed input
  ## dataset for simulation (this is slightly easier than rearranging the dataset).
  eta_draws <- eta_draws[obs_to_sim, ] |>
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
