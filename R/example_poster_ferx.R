# FeRx version of the PAGE-poster figure with a real-world dataset,
# "Visual Predictive Checks for Real-World Data using Propensity-Score Matching"
# (Keizer, Bergstrand, Hughes).
#
# This is the FeRx counterpart of `example_poster_nonmem.R`. The whole point:
# in FeRx the propensity-score-matched VPC (pmVPC) is BUILT IN. The entire
# manual workflow from the NONMEM script -- MAXEVAL=0 to get ETAs, draw from
# omega, MatchIt propensity matching in a loop, rebuild the simulation dataset,
# resimulate in NONMEM -- collapses into a single argument: `match = TRUE`.
#
# Install FeRx (R wrapper) -- needs ferx-core >= the commit that adds
# `simulate_with_options` (ferx-core #288), and the `vpc` package for plotting:
#
#   remotes::install_github("FeRx-NLME/ferx-r")
#   install.packages("vpc")        # or remotes::install_github("ronkeizer/vpc")
#
# As in the NONMEM script we walk through Case 1 (interval-adaptation), which
# induces a bias that the regular and prediction-corrected VPCs cannot remove
# but the pmVPC can.

library(ferx)
library(vpc)
library(ggplot2)

## setwd("<path to the folder where you put this repo>")

## The model and the observed real-world dataset.
## `run1.ferx` is the FeRx translation of `nm/run1.mod` (1-cmt IV, CL/V, two
## ETAs, block omega, proportional error). Parameters are fixed at the poster
## values, mirroring the NONMEM `MAXEVAL=0 POSTHOC` run -- no estimation, the
## posthoc ETAs are computed internally.
model <- "ferx/run1.ferx"
data  <- "data/obs_nm.csv"   # NONMEM-format CSV: ID, TIME, AMT, EVID, DV, CMT

## Observed observations, for the `obs` side of every VPC below.
obs <- read.csv(data) |>
  subset(EVID == 0)

# -----------------------------------------------------------------------------
# 1. Regular VPC -- shows large (but unwarranted) bias.
#
# `ferx_simulate()` draws each subject's eta independently of its design, so the
# design<->eta association created by interval-adaptation is lost. With fixed
# parameters we don't need a fit object; the model file's values are used and
# the posthoc step is implicit.
# -----------------------------------------------------------------------------
sim_reg <- ferx_simulate(model, data, n_sim = 200)
vpc(
  sim = sim_reg,
  obs = obs,
  show = list(obs_dv = TRUE),
  obs_cols = list(dv = "DV",     idv = "TIME", id = "ID"),
  sim_cols = list(dv = "DV_SIM", idv = "TIME", id = "ID", sim = "SIM")
) +
  scale_y_log10()

# -----------------------------------------------------------------------------
# 2. Propensity-score-matched VPC (pmVPC) -- corrects the bias.
#
# `match = TRUE` is the entire correction. Per replicate, FeRx draws a pool of
# etas, optimally matches them (Mahalanobis, under the model omega) to the
# subjects' fitted (posthoc) etas, and simulates each subject's OWN observed
# design with its matched draw -- so a high-clearance draw lands on a subject
# whose adaptive design (e.g. longer interval) reflects high clearance.
#
# Each simulated subject keeps its observed ID and design, so the VPC columns
# are identical to the regular VPC above -- no OID/ITER bookkeeping needed.
# -----------------------------------------------------------------------------
sim_pm <- ferx_simulate(model, data, n_sim = 200, match = TRUE)
vpc(
  sim = sim_pm,
  obs = obs,
  show = list(obs_dv = TRUE),
  obs_cols = list(dv = "DV",     idv = "TIME", id = "ID"),
  sim_cols = list(dv = "DV_SIM", idv = "TIME", id = "ID", sim = "SIM")
) +
  scale_y_log10()
  