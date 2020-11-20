## Load needed Libraries
library(RColorBrewer)
library(tidyverse)
library(gtools)
library(rstan)
library(reshape2)
library(MASS)
library(extrafont)
library(gtable)

#####################################################################################################
#####################################################################################################
### Read in DATA used by Stan to estimates ocean distributions.
### ---- see Code_ReadMe for definitions of all objects contained in 'stan_data.Rdata'
#####################################################################################################
#####################################################################################################

load("./Data/stan_data.RData") # R list that is used by the ".stan" file to run model 
load("./Data/Priors.RData")    # Priors for all parameters.  Included in stan_data.RData as well.

# Observations for all CWT recoveries.  
#Each row is an observation for a single season-ocean region- release group 
dat.all <- read.csv("./Data/dat_all.csv") 

# Observations for all CWT recoveries where CWT are greater than 0.  
# Each row is an observation for a single season-ocean region- release group 
dat.pos <- read.csv("./Data/dat_pos.csv") # Observations for  CWT recoveries

# This is the master list of 1400 release group used in the model.
REL <- read.csv("./Data/REL.csv")
AGE <- read.csv("./Data/AGE.csv")

COV <- read.csv("./Data/COV.csv")
spawn_loc <- read.csv("./Data/spawn_locations.csv")
ocean_temp_avg <- read.csv("./Data/Ocean_Temp_C.csv")
YEARS.RELEASE <- read.csv("./Data/Years_Release.csv")
YEARS.RECOVER <- read.csv("./Data/Years_Recover.csv")
YEARS.BROOD   <- read.csv("./Data/Years_Brood.csv")
knot.loc.sum.fall <- read.csv("./Data/knot_locations_summer_fall.csv")
knot.loc.wint.spr <- read.csv("./Data/knot_locations_winter_spring.csv")


# parameters to monitor during Model run.
stan_pars = c(
  "log_q_troll_pos","log_q_troll_start","log_q_troll_slope",
  "log_q_treaty_pos","log_q_treaty_start","log_q_treaty_slope",
  "log_q_rec_pos","log_q_rec_start","log_q_rec_slope",
  "log_q_rec_can_pos","log_q_rec_can_start","log_q_rec_can_slope",
  "log_q_rec_can_irec_pos",  "log_q_rec_can_irec_start",
  "log_q_hake_ashop_start", "log_q_hake_ashop_pos",
  "log_q_hake_shoreside_start","log_q_hake_shoreside_pos",
  "q_int",
  "epsilon",
  "sigma_cv",
  "sigma_cv_hake",
  "spawn_smooth",
  "prob_age_year",
  "beta_vuln",
  "beta_vuln_hake",
  "alpha_pay",
  "beta_pay",
  "rel_year_all",
  "log_rel_year_mu", 
  "log_rel_year_sigma",
  "origin_sea_int",
  "origin_mat",
  "prop_D",
  "D",
  "log_N_all",
  "F_rec",
  "F_troll",
  "F_hake_ashop",
  "log_F_rec_mean",
  "log_F_troll_mean",
  "F_rec_sigma",
  "F_troll_sigma",
  "F_troll_array",
  "F_treaty_array",
  "F_rec_array",
  "F_hake_ashop_array",
  "F_hake_shoreside_array",
  "cum_M2_temp",
  "theta_space",
  "w_star_sf",
  "w_star_ws",
  "w_star_salish"
 )

### Initial Values
  stan_init_f1 <- function(n.chain,N_loc_spawn,
                           N_knot_sf,N_knot_ws,
                           N_f_rec_idx_param,N_f_troll_idx_param,
                           N_troll_idx, N_rec_us_idx, N_sigma_cv_idx){ 
        A <- list()
  for(i in 1:n.chain){
    A[[i]] <- list(
      log_q_troll_start  = rnorm(N_troll_idx,-10,0.5),
      log_q_troll_slope  = rnorm(1,1,0.2),
      log_q_treaty_start  = rnorm(1,-10,0.5),
      log_q_treaty_slope  = rnorm(1,1,0.1),
      log_q_rec_start    = rnorm(N_rec_us_idx,-13,0.5),
      log_q_rec_slope    = rnorm(1,1,0.1),
      log_q_rec_can_start    = rnorm(1,-14,0.5),
      log_q_rec_can_slope    = rnorm(1,1,0.1),
      log_q_rec_can_irec_start    = rnorm(1,-14,0.5),
      log_q_hake_ashop_start    = rnorm(1,-14,0.5),
      log_q_hake_shoreside_start    = rnorm(1,-14,0.5),
      q_int = rnorm(1,-3,0.5),
      sigma_cv = runif(2,0.8,1.2),
      sigma_cv_hake = runif(1,0.8,1.2),
      spawn_smooth = runif(1,3,5),
      beta_vuln = runif(2,0.5,1.0),
      beta_vuln_hake = c(runif(1,0.1,0.2),runif(1,-0.05,0)),
      log_rel_year_mu = rnorm(1,0.5,0.2),
      log_rel_year_sigma = runif(1,0.2,0.5),
      log_F_rec_mean = rnorm(1,-4,0.5),
      log_F_troll_mean = rnorm(1,-4,0.5),
      F_troll = runif(N_f_troll_idx_param,0.001,0.04),
      F_rec = runif(N_f_rec_idx_param,0.001,0.04),
      log_beta_pay = rnorm(N_loc_spawn,0,0.5),
      alpha_pay_mean = rnorm(1,-5,0.5),
      alpha_pay_sd = runif(1,0.1,2),
      w_star_sf = array(runif(N_loc_spawn*2*N_knot_sf,-0.1,0.1),dim=c(2,N_loc_spawn,N_knot_sf)),
      w_star_ws = array(runif(N_loc_spawn*N_knot_ws,-0.1,0.1),dim=c(N_loc_spawn,N_knot_ws)),
      w_star_salish = array(runif(N_loc_spawn*2*3,-0.1,0.1),dim=c(3,N_loc_spawn,3))
      )
  }  
  return(A)
  }

  
# Initiate rstan options for using multiple cores.
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

###############################################################################
###############################################################################
# DEFINE STATISTICAL MODEL OF INTEREST
Warm        = 10
Iter        = 10
N_CHAIN     = 1
Treedepth   = 9
Adapt_delta = 0.5

SAMP.FILE = "./Output/FINAL.csv"

###############################################################################
###############################################################################
NAME = "Ocean_Climate_Model"

stanMod = stan(file = "./Code/Chinook_Salmon_Distribution.stan",
               data = stan_data, 
               verbose = FALSE, chains = N_CHAIN, thin = 1, 
               warmup = Warm, iter = Warm + Iter, 
               control = list(max_treedepth=Treedepth,
                              adapt_delta=Adapt_delta,
                              #adapt_init_buffer=75,
                              stepsize = 0.01,
                              metric="diag_e"),
               pars = stan_pars,
               boost_lib = NULL,
               sample_file = SAMP.FILE,
               init =   stan_init_f1(n.chain=N_CHAIN,
                                     N_loc_spawn = stan_data$N_loc_spawn,
                                     N_knot_sf = stan_data$N_knot_sf,
                                     N_knot_ws = stan_data$N_knot_ws,
                                     N_f_rec_idx_param = stan_data$N_f_rec_idx_param,
                                     N_f_troll_idx_param = stan_data$N_f_troll_idx_param,
                                     N_troll_idx=stan_data$N_troll_idx,
                                     N_rec_us_idx=stan_data$N_rec_us_idx, 
                                     N_sigma_cv_idx=stan_data$N_sigma_cv_idx),
               init_r = 1
                                   ) 

# extract information from the fitted object.
pars <- rstan::extract(stanMod, permuted = T)
samp_params <- get_sampler_params(stanMod)

base_params <- c(#"tau_process",
                 #"tau_process_prod",
                  "log_q_rec_start","log_q_troll_start","log_q_rec_can_start","log_q_treaty_start",
                 "log_q_rec_slope","log_q_troll_slope","log_q_rec_can_slope","log_q_treaty_slope",
                 "log_q_rec_can_irec_start","log_q_hake_ashop_start","log_q_hake_shoreside_start",
                 "q_int",
                 "sigma_cv","sigma_cv_hake",
                 "log_rel_year_mu", "log_rel_year_sigma",
                 "beta_vuln","beta_vuln_hake",
                 "log_F_rec_mean","log_F_troll_mean","F_rec_sigma","F_troll_sigma"
                 )

# Model Summaries.
stanMod_summary <- summary(stanMod)$summary
summary(stanMod,pars=base_params)$summary

# Write fitted model to file along with all data and options used in the estimation. 
Output <- list(stanMod=stanMod,stanMod_summary=stanMod_summary,
               converge=samp_params, 
               pars=pars,
               raw.dat.all=dat.all,
               raw.dat.pos=dat.pos,
               N_CHAIN=N_CHAIN,
               PRIORS = PRIORS,
               cum_M2_fixed = stan_data$cum_M2_fixed, 
               #origin_vec=origin_vec,
               diri_constant = stan_data$diri_constant,
               ashop_year_break = stan_data$ashop_year_break,
               stan_pars = stan_pars,
               AGE = AGE,
               COV = COV, 
               age_month_cal=stan_data$age_month_cal,
               ocean_temp_dev = stan_data$ocean_temp_dev, 
               ocean_temp = stan_data$ocean.temp, 
               ocean_temp_avg = ocean_temp_avg,
               river_entry= stan_data$river_entry,
               shaker_mort = stan_data$shaker_mort,
               knot.loc.sum.fall = knot.loc.sum.fall, 
               knot.loc.wint.spr = knot.loc.wint.spr,
               q_year_vec = stan_data$q_year_vec,
               phi_space_fix = stan_data$phi_space_fix,
               vuln_age = stan_data$vuln_age,
               vuln_age_trawl = stan_data$vuln_age_trawl, 
               vuln_fixed = stan_data$vuln_fixed,
               vuln_troll_mat = stan_data$vuln_troll_mat, 
               vuln_treaty_mat = stan_data$vuln_treaty_mat,  
               vuln_rec_mat = stan_data$vuln_rec_mat,
               spawn_time_fraction = REL$spawn_time_fraction,
               temperature_season_idx = stan_data$temperature_season_idx,
               N_years_recover = stan_data$N_years_recover,  
               N_years_release = stan_data$N_years_release, 
               N_month = stan_data$N_month, 
               N_time_mod = stan_data$N_time_mod,
               spawn_loc=spawn_loc,
               #TREATY= TREATY,SEASON=SEASON, 
               NAME = NAME, 
               REL=REL,
               loc_18 = "TWO_OR", # This is name of the spatial grouping structure in the model.
               # MONTH.STRUCTURE=MONTH.STRUCTURE, 
               first.season = "month.spring", 
               last.season= "month.fall",
               YEARS.RELEASE=YEARS.RELEASE,
               YEARS.RECOVER=YEARS.RECOVER,
               YEARS.BROOD=YEARS.BROOD)

save(Output,file=paste("./Output/",NAME,".RData",sep=""))

