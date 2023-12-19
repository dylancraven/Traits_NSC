################################
# Hierarchical Bayesian Model ##
# Carbohydrates ################
################################
# leaves
###############################

require(brms)
require(ape)
require(coda)
require(tidyverse)
# require(reshape2)
# require(stringr)
require(phytools)
require(pez)
require(ape)
require(broom)
require(psych)
require(car)

########
# Data #
########

load(file="NSC_trts_fastslow_clean.RData") # no imputed data

#############
# NSC leaves #
#############

leaves<- trees_nsc_trts %>% 
          mutate(species=phylo) %>% 
          filter(.,Organ=="leaves") %>%
          filter(.,Season=="S2")

# leaf_phy <- drop.tip(phy, setdiff(phy$tip.label, leaves$phylo))
# 
# #scale co-variance matrix
# 
# nspp<-length(unique(leaf_phy$tip.label))
# Vphy <- vcv(leaf_phy)
# Vphy <- Vphy/(det(Vphy)^(1/nspp))

#############
# quick QC ##
#############

length(unique(leaves$phylo)) #  57
# length(unique(leaf_phy$tip.label)) #  57

############

#leaves$H<-scale(leaves$H,center=TRUE,scale=TRUE)
leaves$Biome<-as.factor(leaves$Biome)
leaves$Biome<-droplevels(leaves$Biome)

leaves$Leaf.habit<-as.factor(leaves$Leaf.habit)
leaves$Leaf.habit<-droplevels(leaves$Leaf.habit)

leaves$FS1s<-scale(leaves$FS.1,center=TRUE,scale=TRUE)
leaves$FS2s<-scale(leaves$FS.2,center=TRUE,scale=TRUE)

#########################################################################
# test for correlations among predictor variables & variance inflation  #
#########################################################################

leaves<-ungroup(leaves)
stemNSC_c<-select(leaves, H, FS.1, FS.2)

corr.test(stemNSC_c,method="pearson") # all below 0.7

#############
# Fit model #
#############

# 1. compare distributions: gaussian vs log-normal

# Gaussian model

prior <- get_prior(Concen ~ Biome+Leaf.habit+FS1s+FS2s+
                     Biome:FS1s+Biome:FS2s+
                     Leaf.habit:FS1s+Leaf.habit:FS2s+
                     (1+H|species),
                   data = leaves, 
                   family = gaussian())

# fit model

leavesNSC_bm <- brm(Concen ~Biome+Leaf.habit+FS1s+FS2s+
                      Biome:FS1s+Biome:FS2s+
                      Leaf.habit:FS1s+Leaf.habit:FS2s+
                      (1+H|species),
                    data = leaves, 
                    family = gaussian(), 
                    prior = prior,
                    iter = 6000, warmup = 1500,
                    sample_prior = TRUE, chains = 4, cores =4,
                    control=list(adapt_delta=0.99))

# Lognormal model

prior2 <- get_prior(Concen ~Biome+Leaf.habit+FS1s+FS2s+
                      Biome:FS1s+Biome:FS2s+
                      Leaf.habit:FS1s+Leaf.habit:FS2s+
                      (1+H|species),
                    data = leaves,
                    family = lognormal())

# fit model

leavesNSC_bm2 <- brm(Concen ~Biome+Leaf.habit+FS1s+FS2s+
                       Biome:FS1s+Biome:FS2s+
                       Leaf.habit:FS1s+Leaf.habit:FS2s+
                       (1+H|species),
                     data = leaves,
                     family = lognormal(), 
                     prior = prior2,
                     iter = 6000, warmup = 1500,
                     sample_prior = TRUE, chains = 4, cores =4,
                     control=list(adapt_delta=0.99))

# compare with k-fold 

kfold1 <- kfold(leavesNSC_bm, chains = 4, cores = 4)
kfold2 <- kfold(leavesNSC_bm2, chains = 4, cores = 4)
gaus_v_log<-loo_compare(kfold1, kfold2) # keep gaussian

# Model check (convergence)

## posterior prediction checks

pp_check(leavesNSC_bm,ndraws=100)

np <- nuts_params(leavesNSC_bm)
str(np)
# extract the number of divergence transitions
sum(subset(np, Parameter == "divergent__")$Value) #

summary(rhat(leavesNSC_bm)) 

# Model evaluation

summary(leavesNSC_bm,waic=TRUE, loo=TRUE)

plot(conditional_effects(leavesNSC_bm), points = FALSE) 

# remove non-significant interactions

prior3 <- get_prior(Concen ~ Biome+Leaf.habit+FS1s+FS2s+
                      Biome:FS1s+
                      (1+H|species),
                    data = leaves, 
                    family = gaussian())

# fit model

leavesNSC_bm3 <- brm(Concen ~Biome+Leaf.habit+FS1s+FS2s+
                       Biome:FS1s+
                       (1+H|species),
                     data = leaves, 
                     family = gaussian(), 
                     prior = prior3,
                     iter = 6000, warmup = 1500,
                     sample_prior = TRUE, chains = 4, cores =4,
                     control=list(adapt_delta=0.99))

# kfold3 <- kfold(leavesNSC_bm3, chains = 4, cores = 4)
# only_Biome<-loo_compare(kfold1, kfold3) # keep gaussian

pp_check(leavesNSC_bm3,ndraws = 100)

np <- nuts_params(leavesNSC_bm3)
str(np)
# extract the number of divergence transitions
sum(subset(np, Parameter == "divergent__")$Value) #

summary(rhat(leavesNSC_bm3)) 

# Model evaluation

summary(leavesNSC_bm3,waic=TRUE, loo=TRUE)

plot(conditional_effects(leavesNSC_bm3), points = FALSE) 

#
save(leavesNSC_bm, leavesNSC_bm3, gaus_v_log,
     file="leaves_NSC_NOPHYLO_model_dec2023.RData")