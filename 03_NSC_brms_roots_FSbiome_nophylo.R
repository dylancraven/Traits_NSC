################################
# Hierarchical Bayesian Model ##
# Carbohydrates ################
################################
# branches - no phylo
###############################

require(brms)
require(ape)
require(coda)
require(tidyverse)
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

trees_nsc_trts$species<-trees_nsc_trts$phylo

#############
# NSC roots #
#############

branches<- filter(trees_nsc_trts, Organ=="roots") %>%
           filter(., Season=="S2")

# branch_phy <- drop.tip(phy2, setdiff(phy2$tip.label, branches$phylo))

#scale co-variance matrix

# nspp<-length(unique(branch_phy$tip.label))
# Vphy <- vcv(branch_phy)
# Vphy <- Vphy/(det(Vphy)^(1/nspp))

#############
# quick QC ##
#############

length(unique(branches$phylo)) #  59

############

#branches$H<-scale(branches$H,center=TRUE,scale=TRUE)
branches$Biome<-as.factor(branches$Biome)
branches$Biome<-droplevels(branches$Biome)

branches$Leaf.habit<-as.factor(branches$Leaf.habit)
branches$Leaf.habit<-droplevels(branches$Leaf.habit)

branches$FS1s<-scale(branches$FS.1,center=TRUE,scale=TRUE)
branches$FS2s<-scale(branches$FS.2,center=TRUE,scale=TRUE)

#########################################################################
# test for correlations among predictor variables & variance inflation  #
#########################################################################

branches<-ungroup(branches)
stemNSC_c<-select(branches, H, FS.1, FS.2)

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
                   data = branches, 
                   family = gaussian())

# fit model

rootsNSC_bm <- brm(Concen ~Biome+Leaf.habit+FS1s+FS2s+
                        Biome:FS1s+Biome:FS2s+
                        Leaf.habit:FS1s+Leaf.habit:FS2s+
                        (1+H|species),
                      data = branches, 
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
                    data = branches, 
                    family = lognormal())

# fit model

rootsNSC_bm2 <- brm(Concen ~Biome+Leaf.habit+FS1s+FS2s+
                         Biome:FS1s+Biome:FS2s+
                         Leaf.habit:FS1s+Leaf.habit:FS2s+
                         (1+H|species),
                       data = branches, 
                       family = lognormal(),
                       prior = prior2,
                       iter = 6000, warmup = 1500,
                       sample_prior = TRUE, chains = 4, cores =4,
                       control=list(adapt_delta=0.99))

kfold1 <- kfold(rootsNSC_bm, chains = 4, cores = 4)
kfold2 <- kfold(rootsNSC_bm2, chains = 4, cores = 4)
gaus_v_log<-loo_compare(kfold1, kfold2) # keep gauss 

# Model check (convergence)

## posterior prediction checks

pp_check(rootsNSC_bm,ndraws=100)

np <- nuts_params(rootsNSC_bm)
str(np)
# extract the number of divergence transitions
sum(subset(np, Parameter == "divergent__")$Value) #

summary(rhat(rootsNSC_bm)) 

# Model evaluation

summary(rootsNSC_bm,waic=TRUE, loo=TRUE)

plot(rootsNSC_bm)

plot(conditional_effects(rootsNSC_bm), points = FALSE) 

# no interactions

prior3 <- get_prior(Concen ~Biome+Leaf.habit+FS1s+FS2s+
                      (1+H|species),
                    data = branches, 
                    family = gaussian())

# fit model

rootsNSC_bm3 <- brm(Concen ~Biome+Leaf.habit+FS1s+FS2s+
                         (1+H|species),
                       data = branches, 
                       family = gaussian(),
                       prior = prior3,
                       iter = 6000, warmup = 1500,
                       sample_prior = TRUE, chains = 4, cores =4,
                       control=list(adapt_delta=0.99))

pp_check(rootsNSC_bm3,ndraws=100)

np <- nuts_params(rootsNSC_bm3)
str(np)
# extract the number of divergence transitions
sum(subset(np, Parameter == "divergent__")$Value) #

summary(rhat(rootsNSC_bm3)) 

# Model evaluation

summary(rootsNSC_bm3,waic=TRUE, loo=TRUE)

plot(rootsNSC_bm3)

plot(conditional_effects(rootsNSC_bm3), points = FALSE) 

save(rootsNSC_bm3, rootsNSC_bm, gaus_v_log,
     file="roots_NSC_NOPHYLO_model_dec2023.RData")