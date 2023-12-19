################################
# Hierarchical Bayesian Model ##
# Carbohydrates ################
################################
# roots
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
# NSC roots #
#############

roots<- trees_nsc_trts %>% 
        mutate(species=phylo) %>% 
        filter(.,Organ=="roots") %>%
        filter(.,Season=="S2")

root_phy <- drop.tip(phy, setdiff(phy$tip.label, roots$phylo))

#scale co-variance matrix

nspp<-length(unique(root_phy$tip.label))
Vphy <- vcv(root_phy)
Vphy <- Vphy/(det(Vphy)^(1/nspp))

#############
# quick QC ##
#############

length(unique(roots$phylo)) #  59
length(unique(root_phy$tip.label)) #  59

############

roots$Biome<-as.factor(roots$Biome)
roots$Biome<-droplevels(roots$Biome)

roots$Leaf.habit<-as.factor(roots$Leaf.habit)
roots$Leaf.habit<-droplevels(roots$Leaf.habit)

roots$FS1s<-scale(roots$FS.1,center=TRUE,scale=TRUE)
roots$FS2s<-scale(roots$FS.2,center=TRUE,scale=TRUE)

#########################################################################
# test for correlations among predictor variables & variance inflation  #
#########################################################################

roots<-ungroup(roots)
rootNSC_c<-select(roots, H, FS.1, FS.2)

p_out<-corr.test(rootNSC_c,method="pearson") # all below 0.7

#############
# Fit model #
#############

# 1. compare distributions: gaussian vs log-normal

# Gaussian model

prior <- get_prior(Concen ~Biome+Leaf.habit+FS1s+FS2s+
                     Biome:FS1s+Biome:FS2s+
                     Leaf.habit:FS1s+Leaf.habit:FS2s+
                     (1+H|species)+(1|gr(phylo, cov=Vphy)),
                   data = roots,
                   data2=list(Vphy=Vphy),
                   family = gaussian())

# fit model

rootsNSC_bm <- brm(Concen ~Biome+Leaf.habit+FS1s+FS2s+
                     Biome:FS1s+Biome:FS2s+
                     Leaf.habit:FS1s+Leaf.habit:FS2s+
                     (1+H|species)+(1|gr(phylo, cov=Vphy)),
                   data = roots,
                   data2=list(Vphy=Vphy),
                   family = gaussian(), 
                   prior = prior,
                   iter = 6000, warmup = 1500,
                   sample_prior = TRUE, chains = 4, cores =4,
                   control=list(adapt_delta=0.99, max_treedepth = 12))

# Lognormal model

prior2 <- get_prior(Concen ~Biome+Leaf.habit+FS1s+FS2s+
                      Biome:FS1s+Biome:FS2s+
                      Leaf.habit:FS1s+Leaf.habit:FS2s+
                      (1+H|species)+(1|gr(phylo, cov=Vphy)),
                    data = roots,
                    data2=list(Vphy=Vphy),
                    family = lognormal())

# fit model

rootsNSC_bm2 <- brm(Concen ~Biome+Leaf.habit+FS1s+FS2s+
                      Biome:FS1s+Biome:FS2s+
                      Leaf.habit:FS1s+Leaf.habit:FS2s+
                      (1+H|species)+(1|gr(phylo, cov=Vphy)),
                    data = roots,
                    data2=list(Vphy=Vphy),
                    family = lognormal(),
                    prior = prior2,
                    iter = 6000, warmup = 1500,
                    sample_prior = TRUE, chains = 4, cores =4,
                    control=list(adapt_delta=0.99, max_treedepth = 12))

# compare with k-fold 

kfold1 <- kfold(rootsNSC_bm, chains = 4, cores = 4)
kfold2 <- kfold(rootsNSC_bm2, chains = 4, cores = 4)
gaus_v_log<-loo_compare(kfold1, kfold2) # keep gauss (similar results)

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


# model with no interactions

prior3 <- get_prior(Concen ~Biome+Leaf.habit+FS1s+FS2s+
                     (1+H|species)+(1|gr(phylo, cov=Vphy)),
                   data = roots,
                   data2=list(Vphy=Vphy),
                   family = gaussian())

# fit model

rootsNSC_bm3 <- brm(Concen ~Biome+Leaf.habit+FS1s+FS2s+
                     (1+H|species)+(1|gr(phylo, cov=Vphy)),
                   data = roots,
                   data2=list(Vphy=Vphy),
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

plot(conditional_effects(rootsNSC_bm3), points = FALSE) 

save(rootsNSC_bm,rootsNSC_bm3, gaus_v_log,
     file="Roots_NSC_model_dec2023.RData")