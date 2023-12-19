################################
# Hierarchical Bayesian Model ##
# Carbohydrates ################
################################
# stems
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

#############
# NSC roots #
#############

stems<- trees_nsc_trts %>% 
        mutate(species=phylo) %>% 
        filter(.,Organ=="stems") %>%
        filter(.,Season=="S2")

stem_phy <- drop.tip(phy, setdiff(phy$tip.label, stems$phylo))

#scale co-variance matrix

nspp<-length(unique(stem_phy$tip.label))
Vphy <- vcv(stem_phy)
Vphy <- Vphy/(det(Vphy)^(1/nspp))

#############
# quick QC ##
#############

length(unique(stems$phylo)) #  53
length(unique(stem_phy$tip.label)) #  53

############

#stems$H<-scale(stems$H,center=TRUE,scale=TRUE)
stems$Biome<-as.factor(stems$Biome)
stems$Biome<-droplevels(stems$Biome)

stems$Leaf.habit<-as.factor(stems$Leaf.habit)
stems$Leaf.habit<-droplevels(stems$Leaf.habit)

stems$FS1s<-scale(stems$FS.1,center=TRUE,scale=TRUE)
stems$FS2s<-scale(stems$FS.2,center=TRUE,scale=TRUE)

#########################################################################
# test for correlations among predictor variables & variance inflation  #
#########################################################################

stems<-ungroup(stems)
stemNSC_c<-select(stems, H, FS.1, FS.2)

corr.test(stemNSC_c,method="pearson") # all below 0.7

#############
# Fit model #
#############

# 1. compare distributions: gaussian vs log-normal

# Gaussian model

prior <- get_prior(Concen ~ Biome+Leaf.habit+FS1s+FS2s+
                     Biome:FS1s+Biome:FS2s+
                     Leaf.habit:FS1s+Leaf.habit:FS2s+
                     (1+H|species)+(1|gr(phylo, cov=Vphy)),
                   data = stems,
                   data2=list(Vphy=Vphy),
                   family = gaussian())

# fit model

stemsNSC_bm <- brm(Concen ~Biome+Leaf.habit+FS1s+FS2s+
                     Biome:FS1s+Biome:FS2s+
                     Leaf.habit:FS1s+Leaf.habit:FS2s+
                     (1+H|species)+(1|gr(phylo, cov=Vphy)),
                   data = stems,
                   data2=list(Vphy=Vphy),
                   family = gaussian(), 
                   prior = prior,
                   iter = 6000, warmup = 1500,
                   sample_prior = TRUE, chains = 4, cores =4,
                   control=list(adapt_delta=0.99))

# Lognormal model

prior2 <- get_prior(Concen ~Biome+Leaf.habit+FS1s+FS2s+
                      Biome:FS1s+Biome:FS2s+
                      Leaf.habit:FS1s+Leaf.habit:FS2s+
                      (1+H|species)+(1|gr(phylo, cov=Vphy)),
                    data = stems,
                    data2=list(Vphy=Vphy),
                    family = lognormal())

# fit model

stemsNSC_bm2 <- brm(Concen ~Biome+Leaf.habit+FS1s+FS2s+
                      Biome:FS1s+Biome:FS2s+
                      Leaf.habit:FS1s+Leaf.habit:FS2s+
                      (1+H|species)+(1|gr(phylo, cov=Vphy)),
                    data = stems,
                    data2=list(Vphy=Vphy),
                    family = lognormal(), 
                    prior = prior2,
                    iter = 6000, warmup = 1500,
                    sample_prior = TRUE, chains = 4, cores =4,
                    control=list(adapt_delta=0.99))

# compare with k-fold 

kfold1 <- kfold(stemsNSC_bm, chains = 4, cores = 4)
kfold2 <- kfold(stemsNSC_bm2, chains = 4, cores = 4)
gaus_v_log<-loo_compare(kfold1, kfold2) # keep log 

# Model check (convergence)

## posterior prediction checks

pp_check(stemsNSC_bm2,ndraws=100)

np <- nuts_params(stemsNSC_bm2)
str(np)
# extract the number of divergence transitions
sum(subset(np, Parameter == "divergent__")$Value) #

summary(rhat(stemsNSC_bm2)) 

# Model evaluation

summary(stemsNSC_bm2,waic=TRUE, loo=TRUE)

plot(stemsNSC_bm2)

plot(conditional_effects(stemsNSC_bm2), points = FALSE) 

# model with no interactions

prior3 <- get_prior(Concen ~Biome+Leaf.habit+FS1s+FS2s+
                      (1+H|species)+(1|gr(phylo, cov=Vphy)),
                    data = stems,
                    data2=list(Vphy=Vphy),
                    family = lognormal())

# fit model

stemsNSC_bm3 <- brm(Concen ~Biome+Leaf.habit+FS1s+FS2s+
                      (1+H|species)+(1|gr(phylo, cov=Vphy)),
                    data = stems,
                    data2=list(Vphy=Vphy),
                    family = lognormal(), 
                    prior = prior3,
                    iter = 6000, warmup = 1500,
                    sample_prior = TRUE, chains = 4, cores =4,
                    control=list(adapt_delta=0.99))

# model check

pp_check(stemsNSC_bm3, ndraws=100)

summary(rhat(stemsNSC_bm3)) 

# Model evaluation (short version)

summary(stemsNSC_bm3,waic=TRUE, loo=TRUE)

#

save(stemsNSC_bm3, stemsNSC_bm2, gaus_v_log,
     file="Stems_NSC_model_dec2023.RData")