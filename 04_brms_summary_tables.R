require(knitr,quietly = TRUE)
require(tidyverse,quietly = TRUE)
require(tibble)
require(brms)

load(file="branches_NSC_model_dec2023.RData")
load(file="leaves_NSC_model_dec2023.RData")
load(file="Stems_NSC_model_dec2023.RData")
load(file="Roots_NSC_model_dec2023.RData")

# Roots
root_summ<-data.frame(summary(rootsNSC_bm3)$fixed)
root_summ<-rownames_to_column(root_summ,var="Fixed_effect")
root_summ<-select(root_summ, Fixed_effect,Estimate, l95=l.95..CI,u95= u.95..CI, ESS= Bulk_ESS)
root_summ<-root_summ%>%
  mutate(Estimate=round(Estimate, digits=2),
         l95=round(l95,digits=2),
         u95=round(u95,digits=2),
         Fixed_effect = fct_recode(Fixed_effect, "Intercept"="Intercept", "Biome LTF"="BiomeLTF","Biome UMF"="BiomeUMF",
                                                 "Leaf habit (Evergreen)"="Leaf.habitEvergreen",
                                                 "Fast-slow PC1"="FS1s", "Fast-slow PC2"="FS2s"))

root_summ$Estimate<-paste(root_summ$Estimate," (",root_summ$l95," - ",root_summ$u95,")",sep="")
root_summ<-select(root_summ, "Fixed effect"=Fixed_effect, Estimate, "Effective sample size"=ESS)

brms_root<-kable(root_summ, format = "markdown",caption = "Table S1. Summary of phylogenetic hierarchical models that examine variation in NSC concentrations of roots.", align="c", longtable=TRUE)

# Stems
stem_summ<-data.frame(summary(stemsNSC_bm3)$fixed)
stem_summ<-rownames_to_column(stem_summ,var="Fixed_effect")
stem_summ<-select(stem_summ, Fixed_effect,Estimate, l95=l.95..CI,u95= u.95..CI, ESS= Bulk_ESS)

stem_summ<-stem_summ%>%
  mutate(Estimate=round(Estimate, digits=2),
         l95=round(l95,digits=2),
         u95=round(u95,digits=2),
         Fixed_effect = fct_recode(Fixed_effect, "Intercept"="Intercept", "Biome LTF"="BiomeLTF","Biome UMF"="BiomeUMF",
                                   "Leaf habit (Evergreen)"="Leaf.habitEvergreen",
                                   "Fast-slow PC1"="FS1s", "Fast-slow PC2"="FS2s"))

stem_summ$Estimate<-paste(stem_summ$Estimate," (",stem_summ$l95," - ",stem_summ$u95,")",sep="")
stem_summ<-select(stem_summ, "Fixed effect"=Fixed_effect, Estimate, "Effective sample size"=ESS)

brms_stem<-kable(stem_summ, format = "markdown",caption = "Table S2. Summary of phylogenetic hierarchical models that examine variation in NSC concentrations of stems.", align="c", longtable=TRUE)

# Branches

branches_summ<-data.frame(summary(branchesNSC_bm3)$fixed)
branches_summ<-rownames_to_column(branches_summ,var="Fixed_effect")
branches_summ<-select(branches_summ, Fixed_effect,Estimate, l95=l.95..CI,u95= u.95..CI, ESS= Bulk_ESS)

branches_summ<-branches_summ%>%
  mutate(Estimate=round(Estimate, digits=2),
         l95=round(l95,digits=2),
         u95=round(u95,digits=2),
         Fixed_effect = fct_recode(Fixed_effect, "Intercept"="Intercept", "Biome LTF"="BiomeLTF","Biome UMF"="BiomeUMF",
                                   "Leaf habit (Evergreen)"="Leaf.habitEvergreen",
                                   "Fast-slow PC1"="FS1s", "Fast-slow PC2"="FS2s"))

branches_summ$Estimate<-paste(branches_summ$Estimate," (",branches_summ$l95," - ",branches_summ$u95,")",sep="")
branches_summ<-select(branches_summ, "Fixed effect"=Fixed_effect, Estimate, "Effective sample size"=ESS)

brms_branch<-kable(branches_summ, format = "markdown",caption = "Table S3. Summary of phylogenetic hierarchical models that examine variation in NSC concentrations of stems", align="c", longtable=TRUE)

# Leaves

leaves_summ<-data.frame(summary(leavesNSC_bm3)$fixed)
leaves_summ<-rownames_to_column(leaves_summ,var="Fixed_effect")
leaves_summ<-select(leaves_summ, Fixed_effect,Estimate, l95=l.95..CI,u95= u.95..CI, ESS= Bulk_ESS)
leaves_summ<-leaves_summ%>%
  mutate(Estimate=round(Estimate, digits=2),
         l95=round(l95,digits=2),
         u95=round(u95,digits=2),
         Fixed_effect = fct_recode(Fixed_effect, "Intercept"="Intercept", "Biome LTF"="BiomeLTF","Biome UMF"="BiomeUMF",
                                   "Leaf habit (Evergreen)"="Leaf.habitEvergreen",
                                   "Fast-slow PC1"="FS1s", "Fast-slow PC2"="FS2s",
                                   "LTF x Fast-slow PC1"="BiomeLTF:FS1s","UMF x Fast-slow PC1"="BiomeUMF:FS1s"))

leaves_summ$Estimate<-paste(leaves_summ$Estimate," (",leaves_summ$l95," - ",leaves_summ$u95,")",sep="")
leaves_summ<-select(leaves_summ, "Fixed effect"=Fixed_effect, Estimate, "Effective sample size"=ESS)

brms_leaf<-kable(leaves_summ, format = "markdown",caption = "Table S4. Summary of phylogenetic hierarchical models that examine variation in NSC concentrtions of leaves", align="c", longtable=TRUE)

save(brms_leaf, brms_branch, brms_stem, brms_root, file="Tables_brms_model_summary.RData")
