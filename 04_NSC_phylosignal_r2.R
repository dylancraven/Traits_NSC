####################
# NSC: main models #
# 1) bayes r2      #
# 2) phylo. signal #
####################

require(brms)
require(knitr)
require(ggplot2)
require(tidyverse)

##########
# models #
##########

load(file="branches_NSC_model_dec2023.RData")
load(file="leaves_NSC_model_dec2023.RData")
load(file="Roots_NSC_model_dec2023.RData")
load(file="Stems_NSC_model_dec2023.RData")

############################
# 1. Phylogenetic signal   #
############################

# Roots
hyp1 <- paste(
  "sd_phylo__Intercept^2 /", 
  "(sd_phylo__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0")

(hyp1 <- hypothesis(rootsNSC_bm3, hyp1, class = NULL)) # test hypothesis

rootss<-c(hyp1$hypothesis$Estimate, hyp1$hypothesis$CI.Lower,hyp1$hypothesis$CI.Upper, "Roots")

# Stems

hyp2 <- paste(
  "sd_phylo__Intercept^2 /", 
  "(sd_phylo__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0")

(hyp2 <- hypothesis(stemsNSC_bm3, hyp2, class = NULL))
stemss<-c(hyp2$hypothesis$Estimate, hyp2$hypothesis$CI.Lower,hyp2$hypothesis$CI.Upper, "Stems")

# Branches

hyp3 <- paste(
  "sd_phylo__Intercept^2 /", 
  "(sd_phylo__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0")
(hyp3 <- hypothesis(branchesNSC_bm3, hyp3, class = NULL))
branches<-c(hyp3$hypothesis$Estimate, hyp3$hypothesis$CI.Lower,hyp3$hypothesis$CI.Upper, "Branches")

# Leaves

hyp4 <- paste(
  "sd_phylo__Intercept^2 /", 
  "(sd_phylo__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0")
(hyp4 <- hypothesis(leavesNSC_bm3, hyp4, class = NULL))
leaves<-c(hyp4$hypothesis$Estimate, hyp4$hypothesis$CI.Lower,hyp4$hypothesis$CI.Upper, "Leaves")

# combine results

phylo_sig<-rbind.data.frame(rootss, stemss,branches,leaves, stringsAsFactors = FALSE)
colnames(phylo_sig)[1]<-"lambda"
colnames(phylo_sig)[2]<-"l-95%"
colnames(phylo_sig)[3]<-"u-95%"
colnames(phylo_sig)[4]<-"Organ"
phylo_sig$lambda<-as.numeric(phylo_sig$lambda)
phylo_sig$'u-95%'<-as.numeric(phylo_sig$'u-95%')
phylo_sig$'l-95%'<-as.numeric(phylo_sig$'l-95%')

write.table(phylo_sig, "NSC_brms_phylosig.csv",sep=",",row.names = FALSE)



# phylo_sigg<-knitr::kable(phylo_sig, format="markdown",digits=c(2,2,2,NA), align="c",
#                          caption="Phylogenetic signal of hierarchical Bayesian models examining seasonal variation 
#                          in the relationship between NSC and functional traits")



###############
# 2. Bayes r2 #
###############

r_R2<-data.frame(bayes_R2(rootsNSC_bm3))
r_R2$Organ<-"Roots"

s_R2<-data.frame(bayes_R2(stemsNSC_bm3))
s_R2$Organ<-"Stems"

b_R2<-data.frame(bayes_R2(branchesNSC_bm3))
b_R2$Organ<-"Branches"

l_R2<-data.frame(bayes_R2(leavesNSC_bm3))
l_R2$Organ<-"Leaves"

R222<-rbind.data.frame(r_R2, s_R2, b_R2, l_R2)
R222<-dplyr::select(R222, Organ, "Bayesian R-sq."=Estimate, "l-95%"="Q2.5","u-95%"="Q97.5")
write.table(R222, "NSC_R2.csv",sep=",",row.names = FALSE)

##############
# r hat ######
##############

r_rhat<-data.frame(rhat=rhat(rootsNSC_bm3))
r_rhat$Organ<-"Roots"

s_rhat<-data.frame(rhat=rhat(stemsNSC_bm3))
s_rhat$Organ<-"Stems"

b_rhat<-data.frame(rhat=rhat(branchesNSC_bm3))
b_rhat$Organ<-"Branches"

l_rhat<-data.frame(rhat=rhat(leavesNSC_bm3))
l_rhat$Organ<-"Leaves"

all_rhat<-rbind.data.frame(r_rhat, s_rhat, b_rhat, l_rhat)

all_rhat<-  all_rhat %>% 
            group_by(Organ)%>%
            do(data.frame(rbind(Hmisc::smean.cl.normal(.$rhat))))

write.table(all_rhat, "NSC_rhat.csv",sep=",",row.names = FALSE)
