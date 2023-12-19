#########################
# Model fit / assessment#
#########################

require(brms)
require(tidyverse)
require(ggplot2)

load(file="branches_NSC_model_dec2023.RData")
load(file="leaves_NSC_model_dec2023.RData")
load(file="Stems_NSC_model_dec2023.RData")
load(file="Roots_NSC_model_dec2023.RData")

# 1. Posterior predictive distribution

root_pp<-pp_check(rootsNSC_bm3, ndraws = 1000) +
          labs(y="Density",x="Root NSC")+
         theme(axis.text.x=element_text(face="bold",size=6),
               axis.text.y=element_text(face="bold",size=6),
               axis.title.y=element_text(face="bold",size=7),
               axis.title.x=element_text(face="bold",size=7),
               legend.text = element_text(face="bold",size=6),
               legend.position="none",
               panel.grid.major = element_blank(), panel.grid.minor = element_blank())

stem_pp<-pp_check(stemsNSC_bm3, ndraws = 1000) +
        labs(x="Stem NSC",y="Density")+
         scale_x_continuous(trans="log",breaks=c(0.1,1,5, 10, 20,60))+
         theme(axis.text.x=element_text(face="bold",size=6),
               axis.text.y=element_text(face="bold",size=6),
               axis.title.y=element_text(face="bold",size=7),
               axis.title.x=element_text(face="bold",size=7),
         legend.text = element_text(face="bold",size=6),
         legend.position="none",
         panel.grid.major = element_blank(), panel.grid.minor = element_blank())

branch_pp<-pp_check(branchesNSC_bm3, ndraws = 1000) +
          scale_x_continuous(trans="log",breaks=c(0.1,1,5, 10, 20,60))+
          scale_y_continuous(limit=c(0,1.2),breaks=c(0, 0.2, 0.4,0.6,0.8,1.0))+
          labs(x="Branch NSC",y="Density")+
           theme(axis.text.x=element_text(face="bold",size=6),
                 axis.text.y=element_text(face="bold",size=6),
                 axis.title.y=element_text(face="bold",size=7),
                 axis.title.x=element_text(face="bold",size=7),
           legend.text = element_text(face="bold",size=6),
           legend.position="none",
           panel.grid.major = element_blank(), panel.grid.minor = element_blank())

leaf_pp<-pp_check(leavesNSC_bm3, ndraws = 1000)+
  
         scale_y_continuous(limit=c(0,0.18),breaks=c(0, 0.04,0.08,0.12,.16))+
  
       labs(x="Leaf NSC",y="Density")+
         theme(axis.text.x=element_text(face="bold",size=6),
               axis.text.y=element_text(face="bold",size=6),
               axis.title.y=element_text(face="bold",size=7),
               axis.title.x=element_text(face="bold",size=7),
         legend.text = element_text(face="bold",size=6),
         legend.position="none",
         panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#### together

pp_all<-cowplot::plot_grid(root_pp, 
                  stem_pp+theme(axis.title.y=element_blank()), 
                  branch_pp, 
                  leaf_pp+theme(axis.title.y=element_blank()), 
                  labels=c("a)","b)","c)","d)"), label_size=6, align="hv")

ggsave(pp_all,filename="PP_check_mainfigs.png", 
    units="in", 
    width=7, 
    height=5, 
    pointsize=2, 
    dpi=500)