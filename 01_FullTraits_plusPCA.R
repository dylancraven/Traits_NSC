###########################
#  pca of all traits      #
#  complete data set      #
#  no imputation          #
###########################

# Dic. 2023

# require(reshape2)
# require(stringr)
require(tidyverse)
# require(ggplot2)
require(FactoMineR)
require(phytools)
require(ggrepel)

########
# Data #
########

load(file="NSC_trts_fastslow_original.RData")

trees_trts<-trees_nscc%>%
            ungroup(.)%>%
            select(.,-species,-Organ,-Concen)%>%
            distinct(.)

trees_trts$Biome<-droplevels(trees_trts$Biome)
trees_trts$ID<-droplevels(trees_trts$ID)

#---------------------------#
# create full trait data    #
# --------------------------#

trts_full<-dplyr::select(trees_trts, -Aarea) # or remove Amax
trts_full<-data.frame(trts_full)

trts_full<-trts_full %>% 
          complete(.) %>% 
          drop_na(.)

trts_full$Biome<-droplevels(trts_full$Biome)
trts_full$Leaf.habit<-droplevels(trts_full$Leaf.habit)
trts_full$ID<-droplevels(trts_full$ID)
trts_full$ID<-as.character(trts_full$ID)

par(mfrow=c(4, 3),mar=c(1, 1, 1, 1))
qqnorm(trts_full$Amass)
qqnorm(trts_full$SLA)
qqnorm(trts_full$LT)
qqnorm(trts_full$leaf_area)
qqnorm(trts_full$LeafCa)
qqnorm(trts_full$LeafMg)
qqnorm(trts_full$LeafP)
qqnorm(trts_full$LeafN)
qqnorm(trts_full$SD)
qqnorm(trts_full$BD)

# select trts and make matrix

trts_spp_m<-trts_full %>% 
            mutate_at(c("LDMC", "SD", "BD","LeafN","LeafP","LeafMg"), ~scale(., center = TRUE, scale=TRUE)) %>% 
            mutate_at(c("SLA","Amass","leaf_area","LeafCa","LT"), ~scale(log(.), center = TRUE, scale=TRUE)) %>% 
            column_to_rownames(.,var="ID") %>% 
            dplyr::select(., leaf_area, LT, SLA, LeafCa, Amass, LDMC, SD, BD, LeafN, LeafP, LeafMg )

trts_spp_mm<-as.matrix(trts_spp_m)

#############
#  pca      #
#############

nsc_pca<-PCA(trts_spp_mm, scale.unit = T, ncp=5)


# Variables (the 10 first)
#   Dim.1    ctr   cos2                  Dim.2    ctr   cos2    Dim.3    ctr   cos2  
#   leaf_area |  0.305  2.468  0.093 | -0.383  5.768  0.147 | -0.111  0.921  0.012 |
#   LT        | -0.476  6.014  0.227 | -0.717 20.209  0.514 | -0.352  9.340  0.124 |
#   SLA       |  0.712 13.436  0.507 |  0.601 14.210  0.362 |  0.130  1.280  0.017 |
#   LeafCa    |  0.347  3.185  0.120 | -0.497  9.699  0.247 |  0.661 32.860  0.437 |
#   Amass     |  0.739 14.485  0.546 |  0.516 10.466  0.266 |  0.172  2.231  0.030 |
#   LDMC      | -0.569  8.596  0.324 |  0.466  8.542  0.217 |  0.460 15.883  0.211 |
#   SD        | -0.615 10.039  0.379 |  0.308  3.738  0.095 | -0.137  1.410  0.019 |
#   BD        | -0.636 10.732  0.405 |  0.505 10.012  0.255 |  0.003  0.001  0.000 |
#   LeafN     |  0.783 16.241  0.612 |  0.106  0.438  0.011 | -0.215  3.464  0.046 |
#   LeafP     |  0.729 14.096  0.532 | -0.085  0.283  0.007 | -0.434 14.186  0.189 |


spp_scrs<-data.frame(nsc_pca$ind$coord[,1:2], ID=rownames(nsc_pca$ind$coord))
trt_loadings<-sweep(nsc_pca$var$coord,2,sqrt(nsc_pca$eig[1:ncol(nsc_pca$var$coord),1]),FUN="/")

# add biome to spp scores

spp_biome<- distinct(select(ungroup(trts_full),ID, Biome, Leaf.habit))
spp_scrs<-left_join(spp_scrs, spp_biome, by="ID")
spp_scrs$Biome<-droplevels(spp_scrs$Biome)

spp_scrs$Biome<-as.factor(spp_scrs$Biome)
spp_scrs$Biome<-factor(spp_scrs$Biome,levels=c("LTF","UMF","DTF"))
spp_scrs$Leaf.habit<-as.factor(spp_scrs$Leaf.habit)
spp_scrs$Leaf.habit<-factor(spp_scrs$Leaf.habit,levels=c("Deciduous","Evergreen"))

trts_scrs<-data.frame(trt_loadings[,1:2],Traits=rownames(trt_loadings))
trts_scrs$Traits<-str_remove(trts_scrs$Traits, "value.cons_")

trts_scrs$Traits<-ifelse(trts_scrs$Traits=="leaf_area","Leaf area",trts_scrs$Traits)
trts_scrs$Traits<-ifelse(trts_scrs$Traits=="Amass","Amass",trts_scrs$Traits)
trts_scrs$Traits<-ifelse(trts_scrs$Traits=="LeafMg","Leaf Mg",trts_scrs$Traits)
trts_scrs$Traits<-ifelse(trts_scrs$Traits=="LeafCa","Leaf Ca",trts_scrs$Traits)
trts_scrs$Traits<-ifelse(trts_scrs$Traits=="LeafP","Leaf P",trts_scrs$Traits)
trts_scrs$Traits<-ifelse(trts_scrs$Traits=="LeafN","Leaf N",trts_scrs$Traits)


# make pretty PCA

normal_pca <- ggplot(spp_scrs) +
              geom_vline(xintercept=0, color="black")+
              geom_hline(yintercept=0, color="black")+
              geom_point(mapping = aes(x = Dim.1, y = Dim.2, group=Biome, colour=Biome, shape=Leaf.habit),
                         alpha=0.7) +
              coord_fixed() + ## need aspect ratio of 1!
              scale_colour_manual(name="Biome",values=c("LTF"="#ca0020","UMF"="#f4a582","DTF"="#0571b0"))+
              scale_shape_manual(name="Leaf habit",values=c("Evergreen"=21,"Deciduous"=24))+
              
              geom_segment(data = trts_scrs,
                           aes(x = 0, xend = Dim.1*10, y = 0, yend = Dim.2*10),
                           arrow = arrow(length = unit(0.5, "cm")),colour = "grey") +
              
              geom_text_repel(data = trts_scrs, aes(x = Dim.1*10, y = Dim.2*10, label = Traits),
                             size = 3, hjust=0.9)+
              #vjust=c(-0.2,-0.2,-0.2,1,-0.2,-0.2,-0.2))+
              
              xlab("PC1 (34.3%)") + ylab("PC2 (23.1%)")+
              theme_bw()+
            theme( axis.title.x=element_text(colour="black",face="bold",size=7),
                   axis.title.y=element_text(colour="black",face="bold",size=7),
                   axis.text.x=element_text(colour=c("black"),face="bold",size=6),
                   axis.text.y=element_text(colour=c("black"),face="bold",size=6),
                   text=element_text(colour="black",face="bold",size=7),
                   legend.position="right", legend.direction="vertical",
                   legend.key = element_rect(fill="transparent"),
                   legend.key.size =  unit(1,"line"),
                   legend.title=element_text(size=6, color="black",face="bold"),legend.text=element_text(size=6, color="black"),
                   panel.background =element_rect(fill="transparent",colour="black"),
                   panel.grid.minor=element_blank())

ggsave(normal_pca, filename="PCA_Allspp_11trts_nophylo.png", 
        # type="cairo",
        units="cm", 
        width=25, 
        height=25, 
        dpi=1000)

#--------------------------#
# data for analysis -------#
#--------------------------#

nsc_out<-spp_scrs%>%
          select(., ID, FS.1=Dim.1, FS.2=Dim.2 )%>%
          left_join(trts_full,., by="ID")

phy<-drop.tip(phy, setdiff(phy$tip.label, trees_nsc_trts$phylo))

save(trees_nsc_trts,phy, file="NSC_trts_fastslow_clean.RData") # 