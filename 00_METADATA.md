## Table of Contents  
1. **Data**  
  + *NSC_trts_fastslow_original.RData*: This RData object includes a phylogenetic tree build for species included in our study (based on Smith & Brown 2018) and the _original_ trait and NSC data for 81 species.  
+ *phy* : phylogeny with 81 tree species.  
+ *trees_nscc*: data frame with trait and NSC data.  
_ID_: unique numeric identified of each tree in the data set.
_species_: scientific name of each tree.  
_phylo_: scientific name of each tree, with genus and species separated by an underscore.  
_Organ_: tree organ measured for NSC concentrations.  
_Concen_: NSC concentrations. 
_Biome_: Biome of sampled trees.  
_H_: Individual tree height.  
_Leaf.habit_: Leaf phenology (evergreen, deciduous).   
_leaf_area_: Leaf area.  
_LT_: Leaf thickness.  
_SLA_: Specific leaf area.    
_LDMC_: Leaf dry matter content.  
_Aarea_: Area-based bet photosynthetic capacity.  
_Amass_: Mass-based Net photosynethetic capacity.
_SD_: Wood density of stems.  
_BD_: Wood density of branches.
_LeafN_: Foliar N Content.  
_LeafP_: Foliar P Content.  
_LeafMg_: Foliar Mg Content.  
_LeafCa_: Foliar Ca Content.  
+ *NSC_trts_fastslow_clean.RData*: This RData object includes a phylogenetic tree build for species included in our study (based on Smith & Brown 2018) and the trait and NSC data for 61 species, including species scores along 'fast-slow' plant economics spectrum.    
+ _phy_ : phylogeny with 61 tree species.  
+ *trees_nsc_trts*: data frame with trait and NSC data.  
_ID_: unique numeric identified of each tree in the data set.
_Family_: species' Family.
      _Genus_: species' Genus.
_Species_: species' epithet.  
      _phylo_: scientific name of each tree, with genus and species separated by an underscore.  
      _Organ_: tree organ measured for NSC concentrations.  
      _Season_: Season when NSC concentrations were measured (S1 = early growing season; S2 = late growing season).  
      _Concen_: NSC concentrations. 
      _Biome_: Biome of sampled trees.  
      _H_: Individual tree height.  
      _Leaf.habit_: Leaf phenology (evergreen, deciduous).   
      _leaf_area_: Leaf area.  
      _LT_: Leaf thickness.  
      _SLA_: Specific leaf area.    
      _LDMC_: Leaf dry matter content.    
      _Amass_: Mass-based Net photosynethetic capacity.
      _SD_: Wood density of stems.  
      _BD_: Wood density of branches.
      _LeafN_: Foliar N Content.  
      _LeafP_: Foliar P Content.  
      _LeafMg_: Foliar Mg Content.  
      _LeafCa_: Foliar Ca Content.  
      _FS.1_: species' scores from PCA of 11 plant functional traits along the first dimension. 
_FS.2_: species' scores from PCA of 11 plant functional traits along the second dimension. 
2. **Code**    

+ 01_FullTraits_plusPCA.R: Code for PCA and calculating species' positions along fast-slow dimensions.  
+ 02_NSC_brms_roots_FSbiome.R: Code for fitting phylogenetic hierarchical Bayesian models for root NSC concentrations.   
+ 02_NSC_brms_stems_FSbiome.R:  Code for fitting phylogenetic hierarchical Bayesian models for stems NSC concentrations.   
+ 02_NSC_brms_branches_FSbiome.R:  Code for fitting phylogenetic hierarchical Bayesian models for branch NSC concentrations.   
+ 02_NSC_brms_leaves_FSbiome.R:  Code for fitting phylogenetic hierarchical Bayesian models for leaf NSC concentrations.   
+ 03_NSC_brms_roots_FSbiome_nophylo.R: Code for fitting non-phylogenetic hierarchical Bayesian models for root NSC concentrations.   
+ 03_NSC_brms_stems_FSbiome_nophylo.R:  Code for fitting non-phylogenetic hierarchical Bayesian models for stems NSC concentrations.   
+ 03_NSC_brms_branches_FSbiome_nophylo.R:  Code for fitting non-phylogenetic hierarchical Bayesian models for branch NSC concentrations.   
+ 03_NSC_brms_leaves_FSbiome_nophylo.R:  Code for fitting non-phylogenetic hierarchical Bayesian models for leaf NSC concentrations.   
+ 04_brms_summary_tables.R: Code for summarizing phylogenetic hierarchical Bayesian models.  
+ 04_brms_model_ppcheck.R: Code for Posterior (or prior) predictive checks.    
+ 04_NSC_phylosignal_r2.R: Code for estimating Bayesian R^2^, Rhat, and Pagel's lambda.    
