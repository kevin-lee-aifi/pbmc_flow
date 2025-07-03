################################################################################
################################################################################
###
### Date 05/15/2025
### Authors: Samir & Aishwarya
###
################################################################################
################################################################################
setwd('~/Google Drive/Shared drives/Imm - Projects/PID-00002 Oncology/PLAN-00051 FH1 Longitudinal MM (Landscape)/Samir MM/')
rm(list=ls())

### load libraries 
require(data.table)
require(ggplot2)
require(dplyr)
library(ggplot2)
library(ggalluvial)
library(dplyr)
library(RColorBrewer)
library(gridExtra)
library(rlang)
library(dplyr)
require(ggpubr)
source('scripts/compositional analyses/helper_functions.R')

### load data 
pbmc_files = fread('data/pbmc_l1_counts.csv')

### Remove Dara Patients
pbmc_files = pbmc_files[Dara=='non-dara']

### rename metadata for ease of reference
pbmc_files$subject =pbmc_files$subject.subjectGuid
pbmc_files$visit = pbmc_files$sample.visitDetails
pbmc_files$Cell = pbmc_files$tidy.aifi_l1

### rename factor names for analysis 
pbmc_files[visit=='MM Pre-Treatment']$visit <-'Pre-Trt'
pbmc_files[visit=='MM End Induction 1st Draw']$visit <-'End-Ind'
pbmc_files[visit=='MM Post Transplant 60 Days']$visit <-'ASCT 60d'
pbmc_files[visit=='MM Post Transplant 1 year']$visit <-'ASCT 1y'
pbmc_files[visit=='MM Post Transplant 2 year']$visit <-'ASCT 2y'

pbmc_files$visit = factor(pbmc_files$visit,
                          levels = c('Pre-Trt','End-Ind',
                                     'ASCT 60d', 'ASCT 1y','ASCT 2y'))

pbmc_files$cell_type_frac_total = pbmc_files$cell_type_frac

################################################################################

# Load required library
plotlist = lapply(unique(pbmc_files$Cell),
                  function(x){
                    cMat = pbmc_files[Cell==x & visit %in% c('Pre-Trt','End-Ind')]
                    
                    pairwise_data <- subset(cMat, visit %in% c('Pre-Trt','End-Ind'))
                    pairwise_data = dcast(pairwise_data, subject ~ visit, value.var = 'cell_type_frac_total')
                    pairwise_data = pairwise_data[ !is.na(pairwise_data[,2][[1]]) & 
                                                     !is.na(pairwise_data[,3][[1]])]   
                    pairwise_data <- melt(pairwise_data)
                    
                    ggplot(pairwise_data, 
                           aes(x = variable, y = value)) +
                      geom_boxplot(fill = "skyblue", color = "black") +geom_point(alpha=0.3)+
                      geom_line(data=pairwise_data,
                                aes(x=variable,
                                    y=value,
                                    group=subject),
                                alpha=0.3)+
                      labs(title = x,
                           x = "",
                           y = "") +
                      theme_minimal() +
                      theme(text=element_text(size=14),
                            axis.text.x = element_text(angle=90))+
                      ggpubr::stat_compare_means(paired = T)
                  }
)
fig=ggpubr::ggarrange(plotlist = plotlist, common.legend = TRUE,ncol = 5,nrow=2)
fig=annotate_figure(fig,
                    top=text_grob('Relative abundances in PBMCs', size=20))

annotate_figure(fig,
                #top=text_grob('Relative Abundances of Cells Pre- & post-VRD Therapy',size = 20),
                bottom=text_grob('Clinical Visits',size = 20),
                left = text_grob('% of Peripheral Blood',size = 20, rot = 90)
)

################################################################################
################################################################################

# Load required library
plotlist = lapply(unique(pbmc_files$Cell),
                  function(x){
                    print(x)
                    testCelltype_paired(pbmc_files, x, 
                                            pairsOfTimes=list(
                                              c('Pre-Trt','End-Ind'),
                                              c('End-Ind','ASCT 60d'),
                                              c('ASCT 60d','ASCT 1y')
                                            ),
                                            size=4)$Plot
                  }
)
fig=ggpubr::ggarrange(plotlist = plotlist, common.legend = TRUE,ncol = 5,nrow=2, 
                      label.x = unique(pbmc_files$visit))
fig=annotate_figure(fig,
                    top=text_grob('Relative abundances in PBMCs', size=20))

annotate_figure(fig,
                #top=text_grob('Relative Abundances of Cells Pre- & post-VRD Therapy',size = 20),
                bottom=text_grob('Clinical Visits',size = 20),
                left = text_grob('% of Peripheral Blood',size = 20, rot = 90)
)



# Load required library
res_list = lapply(unique(pbmc_files$Cell),
                  function(x){
                    print(x)
                    testCelltype_paired(pbmc_files, x, 
                                        pairsOfTimes=list(
                                          c('Pre-Trt','End-Ind'),
                                          c('End-Ind','ASCT 60d'),
                                          c('ASCT 60d','ASCT 1y')
                                        ),
                                        size=4)$Res
                  }
)

res_list <- rbindlist(res_list)
#res_list$Log2FC <- log2(res_list$Time1Val)-log2(res_list$Time2Val)

res_list$Time_comparison = paste(res_list$Time1, res_list$Time2,sep=' vs. ')
res_list$Time_comparison = factor(res_list$Time_comparison,
                                  levels=c('Pre-Trt vs. End-Ind',
                                           'End-Ind vs. ASCT 60d',
                                           'ASCT 60d vs. ASCT 1y'))

ggplot(res_list,
       aes(x=Log2FC,
           y=-log10(pval),
           label=Celltype))+geom_point()+
  ggrepel::geom_text_repel(size=6)+
  theme_linedraw()+facet_wrap(~Time_comparison)+
  theme(text=element_text(size=28))+
  ylab('-log10(P value)')+
  xlab('Log2 FC')+
  ggtitle('PBMC scRNA Compositions')
