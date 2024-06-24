### figureInference-6.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Oct  9 2023 (11:38) 
## Version: 
## Last-Updated: feb 17 2024 (14:54) 
##           By: Brice Ozenne
##     Update #: 27
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(ggplot2)
library(ggpubr)
library(data.table)
library(cowplot)
library(grid)
library(gridExtra) 
library(patchwork)                        

## * generate figure 6
allResS <- readRDS("results/aggregated-power.rds")
allResS[, approach.legend := factor(method, levels = c("Ustat1","Ustat2","Ustat1-trans","Ustat2-trans","boot-basic","boot-perc","boot-stud","perm-perc","perm-stud"),
                                    labels = c(rep("U-statistic",4),rep("Bootstrap",3),rep("Permutation",2)))]
allResS[, method.legend := factor(method, levels = c("Ustat1","Ustat2","Ustat1-trans","Ustat2-trans","boot-basic","boot-perc","boot-stud","perm-perc","perm-stud"),
                                  c("U-statistic1","U-statistic","U-statistic1 with transformation","U-statistic with transformation","Basic bootstrap","Percentile bootstrap","Studentized bootstrap","Permutation","Studentized permutation"))]
allResS.NTBc <- allResS[outcome == "continuous" & statistic == "netBenefit" & method.legend %in% c("U-statistic1","U-statistic1 with transformation")==FALSE]
allResS.WRc <- allResS[outcome == "continuous" & statistic == "winRatio" & method.legend %in% c("U-statistic1","U-statistic1 with transformation")==FALSE]

convertion <- function(txt){switch(txt,
    "U-statistic1" = "No transformation",
    "U-statistic" = "No transformation",
    "U-statistic1 with transformation" = "With transformation",
    "U-statistic with transformation" = "With transformation",
    "Basic bootstrap" = "Basic\n with transformation",
    "Percentile bootstrap" = "Percentile",
    "Studentized bootstrap" = "Studentized\n with transformation",
    "Permutation" = "Percentile",
    "Studentized permutation" = "Studentized\n with transformation")
}


## ** type 1 error
allResS.NTBc[mu==0 & n==10 & method == "boot-basic", 100*power]
## [1] 0.496

theme_set(theme_bw())
figure6.A1 <- ggplot(allResS.NTBc[mu==0],
                     aes(x = n, y = power, group = method.legend, shape = method.legend, linetype = method.legend))
figure6.A1 <- figure6.A1 + geom_hline(yintercept = 0.05, linewidth = 1.5, color = "darkgrey")
figure6.A1 <- figure6.A1 + geom_line(linewidth = 1) + geom_point(size = 3)
figure6.A1 <- figure6.A1 + scale_linetype_manual(values = c(1,3,1:3,1,3), labels = sapply(levels(allResS$method.legend),convertion)) 
figure6.A1 <- figure6.A1 + scale_shape_manual(values = c(16,8,c(16:17,8),16,8), labels = sapply(levels(allResS$method.legend),convertion)) 
figure6.A1 <- figure6.A1 + coord_cartesian(y = c(0.025,0.09), x = c(0,200))
figure6.A1 <- figure6.A1 + scale_y_continuous(labels = paste0(seq(3,10,by=2),"%"), breaks = seq(3,10,by=2)/100)
figure6.A1 <- figure6.A1 + scale_x_continuous(breaks = seq(0,200,by=50))
figure6.A1 <- figure6.A1 + guides(linetype = guide_legend(nrow = 4, byrow = TRUE))
figure6.A1 <- figure6.A1 + labs(x = "Sample size (per group)", y = "Type 1 error", shape = "", linetype = "")
figure6.A1 <- figure6.A1 + facet_grid( ~ approach.legend)
figure6.A1 <- figure6.A1 + guides(shape = guide_legend(override.aes = list(size = 3), nrow = 3, byrow = FALSE),
                                  linetype = guide_legend(nrow = 3, byrow = FALSE))
figure6.A1 <- figure6.A1 + theme(text = element_text(size=15), 
                                 axis.line = element_line(linewidth = 1.25),
                                 axis.ticks = element_line(linewidth = 2),
                                 axis.ticks.length=unit(.25, "cm"),
                                 legend.text=element_text(size=9),
                                 legend.position="bottom",
                                 legend.key.width = unit(4,"line"))
figure6.A1 <- figure6.A1 + ggtitle("A") + xlab("")
figure6.A1

## figure6.A1 %+% allResS[mu==0 & outcome == "categorical" & statistic == "netBenefit"]
## figure6.A1 %+% allResS[mu==0 & outcome == "continuous" & statistic == "winRatio"]

## ** coverage
figure6.A2 <- ggplot(allResS.NTBc[mu==1],
                     aes(x = n, y = coverage, group = method.legend, shape = method.legend, linetype = method.legend))
figure6.A2 <- figure6.A2 + geom_hline(yintercept = 0.95, linewidth = 1.5, color = "darkgrey")
figure6.A2 <- figure6.A2 + geom_line(linewidth = 1) + geom_point(size = 3)
figure6.A2 <- figure6.A2 + scale_linetype_manual(values = c(1,3,1:3,1,3)) 
figure6.A2 <- figure6.A2 + scale_shape_manual(values = c(16,8,c(16:17,8),16,8)) 
figure6.A2 <- figure6.A2 + coord_cartesian(y = c(0.9,1), x = c(0,200))
figure6.A2 <- figure6.A2 + scale_y_continuous(labels = paste0(seq(90,100,by=2.5),"%"), breaks = seq(90,100,by=2.5)/100)
figure6.A2 <- figure6.A2 + scale_x_continuous(breaks = seq(0,200,by=50))
figure6.A2 <- figure6.A2 + guides(linetype = guide_legend(nrow = 4, byrow = TRUE))
figure6.A2 <- figure6.A2 + labs(x = "Sample size (per group)", y = "Coverage", shape = "", linetype = "")
figure6.A2 <- figure6.A2 + facet_grid( ~ approach.legend)
figure6.A2 <- figure6.A2 + guides(shape = guide_legend(override.aes = list(size = 3), nrow = 3, byrow = FALSE),
                                  linetype = guide_legend(nrow = 3, byrow = FALSE))
figure6.A2 <- figure6.A2 + theme(text = element_text(size=15), 
                                 axis.line = element_line(linewidth = 1.25),
                                 axis.ticks = element_line(linewidth = 2),
                                 axis.ticks.length=unit(.25, "cm"),
                                 legend.text=element_text(size=9),
                                 legend.position="bottom",
                                 legend.key.width = unit(4,"line"))
figure6.A2 <- figure6.A2 + ggtitle("B")
figure6.A2

## figure6.A2 %+% allResS[mu==0 & outcome == "categorical" & statistic == "netBenefit"]
## figure6.A2 %+% allResS[mu==0 & outcome == "continuous" & statistic == "winRatio"]

## ** merge done manually via inkscape
legend.U <- get_legend(figure6.A1 %+% allResS.NTBc[mu==0 & method.legend %in% c("U-statistic","U-statistic with transformation")])
legend.Boot <- get_legend(figure6.A1 %+% allResS.NTBc[mu==0 & method.legend %in% c("Basic bootstrap","Percentile bootstrap","Studentized bootstrap")] + scale_linetype_manual(values = c(1:3), labels = sapply(levels(allResS$method.legend),convertion)) + scale_shape_manual(values = c(16:17,8), labels = sapply(levels(allResS$method.legend),convertion)))
legend.Perm <- get_legend(figure6.A1 %+% allResS.NTBc[mu==0 & method.legend %in% c("Permutation","Studentized permutation")])

design <- "AAA
           BBB
           CDE"
figure6.A1none <- figure6.A1 + guides(shape = "none", linetype = "none")
figure6.A2none <- figure6.A2 + guides(shape = "none", linetype = "none")

figure6 <- figure6.A1none + figure6.A2none + ggdraw(legend.U) + ggdraw(legend.Boot) + ggdraw(legend.Perm) + plot_layout(design = design, heights = c(1,1,0.5))
figure6

pdf(file = "figures/fig_inference_simulation-Type1coverage.pdf", width = 9, height = 8)
print(figure6)
dev.off()

svg(file = "figures/fig_inference_simulation-Type1coverage.svg", width = 9, height = 8)
print(figure6)
dev.off()

##
figure6bis.A1none <- figure6.A1 %+% allResS.WRc[mu==0] + guides(shape = "none", linetype = "none")
figure6bis.A2none <- figure6.A2 %+% allResS.WRc[mu==1] + guides(shape = "none", linetype = "none")

figure6bis <- figure6bis.A1none + figure6bis.A2none + ggdraw(legend.U) + ggdraw(legend.Boot) + ggdraw(legend.Perm) + plot_layout(design = design, heights = c(1,1,0.5))
figure6bis










##----------------------------------------------------------------------
### figureInference-6.R ends here

