### figureInference-6.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Oct  9 2023 (11:38) 
## Version: 
## Last-Updated: feb  9 2024 (19:35) 
##           By: Brice Ozenne
##     Update #: 23
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
allResS.tempo <- readRDS("results/aggregated-power.rds")[outcome == "continuous" & statistic == "netBenefit"]
allResS.tempo[approach.legend == "Asymptotic", approach.legend := "U-statistic"]
allResS.tempo[, approach.legend := factor(approach.legend, levels = c("U-statistic","Bootstrap","Permutation"))]
allResS.tempo[, method.legend := gsub("Asymptotic","U-statistic",method.legend)]
allResS.tempo[, method.legend := factor(method.legend,
                                        c("U-statistic","U-statistic with transformation","Basic bootstrap","Percentile bootstrap","Studentized bootstrap","Permutation","Studentized permutation"))]
convertion <- function(txt){switch(txt,
    "U-statistic" = "No transformation",
    "U-statistic with transformation" = "With transformation",
    "Basic bootstrap" = "Basic\n with transformation",
    "Percentile bootstrap" = "Percentile",
    "Studentized bootstrap" = "Studentized\n with transformation",
    "Permutation" = "Percentile",
    "Studentized permutation" = "Studentized\n with transformation")
}


## ** type 1 error
figure6.A1 <- ggplot(allResS.tempo[mu==0],
                     aes(x = n, y = power, group = method.legend, shape = method.legend, linetype = method.legend))
figure6.A1 <- figure6.A1 + geom_hline(yintercept = 0.05, linewidth = 1.5, color = "darkgrey")
figure6.A1 <- figure6.A1 + geom_line(linewidth = 1) + geom_point(size = 3)
figure6.A1 <- figure6.A1 + scale_linetype_manual(values = c(1,3,1:3,1,3), labels = sapply(levels(allResS.tempo$method.legend),convertion)) 
figure6.A1 <- figure6.A1 + scale_shape_manual(values = c(16,8,c(16:17,8),16,8), labels = sapply(levels(allResS.tempo$method.legend),convertion)) 
figure6.A1 <- figure6.A1 + coord_cartesian(y = c(0.03,0.09), x = c(0,200))
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

## ** coverage
figure6.A2 <- ggplot(allResS.tempo[mu==1],
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

## ** merge done manually via inkscape
legend.U <- get_legend(figure6.A1 %+% allResS.tempo[mu==0 & method.legend %in% c("U-statistic","U-statistic with transformation")])
legend.Boot <- get_legend(figure6.A1 %+% allResS.tempo[mu==0 & method.legend %in% c("Basic bootstrap","Percentile bootstrap","Studentized bootstrap")] + scale_linetype_manual(values = c(1:3), labels = sapply(levels(allResS.tempo$method.legend),convertion)) + scale_shape_manual(values = c(16:17,8), labels = sapply(levels(allResS.tempo$method.legend),convertion)))
legend.Perm <- get_legend(figure6.A1 %+% allResS.tempo[mu==0 & method.legend %in% c("Permutation","Studentized permutation")])

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
##----------------------------------------------------------------------
### figureInference-6.R ends here

