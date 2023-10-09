### figureInference-4.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Oct  9 2023 (11:38) 
## Version: 
## Last-Updated: Oct  9 2023 (12:17) 
##           By: Brice Ozenne
##     Update #: 11
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

## * generate figure 4
allResS.tempo <- readRDS("results/aggregated-power.rds")
convertion <- c("Ustat" = "Asymptotic",
                "Ustat-trans" = "Asymptotic with transformation",
                "boot-perc" = "Percentile bootstrap",
                "boot-stud" = "Studentized bootstrap",
                "perm-perc" = "Permutation",
                "perm-stud" = "Studentized permutation")


## ** type 1 error
figure4.A <- ggplot(allResS.tempo[mu==0], aes(x = n, y = power,
                                             group = method.legend, shape = method.legend, color = method.legend, linetype = method.legend))
figure4.A <- figure4.A + geom_hline(yintercept = 0.05, size = 2, color = "lightgray")
figure4.A <- figure4.A + geom_line(linewidth = 1.15) + geom_point(size = 3)
figure4.A <- figure4.A + facet_grid(mu.legend ~ statistic.legend)
figure4.A <- figure4.A + scale_linetype_manual(values = c(1,1,2,2,3,3), breaks = convertion) 
figure4.A <- figure4.A + scale_shape_manual(values = rep(c(8,15),3), breaks = convertion) 
figure4.A <- figure4.A + scale_color_manual(values = rep(c("darkgray","black"),3), breaks = convertion) 
figure4.A <- figure4.A + coord_cartesian(y = c(0.04,0.1))
figure4.A <- figure4.A + guides(color = guide_legend(nrow = 6, byrow = TRUE))
figure4.A <- figure4.A + labs(x = "Sample size (per group)", y = "Type 1 error", color = "Inferential method\n", shape = "Inferential method\n", linetype = "Inferential method\n")
figure4.A <- figure4.A + theme(text = element_text(size=15), 
                             axis.line = element_line(linewidth = 1.25),
                             axis.ticks = element_line(linewidth = 2),
                             axis.ticks.length=unit(.25, "cm"),
                             legend.text=element_text(size=11),
                             legend.key.width = unit(4,"line"))
figure4.A


## ** coverage
figure4.B <- ggplot(allResS.tempo[allResS.tempo$method.legend %in% convertion[1:4]],
                 aes(x = n, y = coverage,
                     group = method.legend, shape = method.legend, color = method.legend, linetype = method.legend))
figure4.B <- figure4.B + geom_hline(yintercept = 0.95, linewidth = 2, color = "lightgray")
figure4.B <- figure4.B + geom_line(linewidth = 1.15) + geom_point(size = 3)
figure4.B <- figure4.B + facet_grid(mu.legend ~ statistic.legend) 
figure4.B <- figure4.B + scale_linetype_manual(values = c(1,1,2,2), breaks = convertion[1:4]) 
figure4.B <- figure4.B + scale_shape_manual(values = c(8,15,8,15), breaks = convertion[1:4]) 
figure4.B <- figure4.B + scale_color_manual(values = rep(c("darkgray","black"),2), breaks = convertion[1:4]) 
figure4.B <- figure4.B + coord_cartesian(y = c(0.9,0.96))
figure4.B <- figure4.B + labs(x = "Sample size (per group)", y = "Coverage", color = "Inferential method\n", shape = "Inferential method\n", linetype = "Inferential method\n")
figure4.B <- figure4.B + guides(shape = guide_legend(nrow = 2, byrow = TRUE))
figure4.B <- figure4.B + guides(color = guide_legend(nrow = 2, byrow = TRUE))
figure4.B <- figure4.B + guides(linetype = guide_legend(nrow = 2, byrow = TRUE))
figure4.B <- figure4.B + theme(text = element_text(size=15), 
                         axis.line = element_line(linewidth = 1.25),
                         axis.ticks = element_line(linewidth = 2),
                         axis.ticks.length=unit(.25, "cm"),
                         legend.key.width = unit(4,"line"),
                         legend.position = "bottom")
figure4.B

## ** merge
figure4 <- ggarrange(figure4.A + guides(color = guide_legend(nrow = 3, byrow = TRUE), linetype = guide_legend(nrow = 3, byrow = TRUE), shape = guide_legend(nrow = 3, byrow = TRUE)) + ggtitle("A"),
                     figure4.B + ggtitle("B"), common.legend = TRUE, legend = "bottom", ncol = 1,
                     heights = c(3, 4, 4), nrow = 2, align = "v")
print(figure4)

ggsave(figure4, filename = "figures/fig_inference_simulation-type1.pdf", width = 9, height = 8,
       device = cairo_pdf)


##----------------------------------------------------------------------
### figureInference-4.R ends here
