### figureInference-8.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt  9 2023 (14:57) 
## Version: 
## Last-Updated: feb  9 2024 (15:11) 
##           By: Brice Ozenne
##     Update #: 8
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(ggplot2)

## * generate figure 6
allResS.tempo2 <- readRDS("results/aggregated-FWER.rds")
allResSL.tempo2 <- data.table::melt(allResS.tempo2[,.(n,mu,power,power.band,power.bonf)], id.vars = c("n","mu"))
allResSL.tempo2[, mu.legend := paste0("\u0394\u03bc=",mu)]
allResSL.tempo2[, variable.legend := factor(variable,
                                            levels = c("power","power.band","power.bonf"),
                                            labels = c("None","Dunnett","Bonferroni"))]

figure6 <- ggplot(allResSL.tempo2,aes(x = n, y = value))
figure6 <- figure6 + geom_hline(data = data.frame(x =  unique(allResSL.tempo2$n), y = 0.05, mu.legend = paste0("\u0394\u03bc=",0)), aes(yintercept = y) ,color = "gray", linewidth = 1.5)
figure6 <- figure6 + geom_line(aes(group = variable.legend, linetype = variable.legend), linewidth = 1.15) + geom_point(aes(group = variable, shape = variable.legend), size = 3)
figure6 <- figure6 + facet_wrap(~mu.legend)
figure6 <- figure6 + labs(x = "Sample size (per group)", y = "Rejection rate", linetype = "Multiplicity adjusment", shape = "Multiplicity adjusment")
figure6 <- figure6 + theme(text = element_text(size=15), 
                           axis.line = element_line(linewidth = 1.25),
                           axis.ticks = element_line(linewidth = 2),
                           axis.ticks.length=unit(.25, "cm"),
                           legend.text=element_text(size=11),
                           legend.key.width = unit(4,"line"),
                           legend.position = "bottom")
figure6

ggsave(figure6, filename = "figures/fig_inference_simulation-FWER.pdf", width = 9, height = 6,
       device = cairo_pdf)


##----------------------------------------------------------------------
### figureInference-8.R ends here
