### figureInference-5.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Oct  9 2023 (12:30) 
## Version: 
## Last-Updated: Oct  9 2023 (12:35) 
##           By: Brice Ozenne
##     Update #: 5
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(ggplot2)
library(ggforce)
library(data.table)
library(ggpubr)
library(mvtnorm)
library(ggpattern) ## remotes::install_github("coolbutuseless/ggpattern")
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

## * generate figure 5

## ** rejection boundaries
grid1D <- seq(-3,3,0.025)
rejection1D <- c(qnorm(0.025), qnorm(0.975))

rejection2D.2uni <- data.table(xmax = qnorm(0.975, mean = 0, sd = 1),
                               ymax = qnorm(0.975, mean = 0, sd = 1),
                               xmin = qnorm(0.025, mean = 0, sd = 1),
                               ymin = qnorm(0.025, mean = 0, sd = 1))

qq <- qmvnorm(0.95, mean = c(0,0), sigma = diag(1,2), tail = "both")
rejection2D.2uniadj <- data.table(
  xmax = qq$quantile, ymax = qq$quantile,
  xmin = -qq$quantile, ymin = -qq$quantile
)

rejection2D.chisq <- sqrt(qchisq(0.95, df = 2))

## ** Multivariate Wald test
mismatch <- sqrt(rejection2D.chisq^2 - qq$quantile^2)
mismatch2 <- -(0.7*qq$quantile+1.3*rejection2D.chisq)/2
dtFisher.areaTrunCi <- data.table(x = c(-5, -5, -qq$quantile, -qq$quantile, mismatch2, -rejection2D.chisq, mismatch2, -qq$quantile, -qq$quantile),
                                  y = c(-qq$quantile,  qq$quantile, qq$quantile,     mismatch, mismatch/2,                 0, -mismatch/2,   -mismatch, -qq$quantile))
dtFisher.areaTrunCiC <- data.table(x = c(-qq$quantile, mismatch2, -rejection2D.chisq, mismatch2, -qq$quantile),
                                   y = c(mismatch, mismatch/2,                 0, -mismatch/2,   -mismatch))
dtFisher.areaTrunCiC2 <- data.table(x = c(-qq$quantile, -qq$quantile, -mismatch, -(0.4*qq$quantile+0.5*mismatch), -(0.6*qq$quantile+0.5*mismatch)),
                                    y = c(-mismatch, -qq$quantile, -qq$quantile, -(0.5*qq$quantile+0.9*mismatch), -(qq$quantile+mismatch)/2))
dtFisher.areaRect <- data.table(xmin = -5, xmax = -qq$quantile, ymin = -5, ymax = -qq$quantile)
dtFisher.label12 <- data.frame(x = c(-3.35,3.35,-3.35,3.35), y = c(-3.25,-3.25,3.25,3.25), label = as.character(expression(reject~H[0]^1~and~H[0]^2)))
dtFisher.label1 <- data.frame(x = c(-3.35,3.35,0,0), y = c(0,0,3.25,-3.25),
                              label = c(rep(as.character(expression("reject H"[0]^1)),2),rep(as.character(expression("reject H"[0]^2)),2)))


figure5.A <- ggplot() + labs(x=expression(hat(Delta)[1]/sigma[hat(Delta)[1]]), y=expression(hat(Delta)[2]/sigma[hat(Delta)[2]]))
figure5.A <- figure5.A + geom_text(data = data.frame(x=0,y=0,label="No rejection"), aes(x=x,y=y,label=label), size = 3.5)
figure5.A <- figure5.A + geom_text(data = dtFisher.label12, aes(x=x,y=y,label=label), size = 3.5, parse = TRUE)
figure5.A <- figure5.A + geom_text(data = dtFisher.label1, aes(x=x,y=y,label=label), size = 3.5, parse = TRUE)
figure5.A <- figure5.A + geom_rect(data = dtFisher.areaRect, aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax), fill = "red", alpha = 0.4)
figure5.A <- figure5.A + geom_rect(data = dtFisher.areaRect, aes(xmin = xmin, ymin = -ymin, xmax = xmax, ymax = -ymax), fill = "red", alpha = 0.4)
figure5.A <- figure5.A + geom_rect(data = dtFisher.areaRect, aes(xmin = -xmin, ymin = ymin, xmax = -xmax, ymax = ymax), fill = "red", alpha = 0.4)
figure5.A <- figure5.A + geom_rect(data = dtFisher.areaRect, aes(xmin = -xmin, ymin = -ymin, xmax = -xmax, ymax = -ymax), fill = "red", alpha = 0.4)
figure5.A <- figure5.A + geom_polygon_pattern(data = dtFisher.areaTrunCiC2, aes(x=x,y=y), fill = NA, pattern_density = 0.1, pattern_spacing = 0.01, pattern = "circle")
figure5.A <- figure5.A + geom_polygon_pattern(data = dtFisher.areaTrunCiC2, aes(x=x,y=-y), fill = NA, pattern_density = 0.1, pattern_spacing = 0.01, pattern = "circle")
figure5.A <- figure5.A + geom_polygon_pattern(data = dtFisher.areaTrunCiC2, aes(x=-x,y=y), fill = NA, pattern_density = 0.1, pattern_spacing = 0.01, pattern = "circle")
figure5.A <- figure5.A + geom_polygon_pattern(data = dtFisher.areaTrunCiC2, aes(x=-x,y=-y), fill = NA, pattern_density = 0.1, pattern_spacing = 0.01, pattern = "circle")
figure5.A <- figure5.A + geom_polygon_pattern(data = dtFisher.areaTrunCiC, aes(x=x,y=y), fill = NA, pattern_density = 0.1, pattern_spacing = 0.01, pattern = "stripe")
figure5.A <- figure5.A + geom_polygon_pattern(data = dtFisher.areaTrunCiC, aes(x=-x,y=y), fill = NA, pattern_density = 0.1, pattern_spacing = 0.01, pattern = "stripe")
figure5.A <- figure5.A + geom_polygon_pattern(data = dtFisher.areaTrunCiC, aes(x=y,y=x), fill = NA, pattern_density = 0.1, pattern_spacing = 0.01, pattern = "stripe")
figure5.A <- figure5.A + geom_polygon_pattern(data = dtFisher.areaTrunCiC, aes(x=y,y=-x), fill = NA, pattern_density = 0.1, pattern_spacing = 0.01, pattern = "stripe")
figure5.A <- figure5.A + geom_polygon(data = dtFisher.areaTrunCi, aes(x=x,y=y), fill = "red", alpha = 0.2)
figure5.A <- figure5.A + geom_polygon(data = dtFisher.areaTrunCi, aes(x=-x,y=y), fill = "red", alpha = 0.2)
figure5.A <- figure5.A + geom_polygon(data = dtFisher.areaTrunCi, aes(x=y,y=x), fill = "red", alpha = 0.2)
figure5.A <- figure5.A + geom_polygon(data = dtFisher.areaTrunCi, aes(x=y,y=-x), fill = "red", alpha = 0.2)
figure5.A <- figure5.A + geom_rect(data = rejection2D.2uni,
                                 aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                                     color = "Gaussian (\u03B1=0.05)", linetype = "Gaussian (\u03B1=0.05)"),
                                 size = 2,
                                 fill = NA)
figure5.A <- figure5.A + geom_rect(data = rejection2D.2uniadj,
                                 aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                                     color = "Gaussian (\u03B1=0.025)", linetype = "Gaussian (\u03B1=0.025)"),
                                 size = 2,
                                 fill = NA)
figure5.A <- figure5.A + geom_circle(aes(x0=0, y0=0, r = rejection2D.chisq,
                                       color = "Chi-squared (df=2)", linetype = "Chi-squared (df=2)"),
                                   size = 2)
figure5.A <- figure5.A + scale_linetype_manual(values = c(1,2,3))
figure5.A <- figure5.A + labs(linetype = "critical quantile", color = "critical quantile")
figure5.A <- figure5.A + theme_light()
figure5.A <- figure5.A + theme(text = element_text(size=14), 
                             axis.line = element_line(linewidth = 1.25),
                             axis.ticks = element_line(linewidth = 2),
                             axis.ticks.length=unit(.25, "cm"),
                             legend.key.size = unit(1.5,"line"))
figure5.A <- figure5.A + coord_cartesian(xlim = c(-4,4), ylim = c(-4,4))
figure5.A <- figure5.A + ggtitle("(A) Multivariate Wald test")
figure5.A

## ** Sequential testing
rejection2D.gate <- data.table(
  xmax = c(5,qnorm(0.025, mean = 0, sd = 1)), ymax = qnorm(0.975, mean = 0, sd = 1),
  xmin = c(qnorm(0.975, mean = 0, sd = 1), -5), ymin = qnorm(0.025, mean = 0, sd = 1)
)
rejection2D.gate2 <- rbind(data.table(
  ymax = c(5,qnorm(0.025, mean = 0, sd = 1)), xmax = qnorm(0.025, mean = 0, sd = 1),
  ymin = c(qnorm(0.975, mean = 0, sd = 1), -5), xmin = -5
),
data.table(
  ymax = c(5,qnorm(0.025, mean = 0, sd = 1)), xmax = 5,
  ymin = c(qnorm(0.975, mean = 0, sd = 1), -5), xmin = qnorm(0.975, mean = 0, sd = 1)
))
dtGate.label1 <- data.frame(x = c(-3,3),
                            y = 0,
                            label = as.character(expression("reject H"[0]^1)))
dtGate.label12 <- data.frame(x = c(-3,3,-3,3),
                             y = c(-3,-3,3,3),
                             label = as.character(expression(reject~H[0]^1~and~H[0]^2)))

figure5.B <- ggplot() + labs(x=expression(hat(Delta)[1]/sigma[hat(Delta)[1]]), y=expression(hat(Delta)[2]/sigma[hat(Delta)[2]]))
figure5.B <- figure5.B + geom_vline(data = data.table(xintercept = c(qnorm(0.025, mean = 0, sd = 1),qnorm(0.975, mean = 0, sd = 1))),
                              aes(xintercept = xintercept), linewidth = 2, linetype = 3, color = gg_color_hue(3)[3])
figure5.B <- figure5.B + geom_rect(data = rejection2D.gate,
                             aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                             size = 2,
                             fill = "red", alpha = 0.2)
figure5.B <- figure5.B + geom_rect(data = rejection2D.gate2,
                             aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                             size = 2,
                             fill = "red", alpha = 0.4)
figure5.B <- figure5.B + geom_text(data = data.frame(x=0,y=0,label="No rejection"), aes(x=x,y=y,label=label), size = 3.5)
figure5.B <- figure5.B + geom_text(data = dtGate.label1, aes(x=x,y=y,label=label), size = 3.5, parse = TRUE)
figure5.B <- figure5.B + geom_text(data = dtGate.label12, aes(x=x,y=y,label=label), size = 3.5, parse = TRUE)
figure5.B <- figure5.B + coord_cartesian(xlim = c(-4,4), ylim = c(-4,4))
figure5.B <- figure5.B + theme_light()
figure5.B <- figure5.B + ggtitle("(B) Sequential testing")
figure5.B <- figure5.B + theme(text = element_text(size=14), 
                         axis.line = element_line(linewidth = 1.25),
                         axis.ticks = element_line(linewidth = 2),
                         axis.ticks.length=unit(.25, "cm"),
                         legend.key.size = unit(1.5,"line"))
figure5.B

## ** Selective testing
dtMax.areaRect1 <- data.table(xmin = -qq$quantile, xmax = qq$quantile, ymin = -5, ymax = -qq$quantile)
dtMax.areaRect2 <- data.table(xmin = -5, xmax = -qq$quantile, ymin = -5, ymax = -qq$quantile)
dtMax.label12 <- data.frame(x = c(-3.35,3.35,-3.35,3.35), y = c(-3.25,-3.25,3.25,3.25), label = as.character(expression(reject~H[0]^1~and~H[0]^2)))
dtMax.label1 <- data.frame(x = c(-3.35,3.35,0,0), y = c(0,0,3.25,-3.25),
                           label = c(rep(as.character(expression("reject H"[0]^1)),2),rep(as.character(expression("reject H"[0]^2)),2)))
figure5.C <- ggplot() + labs(x=expression(hat(Delta)[1]/sigma[hat(Delta)[1]]), y=expression(hat(Delta)[2]/sigma[hat(Delta)[2]]))
figure5.C <- figure5.C + geom_rect(data = data.table(xmin = -qq$quantile, xmax = qq$quantile, ymin = -qq$quantile, ymax = qq$quantile),
                           aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax), fill=NA, size = 2, linetype = 2, color = gg_color_hue(3)[2])
figure5.C <- figure5.C + geom_rect(data = dtMax.areaRect1, aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax), fill = "red", alpha = 0.2)
figure5.C <- figure5.C + geom_rect(data = dtMax.areaRect1, aes(xmin = xmin, ymin = -ymin, xmax = xmax, ymax = -ymax), fill = "red", alpha = 0.2)
figure5.C <- figure5.C + geom_rect(data = dtMax.areaRect1, aes(xmin = ymin, ymin = xmin, xmax = ymax, ymax = xmax), fill = "red", alpha = 0.2)
figure5.C <- figure5.C + geom_rect(data = dtMax.areaRect1, aes(xmin = -ymin, ymin = xmin, xmax = -ymax, ymax = xmax), fill = "red", alpha = 0.2)
figure5.C <- figure5.C + geom_rect(data = dtMax.areaRect2, aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax), fill = "red", alpha = 0.4)
figure5.C <- figure5.C + geom_rect(data = dtMax.areaRect2, aes(xmin = xmin, ymin = -ymin, xmax = xmax, ymax = -ymax), fill = "red", alpha = 0.4)
figure5.C <- figure5.C + geom_rect(data = dtMax.areaRect2, aes(xmin = -xmin, ymin = ymin, xmax = -xmax, ymax = ymax), fill = "red", alpha = 0.4)
figure5.C <- figure5.C + geom_rect(data = dtMax.areaRect2, aes(xmin = -xmin, ymin = -ymin, xmax = -xmax, ymax = -ymax), fill = "red", alpha = 0.4)
figure5.C <- figure5.C + geom_text(data = data.frame(x=0,y=0,label="No rejection"), aes(x=x,y=y,label=label), size = 3.5)
figure5.C <- figure5.C + geom_text(data = dtMax.label12, aes(x=x,y=y,label=label), size = 3.5, parse = TRUE)
figure5.C <- figure5.C + geom_text(data = dtMax.label1, aes(x=x,y=y,label=label), size = 3.5, parse = TRUE)
figure5.C <- figure5.C + labs(linetype = "critical quantile", color = "critical quantile")
figure5.C <- figure5.C + theme_light()
figure5.C <- figure5.C + theme(text = element_text(size=14), 
                       axis.line = element_line(linewidth = 1.25),
                       axis.ticks = element_line(linewidth = 2),
                       axis.ticks.length=unit(.25, "cm"),
                       legend.key.size = unit(1.5,"line"))
figure5.C <- figure5.C + coord_cartesian(xlim = c(-4,4), ylim = c(-4,4))
figure5.C <- figure5.C + ggtitle("(C) Selective testing")
figure5.C

## ** Merge
pdf("figures/fig_inference_consonanceCoherence.pdf", width = 12, height = 10)
ggarrange(ggarrange(figure5.A + guides(color="none",linetype="none"),
                    as_ggplot(cowplot::get_legend(figure5.A)),ncol=2),
          ggarrange(figure5.B,figure5.C,ncol=2),nrow=2)
dev.off()

##----------------------------------------------------------------------
### figureInference-5.R ends here
