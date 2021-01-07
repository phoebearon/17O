
########################### Supplement 5. R script that calculates and plots oxygen and hydrogen isotope variation during mixing of water bodies ###########################
#
# Triple oxygen isotopes in the water cycle
# Author: Phoebe Aron
# Contact Phoebe Aron (paron@umich.edu) or Naomi Levin (nelevin@umich.edu) with questions and comments
# 
#
#   This scripts calculates d18O, dp18O (dp = delta prime), d17O, dp17O, d2H, d-excess, and Dp17O as two water bodies mix. Because Dp17O is defined 
#     with logarithmic dp notation, the mixing response in non linear. d18O and d-excess are defined with normal delta notation and mixing of these
#     results in a linear function of the fraction of each mixed water.
#
#   This script shows how to do the mixing calculation and generates plots of d18O, d-excess, and Dp17O as a function of the mixing fraction. 



## --------------------- User Setup --------------------- #


# required packages
library(ggplot2)
library(ggpubr)

# define isotope values of Water 1 (SMOW)
water1_d18O <- 0
water1_d17O <- 0
water1_d2H <- 0 

# define isotope values of Water 2 (SLAP)
water2_d18O <- -55.5
water2_d17O <- -29.6968
water2_d2H <- -428 

# define mixing fractions
mixing_f <- seq(0,1,0.1)

# define triple oxygen isotope reference slope 
lambda_ref_triple.oxygen <- 0.528 # Luz and Barkan, 2010

# define oxygen hydrogen reference slope
lambda_ref_oxygen.hydrogen <- 8 # Craig, 1961

## --------------------- Mix Waters --------------------- #

# mixed delta values
mix_d18O <- water1_d18O*mixing_f + water2_d18O*(1-mixing_f)
mix_d17O <- water1_d17O*mixing_f + water2_d17O*(1-mixing_f)
mix_d2H <- water1_d2H*mixing_f + water2_d2H*(1-mixing_f)

# convert to delta prime for Dp17O calculation
mix_dp18O <- log(mix_d18O/1000 + 1) * 1000
mix_dp17O <- log(mix_d17O/1000 + 1) * 1000
mix_dp2H <- log(mix_d2H/1000 + 1) * 1000

# mixed Dp17O
mix_Dp17O <- mix_dp17O - lambda_ref_triple.oxygen*mix_dp18O

# mixed d-excess 
mix.d_excess <- mix_d2H - lambda_ref_oxygen.hydrogen*mix_d18O

## --------------------- Plots --------------------- #

#define theme
supplement_5.theme <- theme(axis.text.x = element_text(margin = margin(t = 0.3, unit = "cm")),
                                axis.text.y = element_text(margin = margin(r = 0.3, unit = "cm")),
                                axis.ticks.length=unit(-0.15, "cm"),
                                axis.ticks = element_line(colour = "black"),
                                text = element_text(size = 10),
                                axis.title = element_text(size = 10))
# d18O
mix_d18O_plot <- ggplot() + 
  geom_point(data=data.frame(mixing_f,mix_d18O), aes(x=mixing_f, y=mix_d18O)) +
  theme_bw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) + supplement_5.theme +
  labs(x="Fraction of Water 1", y=expression("Mixture "*delta^"18"*"O (\u2030)")) +
  ggtitle(expression(delta^"18"*"O Mixing Response"))
plot(mix_d18O_plot)

# d-excess
mix_dxs_plot <- ggplot() + 
  geom_point(data=data.frame(mixing_f,mix.d_excess), aes(x=mixing_f, y=mix.d_excess)) +
  theme_bw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) + supplement_5.theme +
  labs(x="Fraction of Water 1", y="Mixture d-excess (\u2030)") +
  ggtitle("d-excess Mixing Response")
plot(mix_dxs_plot)

# Dp17O
mix_Dp17O_plot <- ggplot() + 
  geom_point(data=data.frame(mixing_f,mix_Dp17O), aes(x=mixing_f, y=mix_Dp17O*1000)) +
  theme_bw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) + supplement_5.theme +
  labs(x="Fraction of Water 1", y=expression(Delta*"'"^"17"*"O (per meg)")) +
  ggtitle(expression(Delta*"'"^"17"*"O Mixing Response"))
plot(mix_Dp17O_plot)


## compile plots

mixing_plots <- ggarrange(mix_d18O_plot,mix_dxs_plot,mix_Dp17O_plot,nrow=1,ncol=3,align="hv")
plot(mixing_plots)

