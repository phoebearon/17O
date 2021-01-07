
########################### Supplement 6. R script that calculates and plots oxygen and hydrogen isotope fractionation during pan evaporation ###########################
# 
# Triple oxygen isotopes in the water cycle
# Author: Phoebe Aron
# Contact Phoebe Aron (paron@umich.edu) or Naomi Levin (nelevin@umich.edu) with questions and comments
# 

#   This script calculates d18O, d17O, Dp17O, d2H, d-excess of water in isolated water bodies (pan evaporation). 
#   These calculations are based on Equations 11-13 in Aron et al. 2020.

#   This script generates two plots, dp18O vs. Dp17O and d18O vs d-excess, during pan evaporation.


## --------------------- User Setup --------------------- #

# required packages
  library(ggplot2)
  library(ggpubr)

# environmental conditions 
  temp <- 15 # air temperature (deg C)
  rh <- 0.7 # define humidity above the ocean, normalized to the ocean temperature (%)

# % liquid remaining   
  remaining_f <- seq(0.1,1,.1) # % of evaporating water body that remains as liquid (assuming pan evaporation)

# constants
  w <- 0.2 # turbulence factor (1 - pure diffusion, 0 - pure turbulence)
  lambda_ref <- 0.528   # lambda_ref. Luz and Barkan, 2010

# initial liquid water isotope values (delta values, not delta prime!)
  ini_l_18O <- -6.208 # initial liquid d18O (per mil), DI from North University Building, Ann Arbor, MI
  ini_l_17O <- -3.270 # initial liquid d17O (per mil), DI from North University Building, Ann Arbor, MI
  ini_l_2H <- -46.2 # initial liquid d2H (per mil), DI from North University Building, Ann Arbor, MI

  
## --------------------- Define Constants --------------------- #

# R values from IAEA reference sheet 2006
  R_VSMOW_18O <- 0.0020052 # ratio
  R_VSMOW_17O <- 0.0003799 # ratio
  R_VSMOW_d2H <- 0.00015576 # ratio

  theta_eq <- 0.529 # triple oxygen isotope exponent for equilibrium liquid-vapor fractionation. Barkan and Luz, 2005
  theta_diff <- 0.5185 # triple oxygen isotope exponent for diffusion of water vapor through air. Barkan and Luz, 2007
  diffratio_18 <- 1/0.9723 # oxygen D/D* for diffusion of water vapor through air.  Merlivat, 1978
  diffratio_2 <- 1/0.9755 # hydrogen D/D* for diffusion of water vapor through air. Merlivat, 1978

## --------------------- Calculate Fractionation Factors --------------------- #

#temperature dependent equilibrium fractionation factors, following Majoube, 1971 
  alpha18_eq <- exp(-2.0667 * 10^-3 - 0.4156/(temp+273.15) + (1.137*10^3)/(temp+273.15)^2) # Equilibrium  18O/16O liquid-vapor fractionation factor
  alpha17_eq <- alpha18_eq^theta_eq # Equilibrium 17O/16O liquid-vapor fractionation factor
  alpha2_eq <- exp(52.612*10^-3 - 76.248/(temp+273.15) + (24.844*10^3)/(temp+273.15)^2) # Equilibrium 2H/1H liquid-vapor fractionation factor

# fractionation factors for water vapor diffusing through air, following Passey and Ji, 2019
# Equation calculates alpha diff based on 'w', the ratio of non-fractionating turbulent water vapor transport to fractionating molecular diffusion transport.
  alpha18_diff <- w*diffratio_18+(1-w) 
  alpha17_diff <- alpha18_diff^theta_diff 
  alpha2_diff <- w*diffratio_2+(1-w) 


## --------------------- Calculate Isotopic Compositions as Water Evaporates --------------------- #

# assume that the vapor and unevaporated initial liquid are in isotopic equilibrium!! this is a simplifying assumption, but may not reflect natural conditions. 
 
# d18O
  ini_l_18O_R <- ((ini_l_18O/1000) +1)*R_VSMOW_18O # liquid d18O R (necessary to calculate R_w, Equation 11)
  ini_v_18O <- ((ini_l_18O+1000)/alpha18_eq)-1000 # vapor d18O, per mil (note that this vapor is in isotopic equilibrium with the initial unevaporated water)
  ini_v_18O_R <- ((ini_v_18O/1000) +1)*R_VSMOW_18O # vapor d18O R (necessary to calculate R_wss, Equation 13)

# d17O
  ini_l_17O_R <- ((ini_l_17O/1000) +1)*R_VSMOW_17O # liquid d17O R
  ini_v_17O <- ((ini_l_17O+1000)/alpha17_eq)-1000  # vapor d17O, per mil
  ini_v_17O_R <- ((ini_v_17O/1000) +1)*R_VSMOW_17O # vapor d17O R

# d2H
  ini_l_2H_R <- ((ini_l_2H/1000) +1)*R_VSMOW_d2H # liquid d2H R
  ini_v_2H <- ((ini_l_2H+1000)/alpha2_eq)-1000 # vapor d2H, per mil
  ini_v_2H_R <- ((ini_v_2H/1000) +1)*R_VSMOW_d2H # vapor d2H R


# u exponent, Equation 12 
  u18 <- (1 - alpha18_eq*alpha18_diff*(1 - rh))/(alpha18_eq*alpha18_diff*(1-rh))
  u17 <- (1 - alpha17_eq*alpha17_diff*(1 - rh))/(alpha17_eq*alpha17_diff*(1-rh))
  u2 <- (1 - alpha2_eq*alpha2_diff*(1 - rh))/(alpha2_eq*alpha2_diff*(1-rh))

# R_wss, Equation 13 
  R18_wss <- (alpha18_eq*rh*ini_v_18O_R)/(1-alpha18_eq*alpha18_diff*(1-rh))
  R17_wss <- (alpha17_eq*rh*ini_v_17O_R)/(1-alpha17_eq*alpha17_diff*(1-rh))
  R2_wss <- (alpha2_eq*rh*ini_v_2H_R)/(1-alpha2_eq*alpha2_diff*(1-rh))

# R_w, Equation 11
  R18_w <- remaining_f^u18*(ini_l_18O_R - R18_wss) + R18_wss
  R17_w <- remaining_f^u17*(ini_l_17O_R - R17_wss) + R17_wss
  R2_w <- remaining_f^u2*(ini_l_2H_R - R2_wss) + R2_wss

# convert R values of the remaining liquid to delta notation
  remaining_l_18O <- (R18_w/R_VSMOW_18O-1)*1000 # remaining liquid d18O, per mil
  remaining_l_17O <- (R17_w/R_VSMOW_17O-1)*1000 # remaining liquid d17O, per mil
  remaining_l_Dp17O <- (log(remaining_l_17O/1000+1) - lambda_ref*log(remaining_l_18O/1000+1))*10^6 # remaining liquid Dp17O, per meg
  remaining_l_2H <- (R2_w/R_VSMOW_d2H-1)*1000 # remaining liquid d2H, per mil
  remaining_l_dxs <- remaining_l_2H - 8*remaining_l_18O # remaining liquid d-excess, per mil


## --------------------- Plotting Section --------------------- #


#define theme
  supplement_6.theme <- theme(axis.text.x = element_text(margin = margin(t = 0.3, unit = "cm")),
                                axis.text.y = element_text(margin = margin(r = 0.3, unit = "cm")),
                                axis.ticks.length=unit(-0.15, "cm"),
                                axis.ticks = element_line(colour = "black"),
                                text = element_text(size = 10),
                                axis.title = element_text(size = 10))

dp18O_Dp17O_evap <- ggplot() +
  geom_point(data=data.frame(remaining_l_18O,remaining_l_Dp17O),aes(x=log(remaining_l_18O/1000+1)*1000,y=remaining_l_Dp17O)) +
  theme_bw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) + supplement_6.theme +
  labs(x=expression(delta*"'"^"18"*"O (\u2030)"), y=expression(Delta*"'"^"17"*"O (per meg)"))
plot(dp18O_Dp17O_evap)

d18O_dxs_evap <- ggplot() +
  geom_point(data=data.frame(remaining_l_18O,remaining_l_dxs),aes(x=remaining_l_18O,y=remaining_l_dxs)) +
  theme_bw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) + supplement_6.theme +
  labs(x=expression(delta^"18"*"O (\u2030)"), y="d-excess (\u2030)")
plot(d18O_dxs_evap)


## compile plots

S6_plots <- ggarrange(dp18O_Dp17O_evap,d18O_dxs_evap,nrow=1,ncol=2,align="hv")
plot(S6_plots)





