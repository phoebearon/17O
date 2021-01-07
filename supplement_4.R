
########################### Supplement 4. R script that calculates and plots oxygen and hydrogen isotope fractionation of waters in a closed system ###########################
# 
# Triple oxygen isotopes in the water cycle
# Author: Phoebe Aron
# Contact Phoebe Aron (paron@umich.edu) or Naomi Levin (nelevin@umich.edu) with questions and comments
# 
# This script steps through d18O, dp18O (dp = delta prime), d17O, dp17O, d2H, Dp17O, and d-excess variation as water fractionates in a closed system. 
#   
#   The script is broken into 3 sections: 
#       Section 1 - user setup, constants, and calculations of equilibrium and kinetic fractionation factors
#       Section 2 - isotope calculations as water evaporates from the ocean, condenses, and evaporates from an isolated water body (ie, pan evaporation)
#                   this is similar to the isotopic variation show in Figures 9 and 10 in Aron et al. 2020.
#       Section 3 - rayleigh distillation. this section starts with evaporated vapor from section 2 and calculates isotopic variation during Rayleigh distillation. 
#
#   Sections 2 and 3 are followed by separate plotting sections
#       Section 2 returns 4 plots (dp18O vs dp17O, dp18O vs Dp17O, d18O vs d-excess, and d18O vs d2H) of ocean water, equilibrium vapor, 
#         diffused vapor, condensation, evaporated vapor, and residual evaporated liquid (via pan evaporation).
#       Section 3 returns 3 plots that show liquid and vapor d18O, Dp17O, and d-excess variation during Rayleigh distillation
#  
#

# required packages
library(ggplot2)
library(ggpubr)

##################################################################
#                         SECTION 1                              #
##################################################################


## --------------------- User Setup --------------------- #

ocean_temp <- 10 # define ocean temperature (deg C). 
atm_temp <- 20 # define atmosphere (atm) temperature (deg C)
sfc_temp <- 16 # define surface (sfc) temperature of the evaporating isolated surface water body (deg C). 
rh <- 0.6 # define humidity above the ocean humidity, normalized to the ocean temperature (%)

rayleigh_f <- seq(0.1,1,0.1) # % vapor remaining after Rayleigh distillation 
remaining_f <- 0.9 # % of evaporating water body that remains as liquid (assuming pan evaporation)

w <- 0.2 # turbulence factor (1 - pure diffusion, 0 - pure turbulence)

lambda_ref <- 0.528 # define lambda_ref. Luz and Barkan, 2010
#lambda_ref <- 0.5268 #lambda_ref from Equation 8

## Starting conditions in the ocean
ocean_d18O_mil <- 0 # starting ocean water, per mil
ocean_D17O_meg <- 0 # starting ocean water, per meg
ocean_d17O_mil <- (exp(ocean_D17O_meg/10^6+lambda_ref*log(ocean_d18O_mil/1000+1))-1)*1000 # starting ocean water, per mil
ocean_d2H_mil <- 0  # starting ocean water, per mil
ocean_dxs <- 0 # starting ocean water, per mil

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
# ocean_temp = initial temperature of the ocean
ocean_alpha18_eq <- exp(-2.0667 * 10^-3 - 0.4156/(ocean_temp+273.15) + (1.137*10^3)/(ocean_temp+273.15)^2) # Equilibrium  18O/16O liquid-vapor fractionation factor
ocean_alpha17_eq <- ocean_alpha18_eq^theta_eq # Equilibrium 17O/16O liquid-vapor fractionation factor
ocean_alpha2_eq <- exp(52.612*10^-3 - 76.248/(ocean_temp+273.15) + (24.844*10^3)/(ocean_temp+273.15)^2) # Equilibrium 2H/1H liquid-vapor fractionation factor

# atm_temp = temperature during initial condensation and subsequent rayleigh distillation
atm_alpha18_eq <- exp(-2.0667 * 10^-3 - 0.4156/(atm_temp+273.15) + (1.137*10^3)/(atm_temp+273.15)^2) # Equilibrium 18O/16O liquid-vapor fractionation factor; Majoube 1971
atm_alpha17_eq <- atm_alpha18_eq^theta_eq # Equilibrium 17O/16O liquid-vapor fractionation factor
atm_alpha2_eq <- exp(52.612*10^-3 - 76.248/(atm_temp+273.15) + (24.844*10^3)/(atm_temp+273.15)^2) # Equilibrium 2H/1H liquid-vapor fractionation factor; Majoube 1971

# sfc_temp = temperature of the evaporating surface water body
sfc_alpha18_eq <- exp(-2.0667 * 10^-3 - 0.4156/(sfc_temp+273.15) + (1.137*10^3)/(sfc_temp+273.15)^2) # Equilibrium 18O/16O liquid-vapor fractionation factor; Majoube 1971
sfc_alpha17_eq <- sfc_alpha18_eq^theta_eq # Equilibrium 17O/16O liquid-vapor fractionation factor
sfc_alpha2_eq <- exp(52.612*10^-3 - 76.248/(sfc_temp+273.15) + (24.844*10^3)/(sfc_temp+273.15)^2) # Equilibrium 2H/1H liquid-vapor fractionation factor; Majoube 1971

# fractionation factors for water vapor diffusing through air, following Passey and Ji, 2019
# Equation calculates alpha diff based on 'w', the ratio of non-fractionating turbulent water vapor transport to fractionating molecular diffusion transport.
alpha18_diff <- w*diffratio_18+(1-w) 
alpha17_diff <- alpha18_diff^theta_diff 
alpha2_diff <- w*diffratio_2+(1-w) 

##################################################################
#                         SECTION 2                              #
##################################################################

## --------------------- Calculate isotopic compositions during ocean evaporation, condensation, surface evaporation --------------------- #

  
## 1. Evaporation from ocean (assume an infinite water source)
  # Redefine values in per mil to delta notation and as absolute R-values
  ocean_d18O_rat <- ocean_d18O_mil/1000 # initial water, delta ratio
  ocean_d18O_R <- (ocean_d18O_rat/1000+1)*R_VSMOW_18O # initial water, 18O/16O
  
  ocean_d17O_rat <- ocean_d17O_mil/1000 # initial water, delta ratio
  ocean_d17O_R <- (ocean_d17O_rat/1000+1)*R_VSMOW_17O # initial water, 17O/16O
  
  ocean_D17O_rat <- log(ocean_d17O_rat+1) - lambda_ref*log(ocean_d18O_rat+1) # initial water, delta ratio
  
  ocean_d2H_rat <- ocean_d2H_mil/1000 # initial water, delta ratio
  ocean_d2H_R <- (ocean_d2H_rat/1000+1)*R_VSMOW_d2H # initial water, 18O/16O
  
  # first, equilibrium vapor (eq = equilibrium, v = vapor)
  eq_v_d18O <- (ocean_d18O_R/R_VSMOW_18O*(1/(ocean_alpha18_eq))-1)*1000 # equilibrium vapor d18O, per mil
  eq_v_d17O <- (ocean_d17O_R/R_VSMOW_17O*(1/(ocean_alpha17_eq))-1)*1000 # equilibrium vapor d17O, per mil
  eq_v_Dp17O <- (log(eq_v_d17O/1000+1) - lambda_ref*log(eq_v_d18O/1000+1))*10^6 # equilibrium vapor Dp17O, per meg
  eq_v_d2H <- (ocean_d2H_R/R_VSMOW_d2H*(1/(ocean_alpha2_eq))-1)*1000 # equilibrium vapor d2H, per mil
  eq_v_dxs <- eq_v_d2H - 8*eq_v_d18O # equilibrium vapor d-excess, per mil
  
  #second, diffused atmospheric vapor (diff = diffused, v = vapor)
  diff_v_d18O <- (ocean_d18O_R/R_VSMOW_18O*(1/(ocean_alpha18_eq*(alpha18_diff*(1-rh)+rh)))-1)*1000 # diffused vapor d18O, per mil
  diff_v_d17O <- (ocean_d17O_R/R_VSMOW_17O*(1/(ocean_alpha17_eq*(alpha17_diff*(1-rh)+rh)))-1)*1000 # diffused vapor d17O, per mil
  diff_v_Dp17O <- (log(diff_v_d17O/1000+1) - lambda_ref*log(diff_v_d18O/1000+1))*10^6 # diffused vapor Dp17O, per meg
  diff_v_d2H <- (ocean_d2H_R/R_VSMOW_d2H*(1/(ocean_alpha2_eq*(alpha2_diff*(1-rh)+rh)))-1)*1000 # diffused vapor d2H, per mil
  diff_v_dxs <- diff_v_d2H - 8*diff_v_d18O # diffused vapor d-excess, per mil


## 2. Equilibrium condensation from diffused vapor 
  
  # to condense, convert vapor per mil delta values to isotope ratios and then absolute R values, then use the equilibrium 
    # fractionation factor to calculate the R value of liquid, then convert the liquid R value back to delta and per mil notation.
    # cond = condensation
  diff_v_d18O_rat <- diff_v_d18O/1000 # vapor ratio
  diff_v_d18O_R <- (diff_v_d18O_rat + 1)*R_VSMOW_18O # vapor R
  cond_l_d18O_R <- diff_v_d18O_R*atm_alpha18_eq # liquid R
  cond_l_d18O_rat <- cond_l_d18O_R/R_VSMOW_18O - 1 # liquid ratio
  cond_l_d18O_mil <- cond_l_d18O_rat*1000 # liquid per mil
  
  diff_v_d17O_rat <- diff_v_d17O/1000 # vapor ratio
  diff_v_d17O_R <- (diff_v_d17O_rat + 1)*R_VSMOW_17O # vapor R
  cond_l_d17O_R <- diff_v_d17O_R*atm_alpha17_eq # liquid R
  cond_l_d17O_rat <- cond_l_d17O_R/R_VSMOW_17O -1 # liquid ratio
  cond_l_d17O_mil <- cond_l_d17O_rat*1000 # liquid per mil
  
  cond_l_Dp17O <- (log(cond_l_d17O_rat+1) - lambda_ref*log(cond_l_d18O_rat+1))*10^6 # liquid Dp17O, per meg
  
  diff_v_d2H_rat <- diff_v_d2H/1000 # vapor ratio
  diff_v_d2H_R <- (diff_v_d2H_rat + 1)*R_VSMOW_d2H # vapor R
  cond_l_d2H_R <- diff_v_d2H_R*atm_alpha2_eq # liquid R
  cond_l_d2H_rat <- cond_l_d2H_R/R_VSMOW_d2H -1 # liquid ratio
  cond_l_d2H_mil <- cond_l_d2H_rat*1000 # liquid per mil
  
  cond_l_dxs <- cond_l_d2H_mil-8*cond_l_d18O_mil # liquid d-excess, per mil


## 3. Surface evaporation (assume pan evaporation -- an isolated water body evaporating to dryness). Equations are referenced from Aron et al. 2020.
  
  # u exponent, Equation 12
  u18 <- (1 - sfc_alpha18_eq*alpha18_diff*(1 - rh))/(sfc_alpha18_eq*alpha18_diff*(1-rh))
  u17 <- (1 - sfc_alpha17_eq*alpha17_diff*(1 - rh))/(sfc_alpha17_eq*alpha17_diff*(1-rh))
  u2 <- (1 - sfc_alpha2_eq*alpha2_diff*(1 - rh))/(sfc_alpha2_eq*alpha2_diff*(1-rh))
  
  # R_wss, Equation 13 
  R18_wss <- (sfc_alpha18_eq*rh*diff_v_d18O_R)/(1-sfc_alpha18_eq*alpha18_diff*(1-rh))
  R17_wss <- (sfc_alpha17_eq*rh*diff_v_d17O_R)/(1-sfc_alpha17_eq*alpha17_diff*(1-rh))
  R2_wss <- (sfc_alpha2_eq*rh*diff_v_d2H_R)/(1-sfc_alpha2_eq*alpha2_diff*(1-rh))
  
  # R_w, Equation 11
  R18_w <- remaining_f^u18*(cond_l_d18O_R - R18_wss) + R18_wss
  R17_w <- remaining_f^u17*(cond_l_d17O_R - R17_wss) + R17_wss
  R2_w <- remaining_f^u2*(cond_l_d2H_R - R2_wss) + R2_wss
  
  
  # convert R values to delta prime for the remaining liquid
  remaining_l_d18O <- (R18_w/R_VSMOW_18O-1)*1000 # remaining liquid d18O, per mil
  remaining_l_d17O <- (R17_w/R_VSMOW_17O-1)*1000 # remaining liquid d17O, per mil
  remaining_l_Dp17O <- (log(remaining_l_d17O/1000+1) - lambda_ref*log(remaining_l_d18O/1000+1))*10^6 # remaining liquid Dp17O, per meg
  remaining_l_d2H <- (R2_w/R_VSMOW_d2H-1)*1000 # remaining liquid d2H, per mil
  remaining_l_dxs <- remaining_l_d2H - 8*remaining_l_d18O # remaining liquid d-excess, per mil
  
  # evaporated vapor
  evaporated_v_d18O <- (R18_w/R_VSMOW_18O*(1/(sfc_alpha18_eq*(alpha18_diff*(1-rh)+rh)))-1)*1000 # evaporated vapor d18O, per mil
  evaporated_v_d17O <- (R17_w/R_VSMOW_17O*(1/(sfc_alpha17_eq*(alpha17_diff*(1-rh)+rh)))-1)*1000 # evaporated vapor d17O, per mil
  evaporated_v_Dp17O <- (log(evaporated_v_d17O/1000+1) - lambda_ref*log(evaporated_v_d18O/1000+1))*10^6 #evaporated vapor Dp17O, per meg
  evaporated_v_d2H <- (R2_w/R_VSMOW_d2H*(1/(sfc_alpha2_eq*(alpha2_diff*(1-rh)+rh)))-1)*1000 # evaporated vapor d2H, per mil
  evaporated_v_dxs <- evaporated_v_d2H - 8*evaporated_v_d18O # evaporated vapor d-excess, per mil


## --------------------- Plotting Section --------------------- #
  
  # open symbols = vapor 
      # open upward triangle = equilibrium vapor (shape=24)
      # open downward triangle = diffused vapor (shape=25)
      # open diamond = recycled vapor (shape=23)
  # filled symbols = liquid 
      # filled circle = ocean (shape=19)
      # filled square = precipitation (meteoric water) (shape=15)
      # filled diamond = evaporated water (shape=23)
      
  # define theme
  supplement_4.theme <- theme(axis.text.x = element_text(margin = margin(t = 0.3, unit = "cm")),
                                  axis.text.y = element_text(margin = margin(r = 0.3, unit = "cm")),
                                  axis.ticks.length=unit(-0.15, "cm"),
                                  axis.ticks = element_line(colour = "black"),
                                  text = element_text(size = 10),
                                  axis.title = element_text(size = 10))
  
  #combine data to plot
  data.to.plot <- cbind.data.frame(
    c("Ocean","Equilibrium Vapor", "Diffused Vapor","Condensation","Evaporated Vapor","Remaining Liquid"), # water types
    c(ocean_d18O_mil,eq_v_d18O,diff_v_d18O,cond_l_d18O_mil,evaporated_v_d18O,remaining_l_d18O), # d18O
    c(ocean_d17O_mil,eq_v_d17O,diff_v_d17O,cond_l_d17O_mil,evaporated_v_d17O,remaining_l_d17O), # d17O
    c(ocean_D17O_meg,eq_v_Dp17O,diff_v_Dp17O,cond_l_Dp17O,evaporated_v_Dp17O,remaining_l_Dp17O), # Dp17O
    c(ocean_d2H_mil,eq_v_d2H,diff_v_d2H,cond_l_d2H_mil,evaporated_v_d2H,remaining_l_d2H), #d2H
    c(ocean_dxs,eq_v_dxs,diff_v_dxs,cond_l_dxs,evaporated_v_dxs,remaining_l_dxs) # d-excess
  )
  
  names(data.to.plot) <- c("Water_ID","d18O","d17O","Dp17O","d2H","d_excess")
  
  #add delta prime 
  data.to.plot$dp18O <- log(data.to.plot$d18O/1000+1)*1000
  data.to.plot$dp17O <- log(data.to.plot$d17O/1000+1)*1000
  
  ## plots
  dp18O_dp17O_out <- ggplot() +
    geom_abline(slope=lambda_ref, intercept=0,color="black") + #dp18O dp17O reference line
    geom_point(data=data.to.plot, aes(x=dp18O, y=dp17O,shape=as.factor(Water_ID),fill=as.factor(Water_ID))) +
    theme_bw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) + supplement_4.theme +
    scale_shape_manual(values=c(15,25,24,23,19,23),name=NULL) + 
    scale_fill_manual(values=c("black","white","white","white","black","black"),name=NULL) +
    labs(x=expression(delta*"'"^"18"*"O (\u2030)"), y=expression(delta*"'"^"17"*"O (\u2030)")) +
    ggtitle(expression(delta*"'"^"18"*"O vs. "*delta*"'"^"17"*"O"))
  plot(dp18O_dp17O_out)
  
  
  dp18O_D17O_out <- ggplot() +
    geom_point(data=data.to.plot, aes(x=dp18O, y=Dp17O,shape=as.factor(Water_ID),fill=as.factor(Water_ID))) +
    theme_bw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) + supplement_4.theme +
    scale_shape_manual(values=c(15,25,24,23,19,23),name=NULL) + 
    scale_fill_manual(values=c("black","white","white","white","black","black"),name=NULL) +
    labs(x=expression(delta*"'"^"18"*"O (\u2030)"), y=expression(Delta*"'"^"17"*"O (per meg)")) +
    ggtitle(expression(delta*"'"^"18"*"O vs. "*Delta*"'"^"17"*"O"))
  plot(dp18O_D17O_out)
  
  d18O_dxs_out <- ggplot() +
    geom_point(data=data.to.plot, aes(x=d18O, y=d_excess,shape=as.factor(Water_ID),fill=as.factor(Water_ID))) +
    theme_bw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) + supplement_4.theme +
    scale_shape_manual(values=c(15,25,24,23,19,23),name=NULL) + 
    scale_fill_manual(values=c("black","white","white","white","black","black"),name=NULL) +
    labs(x=expression(delta^"18"*"O (\u2030)"), y="d-excess (\u2030)") +
    ggtitle(expression(delta^"18"*"O vs. d-excess"))
  plot(d18O_dxs_out)
  
  d18O_d2H_out <- ggplot() +
    geom_abline(slope=8, intercept=10,color="black") + #gmwl (Craig 1961)
    geom_point(data=data.to.plot, aes(x=d18O, y=d2H,shape=as.factor(Water_ID),fill=as.factor(Water_ID))) +
    theme_bw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) + supplement_4.theme +
    scale_shape_manual(values=c(15,25,24,23,19,23),name=NULL) + 
    scale_fill_manual(values=c("black","white","white","white","black","black"),name=NULL) +
    labs(x=expression(delta^"18"*"O (\u2030)"), y=expression(delta^"2"*"H (\u2030)")) +
    ggtitle(expression(delta^"18"*"O vs. "*delta^"2"*"H"))
  plot(d18O_d2H_out)
  

## compile plots

S4_out <- ggarrange(dp18O_dp17O_out,dp18O_D17O_out,d18O_dxs_out,d18O_d2H_out, nrow=1,ncol=4,align="hv")
plot(S4_out)


##################################################################
#                         SECTION 3                              #
##################################################################

## --------------------- Calculate isotopic compositions during Rayleigh distillation --------------------- #


## this scenario starts with the diffused vapor (v_d18O_diff, v_d17O_diff, v_d2H_diff) from the calculations above. 
## rather than condensing once and then evaporating from the surface, in this scenario the airmass continues to lose vapor via rainout. 

# rayleigh vapor (v = vapor)
rayleigh_v_d18O <- (diff_v_d18O + 1000) * rayleigh_f^(atm_alpha18_eq-1)-1000 # rayleigh distillation vapor d18O, per mil 
rayleigh_v_d17O <- (diff_v_d17O + 1000) * rayleigh_f^(atm_alpha17_eq-1)-1000 # rayleigh distillation vapor d17O, per mil 
rayleigh_v_Dp17O <- (log(rayleigh_v_d17O/1000+1) - lambda_ref*log(rayleigh_v_d18O/1000+1))*10^6 # rayleigh distillation vapor Dp17O, per meg 
rayleigh_v_d2H <- (diff_v_d2H + 1000) * rayleigh_f^(atm_alpha2_eq-1)-1000 # rayleigh distillation vapor d2H, per mil 
rayleigh_v_dxs <- rayleigh_v_d2H - 8*rayleigh_v_d18O # rayleigh distillation vapor d-excess, per mil 

# rayleigh liquid (l = liquid)
rayleigh_l_d18O <- atm_alpha18_eq*(rayleigh_v_d18O+1000)-1000 # rayleigh distillation liquid d18O, per mil
rayleigh_l_d17O <- atm_alpha17_eq*(rayleigh_v_d17O+1000)-1000 # rayleigh distillation liquid d17O, per mil
rayleigh_l_Dp17O <- (log(rayleigh_l_d17O/1000+1) - lambda_ref*log(rayleigh_l_d18O/1000+1))*10^6 # rayleigh distillation liquid Dp17O, per meg
rayleigh_l_d2H <- atm_alpha2_eq*(rayleigh_v_d2H+1000)-1000 # rayleigh distillation liquid d2H, per mil
rayleigh_l_dxs <- rayleigh_l_d2H - 8*rayleigh_l_d18O #rayleigh distillation liquid d-excess, per mil


# combine data to plot 
rayleigh.data <- cbind.data.frame(
  c(rep("Vapor",length(rayleigh_f)),rep("Liquid",length(rayleigh_f))),
  c(rayleigh_f,rayleigh_f),
  c(rayleigh_v_d18O,rayleigh_l_d18O),
  c(rayleigh_v_d17O,rayleigh_l_d17O),
  c(rayleigh_v_Dp17O,rayleigh_l_Dp17O),
  c(rayleigh_v_d2H,rayleigh_l_d2H),
  c(rayleigh_v_dxs,rayleigh_l_dxs)
)

names(rayleigh.data) <- c("Water_Type","Fraction_Remaining","d18O","d17O","Dp17O","d2H","d_excess")

rayleigh.d18O <- ggplot() +
  geom_point(data=rayleigh.data, aes(x=Fraction_Remaining, y=d18O, shape=as.factor(Water_Type))) +
  theme_bw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) + supplement_4.theme +
  scale_shape_manual(values=c(19,21), name=NULL) +
  labs(y=expression(delta^"18"*"O (\u2030)"), x="% Vapor Remaining") 
plot(rayleigh.d18O)

rayleigh.D17O <- ggplot() +
  geom_point(data=rayleigh.data, aes(x=Fraction_Remaining, y=Dp17O, shape=as.factor(Water_Type))) +
  theme_bw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) + supplement_4.theme +
  scale_shape_manual(values=c(19,21), name=NULL) +
  labs(y=expression(Delta*"'"^"17"*"O (per meg)"), x="% Vapor Remaining") 
plot(rayleigh.D17O)

rayleigh.d_excess <- ggplot() +
  geom_point(data=rayleigh.data, aes(x=Fraction_Remaining, y=d_excess, shape=as.factor(Water_Type))) +
  theme_bw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) + supplement_4.theme +
  scale_shape_manual(values=c(19,21), name=NULL) +
  labs(y="d-excess (\u2030)", x="% Vapor Remaining") 
plot(rayleigh.d_excess)

## compile plots

S4_rayleigh_out <- ggarrange(rayleigh.d18O,rayleigh.D17O,rayleigh.d_excess, nrow=1,ncol=3,align="hv")
plot(S4_rayleigh_out)


