#This is code to support the analyses presented in the manuscript Current Biology script
#"Assembly rules of coral reef fish communities along the depth gradient", by Pinhiero, MacDonald, et al. 
# Published in Current Biology, 2023.

#Source codes are available to link from the same repository
source ("FE_metrics.R")
source ("species_to_FE.R")
source ("multidimFD.R")
source ("multidimFbetaD.R")
source ("choosing_transect.R")

## library ----
library("picante")
library("gawdis")
library("hillR")
library("tidyverse")
library("patchwork")
library("reshape2")
library("reshape")
library("cluster")
library("mFD")
library("stringi")
library("lme4")
library("beta")
ZoneBlues <- c("#2171b5", "#6baed6", "#bdd7e7")

#Model numbers (1-37) correspond to models reported in the manuscript.
#A full list of models is available in the supplemental table S1. 
#Model numbers are also referenced in methods and results sections of the manuscript.
#Example code is given below and data is available in the supplemental material of the manuscript.
#Any further information is available from the authors.

#Model 1 is with models 22 - 34 ####

#Models 2 - 5 ####
test_dat <-  read.csv("subsamples_dat.csv") %>% 
  mutate(Zone = factor(Zone, levels = c("Shallow", "Upper", "Lower")))

#'upto' gives an updated starting point should the iteration process be stopped before completion
#This iteration process and rational is detailed in the methods sectoin of the manuscript - refer to detail on minimal sample area (MAS). 
upto <- dim(as.data.frame(outlierTestmSpD_Coastline_scaled_z [1,]) %>% drop_na)[1]+1

# start with upto = 1
#upto <- 1

N <- 10 #1000 used in manuscript
## Species Richness
for(m in upto:N){
  
  
  ####Model 2 - Dist IPR ####
  #run model
   mSpD_IPR_scaled_z <- lmer(Nsp_scaled_z ~  IPR*Zone+(1|Location),  data = test_dat %>% filter(iter == m))

   #in first iteration (m = 1) create matrix to store extracted R2 coefficients
  if(m == 1){ 
    mSpD_R2_IPR_Zone_1= matrix(nrow = 2000, ncol = 1)}
  caca<- as.data.frame(performance::r2(mSpD_IPR_scaled_z)[2])
  mSpD_R2_IPR_Zone_1 [m, 1 ] <- caca[1,1]
  
  #Gather model coefficients
  plotVals_scaled_z <- emmip(mSpD_IPR_scaled_z, Zone~IPR, at=list(IPR = seq(min(test_dat$IPR), max(test_dat$IPR), length = 100)), type = "response", plotit = F)%>%
    mutate(Iter = m)
  
  #in first iteration (m = 1) create matrix to store extracted model coefficients
  if(m == 1){
    plotVals_IPR_Zone_scaled_z_1 <- matrix(NA, nrow = 300000, 8);
    meansmSpD_IPR_Zone_scaled_z_1 <- matrix(NA, nrow = 1000, 3);
    contrmSpD_IPR_Zone_scaled_z_1 <- matrix(NA, nrow = 1000, 3);
    contrmSpD_IPR_Zonepval_scaled_z_1 <- matrix(NA, nrow = 1000, 3);
  }
  
  #Store values for plotting in matrix
  plotVals_IPR_Zone_scaled_z_1 [c((1:300)+(m-1)*300),c(1:8)] <- as.matrix(plotVals_scaled_z)
  
  #Calculate and store pairwise contrasts of model slopes 
  modmeans_scaled_z <-  emtrends(mSpD_IPR_scaled_z, pairwise ~ Zone, var = "IPR")
  meansmSpD_IPR_Zone_scaled_z_1 [m, ] <- as.data.frame(modmeans_scaled_z$emtrends)[,2]
  contrmSpD_IPR_Zone_scaled_z_1 [m, ] <- as.data.frame(modmeans_scaled_z$contrasts)[,2]
  contrmSpD_IPR_Zonepval_scaled_z_1 [m,] <- as.data.frame(modmeans_scaled_z$contrasts)[,6]
  
  #Store model test values from DHARMA
  resTest <- testResiduals(simulateResiduals(mSpD_IPR_scaled_z), plot = F)
  
  #in first iteration (m = 1) create matrix to store extracted model test coefficients
  if(m==1){uniformityTestmSpD_IPR_scaled_z <- matrix(nrow = 2, ncol = 1000);
  dispersionTestmSpD_IPR_scaled_z <- matrix(nrow = 2, ncol = 1000);
  outlierTestmSpD_IPR_scaled_z <- matrix(nrow = 2, ncol = 1000)}
  
  #Store modeltest coefficients
  uniformityTestmSpD_IPR_scaled_z[, m] <- as.vector(c(resTest$uniformity$statistic, resTest$uniformity$p.value))
  dispersionTestmSpD_IPR_scaled_z [, m] <- as.vector(c(resTest$dispersion$statistic, resTest$dispersion$p.value))
  outlierTestmSpD_IPR_scaled_z [, m] <- as.vector(c(resTest$outliers$statistic, resTest$outliers$p.value))
  
  
  #Model 3 - Connectivity and Isolation####
   mSpD_ConIso_scaled_z <- lmer(Nsp_scaled_z ~  ConIso*Zone+(1|Location),  data = test_dat %>% filter(iter == m))
  
  if(m == 1){ 
    mSpD_R2_ConIso_Zone_1= matrix(nrow = 2000, ncol = 1)}
  caca<- as.data.frame(performance::r2(mSpD_ConIso_scaled_z)[2])
  mSpD_R2_ConIso_Zone_1 [m, 1 ] <- caca[1,1]
  
  
  plotVals_scaled_z <- emmip(mSpD_ConIso_scaled_z, Zone~ConIso, at=list(ConIso = seq(min(test_dat$ConIso), max(test_dat$ConIso), length = 100)), type = "response", plotit = F)%>%
    mutate(Iter = m)
  if(m == 1){
    plotVals_ConIso_Zone_scaled_z_1 <- matrix(NA, nrow = 300000, 8);
    meansmSpD_ConIso_Zone_scaled_z_1 <- matrix(NA, nrow = 1000, 3);
    contrmSpD_ConIso_Zone_scaled_z_1 <- matrix(NA, nrow = 1000, 3);
    contrmSpD_ConIso_Zonepval_scaled_z_1 <- matrix(NA, nrow = 1000, 3);
  }
  
  plotVals_ConIso_Zone_scaled_z_1 [c((1:300)+(m-1)*300),c(1:8)] <- as.matrix(plotVals_scaled_z)
  
  modmeans_scaled_z <-  emtrends(mSpD_ConIso_scaled_z, pairwise ~ Zone, var = "ConIso")
  meansmSpD_ConIso_Zone_scaled_z_1 [m, ] <- as.data.frame(modmeans_scaled_z$emtrends)[,2]
  contrmSpD_ConIso_Zone_scaled_z_1 [m, ] <- as.data.frame(modmeans_scaled_z$contrasts)[,2]
  contrmSpD_ConIso_Zonepval_scaled_z_1 [m,] <- as.data.frame(modmeans_scaled_z$contrasts)[,6]
  
  if(m==1){uniformityTestmSpD_ConIso_scaled_z <- matrix(nrow = 2, ncol = 1000);
  dispersionTestmSpD_ConIso_scaled_z <- matrix(nrow = 2, ncol = 1000);
  outlierTestmSpD_ConIso_scaled_z <- matrix(nrow = 2, ncol = 1000)}
  
  resTest <- testResiduals(simulateResiduals(mSpD_ConIso_scaled_z), plot = F)
  
  uniformityTestmSpD_ConIso_scaled_z[, m] <- as.vector(c(resTest$uniformity$statistic, resTest$uniformity$p.value))
  dispersionTestmSpD_ConIso_scaled_z [, m] <- as.vector(c(resTest$dispersion$statistic, resTest$dispersion$p.value))
  outlierTestmSpD_ConIso_scaled_z [, m] <- as.vector(c(resTest$outliers$statistic, resTest$outliers$p.value))
  
  
  #Model 4 - Temperature####
   mSpD_Temp_scaled_z <- lmer(Nsp_scaled_z ~  Temp_z_g_scale*Zone+(1|Location),  data = test_dat %>% filter(iter == m))

  
  if(m == 1){ 
    mSpD_R2_Temp_Zone_1= matrix(nrow = 2000, ncol = 1)}
  caca<- as.data.frame(performance::r2(mSpD_Temp_scaled_z)[2])
  mSpD_R2_Temp_Zone_1 [m, 1 ] <- caca[1,1]
  
  plotVals_scaled_z <- emmip(mSpD_Temp_scaled_z, Zone~Temp_z_g_scale, at=list(Temp_z_g_scale = seq(min(test_dat$Temp_z_g_scale), max(test_dat$Temp_z_g_scale), length = 100)), type = "response", plotit = F)%>%
    mutate(Iter = m)
  if(m == 1){
    plotVals_Temp_Zone_scaled_z_1 <- matrix(NA, nrow = 300000, 8);
    meansmSpD_Temp_Zone_scaled_z_1 <- matrix(NA, nrow = 1000, 3);
    contrmSpD_Temp_Zone_scaled_z_1 <- matrix(NA, nrow = 1000, 3);
    contrmSpD_Temp_Zonepval_scaled_z_1 <- matrix(NA, nrow = 1000, 3);
  }
  
  plotVals_Temp_Zone_scaled_z_1 [c((1:300)+(m-1)*300),c(1:8)] <- as.matrix(plotVals_scaled_z)
  
  modmeans_scaled_z <-  emtrends(mSpD_Temp_scaled_z, pairwise ~ Zone, var = "Temp_z_g_scale")
  meansmSpD_Temp_Zone_scaled_z_1 [m, ] <- as.data.frame(modmeans_scaled_z$emtrends)[,2]
  contrmSpD_Temp_Zone_scaled_z_1 [m, ] <- as.data.frame(modmeans_scaled_z$contrasts)[,2]
  contrmSpD_Temp_Zonepval_scaled_z_1 [m,] <- as.data.frame(modmeans_scaled_z$contrasts)[,6]
  
  if(m==1){uniformityTestmSpD_Temp_scaled_z <- matrix(nrow = 2, ncol = 1000);
  dispersionTestmSpD_Temp_scaled_z <- matrix(nrow = 2, ncol = 1000);
  outlierTestmSpD_Temp_scaled_z <- matrix(nrow = 2, ncol = 1000)}
  
  
  resTest <- testResiduals(simulateResiduals(mSpD_Temp_scaled_z), plot = F)
  
  uniformityTestmSpD_Temp_scaled_z[, m] <- as.vector(c(resTest$uniformity$statistic, resTest$uniformity$p.value))
  dispersionTestmSpD_Temp_scaled_z [, m] <- as.vector(c(resTest$dispersion$statistic, resTest$dispersion$p.value))
  outlierTestmSpD_Temp_scaled_z [, m] <- as.vector(c(resTest$outliers$statistic, resTest$outliers$p.value))
  

  #Model 5 - Coral Reef area####
  mSpD_CRarea_scaled_z <- lmer(Nsp_scaled_z ~  CRarea*Zone+(1|Location),  data = test_dat %>% filter(iter == m))
  
  performance_CRarea_int_scaled_z [(m*2-1): (m*2),6:11]  <- performance_CRarea_int_scaled_z[1:2, 5:10] %>% as.matrix()
  
  if(m == 1){ 
    mSpD_R2_CRarea_Zone_1= matrix(nrow = 2000, ncol = 1)}
  caca<- as.data.frame(performance::r2(mSpD_CRarea_scaled_z)[2])
  mSpD_R2_CRarea_Zone_1 [m, 1 ] <- caca[1,1]
  
  plotVals_scaled_z <- emmip(mSpD_CRarea_scaled_z, Zone~CRarea, at=list(CRarea = seq(min(test_dat$CRarea), max(test_dat$CRarea), length = 100)), type = "response", plotit = F)%>%
    mutate(Iter = m)
  if(m == 1){
    plotVals_CRarea_Zone_scaled_z_1 <- matrix(NA, nrow = 300000, 8);
    meansmSpD_CRarea_Zone_scaled_z_1 <- matrix(NA, nrow = 1000, 3);
    contrmSpD_CRarea_Zone_scaled_z_1 <- matrix(NA, nrow = 1000, 3);
    contrmSpD_CRarea_Zonepval_scaled_z_1 <- matrix(NA, nrow = 1000, 3);
  }
  
  plotVals_CRarea_Zone_scaled_z_1 [c((1:300)+(m-1)*300),c(1:8)] <- as.matrix(plotVals_scaled_z)
  
  modmeans_scaled_z <-  emtrends(mSpD_CRarea_scaled_z, pairwise ~ Zone, var = "CRarea")
  meansmSpD_CRarea_Zone_scaled_z_1 [m, ] <- as.data.frame(modmeans_scaled_z$emtrends)[,2]
  contrmSpD_CRarea_Zone_scaled_z_1 [m, ] <- as.data.frame(modmeans_scaled_z$contrasts)[,2]
  contrmSpD_CRarea_Zonepval_scaled_z_1 [m,] <- as.data.frame(modmeans_scaled_z$contrasts)[,6]
  
  
  if(m==1){uniformityTestmSpD_CRarea_scaled_z <- matrix(nrow = 2, ncol = 1000);
  dispersionTestmSpD_CRarea_scaled_z <- matrix(nrow = 2, ncol = 1000);
  outlierTestmSpD_CRarea_scaled_z <- matrix(nrow = 2, ncol = 1000)}
  
  resTest <- testResiduals(simulateResiduals(mSpD_CRarea_scaled_z), plot = F)
  
  uniformityTestmSpD_CRarea_scaled_z[, m] <- as.vector(c(resTest$uniformity$statistic, resTest$uniformity$p.value))
  dispersionTestmSpD_CRarea_scaled_z [, m] <- as.vector(c(resTest$dispersion$statistic, resTest$dispersion$p.value))
  outlierTestmSpD_CRarea_scaled_z [, m] <- as.vector(c(resTest$outliers$statistic, resTest$outliers$p.value))
  
}

#Example code for plotting scaled effects of distance to IPR ####

#dist IPR _ scaled 
glimpse(meansmSpD_IPR_Zone_scaled_z_1)
meansmSpD_IPR_Zone_scaled_dat <- as.data.frame(t(meansmSpD_IPR_Zone_scaled_z_1))
meansmSpD_IPR_Zone_scaledConfint <- matrix (NA, ncol=7, nrow=3)
colnames(meansmSpD_IPR_Zone_scaledConfint) <- c("mean", "LCL95", "UCL95", "LCL80", "UCL80", "LCL60", "UCL60" )
rownames(meansmSpD_IPR_Zone_scaledConfint)<- c("Shallow", "Upper", "Lower")
for(i in 1:3) {
  meansmSpD_IPR_Zone_scaledConfint[i,1]<- apply(meansmSpD_IPR_Zone_scaled_dat[i,],   MARGIN = 1, FUN = mean, na.rm = TRUE)
  meansmSpD_IPR_Zone_scaledConfint[i,2:3] <- quantile(meansmSpD_IPR_Zone_scaled_dat[i,], c(0.025, 0.975), na.rm = TRUE)
  meansmSpD_IPR_Zone_scaledConfint[i,4:5]<- quantile(meansmSpD_IPR_Zone_scaled_dat[i,], c(0.1, 0.9), na.rm = TRUE)
  meansmSpD_IPR_Zone_scaledConfint[i,6:7] <- quantile(meansmSpD_IPR_Zone_scaled_dat[i,], c(0.2, 0.8), na.rm = TRUE)
}
meansmSpD_IPR_Zone_scaledConfint_1 <- as.data.frame(meansmSpD_IPR_Zone_scaledConfint) %>% 
  mutate(Zone = factor(rownames(meansmSpD_IPR_Zone_scaledConfint), levels <- c("Lower", "Upper", "Shallow")), Metric = "Species richness", Statistic = "trend",
         Predictor = "Distance to coral triangle")

meansmSpD_IPR_Zone_scaledConfint_1

contrastsmSpD_IPR_Zone_scaled_dat_1 <- as.data.frame(contrmSpD_IPR_Zone_scaled_z_1)
contrastsmSpD_IPR_Zone_scaledConfint <- matrix (NA, ncol=7, nrow = 3)
colnames(contrastsmSpD_IPR_Zone_scaledConfint) <- c("mean", "LCL95", "UCL95", "LCL80", "UCL80", "LCL60", "UCL60" )
rownames(contrastsmSpD_IPR_Zone_scaledConfint)<- as.data.frame(modmeans_scaled_z$contrasts)[,1]
for(i in 1:3) {
  contrastsmSpD_IPR_Zone_scaledConfint[i,1]<- apply(as.data.frame(contrastsmSpD_IPR_Zone_scaled_dat_1[,i]),   MARGIN = 2, FUN = mean, na.rm = TRUE)
  contrastsmSpD_IPR_Zone_scaledConfint[i,2:3] <- quantile(contrastsmSpD_IPR_Zone_scaled_dat_1[,i], c(0.025, 0.975), na.rm = TRUE)
  contrastsmSpD_IPR_Zone_scaledConfint[i,4:5]<- quantile(contrastsmSpD_IPR_Zone_scaled_dat_1[,i], c(0.1, 0.9), na.rm = TRUE)
  contrastsmSpD_IPR_Zone_scaledConfint[i,6:7] <- quantile(contrastsmSpD_IPR_Zone_scaled_dat_1[,i], c(0.2, 0.8), na.rm = TRUE)
}

contrastsmSpD_IPR_Zone_scaledConfint_1 <- as.data.frame(contrastsmSpD_IPR_Zone_scaledConfint) %>% 
  mutate(Zone = rownames(contrastsmSpD_IPR_Zone_scaledConfint), Metric = "Species richness", Statistic = "trend contrast",
         Predictor = "Distance to coral triangle")

contrastsmSpD_IPR_Zone_scaledConfint_1 %>% glimpse


#Plot####
dodge <- position_dodge(width=0.3)

pIPReffects_scaled <- ggplot(meansmSpD_IPR_Zone_scaledConfint_1, aes(x = Zone))+
  geom_hline(yintercept = 0, linetype="dashed", color = "blue", alpha = 0.5) +
  geom_linerange( aes( ymin=LCL95 , ymax=UCL95),  size = 1,  alpha = 1, position=dodge)+
  geom_linerange( aes( ymin=LCL80 , ymax=UCL80),  size = 2,  alpha = 1, position=dodge)+
  geom_linerange( aes( ymin=LCL60 , ymax=UCL60),  size = 3,  alpha = 1, position=dodge)+
  geom_point(aes(y = mean), colour = 'red') +
  theme_classic()+
  theme(legend.position = "none")+
  coord_flip()+ #legend.position = "none"
  labs(title="Scaled",
       x = meansmSpD_IPR_Zone_scaledConfint_1$Predictor, y = "Effect") 

pIPRcont_scaled <- ggplot(contrastsmSpD_IPR_Zone_scaledConfint_1, aes(x = Zone))+
  geom_hline(yintercept = 0, linetype="dashed", color = "blue", alpha = 0.5) +
  geom_linerange( aes( ymin=LCL95 , ymax=UCL95),  size = 1,  alpha = 1, position=dodge)+
  geom_linerange( aes( ymin=LCL80 , ymax=UCL80),  size = 2,  alpha = 1, position=dodge)+
  geom_linerange( aes( ymin=LCL60 , ymax=UCL60),  size = 3,  alpha = 1, position=dodge)+
  geom_point(aes(y = mean), colour = 'red') +
  theme_classic()+
  theme(legend.position = "none")+
  coord_flip()+ #legend.position = "none"
  labs(x = "Zone contrast", y = "Difference in trends",
       title='Scaled')

plotVals_IPR_Zone_scaled_dat_1 <- as.data.frame(plotVals_IPR_Zone_scaled_z_1) %>%
  drop_na() 
colnames(plotVals_IPR_Zone_scaled_dat_1) <- colnames(plotVals_scaled_z) 
colnames(plotVals_IPR_Zone_scaled_dat_1)[2] <- "IPR"
plotVals_IPR_Zone_scaled_dat_1 <- plotVals_IPR_Zone_scaled_dat_1 %>% 
  mutate(Zone = factor(Zone, levels <- c("Lower", "Upper", "Shallow")),
         IPR = as.numeric(IPR),
         xvar = as.numeric(xvar),
         yvar = as.numeric(yvar)) 
plotVals_IPR_Zone_scaled_means <- plotVals_IPR_Zone_scaled_dat_1 %>% 
  group_by(Zone, IPR) %>% 
  dplyr::summarise(yvar_mean = mean(yvar))

pIPR_scaled <- ggplot(plotVals_IPR_Zone_scaled_dat_1, aes(IPR, yvar, colour = Zone, group = interaction(Zone, Iter)))+
  geom_line(alpha = 0.05)+
  geom_line(data = plotVals_IPR_Zone_scaled_means, aes(IPR, yvar_mean, group = Zone), alpha = 1, size = 1.5)+
  scale_colour_viridis_d()+
  theme_classic()+
  # ylim(0,120)+
  labs(x = meansmSpD_IPR_Zone_scaledConfint_1$Predictor, 
       y = expression(paste("Scaled local species richness (240", m^-2,")")))

pIPR_scaled_c <- ggplot(plotVals_IPR_Zone_scaled_dat_1, aes(IPR, yvar, colour = Zone, group = interaction(Zone, Iter)))+
  geom_line(alpha = 0.05)+
  scale_colour_viridis_d()+
  theme_classic()+
  # ylim(0,120)+
  labs(x = meansmSpD_IPR_Zone_scaledConfint_1$Predictor, 
       y = expression(paste("Scaled local species richness (240", m^-2,")")))

#view plots
pIPReffects_scaled 
pIPRcont_scaled
pIPR_scaled

#Stitch plots with 'patchwork'
pIPR+pIPR_scaled + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

pIPReffects +pIPReffects_scaled + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

pIPRcont+pIPRcont_scaled+ plot_layout(guides = "collect") & theme(legend.position = 'bottom')



#Models 6 - 9 ####

FRdat <- read.csv("subsamples_dat.csv") %>% 
  mutate(Zone = factor(Zone, levels = c("Shallow", "Upper", "Lower"))) %>% glimpse

FRdatSH <- FRdat %>% dplyr::filter(Zone == "Shallow")
FRdatUP <- FRdat %>% dplyr::filter(Zone == "Upper")
FRdatLW <- FRdat %>% dplyr::filter(Zone == "Lower")

upto <- dim(as.data.frame(contrmSpD_Temp_Zonepval_scaled_z_reduced2[,1]) %>% drop_na)[1]+1

# start with upto = 1
#upto <- 1

N <- 10 #1000 used in manuscript
#start loop####
for(i in upto:N){
  
  #Scaled_z Response variable
  # mSpDr_classic_sh_scaled_z <- lm(Nsp_scaled_z ~ IPR+Temp_z_g_scale+CRarea+ConIso+CenterD,  data = FRdatSH %>% dplyr::filter(iter == i))
  # mSpDr_classic_up_scaled_z <- lm(Nsp_scaled_z ~ IPR+Temp_z_g_scale+CRarea+ConIso+CenterD,  data = FRdatUP %>% dplyr::filter(iter == i))
  # mSpDr_classic_lw_scaled_z <- lm(Nsp_scaled_z ~ IPR+Temp_z_g_scale+CRarea+ConIso+CenterD,  data = FRdatLW %>% dplyr::filter(iter == i))
  
  mSpDr_classic_int_scaled_z_reduced2 <- lmer(Nsp_scaled_z ~  Zone*IPR +   Zone*CRarea +  Zone*ConIso + Zone*Temp_z_g_scale+(1|Location) ,  data = FRdat %>% dplyr::filter(iter == i))
  
  caca<-r2beta(mSpDr_classic_int_scaled_z_reduced2, method = 'sgv',partial = TRUE)
  caca<-caca[order(caca$Effect),]
  if(i == 1){mSpD_R2LR_classic_int_scaled_z_reduced2 <- matrix(ncol = 15, nrow = 1e3)}
  mSpD_R2LR_classic_int_scaled_z_reduced2[i,1:15] <- caca$Rsq
  if(i == 1){values_mod_mSpDLR_classic_int_scaled_z_reduced2 <- matrix(ncol = 15, nrow = 1e3)}
  values_mod_mSpDLR_classic_int_scaled_z_reduced2[i, 1:15] <- fixef (mSpDr_classic_int_scaled_z_reduced2)
  
  # # residual tests_scaled_z ####
  
  # resTest <- testResiduals(simulateResiduals(mSpDr_classic_int_scaled_z_reduced2), plot = F)
  # uniformityTest_classic_int_scaled_z_reduced2 [, i] <- as.vector(c(resTest$uniformity$statistic, resTest$uniformity$p.value))
  # dispersionTest_classic_int_scaled_z_reduced2 [, i] <- as.vector(c(resTest$dispersion$statistic, resTest$dispersion$p.value))
  # outlierTest_classic_int_scaled_z_reduced2 [, i] <- as.vector(c(resTest$outliers$statistic, resTest$outliers$p.value))
  # #
  # resTest <- testResiduals(simulateResiduals(mSpDr_classic_sh_scaled_z), plot = T)
  # uniformityTest_classic_sh_scaled_z [, i] <- as.vector(c(resTest$uniformity$statistic, resTest$uniformity$p.value))
  # dispersionTest_classic_sh_scaled_z [, i] <- as.vector(c(resTest$dispersion$statistic, resTest$dispersion$p.value))
  # outlierTest_classic_sh_scaled_z [, i] <- as.vector(c(resTest$outliers$statistic, resTest$outliers$p.value))
  # #
  # 
  # resTest <- testResiduals(simulateResiduals(mSpDr_classic_up_scaled_z), plot = T)
  # uniformityTest_classic_up_scaled_z [, i] <- as.vector(c(resTest$uniformity$statistic, resTest$uniformity$p.value))
  # dispersionTest_classic_up_scaled_z [, i] <- as.vector(c(resTest$dispersion$statistic, resTest$dispersion$p.value))
  # outlierTest_classic_up_scaled_z [, i] <- as.vector(c(resTest$outliers$statistic, resTest$outliers$p.value))
  # #
  # 
  # resTest <- testResiduals(simulateResiduals(mSpDr_classic_lw_scaled_z), plot = T)
  # uniformityTest_classic_lw_scaled_z [, i] <- as.vector(c(resTest$uniformity$statistic, resTest$uniformity$p.value))
  # dispersionTest_classic_lw_scaled_z [, i] <- as.vector(c(resTest$dispersion$statistic, resTest$dispersion$p.value))
  # outlierTest_classic_lw_scaled_z [, i] <- as.vector(c(resTest$outliers$statistic, resTest$outliers$p.value))
  # #
  # # FRdat_i <- FRdat %>% filter(iter == i)
  # 
  # plotVals_IPR_scaled_z_reduced2 <- emmip(mSpDr_classic_int_scaled_z_reduced2, Zone~IPR, at=list(IPR = seq(min(FRdat_i$IPR), max(FRdat_i$IPR), length = 100)), type = "response", plotit = F)%>%
  #   mutate(Iter = i)
  # if(i == 1){plotVals_IPR_Zone_m_scaled_z_reduced2 <- matrix(ncol = 8, nrow = 300*1e3)}
  # plotVals_IPR_Zone_m_scaled_z_reduced2[1:300+((i-1)*300),1:8] <- as.matrix(plotVals_IPR_scaled_z_reduced2)
  # modmeans_IPR_scaled_z_reduced2 <-  emtrends(mSpDr_classic_int_scaled_z_reduced2, pairwise ~ Zone, var = "IPR")
  # if(i == 1){meansmSpD_IPR_Zone_scaled_z_reduced2 <- matrix(ncol = 3, nrow = 1e3)}
  # meansmSpD_IPR_Zone_scaled_z_reduced2 [i, ] <- as.data.frame(modmeans_IPR_scaled_z_reduced2$emtrends)[,2]
  # sigtest <- 0 >= as.data.frame(modmeans_IPR_scaled_z_reduced2$emtrends)[,5] & 0 <= as.data.frame(modmeans_IPR_scaled_z_reduced2$emtrends)[,6]
  # if(i == 1){sigmSpD_IPR_Zone_scaled_z_reduced2 <- matrix(ncol = 3, nrow = 1e3)}
  # sigmSpD_IPR_Zone_scaled_z_reduced2 [i, ] <- if_else(sigtest == "TRUE", "NS", "Sig")
  # if(i == 1){contrmSpD_IPR_Zonepval_scaled_z_reduced2 <- matrix(ncol = 3, nrow = 1e3)}
  # contrmSpD_IPR_Zonepval_scaled_z_reduced2 [i, ] <- as.data.frame(modmeans_IPR_scaled_z_reduced2$contrasts)[,6]
  # if(i == 1){contrmSpD_IPR_Zone_scaled_z_reduced2 <- matrix(ncol = 3, nrow = 1e3)}
  # contrmSpD_IPR_Zone_scaled_z_reduced2 [i, ] <- as.data.frame(modmeans_IPR_scaled_z_reduced2$contrasts)[,2]
  # 
  # 
  # plotVals_CRarea_scaled_z_reduced2 <- emmip(mSpDr_classic_int_scaled_z_reduced2, Zone~CRarea, at=list(CRarea = seq(min(FRdat_i$CRarea), max(FRdat_i$CRarea), length = 100)), type = "response", plotit = F)%>%
  #   mutate(Iter = i)
  # if(i == 1){plotVals_CRarea_Zone_m_scaled_z_reduced2 <- matrix(ncol = 8, nrow = 300*1e3)}
  # plotVals_CRarea_Zone_m_scaled_z_reduced2[1:300+((i-1)*300),1:8] <- as.matrix(plotVals_CRarea_scaled_z_reduced2)
  # modmeans_CRarea_scaled_z_reduced2 <-  emtrends(mSpDr_classic_int_scaled_z_reduced2, pairwise ~ Zone, var = "CRarea")
  # if(i == 1){meansmSpD_CRarea_Zone_scaled_z_reduced2 <- matrix(ncol = 3, nrow = 1e3)}
  # meansmSpD_CRarea_Zone_scaled_z_reduced2 [i, ] <- as.data.frame(modmeans_CRarea_scaled_z_reduced2$emtrends)[,2]
  # sigtest <- 0 >= as.data.frame(modmeans_CRarea_scaled_z_reduced2$emtrends)[,5] & 0 <= as.data.frame(modmeans_CRarea_scaled_z_reduced2$emtrends)[,6]
  # if(i == 1){sigmSpD_CRarea_Zone_scaled_z_reduced2 <- matrix(ncol = 3, nrow = 1e3)}
  # sigmSpD_CRarea_Zone_scaled_z_reduced2 [i, ] <- if_else(sigtest == "TRUE", "NS", "Sig")
  # if(i == 1){contrmSpD_CRarea_Zonepval_scaled_z_reduced2 <- matrix(ncol = 3, nrow = 1e3)}
  # contrmSpD_CRarea_Zonepval_scaled_z_reduced2 [i, ] <- as.data.frame(modmeans_CRarea_scaled_z_reduced2$contrasts)[,6]
  # if(i == 1){contrmSpD_CRarea_Zone_scaled_z_reduced2 <- matrix(ncol = 3, nrow = 1e3)}
  # contrmSpD_CRarea_Zone_scaled_z_reduced2 [i, ] <- as.data.frame(modmeans_CRarea_scaled_z_reduced2$contrasts)[,2]
  #
  # 
  # plotVals_ConIso_scaled_z_reduced2 <- emmip(mSpDr_classic_int_scaled_z_reduced2, Zone~ConIso, at=list(ConIso = seq(min(FRdat_i$ConIso), max(FRdat_i$ConIso), length = 100)), type = "response", plotit = F)%>%
  #   mutate(Iter = i)
  # if(i == 1){plotVals_ConIso_Zone_m_scaled_z_reduced2 <- matrix(ncol = 8, nrow = 300*1e3)}
  # plotVals_ConIso_Zone_m_scaled_z_reduced2[1:300+((i-1)*300),1:8] <- as.matrix(plotVals_ConIso_scaled_z_reduced2)
  # modmeans_ConIso_scaled_z_reduced2 <-  emtrends(mSpDr_classic_int_scaled_z_reduced2, pairwise ~ Zone, var = "ConIso")
  # if(i == 1){meansmSpD_ConIso_Zone_scaled_z_reduced2 <- matrix(ncol = 3, nrow = 1e3)}
  # meansmSpD_ConIso_Zone_scaled_z_reduced2 [i, ] <- as.data.frame(modmeans_ConIso_scaled_z_reduced2$emtrends)[,2]
  # sigtest <- 0 >= as.data.frame(modmeans_ConIso_scaled_z_reduced2$emtrends)[,5] & 0 <= as.data.frame(modmeans_ConIso_scaled_z_reduced2$emtrends)[,6]
  # if(i == 1){sigmSpD_ConIso_Zone_scaled_z_reduced2 <- matrix(ncol = 3, nrow = 1e3)}
  # sigmSpD_ConIso_Zone_scaled_z_reduced2 [i, ] <- if_else(sigtest == "TRUE", "NS", "Sig")
  # if(i == 1){contrmSpD_ConIso_Zonepval_scaled_z_reduced2 <- matrix(ncol = 3, nrow = 1e3)}
  # contrmSpD_ConIso_Zonepval_scaled_z_reduced2 [i, ] <- as.data.frame(modmeans_ConIso_scaled_z_reduced2$contrasts)[,6]
  # if(i == 1){contrmSpD_ConIso_Zone_scaled_z_reduced2 <- matrix(ncol = 3, nrow = 1e3)}
  # contrmSpD_ConIso_Zone_scaled_z_reduced2 [i, ] <- as.data.frame(modmeans_ConIso_scaled_z_reduced2$contrasts)[,2]
  # 
  # 
  # plotVals_Temp_scaled_z_reduced2 <- emmip(mSpDr_classic_int_scaled_z_reduced2, Zone~Temp_z_g_scale, at=list(Temp_z_g_scale = seq(min(FRdat_i$Temp_z_g_scale), max(FRdat_i$Temp_z_g_scale), length = 100)), type = "response", plotit = F)%>%
  #   mutate(Iter = i)
  # if(i == 1){plotVals_Temp_Zone_m_scaled_z_reduced2 <- matrix(ncol = 8, nrow = 300*1e3)}
  # plotVals_Temp_Zone_m_scaled_z_reduced2[1:300+((i-1)*300),1:8] <- as.matrix(plotVals_Temp_scaled_z_reduced2)
  # modmeans_Temp_scaled_z_reduced2 <-  emtrends(mSpDr_classic_int_scaled_z_reduced2, pairwise ~ Zone, var = "Temp_z_g_scale")
  # if(i == 1){meansmSpD_Temp_Zone_scaled_z_reduced2 <- matrix(ncol = 3, nrow = 1e3)}
  # meansmSpD_Temp_Zone_scaled_z_reduced2 [i, ] <- as.data.frame(modmeans_Temp_scaled_z_reduced2$emtrends)[,2]
  # sigtest <- 0 >= as.data.frame(modmeans_Temp_scaled_z_reduced2$emtrends)[,5] & 0 <= as.data.frame(modmeans_Temp_scaled_z_reduced2$emtrends)[,6]
  # if(i == 1){sigmSpD_Temp_Zone_scaled_z_reduced2 <- matrix(ncol = 3, nrow = 1e3)}
  # sigmSpD_Temp_Zone_scaled_z_reduced2 [i, ] <- if_else(sigtest == "TRUE", "NS", "Sig")
  # if(i == 1){contrmSpD_Temp_Zonepval_scaled_z_reduced2 <- matrix(ncol = 3, nrow = 1e3)}
  # contrmSpD_Temp_Zonepval_scaled_z_reduced2 [i, ] <- as.data.frame(modmeans_Temp_scaled_z_reduced2$contrasts)[,6]
  # if(i == 1){contrmSpD_Temp_Zone_scaled_z_reduced2 <- matrix(ncol = 3, nrow = 1e3)}
  # contrmSpD_Temp_Zone_scaled_z_reduced2 [i, ] <- as.data.frame(modmeans_Temp_scaled_z_reduced2$contrasts)[,2]
  # 
  # 
  
  # if(i == 1){values_mod_mSpDLR_classic_sh_scaled_z <- matrix(ncol = 6, nrow = 1e3)}
  # values_mod_mSpDLR_classic_sh_scaled_z[i,1:6] <- coef (mSpDr_classic_sh_scaled_z)
  # 
  # caca<-r2beta(mSpDr_classic_sh_scaled_z, method = 'sgv', partial = TRUE)
  # caca<-caca[order(caca$Effect),]
  # if(i == 1){mSpD_R2LR_classic_sh_scaled_z <- matrix(ncol = 6, nrow = 1e3)}
  # mSpD_R2LR_classic_sh_scaled_z [i,1:6] <- caca$Rsq
  # if(i == 1){values_mod_mSpDLR_classic_up_scaled_z <- matrix(ncol = 6, nrow = 1e3)}
  # values_mod_mSpDLR_classic_up_scaled_z[i,1:6] <- coef (mSpDr_classic_up_scaled_z)
  # 
  # caca<-r2beta(mSpDr_classic_up_scaled_z, method = 'sgv',partial = TRUE)
  # caca<-caca[order(caca$Effect),]
  # if(i == 1){mSpD_R2LR_classic_up_scaled_z <- matrix(ncol = 6, nrow = 1e3)}
  # mSpD_R2LR_classic_up_scaled_z [i,1:6] <- caca$Rsq
  # if(i == 1){values_mod_mSpDLR_classic_lw_scaled_z <- matrix(ncol = 6, nrow = 1e3)}
  # values_mod_mSpDLR_classic_lw_scaled_z[i,1:6] <- coef(mSpDr_classic_lw_scaled_z)
  # 
  # caca<-r2beta(mSpDr_classic_lw_scaled_z, method = 'sgv',partial = TRUE)
  # caca<-caca[order(caca$Effect),]
  # if(i == 1){mSpD_R2LR_classic_lw_scaled_z <- matrix(ncol = 6, nrow = 1e3)}
  # mSpD_R2LR_classic_lw_scaled_z [i,1:6] <- caca$Rsq
  # if(i == 1){values_mod_mSpDLR_classic_add_scaled_z <- matrix(ncol = 8, nrow = 1e3)}
  # values_mod_mSpDLR_classic_add_scaled_z[i,1:8] <- fixef (mSpDr_classic_add_scaled_z)
  # 
  
  print(paste0("Completed iteration ", i))
}

#Plotting


mSpD_R2LR_classic_sh_dat_scaled_z <- as.data.frame(mSpD_R2LR_classic_sh_scaled_z) %>% drop_na
mSpD_R2LR_classic_sh_confint_scaled_z <- matrix (NA, ncol=7, length(VarsR2))
colnames(mSpD_R2LR_classic_sh_confint_scaled_z) <- c("mean", "LCL95", "UCL95", "LCL80", "UCL80", "LCL60", "UCL60" )
rownames(mSpD_R2LR_classic_sh_confint_scaled_z)<- VarsR2
for(i in 1:length(VarsR2)) {
  
  mSpD_R2LR_classic_sh_confint_scaled_z[i,1]<- mSpD_R2LR_classic_sh_dat_scaled_z[,i] %>% mean
  mSpD_R2LR_classic_sh_confint_scaled_z[i,2:3] <- quantile(mSpD_R2LR_classic_sh_dat_scaled_z[,i], c(0.025, 0.975), na.rm = TRUE)
  mSpD_R2LR_classic_sh_confint_scaled_z[i,4:5]<- quantile(mSpD_R2LR_classic_sh_dat_scaled_z[,i], c(0.1, 0.9), na.rm = TRUE)
  mSpD_R2LR_classic_sh_confint_scaled_z[i,6:7] <- quantile(mSpD_R2LR_classic_sh_dat_scaled_z[,i], c(0.2, 0.8), na.rm = TRUE)
}

mSpD_R2LR_classic_sh_confint_scaled_z <- as.data.frame(mSpD_R2LR_classic_sh_confint_scaled_z)%>% 
  mutate(Predictor = factor(VarsR2))
# %>%
#   arrange(mean)

mSpD_R2LR_classic_sh_confint_scaled_z
mSpD_R2LR_classic_sh_confint_scaled_z <- mSpD_R2LR_classic_sh_confint_scaled_z %>%
  mutate(Predictor = factor(Predictor, levels = as.vector(mSpD_R2LR_classic_sh_confint_scaled_z[,8])))


p_classic_sh_scaled_z_R2 <- ggplot(mSpD_R2LR_classic_sh_confint_scaled_z, aes(x = Predictor))+
  geom_linerange( aes( ymin=LCL95 , ymax=UCL95),  size = 1,  alpha = 1)+
  geom_linerange( aes( ymin=LCL80 , ymax=UCL80),  size = 2,  alpha = 1)+
  geom_linerange( aes( ymin=LCL60 , ymax=UCL60),  size = 3,  alpha = 1)+
  geom_point(aes(y = mean), colour = 'black') +
  scale_colour_viridis_d(direction = -1)+
  theme_classic()+
  theme(legend.position = "none")+
  coord_flip()+ #legend.position = "none"
  ylim(0,1)+
  labs(title="Shallow",
       x ="", y = expression(paste("Semi partial ", R^2)))

#



mSpD_R2LR_classic_up_dat_scaled_z <- as.data.frame(mSpD_R2LR_classic_up_scaled_z) %>% drop_na
mSpD_R2LR_classic_up_confint_scaled_z <- matrix (NA, ncol=7, length(VarsR2))
colnames(mSpD_R2LR_classic_up_confint_scaled_z) <- c("mean", "LCL95", "UCL95", "LCL80", "UCL80", "LCL60", "UCL60" )
rownames(mSpD_R2LR_classic_up_confint_scaled_z)<- VarsR2
for(i in 1:length(VarsR2)) {
  
  mSpD_R2LR_classic_up_confint_scaled_z[i,1]<- mSpD_R2LR_classic_up_dat_scaled_z[,i] %>% mean
  mSpD_R2LR_classic_up_confint_scaled_z[i,2:3] <- quantile(mSpD_R2LR_classic_up_dat_scaled_z[,i], c(0.025, 0.975), na.rm = TRUE)
  mSpD_R2LR_classic_up_confint_scaled_z[i,4:5]<- quantile(mSpD_R2LR_classic_up_dat_scaled_z[,i], c(0.1, 0.9), na.rm = TRUE)
  mSpD_R2LR_classic_up_confint_scaled_z[i,6:7] <- quantile(mSpD_R2LR_classic_up_dat_scaled_z[,i], c(0.2, 0.8), na.rm = TRUE)
}
mSpD_R2LR_classic_up_confint_scaled_z <- as.data.frame(mSpD_R2LR_classic_up_confint_scaled_z)%>% 
  mutate(Predictor = factor(VarsR2))
# %>%
#   arrange(mean)

mSpD_R2LR_classic_up_confint_scaled_z
mSpD_R2LR_classic_up_confint_scaled_z <- mSpD_R2LR_classic_up_confint_scaled_z %>%
  mutate(Predictor = factor(Predictor, levels = as.vector(mSpD_R2LR_classic_up_confint_scaled_z[,8])))


p_classic_up_scaled_z_R2 <- ggplot(mSpD_R2LR_classic_up_confint_scaled_z, aes(x = Predictor))+
  geom_linerange( aes( ymin=LCL95 , ymax=UCL95),  size = 1,  alpha = 1)+
  geom_linerange( aes( ymin=LCL80 , ymax=UCL80),  size = 2,  alpha = 1)+
  geom_linerange( aes( ymin=LCL60 , ymax=UCL60),  size = 3,  alpha = 1)+
  geom_point(aes(y = mean), colour = 'yellow') +
  scale_colour_viridis_d(direction = -1)+
  theme_classic()+
  theme(legend.position = "none")+
  coord_flip()+ #legend.position = "none"
  ylim(0,1)+
  labs(title="Upper",
       x ="", y = expression(paste("Semi partial ", R^2)))


#

LRvarsR2_classic_scaled_z <- c("Distance to center of diversity",  "Geographic isolation", "Coral reef area", "Distance to coral triangle", "Absolute lattitude", "Model",  "Mean temperature within depth zone")

mSpD_R2LR_classic_lw_dat_scaled_z <- as.data.frame(mSpD_R2LR_classic_lw_scaled_z) %>% drop_na
mSpD_R2LR_classic_lw_confint_scaled_z <- matrix (NA, ncol=7, length(VarsR2))
colnames(mSpD_R2LR_classic_lw_confint_scaled_z) <- c("mean", "LCL95", "UCL95", "LCL80", "UCL80", "LCL60", "UCL60" )
rownames(mSpD_R2LR_classic_lw_confint_scaled_z)<- VarsR2
for(i in 1:length(VarsR2)) {
  
  mSpD_R2LR_classic_lw_confint_scaled_z[i,1]<- mSpD_R2LR_classic_lw_dat_scaled_z[,i] %>% mean
  mSpD_R2LR_classic_lw_confint_scaled_z[i,2:3] <- quantile(mSpD_R2LR_classic_lw_dat_scaled_z[,i], c(0.025, 0.975), na.rm = TRUE)
  mSpD_R2LR_classic_lw_confint_scaled_z[i,4:5]<- quantile(mSpD_R2LR_classic_lw_dat_scaled_z[,i], c(0.1, 0.9), na.rm = TRUE)
  mSpD_R2LR_classic_lw_confint_scaled_z[i,6:7] <- quantile(mSpD_R2LR_classic_lw_dat_scaled_z[,i], c(0.2, 0.8), na.rm = TRUE)
}
mSpD_R2LR_classic_lw_confint_scaled_z <- as.data.frame(mSpD_R2LR_classic_lw_confint_scaled_z)%>% 
  mutate(Predictor = factor(VarsR2))
# %>%
#   arrange(mean)

mSpD_R2LR_classic_lw_confint_scaled_z
mSpD_R2LR_classic_lw_confint_scaled_z <- mSpD_R2LR_classic_lw_confint_scaled_z %>%
  mutate(Predictor = factor(Predictor, levels = as.vector(mSpD_R2LR_classic_lw_confint_scaled_z[,8])))


p_classic_lw_scaled_z_R2 <- ggplot(mSpD_R2LR_classic_lw_confint_scaled_z, aes(x = Predictor))+
  geom_linerange( aes( ymin=LCL95 , ymax=UCL95),  size = 1,  alpha = 1)+
  geom_linerange( aes( ymin=LCL80 , ymax=UCL80),  size = 2,  alpha = 1)+
  geom_linerange( aes( ymin=LCL60 , ymax=UCL60),  size = 3,  alpha = 1)+
  geom_point(aes(y = mean), colour = 'yellow') +
  scale_colour_viridis_d(direction = -1)+
  theme_classic()+
  theme(legend.position = "none")+
  coord_flip()+ #legend.position = "none"
  ylim(0,1)+
  labs(title="Lower",
       x ="", y = expression(paste("Semi partial ", R^2)))




(p_classic_sh_scaled_z_R2 +p_classic_up_scaled_z_R2+p_classic_lw_scaled_z_R2)


VarsR2_int_reduced2 <- caca$Effect

ZoneR2_int_reduced2 <- c("Shallow", "Shallow", "Shallow", "Model", "Shallow",
                         "Lower", "Lower","Lower","Lower","Lower",
                         "Upper","Upper","Upper","Upper","Upper")

plot_order <- c("Model", 
                "IPR", "ZoneUpper:IPR", "ZoneLower:IPR",
                "ConIso", "ZoneUpper:ConIso", "ZoneLower:ConIso", 
                "CRarea", "ZoneUpper:CRarea", "ZoneLower:CRarea",
                "Temp_z_g_scale", "ZoneUpper:Temp_z_g_scale", "ZoneLower:Temp_z_g_scale",
                "ZoneUpper", "ZoneLower") %>% rev

mSpD_R2LR_classic_int_dat_scaled_z_reduced2 <- as.data.frame(mSpD_R2LR_classic_int_scaled_z_reduced2) %>% drop_na
mSpD_R2LR_classic_int_confint_scaled_z_reduced2 <- matrix (NA, ncol=7, length(VarsR2_int_reduced2))
colnames(mSpD_R2LR_classic_int_confint_scaled_z_reduced2) <- c("mean", "LCL95", "UCL95", "LCL80", "UCL80", "LCL60", "UCL60" )
rownames(mSpD_R2LR_classic_int_confint_scaled_z_reduced2)<- VarsR2_int_reduced2
for(i in 1:length(VarsR2_int_reduced2)) {
  
  mSpD_R2LR_classic_int_confint_scaled_z_reduced2[i,1]<- mSpD_R2LR_classic_int_dat_scaled_z_reduced2[,i] %>% mean
  mSpD_R2LR_classic_int_confint_scaled_z_reduced2[i,2:3] <- quantile(mSpD_R2LR_classic_int_dat_scaled_z_reduced2[,i], c(0.025, 0.975), na.rm = TRUE)
  mSpD_R2LR_classic_int_confint_scaled_z_reduced2[i,4:5]<- quantile(mSpD_R2LR_classic_int_dat_scaled_z_reduced2[,i], c(0.1, 0.9), na.rm = TRUE)
  mSpD_R2LR_classic_int_confint_scaled_z_reduced2[i,6:7] <- quantile(mSpD_R2LR_classic_int_dat_scaled_z_reduced2[,i], c(0.2, 0.8), na.rm = TRUE)
}

mSpD_R2LR_classic_int_confint_scaled_z_reduced2 <- as.data.frame(mSpD_R2LR_classic_int_confint_scaled_z_reduced2)%>% 
  mutate(Predictor = factor(VarsR2_int_reduced2))
# %>%
#   arrange(mean)
mSpD_R2LR_classic_int_confint_scaled_z_reduced2 <- mSpD_R2LR_classic_int_confint_scaled_z_reduced2 %>%
  mutate(Predictor = factor(Predictor, levels = plot_order),
         Zone = factor(ZoneR2_int_reduced2, levels = c("Model", "Shallow", "Upper", "Lower")))


mSpD_R2LR_classic_int_confint_scaled_z_reduced2$Zone %>% levels()


p_classic_int_scaled_z_reduced2_R2 <- ggplot(mSpD_R2LR_classic_int_confint_scaled_z_reduced2, aes(x = Predictor, colour = Zone))+
  geom_linerange( aes( ymin=LCL95 , ymax=UCL95),  size = 1,  alpha = 1)+
  geom_linerange( aes( ymin=LCL80 , ymax=UCL80),  size = 2,  alpha = 1)+
  geom_linerange( aes( ymin=LCL60 , ymax=UCL60),  size = 3,  alpha = 1)+
  geom_point(aes(y = mean), colour = 'white') +
  scale_colour_viridis_d(direction = -1)+
  theme_classic()+
  theme(legend.position = "none")+
  coord_flip()+ #legend.position = "none"
  ylim(0,1)+
  labs(title="Shallow",
       x ="", y = expression(paste("Semi partial ", R^2)))

##Model 10, Locations * Depth ####
DepthSeqScaled <- (seq(1, 130, length = 260)-dCent)/dScle #values for unscaling data for plotting
DepthSliceScaled <- (c(5, 30, 60, 90, 120)-dCent)/dScle #values for unscaling data for plotting


mDLsp_nb <- glm.nb(Nsp ~ Location*Depth+offset(log(area)), data = LRdat)

caca <- as.data.frame(performance::r2(mDLsp_nb))
mDLsp_nb_R2 [, m] <- caca[1,1]

#print(check_model(mDLsp_nb))
#Plot Values
plotVals <- emmip(mDLsp_nb, Location~Depth, at=list(Depth = DepthSeqScaled), type = 'response', plotit = F)%>%
  mutate(SpRic = yvar, Iter = m)%>%distinct()
mDLsp_info <- plotVals%>%select(Location, Depth)
spLplotVals [,m] <- plotVals$yvar
# Mean values and contrasts at 5m, 30m, 60m, 90m, and 120m
means <- emmip(mDLsp_nb, Location~Depth, at=list(Depth = DepthSliceScaled), plotit = FALSE)
spLmeans [,m] <- means$yvar
#length(means$yvar) 60
ref <-ref_grid(mDLsp_nb, at=list(Depth = DepthSliceScaled))
spLCont <- emmeans(ref, ~  Location|Depth)
con <- as.data.frame(pairs(spLCont))
spLContrasts [,m] <- con$estimate
#length(con$estimate) 330

# resTest <- testResiduals(simulateResiduals(mDLsp_nb), plot = F)#plot = T to see plots
# uniformityTestmSpD_LocDepth [, m] <- as.vector(c(resTest$uniformity$statistic, resTest$uniformity$p.value))
# dispersionTestmSpD_LocDepth [, m] <- as.vector(c(resTest$dispersion$statistic, resTest$dispersion$p.value))
# outlierTestmSpD_LocDepth [, m] <- as.vector(c(resTest$outliers$statistic, resTest$outliers$p.value))



#Models 11 - 21, regional pool relationships ####
#### Model 11 ####
mRegIPRzoneD <- lm(RegionalPool_DepthMatch~Distance_from_the_IPR_km*Zone,regPool_filtered)
summary(mRegIPRzoneD)
confint(mRegIPRzoneD)
emtrends(mRegIPRzoneD, pairwise ~ Zone, var ="Distance_from_the_IPR_km")
performance::r2(mRegIPRzoneD)

#### Model 12 ####
mRegIPRD <- lm(RegionalPool_DepthMatch~Distance_from_the_IPR_km,regPool_filtered)
performance::r2(mRegIPRD)
summary(mRegIPRD)
confint(mRegIPRD)

#### Model 13 ####
mRegIPRzone <- lm(RegionalPool_GenMatch_DepthMatch~Distance_from_the_IPR_km*Zone,regPool_filtered)
summary(mRegIPRzone)
confint(mRegIPRzone)
emtrends(mRegIPRzone, pairwise ~ Zone, var="Distance_from_the_IPR_km", transform="response")
performance::r2(mRegIPRzone)

#### Model 14 ####
mzs<- lm(RegionalPool_GenMatch_DepthMatch~Distance_from_the_IPR_km*Slice, MEANregPoolDepthSlices)
summary(mzs)
emtrends(mzs, pairwise ~ Slice, var="Distance_from_the_IPR_km", transform="response")
performance::r2(mzs)

#### Model 15 ####
mod15 <- lm(MeanSpRic~Distance_from_the_IPR_km*Zone, meanSpRic)
emtrends(mod15, pairwise ~ Zone, var ="Distance_from_the_IPR_km")
performance::r2(mod15)

#### Model 16 ####
mod16 <- lm(Mean~Distance_from_the_IPR_km*Slice, MEANregPoolDepthSlices)
emtrends(mod16, pairwise ~ Slice, var ="Distance_from_the_IPR_km")
performance::r2(mod16)

#### Model 17 ####
modDelta <- lm(MeanDelta~Regional_richness,spRLocRegDeltaMean)
summary(modDelta)
confint(modDelta)
performance::r2(modDelta)

#### Model 18 ####
mod1Sb <- lm(Mean ~ log(RegPoolProv )*Slice, data  = MEANregPoolDepthSlices)
summary(mod1Sb)
confint(mod1Sb)
emtrends(mod1Sb, pairwise ~ Slice, var="RegPoolProv", transform="response")
performance::r2(mod1Sb)

#### Model 19 ####
mod1Sfb <- lm(Mean~ log(RegionalPool_GenMatch_DepthMatch)*Slice, data  = MEANregPoolDepthSlices)
summary(mod1Sfb)
confint(mod1Sfb)
emtrends(mod1Sfb, pairwise ~ Slice, var="RegionalPool_GenMatch_DepthMatch", transform="response")
performance::r2(mod1Sfb)

#### Model 20 ####
mod20 <- lm(LogRatioLoc ~ LogReg*Slice, data  = regPoolDepthSlices)
summary(mod20)
emtrends(mod20, ~ Slice, var="LogReg", transform="response")
mod20.emm <- emmeans(mod20, ~ Slice|LogReg)
mod20trends <-emtrends(mod20, pairwise ~ Slice, var = "LogReg")
mod20Slopes <- as.data.frame(mod20trends$emtrends)
pairs(mod20.emm, simple = "Slice")
eff_size(mod20.emm, sigma = sigma(mod20), edf = df.residual(mod20))
performance::r2(mod20) #0.92
confint(mod20)

#### Model 21 ####
mod21Prov <- lm(LogRatioLocProv ~ LogRegProv*Slice, data  = regPoolDepthSlices)
summary(mod21Prov)
emtrends(mod21Prov, ~ Slice, var="LogRegProv", transform="response")
mod21Prov.emm <- emmeans(mod21Prov, ~ Slice|LogRegProv)
mod21Provtrends <-emtrends(mod21Prov, pairwise ~ Slice, var = "LogRegProv")
mod21ProvSlopes <- as.data.frame(mod21Provtrends$emtrends)
pairs(mod21Prov.emm, simple = "Slice")
eff_size(mod21Prov.emm, sigma = sigma(mod21Prov), edf = df.residual(mod21Prov))
performance::r2(mod21Prov)
confint(mod21)

#Models 22 - 27, Families####

for (m in upto:N) {
  chosen_transec  <- choosing_transec (sub1, msa)
  sub2 <- droplevels (sub1 [sub1$transect_id %in% chosen_transec, ])
  database_red <- droplevels (abund [abund$transect_id %in% sub2$transect_id, ])
  rm(sub2)
  
  
  # Organization abundance data #
  abundance <- ddply (database_red,. (Name,  ID_transect), summarise, abun=sum(Abun))
  abundance <- cast (ID_transect~Name, value = "abun", fun.aggregate = sum, data=abundance)
  rownames (abundance) <- abundance$ID_transect
  abundance <- abundance [,-1]
  abundance <- data.matrix (abundance)
  
  
  # create presence-absence 
  occ <- abundance
  occ[occ > 0] <- 1 
  
  spr <- as.data.frame(rowSums(occ))
  colnames(spr) <- "Nsp"
  spr <- spr %>% mutate(ID_transect = rownames (spr))
  abnd <- as.data.frame(rowSums(abundance))
  colnames(abnd) <- "AbundTot"
  abnd <- abnd %>% mutate(ID_transect = rownames (abnd))
  
  LR <- left_join(abnd, spr)
  
  # MODELS #
  #Create environmental factors
  enviro2 <- abund %>%
    dplyr::group_by(Ocean, Provinces, Location, ID_Locality_Mesophotic_Zone, area, ID_transect) %>%
    dplyr::summarise(Zone = Mesophotic2, 
                     IPR = Distance_from_the_IPR_km, 
                     CenterD = Distance_from_the_center_of_diversity_km, 
                     CRarea = Coral_reef_area_km2, 
                     Depth = Depth_m) %>%
    distinct() %>%
    ungroup() %>%
    mutate(Pos = seq(1:length(ID_Locality_Mesophotic_Zone)))
  
  
  enviroNum <- enviro2 %>%
    select(-Location, -ID_Locality_Mesophotic_Zone, -Ocean, -Provinces, -Zone, -area, -ID_transect, -Pos)
  
  envAttrCenter <- attr(scale(enviroNum, center = TRUE, scale = TRUE), "scaled:center")
  envAttrScale <- attr(scale(enviroNum, center = TRUE, scale = TRUE), "scaled:scale")
  
  enviroNum <- as.data.frame(scale(enviroNum, center = TRUE, scale = TRUE))
  
  dScle <- envAttrScale["Depth"]
  dCent <- envAttrCenter["Depth"]
  DepthSeqScaled <- (seq(1, 130, length = 260)-dCent)/dScle
  DepthSliceScaled <- (c(5, 30, 60, 90, 120)-dCent)/dScle
  
  enviroNum <- enviroNum %>%
    mutate(Pos = seq(1:length(enviro2$ID_Locality_Mesophotic_Zone)))
  
  enviro2 <- enviro2 %>%
    select(Ocean, Provinces, Location, Zone, ID_Locality_Mesophotic_Zone, ID_transect, area, Pos)
  
  enviroNum <- left_join(enviroNum, enviro2, by = c("Pos" = "Pos"))
  
  
  
  LRdat <- left_join(LR, enviroNum)
  
  #Model 1 - Oceans ####
  mDOsp_nb <- glmer.nb(Nsp ~ Ocean*Depth+(1|Location)+offset(log(area)), data = LRdat,
                       control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e5)))
  
  # #R2
  caca<- as.data.frame(performance::r2(mDOsp_nb)[2])
  mDOsp_nb_R2 [, m] <- caca[1,1]
  # 
  # #Plot Values
  plotVals <- emmip(mDOsp_nb, Ocean~Depth, at=list(Depth = DepthSeqScaled), type = "response", plotit = F)%>%
    mutate(Iter = m)
  mDOsp_info <- plotVals%>%select(Ocean,Depth)
  spOplotVals [,m] <- plotVals$yvar
  # Mean values and contrasts at 5m, 30m, 60m, 90m, and 120m
  means <- emmip(mDOsp_nb, Ocean~Depth, at=list(Depth = DepthSliceScaled), type = "response", plotit = FALSE)
  spOmeans [,m] <- means$yvar
  #length(means$yvar) 10
  ref <-ref_grid(mDOsp_nb, at=list(Depth = DepthSliceScaled))
  spOCont <- emmeans(ref, ~  Ocean|Depth, type = "response")
  con <- as.data.frame(contrast(spOCont))%>%filter(contrast == "Atlantic effect")
  spOContrasts [,m] <- con$ratio
  # #length(con$estimate) 10
  # 
  # resTest <- testResiduals(simulateResiduals(mDOsp_nb), plot = F)#plot = T to see plots
  # uniformityTestmSpD_OceanDepth [, m] <- as.vector(c(resTest$uniformity$statistic, resTest$uniformity$p.value))
  # dispersionTestmSpD_OceanDepth [, m] <- as.vector(c(resTest$dispersion$statistic, resTest$dispersion$p.value))
  # outlierTestmSpD_OceanDepth [, m] <- as.vector(c(resTest$outliers$statistic, resTest$outliers$p.value))
  # 
  
  
  
# Species richness by Family 
families <- ddply (database_red,. (ID_transect, Family), summarise,
                   Total_SpR_Family=length(unique(Name)),
                   Total_abd_Family=sum(Abundance))
families_SpR <- cast (Family~ID_transect, value = "Total_SpR_Family", data=families)
families_SpR[is.na(families_SpR)] <- 0
families_SpR <- melt (families_SpR, id=c("Family"))
rownames (families_SpR) <- NULL
colnames (families_SpR)[2] <- "ID_transect"
colnames (families_SpR)[3] <- "Total_SpR_Family"

families_SpR <- left_join(families_SpR, LR)
families_SpR$Prop_Family <- families_SpR$Total_SpR_Family/families_SpR$Nsp
families_SpR <- left_join(families_SpR, LRdat)%>%drop_na()

Labridae <- droplevels (subset (families_SpR, Family=="Labridae"))
Chaetodontidae <- droplevels (subset (families_SpR, Family=="Chaetodontidae"))
Pomacentridae <- droplevels (subset (families_SpR, Family=="Pomacentridae"))
Serranidae <- droplevels (subset (families_SpR, Family=="Serranidae"))
Acanthuridae <- droplevels (subset (families_SpR, Family=="Acanthuridae"))
Scaridae <- droplevels (subset (families_SpR, Family=="Scaridae"))

#Model 22 ####
#Note: this first model has examples of the matrices are built to store model outputs
#copy and change Family names, as appropriate for models 22 - 34

#  #Oceans - Labridae

#Matrices
spOSlopesm_labridae <- matrix (NA, ncol=N, nrow = 2)
spOInterceptsm_labridae <- matrix (NA, ncol=N, nrow = 2)
spOContrastsm_labridae <- matrix (NA, ncol=N, nrow = 1)
dispO_labridae_m <- matrix (NA, ncol=N, nrow = 1)

spOContrasts5m_labridae <- matrix (NA, ncol=N, nrow = 2)
spOContrasts30m_labridae <- matrix (NA, ncol=N, nrow = 2)
spOContrasts60m_labridae <- matrix (NA, ncol=N, nrow = 2)
spOContrasts90m_labridae <- matrix (NA, ncol=N, nrow = 2)
spOContrasts120m_labridae <- matrix (NA, ncol=N, nrow = 2)
spOMean5m_labridae <- matrix (NA, ncol=N, nrow = 2)
spOMean30m_labridae <- matrix (NA, ncol=N, nrow = 2)
spOMean60m_labridae <- matrix (NA, ncol=N, nrow = 2)
spOMean90m_labridae <- matrix (NA, ncol=N, nrow = 2)
spOMean120m_labridae <- matrix (NA, ncol=N, nrow = 2)


hist(Labridae$Nsp_Fam)
# mDOspG_Labridae <- glmer(Total_SpR_Family ~ Ocean*Depth+(1|Location), family = gaussian, data = Labridae)
mDOsp_Labridae <- glmer.nb(Nsp_Fam ~ Depth*Ocean+offset(log(area))+(1|Location),  data = Labridae)
#print(check_model(mDOsp_Labridae))
#R2
caca<- as.data.frame(performance::r2(mDOsp_Labridae)[2])
mDOsp_Labridae_R2 [, m] <- caca[1,1]

#Plot Values
plotVals <- emmip(mDOsp_Labridae, Ocean~Depth, at=list(Depth = DepthSeqScaled), type = "response",plotit = FALSE)%>%
  mutate(spRic = yvar, Iter = m)
spOplotVals_Labridae[,m] <- plotVals$yvar
# Mean values and contrasts at 5m, 30m, 60m, 90m, and 120m
means <- emmip(mDOsp_Labridae, Ocean~Depth, at=list(Depth = DepthSliceScaled), type = "response",plotit = FALSE)
spOmeans_Labridae [,m]<- means$yvar
ref <-ref_grid(mDOsp_Labridae, at=list(Depth = DepthSliceScaled))
spOCont <- emmeans(ref, ~  Ocean|Depth, type = "response")
con <- as.data.frame(contrast(spOCont))%>%filter(contrast == "Atlantic effect")
spOContrasts_Labridae [,m]<- con$ratio
# resTest <- testResiduals(simulateResiduals(mDOsp_Labridae), plot = F)
# uniformityTestmSpD_OceanDepth_Labridae [, m] <- as.vector(c(resTest$uniformity$statistic, resTest$uniformity$p.value))
# dispersionTestmSpD_OceanDepth_Labridae [, m] <- as.vector(c(resTest$dispersion$statistic, resTest$dispersion$p.value))
# outlierTestmSpD_OceanDepth_Labridae [, m] <- as.vector(c(resTest$outliers$statistic, resTest$outliers$p.value))

mDOsp_labridae <- lmer(Total_SpR_Family ~ Ocean*Depth + (1|Location), data = Labridae)
#   ref <- ref_grid(mDOsp_labridae, at=list(Depth = 0))
#   spOIntercepts_labridae <- emmeans(ref, pairwise ~ Ocean)
#   spOIntercepts_labridae <- as.data.frame(spOIntercepts_labridae$emmeans)[,1:2]
#   spOSlopes_labridae <- emtrends(mDOsp_labridae, pairwise ~ Ocean, var = "Depth")
#   spOSlopes_labridae <- as.data.frame(spOSlopes_labridae$emtrends)[,1:2]
#   datspO_labridae <- left_join(spOIntercepts_labridae,spOSlopes_labridae)
#   spOSlopesm_labridae [,m] <- spOSlopes_labridae [,2]
#   spOInterceptsm_labridae [,m] <- spOIntercepts_labridae [,2]
#   
#   # Mean values and contrasts at 5m, 30m, 60m, 90m, and 120m
#   ref <- ref_grid(mDOsp_labridae, at=list(Depth = 5))
#   spOContrasts <- emmeans(ref, pairwise ~ Ocean)
#   test <-  as.data.frame(spOContrasts$contrasts)[,c(2,6)]
#   colnames(test) <- c()
#   spOContrasts5m_labridae [,m] <- t(test)
#   spOMean5m_labridae [,m] <- as.data.frame(spOContrasts$emmeans)[,c(2)]
#   
#   ref1 <- ref_grid(mDOsp_labridae, at=list(Depth = 30))
#   spOContrasts1 <- emmeans(ref1, pairwise ~ Ocean)
#   test <-  as.data.frame(spOContrasts1$contrasts)[,c(2,6)]
#   colnames(test) <- c()
#   spOContrasts30m_labridae [,m] <- t(test)
#   spOMean30m_labridae [,m] <- as.data.frame(spOContrasts1$emmeans)[,c(2)]
#   
#   ref2 <- ref_grid(mDOsp_labridae, at=list(Depth = 60))
#   spOContrasts2 <- emmeans(ref2, pairwise ~ Ocean)
#   test <-  as.data.frame(spOContrasts2$contrasts)[,c(2,6)]
#   colnames(test) <- c()
#   spOContrasts60m_labridae [,m] <- t(test)
#   spOMean60m_labridae [,m] <- as.data.frame(spOContrasts2$emmeans)[,c(2)]
#   
#   ref3 <- ref_grid(mDOsp_labridae, at=list(Depth = 90))
#   spOContrasts3 <- emmeans(ref3, pairwise ~ Ocean)
#   test <-  as.data.frame(spOContrasts3$contrasts)[,c(2,6)]
#   colnames(test) <- c()
#   spOContrasts90m_labridae [,m] <- t(test)
#   spOMean90m_labridae [,m] <- as.data.frame(spOContrasts3$emmeans)[,c(2)]
#   
#   ref4 <- ref_grid(mDOsp_labridae, at=list(Depth = 120))
#   spOContrasts4 <- emmeans(ref4, pairwise ~ Ocean)
#   test <-  as.data.frame(spOContrasts4$contrasts)[,c(2,6)]
#   colnames(test) <- c()
#   spOContrasts120m_labridae [,m] <- t(test)
#   spOMean120m_labridae [,m] <- as.data.frame(spOContrasts4$emmeans)[,c(2)]
#   

#Model 23  ####
#Oceans - Acanthuridae
hist(Acanthuridae$Nsp_Fam)
mDOsp_Acanthuridae <- glmmTMB(Nsp_Fam ~ Ocean*Depth+offset(log(area))+(1|Location),ziformula=~1,family=poisson, data = Acanthuridae)

#plotResiduals(mDOsp_Acanthuridae)
#R2    
caca<- as.data.frame(performance::r2(mDOsp_Acanthuridae)[2])
mDOsp_Acanthuridae_R2 [, m] <- caca[1,1]

#Plot Values     
plotVals <- emmip(mDOsp_Acanthuridae, Ocean~Depth, at=list(Depth = DepthSeqScaled),type = "response", plotit = FALSE)%>%
  mutate(spRic = yvar, Iter = m)
spOplotVals_Acanthuridae [,m]<- plotVals$yvar
# Mean values and contrasts at 5m, 30m, 60m, 90m, and 120m
means <- emmip(mDOsp_Acanthuridae, Ocean~Depth, at=list(Depth = DepthSliceScaled),type = "response", plotit = FALSE)
spOmeans_Acanthuridae [,m]<- means$yvar
ref <-ref_grid(mDOsp_Acanthuridae, at=list(Depth = DepthSliceScaled))
spOCont <- emmeans(ref, ~  Ocean|Depth, type = "response")
con <- as.data.frame(contrast(spOCont))%>%filter(contrast == "Atlantic effect")
spOContrasts_Acanthuridae [,m]<- con$ratio
# resTest <- testResiduals(simulateResiduals(mDOsp_Acanthuridae), plot = F)
# uniformityTestmSpD_OceanDepth_Acanthuridae [, m] <- as.vector(c(resTest$uniformity$statistic, resTest$uniformity$p.value))
# dispersionTestmSpD_OceanDepth_Acanthuridae [, m] <- as.vector(c(resTest$dispersion$statistic, resTest$dispersion$p.value))
# outlierTestmSpD_OceanDepth_Acanthuridae [, m] <- as.vector(c(resTest$outliers$statistic, resTest$outliers$p.value))

#Model 24 ####
#Oceans - Chaetodontidae
hist(Chaetodontidae$Nsp_Fam)
mDOsp_Chaetodontidae <-  glmmTMB(Nsp_Fam ~ Ocean*Depth+offset(log(area))+(1|Location),ziformula=~1,family=poisson, data = Chaetodontidae)

#plotResiduals(mDOsp_Chaetodontidae)
#R2    
caca<- as.data.frame(performance::r2(mDOsp_Chaetodontidae)[2])
mDOsp_Chaetodontidae_R2 [, m] <- caca[1,1]

#Plot Values     
plotVals <- emmip(mDOsp_Chaetodontidae, Ocean~Depth, at=list(Depth = DepthSeqScaled), type = "response",plotit = FALSE)%>%
  mutate(spRic = yvar, Iter = m)
spOplotVals_Chaetodontidae [,m]<- plotVals$yvar
# Mean values and contrasts at 5m, 30m, 60m, 90m, and 120m
means <- emmip(mDOsp_Chaetodontidae, Ocean~Depth, at=list(Depth = DepthSliceScaled), type = "response",plotit = FALSE)
spOmeans_Chaetodontidae [,m]<- means$yvar
ref <-ref_grid(mDOsp_Chaetodontidae, at=list(Depth = DepthSliceScaled))
spOCont <- emmeans(ref, ~  Ocean|Depth, type = "response")
con <- as.data.frame(contrast(spOCont))%>%filter(contrast == "Atlantic effect")
spOContrasts_Chaetodontidae [,m]<- con$ratio
# resTest <- testResiduals(simulateResiduals(mDOsp_Chaetodontidae), plot = F)
# uniformityTestmSpD_OceanDepth_Chaetodontidae [, m] <- as.vector(c(resTest$uniformity$statistic, resTest$uniformity$p.value))
# dispersionTestmSpD_OceanDepth_Chaetodontidae[, m] <- as.vector(c(resTest$dispersion$statistic, resTest$dispersion$p.value))
# outlierTestmSpD_OceanDepth_Chaetodontidae[, m] <- as.vector(c(resTest$outliers$statistic, resTest$outliers$p.value))

#Model 25 #### 
#Oceans - Pomacentridae
hist(Pomacentridae$Nsp_Fam)
mDOsp_Pomacentridae <- glmmTMB(Nsp_Fam ~ Ocean*Depth+offset(log(area))+(1|Location),ziformula=~1,family=poisson, data = Pomacentridae)

#R2    
caca<- as.data.frame(performance::r2(mDOsp_Pomacentridae)[2])
mDOsp_Pomacentridae_R2 [, m] <- caca[1,1]

#Plot Values     
plotVals <- emmip(mDOsp_Pomacentridae, Ocean~Depth, at=list(Depth = DepthSeqScaled), type = "response",plotit = FALSE)%>%
  mutate(spRic = yvar, Iter = m)
spOplotVals_Pomacentridae [,m]<- plotVals$yvar
# Mean values and contrasts at 5m, 30m, 60m, 90m, and 120m
means <- emmip(mDOsp_Pomacentridae, Ocean~Depth, at=list(Depth = DepthSliceScaled), type = "response",plotit = FALSE)
spOmeans_Pomacentridae [,m]<- means$yvar
ref <-ref_grid(mDOsp_Pomacentridae, at=list(Depth = DepthSliceScaled))
spOCont <- emmeans(ref, ~  Ocean|Depth, type = "response")
con <- as.data.frame(contrast(spOCont))%>%filter(contrast == "Atlantic effect")
spOContrasts_Pomacentridae [,m]<- con$ratio
# resTest <- testResiduals(simulateResiduals(mDOsp_Pomacentridae), plot = F)
# uniformityTestmSpD_OceanDepth_Pomacentridae [, m] <- as.vector(c(resTest$uniformity$statistic, resTest$uniformity$p.value))
# dispersionTestmSpD_OceanDepth_Pomacentridae[, m] <- as.vector(c(resTest$dispersion$statistic, resTest$dispersion$p.value))
# outlierTestmSpD_OceanDepth_Pomacentridae[, m] <- as.vector(c(resTest$outliers$statistic, resTest$outliers$p.value))

#Model 26 #### 
#Oceans - Scaridae
hist(Scaridae$Nsp_Fam)
mDOsp_Scaridae <- glmmTMB(Nsp_Fam ~ Ocean*Depth+offset(log(area))+(1|Location),ziformula=~1,family=poisson, data = Scaridae)

#R2    
caca<- as.data.frame(performance::r2(mDOsp_Scaridae)[2])
mDOsp_Scaridae_R2 [, m] <- caca[1,1]

#Plot Values     
plotVals <- emmip(mDOsp_Scaridae, Ocean~Depth, at=list(Depth = DepthSeqScaled), type = "response",plotit = FALSE)%>%
  mutate(spRic = yvar, Iter = m)
spOplotVals_Scaridae [,m]<- plotVals$yvar
# Mean values and contrasts at 5m, 30m, 60m, 90m, and 120m
means <- emmip(mDOsp_Scaridae, Ocean~Depth, at=list(Depth = DepthSliceScaled), type = "response",plotit = FALSE)
spOmeans_Scaridae [,m]<- means$yvar
ref <-ref_grid(mDOsp_Scaridae, at=list(Depth = DepthSliceScaled))
spOCont <- emmeans(ref, ~  Ocean|Depth, type = "response")
con <- as.data.frame(contrast(spOCont))%>%filter(contrast == "Atlantic effect")
spOContrasts_Scaridae [,m]<- con$ratio
# resTest <- testResiduals(simulateResiduals(mDOsp_Scaridae), plot = F)
# uniformityTestmSpD_OceanDepth_Scaridae [, m] <- as.vector(c(resTest$uniformity$statistic, resTest$uniformity$p.value))
# dispersionTestmSpD_OceanDepth_Scaridae[, m] <- as.vector(c(resTest$dispersion$statistic, resTest$dispersion$p.value))
# outlierTestmSpD_OceanDepth_Scaridae[, m] <- as.vector(c(resTest$outliers$statistic, resTest$outliers$p.value))

#Model 27 #### 
#Oceans - Serranidae
hist(Serranidae$Nsp_Fam)
mDOsp_Serranidae <- glmmTMB(Nsp_Fam ~ Ocean*Depth+offset(log(area))+(1|Location),ziformula=~1,family=poisson, data = Serranidae)

#plotResiduals(mDOsp_Serranidae)
#R2    
caca<- as.data.frame(performance::r2(mDOsp_Serranidae)[2])
mDOsp_Serranidae_R2 [, m] <- caca[1,1]
#Plot Values     
plotVals <- emmip(mDOsp_Serranidae, Ocean~Depth, at=list(Depth = DepthSeqScaled), type = "response",plotit = FALSE)%>%
  mutate(spRic = yvar, Iter = m)
spOplotVals_Serranidae [,m]<- plotVals$yvar
# Mean values and contrasts at 5m, 30m, 60m, 90m, and 120m
means <- emmip(mDOsp_Serranidae, Ocean~Depth, at=list(Depth = DepthSliceScaled), type = "response",plotit = FALSE)
spOmeans_Serranidae [,m]<- means$yvar
ref <-ref_grid(mDOsp_Serranidae, at=list(Depth = DepthSliceScaled))
spOCont <- emmeans(ref, ~  Ocean|Depth, type = "response")
con <- as.data.frame(contrast(spOCont))%>%filter(contrast == "Atlantic effect")
spOContrasts_Serranidae [,m]<- con$ratio
# resTest <- testResiduals(simulateResiduals(mDOsp_Serranidae), plot = F)
# uniformityTestmSpD_OceanDepth_Serranidae [, m] <- as.vector(c(resTest$uniformity$statistic, resTest$uniformity$p.value))
# dispersionTestmSpD_OceanDepth_Serranidae[, m] <- as.vector(c(resTest$dispersion$statistic, resTest$dispersion$p.value))
# outlierTestmSpD_OceanDepth_Serranidae[, m] <- as.vector(c(resTest$outliers$statistic, resTest$outliers$p.value))

#Models 28 - 34, trophic guilds ####

## Species richness by Trophic guilds

TG <- ddply (database_red,. (ID_transect, DietLevel), summarise,
             Total_SpR_TG=length(unique(Name)),
             Total_abd_TG=sum(Abundance))
TG_SpR <- cast (DietLevel~ID_transect, value = "Total_SpR_TG", data=TG)
TG_SpR[is.na(TG_SpR)] <- 0
TG_SpR <- melt (TG_SpR, id=c("DietLevel"))
rownames (TG_SpR) <- NULL
colnames (TG_SpR)[2] <- "ID_transect"
colnames (TG_SpR)[3] <- "Total_SpR_TG"

TG_SpR <- left_join(TG_SpR, LR)
TG_SpR$Prop_TG <- TG_SpR$Total_SpR_TG/TG_SpR$Nsp
TG_SpR <- left_join(TG_SpR, LRdat)%>%select(-CMPLX)%>%drop_na()

OM <- droplevels (subset (TG_SpR, DietLevel=="OM"))
HD <- droplevels (subset (TG_SpR, DietLevel=="HD"))
PK <- droplevels (subset (TG_SpR, DietLevel=="PK"))
IM <- droplevels (subset (TG_SpR, DietLevel=="IM"))
IS <- droplevels (subset (TG_SpR, DietLevel=="IS"))
FC <- droplevels (subset (TG_SpR, DietLevel=="FC"))
HM <- droplevels (subset (TG_SpR, DietLevel=="HM"))

#Model 28 ####  
#Oceans - OM
hist(OM$Nsp_TG)
mDOsp_OM <- glmmTMB(Nsp_TG ~ Ocean*Depth+offset(log(area))+(1|Location),ziformula=~1,family=poisson, data = OM)

#R2    
caca<- as.data.frame(performance::r2(mDOsp_OM)[2])
mDOsp_OM_R2 [, m] <- caca[1,1]

#Plot Values     
plotVals <- emmip(mDOsp_OM, Ocean~Depth, at=list(Depth = DepthSeqScaled), type = "response", plotit = F)%>%
  mutate(spRic = yvar, Iter = m)
spOplotVals_OM [,m]<- plotVals$yvar
# Mean values and contrasts at 5m, 30m, 60m, 90m, and 120m
means <- emmip(mDOsp_OM, Ocean~Depth, at=list(Depth = DepthSliceScaled), type = "response", plotit = FALSE)
spOmeans_OM [,m]<- means$yvar
ref <-ref_grid(mDOsp_OM, at=list(Depth = DepthSliceScaled))
spOCont <- emmeans(ref, ~  Ocean|Depth, type = "response")
con <- as.data.frame(contrast(spOCont))%>%filter(contrast == "Atlantic effect")
spOContrasts_OM [,m]<- con$ratio
# resTest <- testResiduals(simulateResiduals(mDOsp_OM), plot = F)
# uniformityTestmSpD_OceanDepth_OM [, m] <- as.vector(c(resTest$uniformity$statistic, resTest$uniformity$p.value))
# dispersionTestmSpD_OceanDepth_OM[, m] <- as.vector(c(resTest$dispersion$statistic, resTest$dispersion$p.value))
# outlierTestmSpD_OceanDepth_OM[, m] <- as.vector(c(resTest$outliers$statistic, resTest$outliers$p.value))

#Model 29 ####   
#Oceans - HD
hist(HD$Nsp_TG) 
mDOsp_HD <- glmmTMB(Nsp_TG ~ Ocean*Depth+offset(log(area))+(1|Location),ziformula=~1,family=poisson, data = HD)

#plotResiduals(mDOsp_HD)
#R2    
caca<- as.data.frame(performance::r2(mDOsp_HD)[2])
mDOsp_HD_R2 [, m] <- caca[1,1]

#Plot Values     
plotVals <- emmip(mDOsp_HD, Ocean~Depth, at=list(Depth = DepthSeqScaled), type = "response", plotit = FALSE)%>%
  mutate(spRic = yvar, Iter = m)
spOplotVals_HD [,m]<- plotVals$yvar
# Mean values and contrasts at 5m, 30m, 60m, 90m, and 120m
means <- emmip(mDOsp_HD, Ocean~Depth, at=list(Depth = DepthSliceScaled), type = "response", plotit = FALSE)
spOmeans_HD [,m]<- means$yvar
ref <-ref_grid(mDOsp_HD, at=list(Depth = DepthSliceScaled))
spOCont <- emmeans(ref, ~  Ocean|Depth, type = "response")
con <- as.data.frame(contrast(spOCont))%>%filter(contrast == "Atlantic effect")
spOContrasts_HD [,m]<- con$ratio
# resTest <- testResiduals(simulateResiduals(mDOsp_HD), plot = F)
# uniformityTestmSpD_OceanDepth_HD [, m] <- as.vector(c(resTest$uniformity$statistic, resTest$uniformity$p.value))
# dispersionTestmSpD_OceanDepth_HD[, m] <- as.vector(c(resTest$dispersion$statistic, resTest$dispersion$p.value))
# outlierTestmSpD_OceanDepth_HD[, m] <- as.vector(c(resTest$outliers$statistic, resTest$outliers$p.value))

#Model 30 ####  
#Oceans - PK
hist(PK$Nsp_TG)
mDOsp_PK <- glmer.nb(Nsp_TG ~ Ocean*Depth+offset(log(area))+(1|Location), data = PK)
#print(check_model(mDOsp_PK))
#plotResiduals(mDOsp_PK)
#R2
caca<- as.data.frame(performance::r2(mDOsp_PK)[2])
mDOsp_PK_R2 [, m] <- caca[1,1]

#Plot Values
plotVals <- emmip(mDOsp_PK, Ocean~Depth, at=list(Depth = DepthSeqScaled), type = "response", plotit = FALSE)%>%
  mutate(spRic = yvar, Iter = m)
spOplotVals_PK [,m]<- plotVals$yvar
# Mean values and contrasts at 5m, 30m, 60m, 90m, and 120m
means <- emmip(mDOsp_PK, Ocean~Depth, at=list(Depth = DepthSliceScaled), type = "response", plotit = FALSE)
spOmeans_PK [,m]<- means$yvar
ref <-ref_grid(mDOsp_PK, at=list(Depth = DepthSliceScaled))
spOCont <- emmeans(ref, ~  Ocean|Depth, type = "response")
con <- as.data.frame(contrast(spOCont))%>%filter(contrast == "Atlantic effect")
spOContrasts_PK [,m]<- con$ratio
# resTest <- testResiduals(simulateResiduals(mDOsp_PK), plot = F)
# uniformityTestmSpD_OceanDepth_PK [, m] <- as.vector(c(resTest$uniformity$statistic, resTest$uniformity$p.value))
# dispersionTestmSpD_OceanDepth_PK[, m] <- as.vector(c(resTest$dispersion$statistic, resTest$dispersion$p.value))
# outlierTestmSpD_OceanDepth_PK[, m] <- as.vector(c(resTest$outliers$statistic, resTest$outliers$p.value))

#Model 31 #### 
#Oceans - IM
hist(IM$Nsp_TG)
mDOsp_IM <- glmer.nb(Nsp_TG ~ Ocean*Depth+offset(log(area))+(1|Location),  data = IM)
#plotResiduals(mDOsp_IM)
#R2
caca<- as.data.frame(performance::r2(mDOsp_IM)[2])
mDOsp_IM_R2 [, m] <- caca[1,1]

#Plot Values
plotVals <- emmip(mDOsp_IM, Ocean~Depth, at=list(Depth = DepthSeqScaled), type = "response", plotit = FALSE)%>%
  mutate(spRic = yvar, Iter = m)
spOplotVals_IM [,m]<- plotVals$yvar
# Mean values and contrasts at 5m, 30m, 60m, 90m, and 120m
means <- emmip(mDOsp_IM, Ocean~Depth, at=list(Depth = DepthSliceScaled), type = "response", plotit = FALSE)
spOmeans_IM[,m]<- means$yvar
ref <-ref_grid(mDOsp_IM, at=list(Depth = DepthSliceScaled))
spOCont <- emmeans(ref, ~  Ocean|Depth, type = "response")
con <- as.data.frame(contrast(spOCont))%>%filter(contrast == "Atlantic effect")
spOContrasts_IM [,m]<- con$ratio
# resTest <- testResiduals(simulateResiduals(mDOsp_IM), plot = F)
# uniformityTestmSpD_OceanDepth_IM [, m] <- as.vector(c(resTest$uniformity$statistic, resTest$uniformity$p.value))
# dispersionTestmSpD_OceanDepth_IM[, m] <- as.vector(c(resTest$dispersion$statistic, resTest$dispersion$p.value))
# outlierTestmSpD_OceanDepth_IM[, m] <- as.vector(c(resTest$outliers$statistic, resTest$outliers$p.value))


#Model 32 ####  
#Oceans - IS
hist(IS$Nsp_TG)
mDOsp_IS <- glmmTMB(Nsp_TG ~ Ocean*Depth+offset(log(area))+(1|Location), ziformula=~1,family=poisson, data = IS)


#R2    
caca<- as.data.frame(performance::r2(mDOsp_IS)[2])
mDOsp_IS_R2 [, m] <- caca[1,1]

#Plot Values     
plotVals <- emmip(mDOsp_IS, Ocean~Depth, at=list(Depth = DepthSeqScaled), type = "response", plotit = F)%>%
  mutate(spRic = yvar, Iter = m)
spOplotVals_IS [,m]<- plotVals$yvar
# Mean values and contrasts at 5m, 30m, 60m, 90m, and 120m
means <- emmip(mDOsp_IS, Ocean~Depth, at=list(Depth = DepthSliceScaled), type = "response", plotit = FALSE)
spOmeans_IS [,m]<- means$yvar
ref <-ref_grid(mDOsp_IS, at=list(Depth = DepthSliceScaled))
spOCont <- emmeans(ref, ~  Ocean|Depth, type = "response")
con <- as.data.frame(contrast(spOCont))%>%filter(contrast == "Atlantic effect")
spOContrasts_IS [,m]<- con$ratio
# resTest <- testResiduals(simulateResiduals(mDOsp_IS), plot = F)
# uniformityTestmSpD_OceanDepth_IS [, m] <- as.vector(c(resTest$uniformity$statistic, resTest$uniformity$p.value))
# dispersionTestmSpD_OceanDepth_IS[, m] <- as.vector(c(resTest$dispersion$statistic, resTest$dispersion$p.value))
# outlierTestmSpD_OceanDepth_IS[, m] <- as.vector(c(resTest$outliers$statistic, resTest$outliers$p.value))

#Model 33 ####  
#Oceans - FC
hist(FC$Nsp_TG)
mDOsp_FC <- glmmTMB(Nsp_TG ~ Ocean*Depth+offset(log(area))+(1|Location),ziformula=~1,family=poisson, data = FC)

#R2    
caca<- as.data.frame(performance::r2(mDOsp_FC)[2])
mDOsp_FC_R2 [, m] <- caca[1,1]

#Plot Values     
plotVals <- emmip(mDOsp_FC, Ocean~Depth, at=list(Depth = DepthSeqScaled), type = "response", plotit = FALSE)%>%
  mutate(spRic = yvar, Iter = m)
spOplotVals_FC [,m]<- plotVals$yvar
# Mean values and contrasts at 5m, 30m, 60m, 90m, and 120m
means <- emmip(mDOsp_FC, Ocean~Depth, at=list(Depth = DepthSliceScaled), type = "response", plotit = FALSE)
spOmeans_FC [,m]<- means$yvar
ref <-ref_grid(mDOsp_FC, at=list(Depth = DepthSliceScaled))
spOCont <- emmeans(ref, ~  Ocean|Depth, type = "response")
con <- as.data.frame(contrast(spOCont))%>%filter(contrast == "Atlantic effect")
spOContrasts_FC [,m]<- con$ratio
# resTest <- testResiduals(simulateResiduals(mDOsp_FC), plot = F)
# uniformityTestmSpD_OceanDepth_FC [, m] <- as.vector(c(resTest$uniformity$statistic, resTest$uniformity$p.value))
# dispersionTestmSpD_OceanDepth_FC[, m] <- as.vector(c(resTest$dispersion$statistic, resTest$dispersion$p.value))
# outlierTestmSpD_OceanDepth_FC[, m] <- as.vector(c(resTest$outliers$statistic, resTest$outliers$p.value))

#Model 34 ####    
#Oceans - HM
hist(HM$Nsp_TG)

mDOsp_HM <- glmmTMB(Nsp_TG ~ Ocean*Depth+offset(log(area))+(1|Location),ziformula=~1,family=poisson, data = HM)


#R2    
caca<- as.data.frame(performance::r2(mDOsp_HM)[2])
mDOsp_HM_R2 [, m] <- caca[1,1]

#Plot Values     
plotVals <- emmip(mDOsp_HM, Ocean~Depth, at=list(Depth = DepthSeqScaled), type = "response", plotit = FALSE)%>%
  mutate(spRic = yvar, Iter = m)
spOplotVals_HM [,m]<- plotVals$yvar
# Mean values and contrasts at 5m, 30m, 60m, 90m, and 120m
means <- emmip(mDOsp_HM, Ocean~Depth, at=list(Depth = DepthSliceScaled), type = "response", plotit = FALSE)
spOmeans_HM [,m]<- means$yvar
ref <-ref_grid(mDOsp_HM, at=list(Depth = DepthSliceScaled))
spOCont <- emmeans(ref, ~  Ocean|Depth, type = "response")
con <- as.data.frame(contrast(spOCont))%>%filter(contrast == "Atlantic effect")
spOContrasts_HM [,m]<- con$ratio

# resTest <- testResiduals(simulateResiduals(mDOsp_HM), plot = F)#plot = T to see plots
# uniformityTestmSpD_OceanDepth_HM [, m] <- as.vector(c(resTest$uniformity$statistic, resTest$uniformity$p.value))
# dispersionTestmSpD_OceanDepth_HM[, m] <- as.vector(c(resTest$dispersion$statistic, resTest$dispersion$p.value))
# outlierTestmSpD_OceanDepth_HM[, m] <- as.vector(c(resTest$outliers$statistic, resTest$outliers$p.value))
}



#Plot species richness by depth and ocean ####
spOplotDat <- Filter(function(x)!all(is.na(x)), as.data.frame(spOplotVals))
spOplotDat <- bind_cols(mDOsp_info, spOplotDat)%>%gather("Iter", "SpRic", -Ocean, -Depth)
ggplot(spOplotDat, aes(Depth*dScle+dCent, SpRic, colour = Ocean, group = interaction(Ocean, Iter)))+
  geom_line(alpha= .3)+
  scale_colour_manual(values = c("royalblue4", "steelblue2"))+ #darkslategrey
  theme_classic()+
  ylab(expression(paste("Species richness .transec", t^-1," (area offset)")))+
  scale_x_continuous(limits = c(0,130), breaks = c(0, 5, 30, 60, 90, 120))+
  xlab("Depth (m)")+
  theme(legend.position = 'none')
#Mean and and contrasts
#Means
OceanSpRicMeans <- as.data.frame(spOmeans)%>%
  mutate(Ocean = rep(c("Atlantic", "Pacific"),times = 5), Depth = rep(c(5, 30, 60, 90, 120), each = 2))%>%
  gather(key = "Iteration", value = "SpRic", -Ocean, -Depth)%>%
  mutate(SpRic = exp(SpRic))%>%
  group_by(Ocean, Depth) %>%
  dplyr::summarise(BioGroup = "Assembelage",
                   Geolevel = "Ocean",
                   Mean = mean(SpRic, na.rm = T), 
                   LCL95 = quantile(SpRic, .025, na.rm = T), UCL95 = quantile(SpRic, .975, na.rm = T), 
                   LCL80 = quantile(SpRic, .1, na.rm = T), UCL80 = quantile(SpRic, .9, na.rm = T), 
                   LCL60 = quantile(SpRic, .2, na.rm = T), UCL60 = quantile(SpRic, .8, na.rm = T))
#Contrasts
OceanSpRicContrasts <- as.data.frame(spOContrasts)%>%
  mutate( Depth = c(5, 30, 60, 90, 120))%>%
  gather(key = "Iteration", value = "estimate",-Depth)%>%
  group_by(Depth) %>%
  dplyr::summarise(BioGroup = "Assembelage",
                   Geolevel = "Ocean",
                   MeanContrast = mean(estimate, na.rm = T), 
                   LCL95 = quantile(estimate, .025, na.rm = T), UCL95 = quantile(estimate, .975, na.rm = T), 
                   LCL80 = quantile(estimate, .1, na.rm = T), UCL80 = quantile(estimate, .9, na.rm = T), 
                   LCL60 = quantile(estimate, .2, na.rm = T), UCL60 = quantile(estimate, .8, na.rm = T))


#Example code for plotting species richness of families  - adapt for trophic groups
spOplotDat_Labridae <- Filter(function(x)!all(is.na(x)), as.data.frame(spOplotVals_Labridae))
spOplotDat_Labridae <- bind_cols(mDOsp_info, spOplotDat_Labridae)%>%gather("Iter", "SpRic", -Ocean, -Depth)%>%
  mutate(BioGroup = "Family", BioLevel = "Labridae")

spOplotDat_Pomacentridae <- Filter(function(x)!all(is.na(x)), as.data.frame(spOplotVals_Pomacentridae))
spOplotDat_Pomacentridae <- bind_cols(mDOsp_info, spOplotDat_Pomacentridae)%>%gather("Iter", "SpRic", -Ocean, -Depth)%>%
  mutate(BioGroup = "Family", BioLevel = "Pomacentridae")

spOplotDat_Acanthuridae <- Filter(function(x)!all(is.na(x)), as.data.frame(spOplotVals_Acanthuridae))
spOplotDat_Acanthuridae <- bind_cols(mDOsp_info, spOplotDat_Acanthuridae)%>%gather("Iter", "SpRic", -Ocean, -Depth)%>%
  mutate(BioGroup = "Family", BioLevel = "Acanthuridae")

spOplotDat_Scaridae <- Filter(function(x)!all(is.na(x)), as.data.frame(spOplotVals_Scaridae))
spOplotDat_Scaridae <- bind_cols(mDOsp_info, spOplotDat_Scaridae)%>%gather("Iter", "SpRic", -Ocean, -Depth)%>%
  mutate(BioGroup = "Family", BioLevel = "Scaridae")

spOplotDat_Serranidae <- Filter(function(x)!all(is.na(x)), as.data.frame(spOplotVals_Serranidae))
spOplotDat_Serranidae <- bind_cols(mDOsp_info, spOplotDat_Serranidae)%>%gather("Iter", "SpRic", -Ocean, -Depth)%>%
  mutate(BioGroup = "Family", BioLevel = "Serranidae")

spOplotDat_Chaetodontidae <- Filter(function(x)!all(is.na(x)), as.data.frame(spOplotVals_Chaetodontidae))
spOplotDat_Chaetodontidae <- bind_cols(mDOsp_info, spOplotDat_Chaetodontidae)%>%gather("Iter", "SpRic", -Ocean, -Depth)%>%
  mutate(BioGroup = "Family", BioLevel = "Chaetodontidae")


spOplotDat_Families <- bind_rows(spOplotDat_Chaetodontidae, spOplotDat_Serranidae, spOplotDat_Scaridae, spOplotDat_Acanthuridae, spOplotDat_Labridae, spOplotDat_Pomacentridae)%>%
  mutate(Ocean = factor(Ocean, levels = c("Pacific", "Atlantic")))

ggplot(spOplotDat_Families, aes(x=(Depth*dScle+dCent), SpRic, colour = BioLevel, group = interaction(BioLevel, Iter)))+
  geom_line(alpha= 0.1)+
  #scale_colour_brewer(palette = "Spectral") +
  scale_colour_manual(name='Family', values=c("#56B4E9", "#E69F00", "#D55E00","#CC79A7", "#F0E442","#0072B2","#009E73"))+
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_classic()+
  facet_wrap(~Ocean)+
  ylim(0,8)+
  scale_x_continuous(limits = c(0,130), breaks = c(0, 5, 30, 60, 90, 120))+
  labs(x ="Depth (m)", y = expression(paste("Species richness .transec", t^-1," (area offset)")))

#example code for means and contrasts in families and trophic groups
#Labridae - adapt for other families and trophic groups
#Means - Labridae
OceanSpRicMeans_Labridae <- as.data.frame(spOmeans_Labridae)%>%
  mutate(Geolevel = rep(c("Atlantic", "Pacific"),times = 5), Depth = rep(c(5, 30, 60, 90, 120), each = 2))%>%
  gather(key = "Iteration", value = "SpRic", -Geolevel, -Depth)%>%
  group_by(Geolevel, Depth) %>%
  dplyr::summarise(BioGroup = "Family", 
                   Biolevel = "Labridae",
                   Test = "Mean",
                   Estimate = mean(SpRic, na.rm = T), 
                   LCL95 = quantile(SpRic, .025, na.rm = T), UCL95 = quantile(SpRic, .975, na.rm = T), 
                   LCL80 = quantile(SpRic, .1, na.rm = T), UCL80 = quantile(SpRic, .9, na.rm = T), 
                   LCL60 = quantile(SpRic, .2, na.rm = T), UCL60 = quantile(SpRic, .8, na.rm = T))%>%
  dplyr::select(Biolevel, Test, Geolevel, Depth, Estimate, LCL95, UCL95)

#Contrasts Labridae
OceanSpRicContrasts_Labridae <- as.data.frame(spOContrasts_Labridae)%>%
  mutate( Depth = c(5, 30, 60, 90, 120))%>%
  gather(key = "Iteration", value = "estimate",-Depth)%>%
  group_by(Depth) %>%
  dplyr::summarise(BioGroup = "Family", 
                   Biolevel = "Labridae",
                   Geolevel = "Atlantic-Pacific",
                   Test = "Contrast ratio",
                   Estimate = mean(estimate, na.rm = T), 
                   LCL95 = quantile(estimate, .025, na.rm = T), UCL95 = quantile(estimate, .975, na.rm = T), 
                   LCL80 = quantile(estimate, .1, na.rm = T), UCL80 = quantile(estimate, .9, na.rm = T), 
                   LCL60 = quantile(estimate, .2, na.rm = T), UCL60 = quantile(estimate, .8, na.rm = T))%>%
  dplyr::select(Biolevel, Test, Geolevel, Depth, Estimate, LCL95, UCL95)

#Example code for R2 values forfamilies and trophic groups
#Labridae - adapt for other families and trophic groups

mDOsp_Labridae_R2_Means <- as.data.frame(t(mDOsp_Labridae_R2))%>%
  mutate(R2 = V1)%>%
  select(R2)%>%drop_na()%>% 
  dplyr::summarise(MeanR2 = mean(R2, na.rm = TRUE), 
                   LCL95R2 = quantile(R2, .025, na.rm = TRUE), UCL95Stat = quantile(R2, .975, na.rm = TRUE), 
                   LCL80R2 = quantile(R2, .1, na.rm = TRUE), UCL80Stat = quantile(R2, .9, na.rm = TRUE), 
                   LCL60R2 = quantile(R2, .2, na.rm = TRUE), UCL60Stat = quantile(R2, .8, na.rm = TRUE))
mDOsp_Labridae_R2_Means


#Models 35 - 37. Functional metrics####

traitCAS <- abund %>%
  dplyr::select(Name, ID_Locality_Mesophotic_Zone, SizeClass, HomeRange, Activity, Schooling, WaterColPosition, DietLevel) %>%
  unique()%>% 
  arrange(Name) %>% 
  remove_rownames() 

traitCAS_sh <- traitCAS %>% 
  filter(grepl("Shallow",ID_Locality_Mesophotic_Zone))%>%
  dplyr::select(-ID_Locality_Mesophotic_Zone)%>%
  distinct()

traitCAS_up <- traitCAS %>% 
  filter(grepl("Upper",ID_Locality_Mesophotic_Zone)) %>% 
  dplyr::select(-ID_Locality_Mesophotic_Zone)%>% 
  distinct() 
# %>% 
#   column_to_rownames("Name") 

traitCAS_lw <- traitCAS %>% 
  filter(grepl("Deeper",ID_Locality_Mesophotic_Zone)) %>% 
  dplyr::select(-ID_Locality_Mesophotic_Zone)%>% 
  distinct() 
# %>% 
# column_to_rownames("Name") 

traitCAS <- traitCAS %>% 
  dplyr::select(-ID_Locality_Mesophotic_Zone)%>% 
  distinct() 
# %>% 
#   column_to_rownames("Name")


traitCAS <- traitCAS %>% 
  mutate(Name = factor(Name),
         SizeClass = factor(SizeClass, order=TRUE,levels=c("S1", "S2", "S3", "S4", "S5", "S6")),
         HomeRange = factor(HomeRange), 
         Activity = factor(Activity),
         Schooling = factor(Schooling, order = T, levels = c("Sol", "Pair","SmallG", "MedG", "LargeG")),
         WaterColPosition = factor(WaterColPosition),
         DietLevel = factor(DietLevel)) %>% glimpse

traitCAS_all <-  traitCAS %>%column_to_rownames ("Name") %>% glimpse

traitCAS_cat <- as.data.frame(colnames(traitCAS_all)) %>% 
  dplyr::mutate(trait_name = colnames(traitCAS_all),
                trait_type = c("O", "N", "N", "O", "N", "N")) %>% 
  dplyr::select(trait_name, trait_type)



commCAS <- abund %>% 
  dplyr::select(ID_Locality_Mesophotic_Zone, Name, Abundance) %>% 
  group_by(Name, ID_Locality_Mesophotic_Zone) %>% 
  summarise(Abundance = as.numeric(sum(Abundance))) %>% 
  mutate(Abundance = if_else(is.na(Abundance), 0, Abundance)) %>% 
  spread(Name, Abundance, fill = 0) 

commCAS_sh <- commCAS %>% 
  filter(grepl("Shallow",ID_Locality_Mesophotic_Zone)) %>% 
  column_to_rownames("ID_Locality_Mesophotic_Zone") %>% 
  dplyr::select(where(~sum(.) != 0))

commCAS_up <- commCAS %>% 
  filter(grepl("Upper",ID_Locality_Mesophotic_Zone)) %>% 
  column_to_rownames("ID_Locality_Mesophotic_Zone")%>% 
  dplyr::select(where(~sum(.) != 0))

commCAS_lw <- commCAS %>% 
  filter(grepl("Deeper",ID_Locality_Mesophotic_Zone)) %>% 
  column_to_rownames("ID_Locality_Mesophotic_Zone")%>% 
  dplyr::select(where(~sum(.) != 0))

commCAS <- commCAS %>%  
  column_to_rownames("ID_Locality_Mesophotic_Zone")




abund$transect_id <- paste (abund$Location, abund$Mesophotic2, abund$ID_transect, sep = "_")
abund$trans <- 1
sub1 <- cast (Location + Mesophotic2 + ID_Locality_Mesophotic_Zone + area + transect_id ~ ., fun.aggregate = sum, value = "trans", data = abund)
sub1 <- droplevels (sub1[,-ncol(sub1)])
sub1$Mesophotic2 <- as.factor(sub1$Mesophotic2) 
# Estimated the minimal sampling area ---------------------------------------
minimal_sampling_area <- tapply (sub1$area, list(sub1$ID_Locality_Mesophotic_Zone), sum)
msa <- min (minimal_sampling_area)



# vol<-matrix(nrow = 3, ncol = N) 


gower_matrix <- daisy (traitCAS_all, metric = c("gower"))


gm <- gower_matrix %>% 
  as.matrix 
colnames(gm) <- traitCAS$Name
rownames(gm) <- colnames(gm) 


gm <- as.dist(gm) 


# function daisy compute all dissimilarities distances
pcoa <- dudi.pco (quasieuclid(gm), scannf = F, nf = 5) # here 6 axes are retained

#Start loop ----

upto = max(as.data.frame(FRV_m)$iter %>% as.numeric(as.character()), na.rm = T)+1
upto  # iteration = 1:1000 # set for null-model construction
N <- 10 #1000 used in manuscript
for (m in upto:N){
  chosen_transec  <- choosing_transec(sub1, msa)
  sub2 <- droplevels (sub1 [sub1$transect_id %in% chosen_transec, ])
  database_red <- droplevels (abund [abund$transect_id %in% sub2$transect_id, ])
  dim(database_red)
  length(levels(database_red$Name))
  
  # # subset traits and communities #
  traits1 <- traitCAS %>% 
    filter(Name %in% database_red$Name) %>%  
    column_to_rownames("Name")
  
  traits1_sh <- traitCAS_sh %>% 
    filter(Name %in% database_red$Name)%>%  
    column_to_rownames("Name")
  
  traits1_up <- traitCAS_up  %>% 
    filter(Name %in% database_red$Name)%>%  
    column_to_rownames("Name")

  traits1_lw <-traitCAS_lw  %>% 
    filter(Name %in% database_red$Name)%>%  
    column_to_rownames("Name")

  comm1 <- commCAS %>% 
    select_if(names(.) %in% database_red$Name) %>% 
    dplyr::select(order(colnames(.)))
  
  comm1_sh <- commCAS_sh %>% 
    select_if(names(.) %in% database_red$Name) %>% 
    dplyr::select(order(colnames(.)))
 
  comm1_up <- commCAS_up %>% 
    select_if(names(.) %in% database_red$Name) %>% 
    dplyr::select(order(colnames(.)))
 
  comm1_lw <- commCAS_lw %>% 
    select_if(names(.) %in% database_red$Name) %>% 
    dplyr::select(order(colnames(.)))
  
  coords <-  as.data.frame(pcoa$li[, 1:4]) %>% 
    mutate(Species = rownames(pcoa$li)) %>% 
    dplyr::filter(Species %in% colnames(comm1)) %>% 
    dplyr::select(-Species) %>% as.matrix()
  
  
  FD_Index <- multidimFD(coords, as.matrix(comm1)) # coordinates of species in the 6 multi-space and the occurrence of species
  FD_Index <- as.data.frame(FD_Index) %>% 
    mutate(iter = m,
           Site = rep(rownames(comm1))) 
  FD_Index%>% dim
  
  if(m == 1){
    FD_Index_m <- matrix(nrow = 34000, ncol = 26)
    
    colnames(FD_Index_m) <- colnames(FD_Index)
  }
  
  
  FD_Index_m [c((1:34)+(m-1)*34),c(1:26)] <- as.matrix(FD_Index)
  
#Functional redundancy 
  
  # create presence-absence for beta FD 
  occ <- as.matrix(comm1)
  occ[occ > 0] <- 1 
  
  FEtraits <- as.data.frame(traits1)
  colnames(FEtraits) <- c("HR","AC","GS","WP","SC","DI")
  FEtraits$HR <- as.factor(FEtraits$HR) 
  FEtraits$AC <- as.factor(FEtraits$AC) 
  FEtraits$GS <- as.factor(FEtraits$GS) 
  FEtraits$WP <- as.factor(FEtraits$WP)
  FEtraits$SC <- as.factor(FEtraits$SC)
  FEtraits$DI <- as.factor(FEtraits$DI)
  
  FEnts <- species_to_FE(FEtraits)
  FRV <- FE_metrics(FEnts, occ) 
  FRV <- FRV %>% as.data.frame() %>% 
    mutate(iter = m,
           Site = rep(rownames(comm1))) %>% 
    as.matrix()
  
  if(m == 1){
    FRV_m <- matrix(nrow = 34000, ncol = 7)
    
    colnames(FRV_m) <- colnames(FRV)
  }
  
  
  FRV_m [c((1:34)+(m-1)*34),c(1:7)] <- FRV
  
  #is_whole function
  is_whole <- function(x) all(floor(x) == x)
  #save every 5 iterations
  if(is_whole(m/5)){
    save.image(paste0("FD", m, "_iter_",Sys.Date(),".RData"))
  }
  print(paste0("Completed, iteration n. ", m))
}

#Matrix to data frame; split location and dpeth zone; numerise####
FD_index <- FD_Index_m %>% as.data.frame %>% drop_na 

FD_index <- FD_index%>% 
  mutate(Location = as.factor(stri_match_first_regex(FD_index$Site, '(\\w+)_(\\w+)')[,2]),
         Zone = as.factor(stri_match_first_regex(FD_index$Site, '(\\w+)_(\\w+)')[,3]),
         FRic = as.numeric(as.character(FRic)),
         FSpe = as.numeric(as.character(FSpe)),
         iter = as.numeric(as.character(iter)),
         FIde_A1 = as.numeric(FIde_A1), 
         FIde_A2 = as.numeric(FIde_A2), 
         Zone = fct_relevel(Zone, "Deeper","Upper", "Shallow"),
         Ocean = as.factor(if_else(Location == c("Bermuda")|
                                     Location == c("Curacao")| 
                                     Location == c("Fernando_de_Noronha")| 
                                     Location == c("St_Pauls_Rocks"), 
                                   "Atlantic",
                                   "Pacific")))  %>% glimpse

FD_index$Zone <- fct_relevel(FD_index$Zone, "Deeper","Upper", "Shallow")#order mesozones by depth

FD_index$FRic %>% summary()

#Fric
for(i in 1:max(FD_index$iter)){
  mFR <- lmer(FRic~Zone*Ocean+(1|Location), data = FD_index %>% filter(iter == i))
  mFR.rg <- ref_grid(mFR) #create emmeans ref grid
  mFR.emm.s <- emmeans(mFR.rg, "Zone", "Ocean") #get emmeans
  effFR <- eff_size(mFR.emm.s, sigma = sigma(mFR), edf = 42) %>% as.data.frame()#store pairwise effects from model
  caca<-performance::r2(mFR)
  
  if(i == 1){
    FR_estimates_m <- matrix(ncol = 7, nrow = 6000)
    FR_mod_r2_m <- matrix(ncol = 1, nrow = 1000)
    FR_contrasts_m <- matrix(ncol = 7, nrow = 6000)
    
    colnames(FR_estimates_m) <- colnames(mFR.emm.s %>% as.data.frame)
    colnames(FR_contrasts_m) <- colnames(effFR)
  }
  
  FR_estimates_m[c((1:6)+(i-1)*6),c(1:7)] <- mFR.emm.s %>% as.data.frame %>% as.matrix
  FR_contrasts_m[c((1:6)+(i-1)*6),c(1:7)]  <- effFR %>% as.matrix
  FR_mod_r2_m [i,] <- as.numeric(caca$R2_marginal)
  
  
  #FSpecialisation  
  
  mFSpec <- lmer(FSpe~Zone*Ocean+(1|Location), data = FD_index %>% filter(iter == i))
  mFSpec.rg <- ref_grid(mFSpec) #create emmeans ref grid
  mFSpec.emm.s <- emmeans(mFSpec.rg, "Zone", "Ocean") #get emmeans
  effFSpec <- eff_size(mFSpec.emm.s, sigma = sigma(mFSpec), edf = 42) %>% as.data.frame()#store pairwise effects from model
  caca<-performance::r2(mFSpec)
  
  if(i == 1){
    FSpec_estimates_m <- matrix(ncol = 7, nrow = 6000)
    FSpec_mod_r2_m <- matrix(ncol = 1, nrow = 1000)
    FSpec_contrasts_m <- matrix(ncol = 7, nrow = 6000)
    
    colnames(FSpec_estimates_m) <- colnames(mFSpec.emm.s %>% as.data.frame)
    colnames(FSpec_contrasts_m) <- colnames(effFSpec)
  }
  
  FSpec_estimates_m[c((1:6)+(i-1)*6),c(1:7)] <- mFSpec.emm.s %>% as.data.frame %>% as.matrix
  FSpec_contrasts_m[c((1:6)+(i-1)*6),c(1:7)]  <- effFSpec %>% as.matrix
  FSpec_mod_r2_m [i,] <- as.numeric(caca$R2_marginal)
}


#calculate and plot quantiles for FUNCTIONAL RICHNESS 
#calculate and plot quantiles for effects

FR_estimates <- FR_estimates_m %>% as.data.frame %>%
  mutate(emmean = as.numeric(as.character(emmean)),
         Zone = fct_relevel(Zone, c("Deeper", "Upper", "Shallow"))) %>% drop_na

FR_estimates_quantiles <- FR_estimates %>% 
  group_by(Zone, Ocean) %>% 
  summarise(value = quantile(emmean, probs = c(0.025, 0.075, 0.2, .8, .925, .975)),
            quants = c("L95", "L85", "L60", "U60", "U85", "U95")) %>% 
  spread(quants, value)

est_Fric <- ggplot(FR_estimates_quantiles, aes(Zone, colour = Ocean))+
  geom_linerange(aes(ymin = L95, ymax = U95), size = .5)+
  geom_linerange(aes(ymin = L85, ymax = U85), size = .75)+
  geom_linerange(aes(ymin = L60, ymax = U60), size = 1)+
  coord_flip()+
  theme_classic()+
  labs(title = "Richness")


#Calculate and plot quantiles for contrasts 
FR_contrasts <- FR_contrasts_m %>% as.data.frame %>% 
  mutate(effect.size = as.numeric(as.character(effect.size)),
         SE = as.numeric(as.character(SE)),
         df = as.numeric(as.character(df)),
         lower.CL = as.numeric(as.character(lower.CL)),
         upper.CL = as.numeric(as.character(upper.CL))) %>% 
  drop_na

FR_contrasts_quantiles <- FR_contrasts %>% 
  group_by(contrast, Ocean) %>% 
  summarise(value = quantile(effect.size, probs = c(0.025, 0.075, 0.2, .8, .925, .975) ),
            quants = c("L95", "L85", "L60", "U60", "U85", "U95")) %>% 
  spread(quants, value)

contr_FRic <- ggplot(FR_contrasts_quantiles, aes(contrast, colour = Ocean))+
  geom_linerange(aes(ymin = L95, ymax = U95), size = .5)+
  geom_linerange(aes(ymin = L85, ymax = U85), size = .75)+
  geom_linerange(aes(ymin = L60, ymax = U60), size = 1)+
  geom_hline(yintercept = 0)+
  coord_flip()+
  theme_classic()+
  labs(title = "Richness")


#calculate and plot quantiles for FUNCTIONAL Specialization 
#calculate and plot quantiles for effects

FSpec_estimates <- FSpec_estimates_m %>% as.data.frame %>%
  mutate(emmean = as.numeric(as.character(emmean)),
         Zone = fct_relevel(Zone, c("Deeper", "Upper", "Shallow"))) %>% drop_na()

FSpec_estimates_quantiles <- FSpec_estimates %>% 
  group_by(Zone, Ocean) %>% 
  summarise(value = quantile(emmean, probs = c(0.025, 0.075, 0.2, .8, .925, .975)),
            quants = c("L95", "L85", "L60", "U60", "U85", "U95")) %>% 
  spread(quants, value)

est_FSpec <- ggplot(FSpec_estimates_quantiles, aes(Zone, colour = Ocean))+
  geom_linerange(aes(ymin = L95, ymax = U95), size = .5)+
  geom_linerange(aes(ymin = L85, ymax = U85), size = .75)+
  geom_linerange(aes(ymin = L60, ymax = U60), size = 1)+
  coord_flip()+
  theme_classic()+
  labs(title = "Specialization")


#Calculate and plot quantiles for contrasts 
FSpec_contrasts <- FSpec_contrasts_m %>% as.data.frame %>% 
  mutate(effect.size = as.numeric(as.character(effect.size)),
         SE = as.numeric(as.character(SE)),
         df = as.numeric(as.character(df)),
         lower.CL = as.numeric(as.character(lower.CL)),
         upper.CL = as.numeric(as.character(upper.CL))) %>% drop_na

FSpec_contrasts_quantiles <- FSpec_contrasts %>% 
  group_by(contrast, Ocean) %>% 
  summarise(value = quantile(effect.size, probs = c(0.025, 0.075, 0.2, .8, .925, .975) ),
            quants = c("L95", "L85", "L60", "U60", "U85", "U95")) %>% 
  spread(quants, value)

contr_FSpec <- ggplot(FSpec_contrasts_quantiles, aes(contrast, colour = Ocean))+
  geom_linerange(aes(ymin = L95, ymax = U95), size = .5)+
  geom_linerange(aes(ymin = L85, ymax = U85), size = .75)+
  geom_linerange(aes(ymin = L60, ymax = U60), size = 1)+
  geom_hline(yintercept = 0)+
  coord_flip()+
  theme_classic()+
  labs(title = "Specialization")


#Functional redundancy

#Matrix to data frame; split location and dpeth zone; numerise
FRV_index <- FRV_m %>% as.data.frame %>% drop_na 

FRV_index <-FRV_index %>%  mutate(Location = as.factor(stri_match_first_regex(FRV_index$Site, '(\\w+)_(\\w+)')[,2]),
                                  Zone = as.factor(stri_match_first_regex(FRV_index$Site, '(\\w+)_(\\w+)')[,3]),
                                  FRed = as.numeric(as.character(F_Redundancy)),
                                  iter = as.numeric(as.character(iter)),
                                  Zone = fct_relevel(Zone, "Deeper","Upper", "Shallow"),
                                  Ocean = as.factor(if_else(Location == c("Bermuda")|
                                                              Location == c("Curacao")| 
                                                              Location == c("Fernando_de_Noronha")| 
                                                              Location == c("St_Pauls_Rocks"), 
                                                            "Atlantic",
                                                            "Pacific"))) %>% #order mesozones by depth) %>% glimpse
  glimpse



for(i in 1:max(FRV_index$iter)){
  mFRed <- lmer(FRed~Zone*Ocean+(1|Location), data = FRV_index %>% filter(iter == i))
  mFRed.rg <- ref_grid(mFRed) #create emmeans ref grid
  mFRed.emm.s <- emmeans(mFRed.rg, "Zone", "Ocean") #get emmeans
  effFRed <- eff_size(mFRed.emm.s, sigma = sigma(mFRed), edf = 42) %>% as.data.frame()#store pairwise effects from model
  caca<-performance::r2(mFRed)
  
  if(i == 1){
    FRed_estimates_m <- matrix(ncol = 7, nrow = 6000)
    FRed_mod_r2_m <- matrix(ncol = 1, nrow = 1000)
    FRed_contrasts_m <- matrix(ncol = 7, nrow = 6000)
    
    colnames(FRed_estimates_m) <- colnames(mFRed.emm.s %>% as.data.frame)
    colnames(FRed_contrasts_m) <- colnames(effFRed)
  }
  
  FRed_estimates_m[c((1:6)+(i-1)*6),c(1:7)] <- mFRed.emm.s %>% as.data.frame %>% as.matrix
  FRed_contrasts_m[c((1:6)+(i-1)*6),c(1:7)]  <- effFRed %>% as.matrix
  FRed_mod_r2_m [i,] <- as.numeric(caca$R2_marginal)
}



FRed_estimates <- FRed_estimates_m %>% as.data.frame %>%
  mutate(emmean = as.numeric(as.character(emmean)),
         Zone = fct_relevel(Zone, c("Deeper", "Upper", "Shallow"))) %>% drop_na %>% glimpse

FRed_estimates_quantiles <- FRed_estimates %>% 
  group_by(Zone, Ocean) %>% 
  summarise(value = quantile(emmean, probs = c(0.025, 0.075, 0.2, .8, .925, .975)),
            quants = c("L95", "L85", "L60", "U60", "U85", "U95")) %>% 
  spread(quants, value)

est_FRed <- ggplot(FRed_estimates_quantiles, aes(Zone, colour = Ocean))+
  geom_linerange(aes(ymin = L95, ymax = U95), size = .5)+
  geom_linerange(aes(ymin = L85, ymax = U85), size = .75)+
  geom_linerange(aes(ymin = L60, ymax = U60), size = 1)+
  coord_flip()+
  theme_classic()+
  labs(title = "Redundancy")


#Calculate and plot quantiles for contrasts 
FRed_contrasts <- FRed_contrasts_m %>% as.data.frame %>% 
  mutate(effect.size = as.numeric(as.character(effect.size)),
         SE = as.numeric(as.character(SE)),
         df = as.numeric(as.character(df)),
         lower.CL = as.numeric(as.character(lower.CL)),
         upper.CL = as.numeric(as.character(upper.CL)))

FRed_contrasts_quantiles <- FRed_contrasts %>% 
  group_by(contrast, Ocean) %>% 
  summarise(value = quantile(effect.size, probs = c(0.025, 0.075, 0.2, .8, .925, .975) ),
            quants = c("L95", "L85", "L60", "U60", "U85", "U95")) %>% 
  spread(quants, value)

contr_FRed <- ggplot(FRed_contrasts_quantiles, aes(contrast, colour = Ocean))+
  geom_linerange(aes(ymin = L95, ymax = U95), size = .5)+
  geom_linerange(aes(ymin = L85, ymax = U85), size = .75)+
  geom_linerange(aes(ymin = L60, ymax = U60), size = 1)+
  geom_hline(yintercept = 0)+
  coord_flip()+
  theme_classic()+
  labs(title = "Redundancy")


#Example of null model construction####
#Null model for local-regional relationships within depth zones
#Hypothesis is that the breakdown in relationship between local and regional species pools in the mesophtic depths 
#is due to the limited number of species in regional pools with affinity to mesophotic depths. 

# The null model randomly selects 1000 sub-samples of n species from each regional pool, where n is the mean number of species observed 
# across all depths within each location. The null test compares the observed mean number of species in mesophotic depths to the mean
# number of species from the null model, standardised by the standard deviation of observations in the null model. 
# If variation in the number of mesophotic species in the local assemblage is not related to the size and proportion of species available 
# to occupy mesophotic depths in the regional species pool we expect non-zero SES values in all, or some of the locations - with either 
# an ubiquitous or non-random pattern along the gradient of regional species pool sizes. 


#Example of null model using data from Philippines

regPool_genMatch_Philippines <- regPool_genMatch %>% filter(Site.name == "Philippines") %>% drop_na

regPool_Philippines <- regPool %>% filter(Site.name == "Philippines") %>% drop_na

localSpRic_Philippines <- localSpRic %>% filter(Location == "Philippines") %>% 
  group_by(Slice) %>% 
  summarise(SpRic = mean(SpRic))

localSpRic_Philippines_all <- localSpRic_Philippines %>%  summarise(SpRic = sum(SpRic) %>% round(0))
localSpRic_Philippines_all <-localSpRic_Philippines_all$SpRic

localSpRic_Philippines_5 <- localSpRic_Philippines %>% filter(Slice == "5m")
localSpRic_Philippines_30 <- localSpRic_Philippines %>% filter(Slice == "30m")
localSpRic_Philippines_60 <- localSpRic_Philippines %>% filter(Slice == "60m")
localSpRic_Philippines_90 <- localSpRic_Philippines %>% filter(Slice == "90m")
localSpRic_Philippines_120 <- localSpRic_Philippines %>% filter(Slice == "120m")


#Randomly sample n spaecies from Philippines pool = to n species for depth slices in Philippines pool----
i = 1

while (i < 1001){
  rand_Philippines_genMatch <- regPool_genMatch_Philippines %>% 
    slice(c(sample.int(nrow(regPool_genMatch_Philippines), localSpRic_Philippines_all))) %>% 
    dplyr::select(Genus, Species, Pool5:Pool120) %>% gather(Zone, Occurs, Pool5:Pool120) %>% group_by(Genus, Species) %>%
    slice_max(as.numeric(Occurs)) %>% sample_n(1) %>% spread(Zone, Occurs, fill = 0)
  
  rand_Philippines <- regPool_Philippines %>% 
    slice(c(sample.int(nrow(regPool_genMatch_Philippines), localSpRic_Philippines_all))) %>% 
    dplyr::select(Genus, Species, Pool5:Pool120) %>% gather(Zone, Occurs, Pool5:Pool120) %>% group_by(Genus, Species) %>%
    slice_max(as.numeric(Occurs)) %>% sample_n(1) %>% spread(Zone, Occurs, fill = 0)
  
  if(i == 1){
    Philippines_5m_m <- matrix(nrow = 1000, ncol = 2)
  }
  print(paste("start", i))
  Philippines_5m_m[i,1] = sum(rand_Philippines$Pool5) 
  Philippines_5m_m[i,2] = sum(rand_Philippines_genMatch$Pool5) 
  
  if(i == 1){
    Philippines_30m_m <- matrix(nrow = 1000, ncol = 2)
  }
  Philippines_30m_m[i,1] = sum(rand_Philippines$Pool30) 
  Philippines_30m_m[i,2] = sum(rand_Philippines_genMatch$Pool30) 
  
  if(i == 1){
    Philippines_60m_m <- matrix(nrow = 1000, ncol = 2)
  }
  Philippines_60m_m[i,1] = sum(rand_Philippines$Pool60) 
  Philippines_60m_m[i,2] = sum(rand_Philippines_genMatch$Pool60) 
  
  if(i == 1){
    Philippines_90m_m <- matrix(nrow = 1000, ncol = 2)
  }
  Philippines_90m_m[i,1] = sum(rand_Philippines$Pool90) 
  Philippines_90m_m[i,2] = sum(rand_Philippines_genMatch$Pool90)  
  
  if(i == 1){
    Philippines_120m_m <- matrix(nrow = 1000, ncol = 2)
  }
  Philippines_120m_m[i,1] = sum(rand_Philippines$Pool120) 
  Philippines_120m_m[i,2] = sum(rand_Philippines_genMatch$Pool120) 
  print(paste("end", i))
  
  i = i+1
}

SES_Philippines_sp_pool_5m <- (localSpRic_Philippines_5$SpRic - mean(Philippines_5m_m[,1]))/sd(Philippines_5m_m[,1])
SES_Philippines_sp_pool_30m <- (localSpRic_Philippines_30$SpRic - mean(Philippines_30m_m[,1]))/sd(Philippines_30m_m[,1])
SES_Philippines_sp_pool_60m <- (localSpRic_Philippines_60$SpRic - mean(Philippines_60m_m[,1]))/sd(Philippines_60m_m[,1])
SES_Philippines_sp_pool_90m <- (localSpRic_Philippines_90$SpRic - mean(Philippines_90m_m[,1]))/sd(Philippines_90m_m[,1])
SES_Philippines_sp_pool_120m <- (localSpRic_Philippines_120$SpRic - mean(Philippines_120m_m[,1]))/sd(Philippines_120m_m[,1])

SES_Philippines_sp_pool <- data.frame(Location = "Philippines",
                                      Slice = c(5, 30, 60, 90, 120),
                                      SES = c(SES_Philippines_sp_pool_5m,
                                              SES_Philippines_sp_pool_30m,
                                              SES_Philippines_sp_pool_60m,
                                              SES_Philippines_sp_pool_90m,
                                              SES_Philippines_sp_pool_120m))

SES_Philippines_sp_pool

SES_Philippines_genMatch_sp_pool_5m <- (localSpRic_Philippines_5$SpRic - mean(Philippines_5m_m[,2]))/sd(Philippines_5m_m[,2])
SES_Philippines_genMatch_sp_pool_30m <- (localSpRic_Philippines_30$SpRic - mean(Philippines_30m_m[,2]))/sd(Philippines_30m_m[,2])
SES_Philippines_genMatch_sp_pool_60m <- (localSpRic_Philippines_60$SpRic - mean(Philippines_60m_m[,2]))/sd(Philippines_60m_m[,2])
SES_Philippines_genMatch_sp_pool_90m <- (localSpRic_Philippines_90$SpRic - mean(Philippines_90m_m[,2]))/sd(Philippines_90m_m[,2])
SES_Philippines_genMatch_sp_pool_120m <- (localSpRic_Philippines_120$SpRic - mean(Philippines_120m_m[,2]))/sd(Philippines_120m_m[,2])

SES_Philippines_genMatch_sp_pool <- data.frame(Location = "Philippines",
                                               Slice = c(5, 30, 60, 90, 120),
                                               SES = c(SES_Philippines_genMatch_sp_pool_5m,
                                                       SES_Philippines_genMatch_sp_pool_30m,
                                                       SES_Philippines_genMatch_sp_pool_60m,
                                                       SES_Philippines_genMatch_sp_pool_90m,
                                                       SES_Philippines_genMatch_sp_pool_120m))

SES_Philippines_genMatch_sp_pool

rand_local_Philippines <- data.frame(Location = "Philippines",
                                     Slice_5 = Philippines_5m_m[,1],
                                     Slice_30 = Philippines_30m_m[,1],
                                     Slice_60 = Philippines_60m_m[,1],
                                     Slice_90 = Philippines_90m_m[,1],
                                     Slice_120 = Philippines_120m_m[,1])

rand_local_Philippines %>% glimpse

rand_local_genMatch_Philippines <- data.frame(Location = "Philippines",
                                              Slice_5 = Philippines_5m_m[,2],
                                              Slice_30 = Philippines_30m_m[,2],
                                              Slice_60 = Philippines_60m_m[,2],
                                              Slice_90 = Philippines_90m_m[,2],
                                              Slice_120 = Philippines_120m_m[,2])

rand_local_genMatch_Philippines %>% glimpse

