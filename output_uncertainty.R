
##################################
### 1. MERGE INPUT-OUTPUT DATA ###
##################################

out=read.csv("unct_final.csv", sep=",", header=F)
dim(out)
colnames(out)<-c("sim","RNI_tot","RNInew_tot",
                 "RP_f1","RP_f2","RP_f3","RP_m1","RP_m2","RP_m3","RP_tot",
                 "RI_f1","RI_f2","RI_f3","RI_m1","RI_m2","RI_m3","RI_tot",
                 "RA_f1","RA_f2","RA_f3","RA_m1","RA_m2","RA_m3","RA_tot",
                 "RNInew_f1","RNInew_f2","RNInew_f3","RNInew_m1","RNInew_m2","RNInew_m3")

input=read.csv("data_unct_hpc.csv", sep=",", header=T)
dim(input)
#merge input and output datasets
result = merge(out,input,by="sim"); dim(result); #View(result)
result_final = result

################################################
### BEST FIT BASELINE                        ###
################################################

summary(result_final$RI_tot); summary(result_final$RP_tot); summary(result_final$RNI_tot); summary(result_final$RA_tot)
quantile(result_final$RR_all,c(0.025,0.5,0.975))
quantile(result_final$cum_all,c(0.025,0.5,0.975))

summary(((prev_all_base/100*N_total_base) - (prev_all_int/100*N_total_int))/(prev_all_base/100*N_total_base))



# 95% CI relative reduction prevalence, best fit ??
quantile(RR_all,c(0.025,0.5,0.975))
quantile(cum_all,c(0.025,0.5,0.975))

##################
### SOME PLOTS ###
##################
##### FIGURE 2
library(Hmisc)
library(corrplot)
cor.data = result_final[,c(17,3,10,24,76:86)]
res = cor(cor.data)
round(res,2)
colnames(res) = c("Incidence reduction","New infections reduction","Prevalence reduction","AIDS-related mortality reduction","ART acute stage","ART chronic stage","ART pre-AIDS stage","ART AIDS stage","VMMC acceptance",
                  "CU low-risk","CU medium-risk","CU high-risk","Proportion diagnosed","PrEP acceptance","HTS acceptance")
rownames(res) = c("Incidence reduction","New infections reduction","Prevalence reduction","AIDS-related mortality reduction","ART acute stage","ART chronic stage","ART pre-AIDS stage","ART AIDS stage","VMMC acceptance",
                  "CU low-risk","CU medium-risk","CU high-risk","Proportion diagnosed","PrEP acceptance","HTS acceptance")

jpeg("G:/My Drive/UHasselt/PhD/HIV/HIV Uganda (Rebecca Nsubuga)/UPDAT DEC 2020/MANUSCRIPT/PLOS ONE/BMCID/revision/Fig2.jpeg", 
     width = 25, height = 20, units = 'cm', res = 1200)
corrplot(res[1:4,], method = 'shade', type="upper", diag=F, tl.col="black", tl.cex=1.5, cl.cex=1, cl.pos = 'b', cl.ratio = 0.4)
dev.off()

jpeg("Fig3new.jpeg", width = 40, height = 20, units = 'cm', res = 300)
par(mfrow=c(1,3))
par(oma = c(1,1,4,1), mar=c(5,5,3,5))
plot(result_final$scaleCOND1, result_final$RI_tot, ylim=c(0,0.9), xlab="Condom use low-risk group", ylab="Relative reduction",
     cex.axis = 1.5, cex.lab = 2, pch = 16)
points(result_final$scaleCOND1, result_final$RP_tot, col="blue", pch = 16)
points(result_final$scaleCOND1, result_final$RA_tot, col="green", pch = 16)
plot(result_final$scaleART2, result_final$RI_tot, ylim=c(0,0.9), xlab="ART uptake chronic stage", ylab="",
     cex.axis = 1.5, cex.lab = 2, pch = 16)
points(result_final$scaleART2, result_final$RP_tot, col="blue", pch = 16)
points(result_final$scaleART2, result_final$RA_tot, col="green", pch = 16)
plot(result_final$pHCT2, result_final$RI_tot, ylim=c(0,0.9), xlab="Proportion of diagnosed individuals", ylab="",
     cex.axis = 1.5, cex.lab = 2, pch = 16)
points(result_final$pHCT2, result_final$RP_tot, col="blue", pch = 16)
points(result_final$pHCT2, result_final$RA_tot, col="green", pch = 16)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('top', legend=c("Incidence","Prevalence","AIDS-related mortality"), pch=c(16,16,16), col=c("black","blue","green"),
       box.lwd=0, box.col="white", cex=3, pt.cex=5, xpd=T, horiz=T,seg.len=1,bty="n")
dev.off()

##### SYNERGY? FIGURE S5
library(visreg)

# Incidence
jpeg("G:/My Drive/UHasselt/PhD/HIV/HIV Uganda (Rebecca Nsubuga)/UPDAT DEC 2020/MANUSCRIPT/PLOS ONE/BMCID/revision/FigS5a1.jpeg", 
     width = 20, height = 20, units = 'cm', res = 300)
par(oma = c(1,1,1,1))
mod <- lm(RI_tot ~ scaleART2*scaleCOND1, data = cor.data)
visreg2d(mod, xvar = 'scaleART2', yvar = 'scaleCOND1', nn = 300, main = 'Incidence reduction', xlab = 'ART chronic stage', ylab = 'Condom use low-risk group',
         cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()
jpeg("G:/My Drive/UHasselt/PhD/HIV/HIV Uganda (Rebecca Nsubuga)/UPDAT DEC 2020/MANUSCRIPT/PLOS ONE/BMCID/revision/FigS5a2.jpeg", 
     width = 20, height = 20, units = 'cm', res = 300)
par(oma = c(1,1,1,1))
mod2 <- lm(RI_tot ~ pHCT2*scaleART2, data = cor.data)
visreg2d(mod2, xvar = 'scaleART2', yvar = 'pHCT2', nn = 300, main = 'Incidence reduction', xlab = 'ART chronic stage', ylab = 'Proportion diagnosed',
         cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()
jpeg("G:/My Drive/UHasselt/PhD/HIV/HIV Uganda (Rebecca Nsubuga)/UPDAT DEC 2020/MANUSCRIPT/PLOS ONE/BMCID/revision/FigS5a3.jpeg", 
     width = 20, height = 20, units = 'cm', res = 300)
par(oma = c(1,1,1,1))
mod3 <- lm(RI_tot ~ pHCT2*scaleCOND1, data = cor.data)
visreg2d(mod3, xvar = 'pHCT2', yvar = 'scaleCOND1', nn = 300, main = 'Incidence reduction', xlab = 'Proportion diagnosed', ylab = 'Condom use low-risk group',
         cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

# Prevalence
jpeg("G:/My Drive/UHasselt/PhD/HIV/HIV Uganda (Rebecca Nsubuga)/UPDAT DEC 2020/MANUSCRIPT/PLOS ONE/BMCID/revision/FigS5b1.jpeg", 
     width = 20, height = 20, units = 'cm', res = 300)
par(oma = c(1,1,1,1))
mod <- lm(RP_tot ~ scaleART2*scaleCOND1, data = cor.data)
visreg2d(mod, xvar = 'scaleART2', yvar = 'scaleCOND1', nn = 300, main = 'Prevalence reduction', xlab = 'ART chronic stage', ylab = 'Condom use low-risk group',
         cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()
jpeg("G:/My Drive/UHasselt/PhD/HIV/HIV Uganda (Rebecca Nsubuga)/UPDAT DEC 2020/MANUSCRIPT/PLOS ONE/BMCID/revision/FigS5b2.jpeg", 
     width = 20, height = 20, units = 'cm', res = 300)
par(oma = c(1,1,1,1))
mod2 <- lm(RP_tot ~ pHCT2*scaleART2, data = cor.data)
visreg2d(mod2, xvar = 'scaleART2', yvar = 'pHCT2', nn = 300, main = 'Prevalence reduction', xlab = 'ART chronic stage', ylab = 'Proportion diagnosed',
         cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()
jpeg("G:/My Drive/UHasselt/PhD/HIV/HIV Uganda (Rebecca Nsubuga)/UPDAT DEC 2020/MANUSCRIPT/PLOS ONE/BMCID/revision/FigS5b3.jpeg", 
     width = 20, height = 20, units = 'cm', res = 300)
par(oma = c(1,1,1,1))
mod3 <- lm(RP_tot ~ pHCT2*scaleCOND1, data = cor.data)
visreg2d(mod3, xvar = 'pHCT2', yvar = 'scaleCOND1', nn = 300, main = 'Prevalence reduction', xlab = 'Proportion diagnosed', ylab = 'Condom use low-risk group',
         cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

# AR-mortality
jpeg("G:/My Drive/UHasselt/PhD/HIV/HIV Uganda (Rebecca Nsubuga)/UPDAT DEC 2020/MANUSCRIPT/PLOS ONE/BMCID/revision/FigS5c1.jpeg", 
     width = 20, height = 20, units = 'cm', res = 300)
par(oma = c(1,1,1,1))
mod <- lm(RA_tot ~ scaleART2*scaleCOND1, data = cor.data)
visreg2d(mod, xvar = 'scaleART2', yvar = 'scaleCOND1', nn = 300, main = 'AIDS-related mortality reduction', xlab = 'ART chronic stage', ylab = 'Condom use low-risk group',
         cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()
jpeg("G:/My Drive/UHasselt/PhD/HIV/HIV Uganda (Rebecca Nsubuga)/UPDAT DEC 2020/MANUSCRIPT/PLOS ONE/BMCID/revision/FigS5c2.jpeg", 
     width = 20, height = 20, units = 'cm', res = 300)
par(oma = c(1,1,1,1))
mod2 <- lm(RA_tot ~ pHCT2*scaleART2, data = cor.data)
visreg2d(mod2, xvar = 'scaleART2', yvar = 'pHCT2', nn = 300, main = 'AIDS-related mortality reduction', xlab = 'ART chronic stage', ylab = 'Proportion diagnosed',
         cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()
jpeg("G:/My Drive/UHasselt/PhD/HIV/HIV Uganda (Rebecca Nsubuga)/UPDAT DEC 2020/MANUSCRIPT/PLOS ONE/BMCID/revision/FigS5c3.jpeg", 
     width = 20, height = 20, units = 'cm', res = 300)
par(oma = c(1,1,1,1))
mod3 <- lm(RA_tot ~ pHCT2*scaleCOND1, data = cor.data)
visreg2d(mod3, xvar = 'pHCT2', yvar = 'scaleCOND1', nn = 300, main = 'AIDS-related mortality reduction', xlab = 'Proportion diagnosed', ylab = 'Condom use low-risk group',
         cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()
