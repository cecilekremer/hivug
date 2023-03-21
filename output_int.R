
##################################
### 1. MERGE INPUT-OUTPUT DATA ###
##################################

out=read.csv("output_int_best.csv", sep=",", header=F)
colnames(out)<-c("sim","RNI_tot","RNInew_tot",
                 "RP_f1","RP_f2","RP_f3","RP_m1","RP_m2","RP_m3","RP_tot",
                 "RI_f1","RI_f2","RI_f3","RI_m1","RI_m2","RI_m3","RI_tot",
                 "RA_f1","RA_f2","RA_f3","RA_m1","RA_m2","RA_m3","RA_tot",
                 "RNInew_f1","RNInew_f2","RNInew_f3","RNInew_m1","RNInew_m2","RNInew_m3")

input=read.csv("input_int_best.csv", sep=",", header=T)
#merge input and output datasets
result = merge(out,input,by="sim"); dim(result); #View(result)
attach(result)
dim(result)
head(result)

################################################
### 2. INVESTIGATE RR AND AVERTED INFECTIONS ###
################################################

summary(result$RP_tot)
summary(result$RI_tot)
summary(result$RA_tot)
summary(result$RNI_tot)
summary(result$RNInew_tot)


## Simulation with highest overall incidence reduction (= sim with highest overall prevalence reduction ?)
max(result$RI_tot)
max(result$RP_tot)
max(result$RA_tot)
max(result$RNInew_tot)
best = result[which(result$RI_tot == max(result$RI_tot)),]
best$sim

# Intervention coverage
coverage_best = best[c(76:86)]
coverage = round(coverage_best,4)
coverage

best$RI_tot; best$RI_f1; best$RI_f2; best$RI_f3; best$RI_m1; best$RI_m2; best$RI_m3
best$RP_tot; best$RP_f1; best$RP_f2; best$RP_f3; best$RP_m1; best$RP_m2; best$RP_m3
best$RA_tot; best$RA_f1; best$RA_f2; best$RA_f3; best$RA_m1; best$RA_m2; best$RA_m3
# best$RNI_tot
best$RNInew_tot; best$RNInew_f1; best$RNInew_f2; best$RNInew_f3; best$RNInew_m1; best$RNInew_m2; best$RNInew_m3

##########################################################

### Run 'model_int.R' with best fitting input parameters 



##########################################################

### FIGURES MANUSCRIPT

## FIGURE 1
jpeg("Fig1.eps", width = 50, height = 20, units = 'cm', res = 1200)
par(mfrow=c(1,2))

# prevalence
par(oma=c(0,1,0,0))
matplot(trial1_FM_BASE[,"time"], trial1_FM_BASE[,c("PRtot")]*100, type = "l", main="(a)", xlab = "Year", xaxt="n", ylab = "Prevalence (%)",
        lwd = 2, col=1, lty=1, cex.axis=1.5, cex.lab = 1.5)
points(seq(1,721,by=48), trial1_FM_BASE[seq(1,721,by=48),c("PRtot")]*100, pch=21, col=1, cex=1.5,lwd=2)
# observed
abline(v=540, lty=3, lwd=2, col="darkgrey")
# scale-up
lines(trial1_FM[540:721,"time"], trial1_FM[540:721,c("PRtot")]*100, lwd = 2, col=2, lty=2)
points(seq(540,721,by=48), trial1_FM[seq(540,721,by=48),c("PRtot")]*100, pch=24, col=2, cex=1.5, lwd=2)
# sens
lines(trial1_FM_SENS4[540:721,"time"], trial1_FM_SENS4[540:721,c("PRtot")]*100, lwd = 2, col=3, lty=3)
points(seq(540,721,by=48), trial1_FM_SENS4[seq(540,721,by=48),c("PRtot")]*100, pch=22, col=3, cex=1.5, lwd=2)

lines(trial1_FM_SENS3[540:721,"time"], trial1_FM_SENS3[540:721,c("PRtot")]*100, lwd = 2, col=4, lty=3)
points(seq(540,721,by=48), trial1_FM_SENS3[seq(540,721,by=48),c("PRtot")]*100, pch=23, col=4, cex=1.5, lwd=2)

lines(trial1_FM_SENS9[540:721,"time"], trial1_FM_SENS9[540:721,c("PRtot")]*100, lwd = 2, col="darkgrey", lty=3)
points(seq(540,721,by=48), trial1_FM_SENS9[seq(540,721,by=48),c("PRtot")]*100, pch=25, col="darkgrey", cex=1.5,lwd=2)

legend(250,7, cex=1.5, legend=c("Baseline","Intervention scale-up","5% yearly ART dropout","No increase in condom use","30% ART initiation chronic stage"),
       lty=c(1,2,3,3,3), lwd=c(2,2,2,2,2), pch=c(21,24,22,23,25), col=c(1,2,3,4,"darkgrey"), box.lwd=0, box.col="white")

axis(1, at = c(1,180,360,540,720), labels = c("1969","1984","1999","2014","2029"), cex.axis = 1.5)

# incidence
matplot(trial1_FM_BASE[,"time"], trial1_FM_BASE[,c("IRtot")], type = "l", main="(b)", xlab = "Year", xaxt="n", ylab = "Incidence per 100 PY",
        lwd = 2, col=1, lty=1, ylim=c(0,3.2), cex.axis = 1.5, cex.lab = 1.5)
points(seq(1,721,by=48), trial1_FM_BASE[seq(1,721,by=48),c("IRtot")], pch=21, col=1, cex=1.5, lwd=2)
abline(v=540, lty=3, col="darkgrey", lwd=2)
# scale-up
lines(trial1_FM[540:721,"time"], trial1_FM[540:721,c("IRtot")], lwd = 2, col=2, lty=2)
points(seq(540,721,by=48), trial1_FM[seq(540,721,by=48),c("IRtot")], pch=24, col=2, cex=1.5, lwd=2)
# sens
lines(trial1_FM_SENS4[540:721,"time"], trial1_FM_SENS4[540:721,c("IRtot")], lwd = 2, col=3, lty=3)
points(seq(540,721,by=48), trial1_FM_SENS4[seq(540,721,by=48),c("IRtot")], pch=22, col=3, cex=1.5, lwd=2)

lines(trial1_FM_SENS3[540:721,"time"], trial1_FM_SENS3[540:721,c("IRtot")], lwd = 2, col=4, lty=3)
points(seq(540,721,by=48), trial1_FM_SENS3[seq(540,721,by=48),c("IRtot")], pch=23, col=4, cex=1.5, lwd=2)

lines(trial1_FM_SENS9[540:721,"time"], trial1_FM_SENS9[540:721,c("IRtot")], lwd = 2, col="darkgrey", lty=3)
points(seq(540,721,by=48), trial1_FM_SENS9[seq(540,721,by=48),c("IRtot")], pch=25, col="darkgrey", cex=1.5, lwd=2)

axis(1, at = c(1,180,360,540,720), labels = c("1969","1984","1999","2014","2029"), cex.axis = 1.5)

dev.off()

##################################################

### Intervention uptake

jpeg('FigS4.jpg', width = 30, height=30, res = 300, units = 'cm')
par(mfrow = c(2,2))
par(mar = c(5.1,5.1,4.1,2.1))
matplot(trial1_FM[,"time"], trial1_FM[,"propART"], type="l",xlim=c(415,722),
        xlab = 'Year', ylab = 'On ART of all HIV infected', xaxt='n', cex.lab = 2, cex.axis = 1.5, main = '(a)'); abline(v=504,lty=2); abline(v=612, lty=3)
axis(1, at = c(420,504,612,720), labels = c("2004","2011","2020","2029"), cex.axis = 1.2)
matplot(trial1_FM[,"time"], trial1_FM[,"ART.aware"], type="l",xlim=c(415,722), xaxt='n',
        xlab = 'Year', ylab = 'On ART of those aware of HIV+ status', cex.lab = 2, cex.axis = 1.5, main = '(b)'); abline(v=504,lty=2); abline(v=612, lty=3)
axis(1, at = c(420,504,612,720), labels = c("2004","2011","2020","2029"), cex.axis = 1.2)
matplot(trial1_FM[,"time"], trial1_FM[,"propCIRC"], type="l", xaxt='n',
        xlab = 'Year', ylab = 'Susceptible males circumcised', cex.lab = 2, cex.axis = 1.5, main = '(c)'); abline(v=612, lty=3)
axis(1, at = c(1,216,420,612,720), labels = c("1969","1986","2004","2020","2029"), cex.axis = 1.2)
matplot(trial1_FM[,"time"], trial1_FM[,"prop.aware"], type="l", xaxt='n',
        xlab = 'Year', ylab = 'Aware of HIV+ status', cex.lab = 2, cex.axis = 1.5, main = '(d)'); abline(v=216,lty=4); abline(v=612, lty=3)
axis(1, at = c(1,216,420,612,720), labels = c("1969","1986","2004","2020","2029"), cex.axis = 1.2)
dev.off()


