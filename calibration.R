
### 1. MERGE INPUT-OUTPUT DATA ###
##################################

out=read.csv("output_v24.csv", sep=",", header=F)
colnames(out)<-c("sim","Ntot",
                 "IRf1","IRf2","IRf3","IRm1","IRm2","IRm3",
                 "PRf1","PRf2","PRf3","PRm1","PRm2","PRm3",
                 "IRtot","PRtot","propART","propCIRC","AMtot","ARTinc","propAware",
                 "GOF_all","GOF_f1","GOF_f2","GOF_f3","GOF_m1","GOF_m2","GOF_m3",
                 "GOF_IR","GOF_ART","GOF_HCT","GOF_circ")
input=read.csv("input_v24.csv", sep=",", header=T)
#merge input and output datasets
result = merge(out,input,by="sim"); dim(result); head(result)

## Create overall weighted GOF estimates ##

attach(result)
result$GOF = GOF_all +
  GOF_f1 +
  GOF_f2 +
  GOF_f3 +
  GOF_m1 +
  GOF_m2 +
  GOF_m3
summary(result$GOF)

#write csv file to use in ARF analysis
write.csv(result,"data_ARFv24.csv",sep=',', col.names=T, row.names=F)

##########################################################
### 2. GOF ESTIMATES: SELECT TOP 1% OF ALL SIMULATIONS ###
##########################################################

# Data from all 100 000 simulations
shyza17 = result
attach(shyza17); dim(shyza17)

## Select top 1% that minimize overall GOF
library(dplyr)
attach(shyza17)
shyza17 = shyza17[order(GOF),]
attach(shyza17)
summary(shyza17$GOF)
best17 = top_n(shyza17,-100,GOF)
dim(best17)

attach(best17); summary(GOF)
# input parameters best fitting model
write.csv(best17[,c(1,33:77)], "input_base.csv", quote=F, sep=",", col.names=T, row.names = F)


## Comparing uniform (input) and output densities
# Plot for manuscript
jpeg("FigS2.jpeg", width = 18, height = 20, units = 'cm', res = 300)

par(mfrow=c(3,2), cex.main=2, cex.lab = 1.2)

fit1=density(input$const3)
plot(fit1,ylim=c(0,11),lty=2, lwd=2, main=expression("x"[casual]), xlab="Range") #plots uniform distribution (input)
d1 <- density(const3) 
lines(d1, lwd=2) #adds the new obtained distribution

fit1=density(input$actsf33_main)
plot(fit1,ylim=c(0,0.55),lty=2, lwd=2, main=expression(eta[main33]), ylab="", xlab="Range") #plots uniform distribution (input)
d1 <- density(actsf33_main) 
lines(d1, lwd=2) #adds the new obtained distribution

fit1=density(input$beta_m)
plot(fit1,ylim=c(0,100000),lty=2, lwd=2,  main=expression(beta[m]), xlab="Range") #plots uniform distribution (input)
d1 <- density(beta_m) 
lines(d1, lwd=2) #adds the new obtained distribution

fit1=density(input$pHCT)
plot(fit1,ylim=c(0,9),lty=2,lwd=2,  main=expression(P[HCT]), ylab="", xlab="Range") #plots uniform distribution (input)
d1 <- density(pHCT) 
lines(d1, lwd=2) #adds the new obtained distribution

fit1=density(input$f_casual)
plot(fit1,ylim=c(0,11),lty=2, lwd=2, main=expression(f[casual]), xlab="Range") #plots uniform distribution (input)
d1 <- density(f_casual) 
lines(d1, lwd=2) #adds the new obtained distribution

dev.off()


#################################
### 3. ACTIVITY REGION FINDER ###
#################################

# can not be run from RStudio -> use R 2.3.1

# Check whether best solution values of the adapted parameters still lie within the interval
attach(shyza17)
summary(GOF)
library(dplyr)
best = top_n(shyza17,-1,GOF); dim(best); summary(best$GOF)
best$tau1

##########################################
### 4. MAXIMAL INFORMATION COEFFICIENT ###
##########################################
library(minerva)

## Association structure in the top 1% solution parameter space: delete parms that are fixed
best = best17[,c("actsf21_main",
"actsf23_main",
"actsf33_main",
"beta_m",
"pHCT",
"tau1",
"const3",
"f_casual"),]
mic_best=mine(x=best, measure="mic_approx")
mic_best2=mic_best$MIC; View(mic_best2)
mic=as.numeric(mic_best2); summary(mic) #mic=1 on the diagonal
# highest MIC value? summary statistics minus the diagonal
max(mic[mic<1])
