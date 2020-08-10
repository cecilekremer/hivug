
##################################
### 1. MERGE INPUT-OUTPUT DATA ###
##################################

out=read.csv("output_v28.csv", sep=",", header=F)
colnames(out)<-c("sim","N_total","prev_all","prev_f1","prev_f2","prev_f3","prev_m1","prev_m2","prev_m3",
                 "GOF_all","GOF_f1","GOF_f2","GOF_f3","GOF_m1","GOF_m2","GOF_m3")
input=read.csv("input_v28.csv", sep=",", header=T)
#merge input and output datasets
result = merge(out,input,by="sim"); dim(result); View(result)

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
write.csv(result,"data_ARF_v28.csv",sep=',', col.names=T, row.names=F)



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
best17 = top_n(shyza17,-1000,GOF)
dim(best17)

attach(best17); summary(GOF)


# combined plot of most influential parameters: 
jpeg("FigB1.jpeg", width = 18, height = 14, units = 'cm', res = 300)

par(mfrow=c(2,2))
fit3=density(input$C_f); plot(fit3,ylim=c(0,4),lty=2, lwd=2, main=expression("C"[f]), xlab="Range"); d3 <- density(C_f); lines(d3,lwd=2)
fit3=density(input$C_m); plot(fit3,ylim=c(0,4),lty=2, lwd=2, main=expression("C"[m]), ylab=" ", xlab="Range"); d3 <- density(C_m); lines(d3,lwd=2)
fit3=density(input$beta_f); plot(fit3,ylim=c(0,600),lty=2, lwd=2, main=expression(beta[f]), xlab="Range"); d3 <- density(beta_f); lines(d3,lwd=2)
fit3=density(input$beta_m); plot(fit3,ylim=c(0,600),lty=2, lwd=2, main=expression(beta[m]),ylab="", xlab="Range"); d3 <- density(beta_m); lines(d3,lwd=2)

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
best$sim

best_parms = round(best[c(18:57)],4)

##########################################
### 4. MAXIMAL INFORMATION COEFFICIENT ###
##########################################
library(minerva)

## Association structure in the top 1% solution parameter space: delete parms that are fixed
best = best17[,c("lambda","C_f","C_m","high_v2","beta_f","beta_m", "cond_m1",
                 "cond_f2","cond_f3","cond_eff","acts_m3_casual","acts_m2_reg",
                 "acts_f1_main",
                 "acts_m2_main","df1","dm2","dm3","N","A_main","A_casual",
                 "nu"),]
mic_best=mine(x=best, measure="mic_approx")
mic_best2=mic_best$MIC; View(mic_best2)
mic=as.numeric(mic_best2); summary(mic) #mic=1 on the diagonal
# highest MIC value? summary statistics minus the diagonal
max(mic[mic<1])
