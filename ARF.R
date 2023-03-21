
#Read in input-output dataset
data = read.csv("data_ARFv24.csv", sep=",", header=T)
attach(data); dim(data); summary(GOF)

#################################
### 3. ACTIVITY REGION FINDER ###
#################################

library(ARF)

## Create binary response variable: 1 if parameter set belongs to top 1%, 0 otherwise

data$top1[rank(data$GOF)<=100]=1
data$top1[rank(data$GOF)>100]=0
attach(data); top1=as.factor(top1)
table(top1)

#check GOF of top 1%
test=data[which(data$top1==1),]
summary(test$GOF)

## Generate an ARF object that contains the info for constructing a tree

attach(data)
# remove fixed parameters !
best.arf2=f.arf(top1 ~ f_casual + actsf21_main
                + actsf23_main + actsf33_main
                + beta_m + tau1 + pHCT
                + const3, data=data, const=25, minn=20)          #const=lambda


## Generate a report:
f.report(best.arf2, file="ARF_v24.pdf")









