library(ggplot2)

#---------USE MASTER-----------#
master$OMC <- 1-master$Cindex50

corr <- rep(NA, 144)
for (i in 1:144){
  corr[i] <- cor(master$qloss50[master$FLC==i], master$OMC[master$FLC==i])
}
hist(corr, breaks=20)

#By Method
master1 <- master[master$Method=="logrank",]
master2 <- master[master$Method=="L1",]
master3 <- master[master$Method=="CIF",]
master4 <- master[master$Method=="RIST",]
master5 <- master[master$Method=="RotSF",]
master6 <- master[master$Method=="ORSF",]

corr1 <- rep(NA, 144)
for (i in 1:144){
  corr1[i] <- cor(master1$qloss50[master1$FLC==i], master1$OMC[master1$FLC==i])
}

corr2 <- rep(NA, 144)
for (i in 1:144){
  corr2[i] <- cor(master2$qloss50[master2$FLC==i], master1$OMC[master2$FLC==i])
}

corr3 <- rep(NA, 144)
for (i in 1:144){
  corr3[i] <- cor(master3$qloss50[master3$FLC==i], master1$OMC[master3$FLC==i])
}

corr4 <- rep(NA, 144)
for (i in 1:144){
  corr4[i] <- cor(master4$qloss50[master4$FLC==i], master1$OMC[master4$FLC==i])
}

corr5 <- rep(NA, 144)
for (i in 1:144){
  corr5[i] <- cor(master5$qloss50[master5$FLC==i], master1$OMC[master5$FLC==i])
}

corr6 <- rep(NA, 144)
for (i in 1:144){
  corr6[i] <- cor(master6$qloss50[master6$FLC==i], master1$OMC[master6$FLC==i])
}

#--------------------#
par(mfrow = c(3, 2))
hist(corr1, breaks=20, main="logrank", xlab="Correlation of log(log(AL)) and E", ylab="Count")
hist(corr2, breaks=20, main="L1", xlab="Correlation of log(log(AL)) and E", ylab="Count")
hist(corr3, breaks=20, main="CIF", xlab="Correlation of log(log(AL)) and E", ylab="Count")
hist(corr4, breaks=20, main="RIST", xlab="Correlation of log(log(AL)) and E", ylab="Count")
hist(corr5, breaks=20, main="RotSF", xlab="Correlation of log(log(AL)) and E", ylab="Count")
hist(corr6, breaks=20, main="ORSF", xlab="Correlation of log(log(AL)) and E", ylab="Count")
#--------------------------------------#
#-------Histograms with ggplot---------#
#--------------------------------------#

# Create a data frame to store the correlation data and method labels
data <- data.frame(
  Correlation = c(corr1, corr2, corr3, corr4, corr5, corr6),
  Method = factor(
    rep(c("logrank", "L1", "CIF", "RIST", "RotSF", "ORSF"), 
        times = c(length(corr1), length(corr2), length(corr3), length(corr4), length(corr5), length(corr6))
    )
  )
)

p <- ggplot(data, aes(x = Correlation, fill = Method)) +
  geom_histogram(position = "dodge", binwidth = 0.1) +
  facet_wrap(~Method, ncol = 3) +
  labs(
    x = "Sample correlation of log(log(AL)) and E",
    y = "Count"
  ) +
  theme_minimal()

print(p)
#----BS Below---#
fullresults <- read.csv("Full_results_merged.csv")
fullresults$aloss50 <- 2*fullresults$qloss50
fullresults$OMC <- 1-fullresults$Cindex50
cor(fullresults$aloss50, fullresults$OMC)
x
plot(log(log(fullresults$aloss50)), fullresults$OMC, xlab="log(log(AL))", ylab="E=1-C")
cor(log(log(fullresults$aloss50)), fullresults$OMC)

#Broken down by Method
plot(log(log(fullresults$aloss50[fullresults$Method=="logrank"])), fullresults$OMC[fullresults$Method=="logrank"])
cor(log(log(fullresults$aloss50[fullresults$Method=="logrank"])), fullresults$OMC[fullresults$Method=="logrank"])

plot(log(log(fullresults$aloss50[fullresults$Method=="L1"])), fullresults$OMC[fullresults$Method=="L1"])
cor(log(log(fullresults$aloss50[fullresults$Method=="L1"])), fullresults$OMC[fullresults$Method=="L1"])

plot(log(log(fullresults$aloss50[fullresults$Method=="CIF"])), fullresults$OMC[fullresults$Method=="CIF"])
cor(log(log(fullresults$aloss50[fullresults$Method=="CIF"])), fullresults$OMC[fullresults$Method=="CIF"])

plot(log(log(fullresults$aloss50[fullresults$Method=="RIST"])), fullresults$OMC[fullresults$Method=="RIST"])
cor(log(log(fullresults$aloss50[fullresults$Method=="RIST"])), fullresults$OMC[fullresults$Method=="RIST"])

plot(log(log(fullresults$aloss50[fullresults$Method=="RotSF"])), fullresults$OMC[fullresults$Method=="RotSF"])
cor(log(log(fullresults$aloss50[fullresults$Method=="RotSF"])), fullresults$OMC[fullresults$Method=="RotSF"])

plot(log(log(fullresults$aloss50[fullresults$Method=="ORSF"])), fullresults$OMC[fullresults$Method=="ORSF"])
cor(log(log(fullresults$aloss50[fullresults$Method=="ORSF"])), fullresults$OMC[fullresults$Method=="ORSF"])

#Broken down by DGM-Y
plot(log(log(fullresults$aloss50[fullresults$DGM.Y=="Weibull"])), fullresults$OMC[fullresults$DGM.Y=="Weibull"])
cor(log(log(fullresults$aloss50[fullresults$DGM.Y=="Weibull"])), fullresults$OMC[fullresults$DGM.Y=="Weibull"])

plot(log(log(fullresults$aloss50[fullresults$DGM.Y=="AFT"])), fullresults$OMC[fullresults$DGM.Y=="AFT"])
cor(log(log(fullresults$aloss50[fullresults$DGM.Y=="AFT"])), fullresults$OMC[fullresults$DGM.Y=="AFT"])

plot(log(log(fullresults$aloss50[fullresults$DGM.Y=="GG"])), fullresults$OMC[fullresults$DGM.Y=="GG"])
cor(log(log(fullresults$aloss50[fullresults$DGM.Y=="GG"])), fullresults$OMC[fullresults$DGM.Y=="GG"])

#Broken down by DGM-X
plot(log(log(fullresults$aloss50[fullresults$DGM.X=="Linear"])), fullresults$OMC[fullresults$DGM.X=="Linear"])
cor(log(log(fullresults$aloss50[fullresults$DGM.X=="Linear"])), fullresults$OMC[fullresults$DGM.X=="Linear"])

plot(log(log(fullresults$aloss50[fullresults$DGM.X=="QuadInt"])), fullresults$OMC[fullresults$DGM.X=="QuadInt"])
cor(log(log(fullresults$aloss50[fullresults$DGM.X=="QuadInt"])), fullresults$OMC[fullresults$DGM.X=="QuadInt"])

plot(log(log(fullresults$aloss50[fullresults$DGM.X=="Sine"])), fullresults$OMC[fullresults$DGM.X=="Sine"])
cor(log(log(fullresults$aloss50[fullresults$DGM.X=="Sine"])), fullresults$OMC[fullresults$DGM.X=="Sine"])

plot(log(log(fullresults$aloss50[fullresults$DGM.X=="Elbow"])), fullresults$OMC[fullresults$DGM.X=="Elbow"])
cor(log(log(fullresults$aloss50[fullresults$DGM.X=="Elbow"])), fullresults$OMC[fullresults$DGM.X=="Elbow"])

#Broken down by a few different FLCs
plot(log(log(fullresults$aloss50[fullresults$FLC==3])), fullresults$OMC[fullresults$FLC==2])

corr <- rep(NA, 144)
for (i in 1:144){
  corr[i] <- cor(log(log(fullresults$aloss50[fullresults$FLC==i])), fullresults$OMC[fullresults$FLC==i])
}
hist(corr, breaks=20)

