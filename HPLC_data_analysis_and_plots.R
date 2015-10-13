# this script is for analyzing HPLC results from my time course experiments
# by Kara Sarver 12.31.14   file is called "hplc.dat"

str(hplc.dat)
head(hplc.dat)
hplc.dat$Block <- as.factor(hplc.dat$Block)
str(hplc.dat)

#models and ANOVAs (RCBD for General Linear Model)
hplcPut.mod <-lm(Adj.Put ~ Treatment + Block, hplc.dat) # = 1.185e-06
anova(hplcPut.mod)

hplcSpd.mod <-lm(Adj.Spd ~ Treatment + Block, hplc.dat)  #0.08517
anova(hplcSpd.mod)

hplcSpm.mod <-lm(Adj.Spm ~ Treatment + Block, hplc.dat) #0.1963
anova(hplcSpm.mod)


############### Putrescine Analysis #################
#testing assumptions
###resids and predicted plot
hplc.dat$Put.resids <-residuals(hplcPut.mod)
hplc.dat$Put.preds <- predict(hplcPut.mod)
hplc.dat$Put.sq_preds <- hplc.dat$Put.preds^2
plot(Put.resids ~ Put.preds, data = hplc.dat, xlab = "Predicted", ylab = "Residual")   # definite funnel effect!

###Shapiro Wilke to test normality of resid
shapiro.test(hplc.dat$Put.resids) #very significant... p = 0.005073 **

###Levene's Test for homogeneity of variances
#library(car)
leveneTest(Adj.Put ~ Block, data = hplc.dat)      # p = 9.7646
leveneTest(Adj.Put ~ Treatment, data = hplc.dat)  # p = 0.08491 *
#both insignificant (0.7647, 0.08491.)

###Tukey 1df test for Non-additivity
hplcPut_1df.mod <- lm(Adj.Put ~ Treatment + Block + Put.sq_preds, hplc.dat)
anova(hplcPut_1df.mod)  #significant effects of sq_preds (0.0107)...

#VIOLATES ALL ASSUMPTIONS TESTS!
boxplot(Adj.Put ~ Treatment, data = hplc.dat, main = "Effect of Treatment on Adjusted Put Amount", xlab = "Treatment", ylab = "Adj Put Amount")
#can see that the smaller means have tiny variances while the high Put treatments have larger ones

# block means
mean.PutBlock <- tapply(hplc.dat$Adj.Put, list(hplc.dat$Block), mean)
mean.PutBlock     #block 4 is much lower than other blocks... why?
Block4 <- subset(hplc.dat, hplc.dat$Block == "4")
Block4$Adj.Put      #the constitutive treatments are relatively low (24.77 and 16.84)

#POWER TRANSFORMATION
meansPut <- aggregate(hplc.dat$Adj.Put, list(hplc.dat$Treatment), mean)
varsPut <- aggregate(hplc.dat$Adj.Put, list(hplc.dat$Treatment), var)
logmeansPut <- log10(meansPut$x)
logvarsPut <- log10(varsPut$x)
power.modPut <- lm(logvarsPut ~ logmeansPut)
summary(power.modPut)   #logmeansPut 2.9259
Putpower = 1-(2.9259/2)
Putpower  # Putpower = -0.46295
hplc.dat$power_Adj.Put <- (hplc.dat$Adj.Put)^(-0.46295)
hplcPut.modpower <-lm(power_Adj.Put ~ Treatment + Block, hplc.dat)
anova(hplcPut.modpower) #(p = 2.048e-13)
hplc.dat$Put.residspower <-residuals(hplcPut.modpower)
hplc.dat$Put.predspower <- predict(hplcPut.modpower)
hplc.dat$Put.sq_predspower <- hplc.dat$Put.predspower^2
plot(Put.residspower ~ Put.predspower, data = hplc.dat, xlab = "Power Transformed Predicted", ylab = "Power Transformed Residual", main = "Put power")
shapiro.test(hplc.dat$Put.residspower) #p=0.71
leveneTest(power_Adj.Put ~ Block, data = hplc.dat) #p=0.5572
leveneTest(power_Adj.Put ~ Treatment, data = hplc.dat) #p=0.9556
hplcPut_1df.modpower <- lm(power_Adj.Put ~ Treatment + Block + Put.sq_predspower, hplc.dat)
anova(hplcPut_1df.modpower) # p = 3.59e-13
library(agricolae)
tukey.Putpower <- HSD.test(lm(power_Adj.Put ~ Block + Treatment, hplc.dat), "Treatment")
tukey.Putpower
Putpowermeans <- aggregate(hplc.dat$power_Adj.Put, list(hplc.dat$Treatment), mean)
Putpowermeans
(1/-0.46295)                                  

Putpowermeans.dt <- c(38.8484,28.84421,10.1412,16.32468,2.785305,2.727532,2.851993,2.451353)

#power put barplot
hplc.dat$power_Adj.Put.dt <- (hplc.dat$power_Adj.Put)^(1/-0.46295)
mean.Put.dt <- tapply(hplc.dat$power_Adj.Put.dt, list(hplc.dat$Treatment), mean)
mean.Put.dt
sd.Put.dt <- tapply(hplc.dat$power_Adj.Put.dt, list(hplc.dat$Treatment), sd)
n.Put.dt <- tapply(hplc.dat$power_Adj.Put.dt, list(hplc.dat$Treatment), length)
se.Put.dt <- sd.Put.dt/(n.Put.dt)**(1/2)
se.Put.dt
names(mean.Put.dt) <- c("Cnst0","Cnst24","Ind_i24","Ind_i48","Ind_n24","Ind_n48","Ind0","Wt0")
Adj.Put.bp = barplot(mean.Put.dt)
mids <- barplot(mean.Put.dt, beside = TRUE, legend = FALSE, ylab = "Adjusted Putrescine Amount (nmol/g FW)", ylim = c(0,60), col=grey(c(0.3,0.3,0.6,0.6,0.9,0.9,0.9,1)))
arrows(mids, mean.Put.dt - se.Put.dt, mids, mean.Put.dt + se.Put.dt, code = 3, angle = 90, length = 0.1) 
letters = c("a","a","b","ab","c","c","c","c")
text(x=Adj.Put.bp, y= mean.Put.dt + se.Put.dt + 5, label=letters)
box()

############### Spermidine Analysis #################
#testing assumptions
###resids and predicted plot
hplc.dat$Spd.resids <-residuals(hplcSpd.mod)
hplc.dat$Spd.preds <- predict(hplcSpd.mod)
hplc.dat$Spd.sq_preds <- hplc.dat$Spd.preds^2
plot(Spd.resids ~ Spd.preds, data = hplc.dat, xlab = "Predicted", ylab = "Residual", main = "Not transformed Spd")

###Shapiro Wilke to test normality of resid
shapiro.test(hplc.dat$Spd.resids) #good, p = 0.9015

###Levene's Test for homogeneity of variances
#library(car)
leveneTest(Adj.Spd ~ Block, data = hplc.dat)  #good, p= 0.8033
leveneTest(Adj.Spd ~ Treatment, data = hplc.dat)   #good, p = 0.412

###Tukey 1df test for Non-additivity
hplcSpd_1df.mod <- lm(Adj.Spd ~ Treatment + Block + Spd.sq_preds, hplc.dat)
anova(hplcSpd_1df.mod)  #p sp preds = 0.22992

##Tukey Test for rankings
tukey.SpdLog <- HSD.test(lm(Adj.Spd ~ Block + Treatment, hplc.dat), "Treatment")
tukey.SpdLog

#barplots for Spd 
mean.Spd <- tapply(hplc.dat$Adj.Spd, list(hplc.dat$Treatment), mean)
mean.Spd
sd.Spd <- tapply(hplc.dat$Adj.Spd, list(hplc.dat$Treatment), sd)
n.Spd <- tapply(hplc.dat$Adj.Spd, list(hplc.dat$Treatment), length)
se.Spd <- sd.Spd/(n.Spd)**(1/2)
se.Spd
names(mean.Spd) <- c("Cnst0","Cnst24","Ind_i24","Ind_i48","Ind_n24","Ind_n48","Ind0","Wt0")
Spd.bp <- barplot(mean.Spd)
mids <- barplot(mean.Spd, beside = TRUE, legend = FALSE, ylab = "Adjusted Spermidine Amount (nmol/g FW)", ylim = c(0,10), col=grey(c(0.3,0.3,0.6,0.6,0.9,0.9,0.9,1)))
arrows(mids, mean.Spd - se.Spd, mids, mean.Spd + se.Spd, code = 3, angle = 90, length = 0.1)
letters.spd = c("a","a","a","a","a","a","a","a")
text(x=Spd.bp, y= mean.Spd+se.Spd+1, label=letters.spd)
box()


############### Spermine Analysis #################
#testing assumptions
###resids and predicted plot
hplc.dat$Spm.resids <-residuals(hplcSpm.mod)
hplc.dat$Spm.preds <- predict(hplcSpm.mod)
hplc.dat$Spm.sq_preds <- hplc.dat$Spm.preds^2
plot(Spm.resids ~ Spm.preds, data = hplc.dat, xlab = "Predicted", ylab = "Residual", main = "Not transformed Spm")   #also a funnel effect for Spm

###Shapiro Wilke to test normality of resid
shapiro.test(hplc.dat$Spm.resids) #good, p = 0.5498

###Levene's Test for homogeneity of variances
#library(car)
leveneTest(Adj.Spm ~ Block, data = hplc.dat)  #good, p = 0.8343
leveneTest(Adj.Spm ~ Treatment, data = hplc.dat)  # SIGNIFICANT! p = 0.001855

###Tukey 1df test for Non-additivity
hplcSpm_1df.mod <- lm(Adj.Spm ~ Treatment + Block + Spm.sq_preds, hplc.dat)
anova(hplcSpm_1df.mod)  #SIGNIFICANT!! p sp preds = 0.00928  

#VIOLATES homogeneity and non-additivity
boxplot(Adj.Spm ~ Treatment, data = hplc.dat, main = "Effect of Treatment on Adjusted Spm Amount", xlab = "Treatment", ylab = "Adj Spm Amount")
#not as clear of a pattern as with Put, WldT0 is way more variable than any of other treatments

# block means
mean.SpmBlock <- tapply(hplc.dat$Adj.Spm, list(hplc.dat$Block), mean)
mean.SpmBlock     #blocks are all pretty similar

#POWER TRANSFORMATION
meansSpm <- aggregate(hplc.dat$Adj.Spm, list(hplc.dat$Treatment), mean)
varsSpm <- aggregate(hplc.dat$Adj.Spm, list(hplc.dat$Treatment), var)
logmeansSpm <- log10(meansSpm$x)
logvarsSpm <- log10(varsSpm$x)
power.modSpm <- lm(logvarsSpm ~ logmeansSpm)
summary(power.modSpm)   #logmeansSpm -2.543
Spmpower = 1-(2.8418/2)
Spmpower   #Spmpower = -0.4209
hplc.dat$power_Adj.Spm <- (hplc.dat$Adj.Spm)^(-0.4209)
hplcSpm.modpower <-lm(power_Adj.Spm ~ Treatment + Block, hplc.dat)
anova(hplcSpm.modpower)
hplc.dat$Spm.residspower <-residuals(hplcSpm.modpower)
hplc.dat$Spm.predspower <- predict(hplcSpm.modpower)
hplc.dat$Spm.sq_predspower <- hplc.dat$Spm.predspower^2
plot(Spm.residspower ~ Spm.predspower, data = hplc.dat, xlab = "Power TransformedPredicted", ylab = "Power Transformed Residual ", main = "Spm power")
shapiro.test(hplc.dat$Spm.residspower) #p=0.1668
leveneTest(power_Adj.Spm ~ Block, data = hplc.dat) #p=0.78
leveneTest(power_Adj.Spm ~ Treatment, data = hplc.dat) #p=0.1156
hplcSpm_1df.modpower <- lm(power_Adj.Spm ~ Treatment + Block + Spm.sq_predspower, hplc.dat)
anova(hplcSpm_1df.modpower) # p=0.2560
tukey.Spmpower <- HSD.test(lm(power_Adj.Spm ~ Block + Treatment, hplc.dat), "Treatment")
tukey.Spmpower
Spmtpowermeans <- aggregate(hplc.dat$power_Adj.Put, list(hplc.dat$Treatment), mean)

#Spm power barplot
###(detransformed)
hplc.dat$dt.spm.means <- (hplc.dat$power_Adj.Spm)^(-1/0.4209)
hplc.dat$dt.CI.Spm.up <- (hplc.dat$power_Adj.Spm +0.5830452)^(-1/0.4209)
hplc.dat$dt.CI.Spm.down <- (hplc.dat$power_Adj.Spm-0.5830452)^(-1/0.4209)
head(hplc.dat)
dt.mean.Spm <- tapply(hplc.dat$dt.spm.means, list(hplc.dat$Treatment), mean)
dt.mean.Spm
dt.CI.Spm.ups <- tapply(hplc.dat$dt.CI.Spm.up, list(hplc.dat$Treatment), mean)
dt.CI.Spm.downs <- tapply(hplc.dat$dt.CI.Spm.down, list(hplc.dat$Treatment), mean)
sd.Spm.dt <- tapply(hplc.dat$dt.spm.means, list(hplc.dat$Treatment), sd)
n.Spm.dt <- tapply(hplc.dat$dt.spm.means, list(hplc.dat$Treatment), length)
se.Spm.dt <- sd.Spm.dt/(n.Spm.dt)**(1/2)
se.Spm.dt

names(dt.mean.Spm) <- c("Cnst0","Cnst24","Ind_i24","Ind_i48","Ind_n24","Ind_n48","Ind0","Wt0")
Spm.bp = barplot(dt.mean.Spm)
mids <- barplot(dt.mean.Spm, beside = TRUE, legend = FALSE, ylab = "Adjusted Spermidine Amount (nmol/g FW)", ylim = c(0,4), col=grey(c(0.3,0.3,0.6,0.6,0.9,0.9,0.9,1)))
arrows(mids, dt.mean.Spm - se.Spm.dt, mids, dt.mean.Spm + se.Spm.dt, code = 3, angle = 90, length = 0.1) 
lettersSpm = c("a","a","a","a","a","a","a","a")
text(x=Spm.bp, y= dt.mean.Spm + se.Spm.dt +0.3, label=lettersSpm)
box()
