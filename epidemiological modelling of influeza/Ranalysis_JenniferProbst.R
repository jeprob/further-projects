library(readr)
library(tidyverse)
library(poweRlaw)

#read in data
Totaldata <- read_csv("C:/Users/probs/Desktop/Compbio/Semesterprojekt/Totaldatacsv.csv")
Totaldata$City <- as.factor(Totaldata$City)
Totaldata$Transportation <- as.factor(Totaldata$Transportation)

#model for infected
fit11 <- nls(Infected ~ a*log(Force*b), data=Totaldata, start=list(a=700, b=12))
summary(fit11)
cor(Totaldata$Infected, predict(fit11))

new_data <- expand.grid(Force=seq(min(Totaldata$Force), max(Totaldata$Force), length=1000))
p1 <- predict(fit11, newdata=new_data, interval="confidence")
n1 <- cbind(new_data, p1)
n1

ggplot(data=Totaldata) +
  geom_point(aes(x=Force, y=Infected, colour=City)) +
  xlab("Force of Infection") +
  ylab("Total Number of Infected People") +
  ggtitle("Force of Infection vs Prevalence") +
  theme_bw()+ 
  geom_smooth(data=n1, mapping=aes(x=Force, y=p1), stat="identity")

#fit for infected with high and low migration:
migrhigh <- filter(Totaldata, Transportation=='high')
migrlow <- filter(Totaldata, Transportation=='low')

fit1ml <- nls(Infected ~ a*log(Force*b), data=migrlow, start=list(a=700, b=12))
summary(fit1ml)
xnew <- expand.grid(Force=seq(min(Totaldata$Force), max(Totaldata$Force), length=1000))
p1ml <- predict(fit1ml, newdata=xnew, interval="confidence")
n1ml <- cbind(new_data, p1ml)

fit1mh <- nls(Infected ~ a*log(Force*b), data=migrhigh, start=list(a=700, b=12))
summary(fit1mh)
xnew <- expand.grid(Force=seq(min(Totaldata$Force), max(Totaldata$Force), length=1000))
p1mh <- predict(fit1mh, newdata=xnew, interval="confidence")
n1mh <- cbind(new_data, p1mh)

ggplot(Totaldata) +
  geom_point(aes(x=Force, y=Infected, colour=Transportation)) +
  xlab("Force of Infection") +
  ylab("Prevalence") +
  ggtitle("Prevalence vs Time of Pandemic") +
  geom_smooth(data=n1ml, mapping=aes(x=Force, y=p1ml), stat="identity", color='deepskyblue') +
  geom_smooth(data=n1mh, mapping=aes(x=Force, y=p1mh), stat="identity", color='salmon') +
  theme_bw()

#model for maxinfected
fit21 <- nls(MaxInfected ~ a*log(Force*b), data=Totaldata, start=list(a=5, b=20))
summary(fit21)
cor(Totaldata$Infected, predict(fit21))

new_data <- expand.grid(Force=seq(min(Totaldata$Force), max(Totaldata$Force), length=1000))
p2 <- predict(fit21, newdata=new_data, interval="confidence")
n2 <- cbind(new_data, p2)
n2

ggplot(Totaldata) +
  geom_point(aes(x=Force, y=MaxInfected, colour=City)) +
  xlab("Force of Infection") +
  ylab("Maximum of Infected People") +
  ggtitle("Force of Infection vs Maximum of Infected People") +
  theme_bw()+
  geom_smooth(data=n2, mapping=aes(x=Force, y=p2), stat="identity")

#model for time of pandemic
fit31 <- nls(Time ~ a* exp(-b*Force), data=Totaldata, start=list(a=10000, b=7))
fit31 <- nls(Time ~ a* Force**b, data=Totaldata, start=list(a=1000, b=0.1))

summary(fit31)
xnew <- expand.grid(Force=seq(min(Totaldata$Force), max(Totaldata$Force), length=1000))
p3 <- predict(fit31, newdata=xnew, interval="confidence")
n3 <- cbind(new_data, p3)

ggplot(Totaldata) +
  geom_point(aes(x=Force, y=Time, colour=City)) +
  xlab("Force of Infection") +
  ylab("Time of Pandemic") +
  ggtitle("Force of Infection vs Time of Pandemic") +
  geom_smooth(data=n3, mapping=aes(x=Force, y=p3), stat="identity") +
  theme_bw()

#fit for time of pandemic with high and low migration::
migrhigh <- filter(Totaldata, Transportation=='high')
migrlow <- filter(Totaldata, Transportation=='low')

fit3ml <- nls(Time ~ a* exp(-b*Force), data=migrlow, start=list(a=10000, b=7))
summary(fit3ml)
xnew <- expand.grid(Force=seq(min(Totaldata$Force), max(Totaldata$Force), length=1000))
p3ml <- predict(fit3ml, newdata=new_data, interval="confidence")
n3ml <- cbind(new_data, p3ml)

fit3mh <- nls(Time ~ a* exp(-b*Force), data=migrhigh, start=list(a=10000, b=7))
summary(fit3mh)
cor(Totaldata$Infected, predict(fit31))
xnew <- expand.grid(Force=seq(min(Totaldata$Force), max(Totaldata$Force), length=1000))
p3mh <- predict(fit3mh, newdata=new_data, interval="confidence")
n3mh <- cbind(new_data, p3mh)

anova(fit3mh, fit3ml)
anova(fit1mh, fit1ml)

ggplot(Totaldata) +
  geom_point(aes(x=Force, y=Time, colour=Transportation)) +
  xlab("Force of Infection") +
  ylab("Time of Pandemic") +
  ggtitle("Force of Infection vs Time of Pandemic") +
  geom_smooth(data=n3ml, mapping=aes(x=Force, y=p3ml), stat="identity", color='deepskyblue') +
  geom_smooth(data=n3mh, mapping=aes(x=Force, y=p3mh), stat="identity", color='salmon') +
  theme_bw()