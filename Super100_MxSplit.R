######################################################################################
# This code is to estimate Mx and split it in small intervals
# Updated by: Jesus-Adrian Alvarez
# Date: 03-12-2020
######################################################################################
# Clear workspace
rm(list = ls())

library(lubridate)
library(tidyverse)
library(data.table)
library(survival)
library(changepoint)
library(lmvar)
library(hrbrthemes)
#
#
# 0. Read the raw data ----------------------------------------------------

setwd("C:/Users/jmartinez/OneDrive - Syddansk Universitet/Super100/R")
#load("HMDData.RData")
IDL <- data.table(read_csv("IDLOfficialFinal.csv"))

IDL$Country[IDL$Country== "DEU"] <- "GER"
#IDL$Country[IDL$Country== "ENW"] <- "GRB"
source("Super100_FUN.R")
#IDL <- subset(IDL, id != "FR1187502211DSC42") # Jean Clament

# 1.1 Calculate ages ------------------------------------------------------

# Date of birth
IDL$DoB <- as.Date(paste(IDL$Db,IDL$Mb, IDL$Yb, sep = "-"), format = "%d-%m-%Y") 
# Date of entry to the DB
IDL$DoE <- as.Date(paste(IDL$De,IDL$Me, IDL$Ye, sep = "-"), format = "%d-%m-%Y") 
# Date when they turned 105 years old
IDL$DoH <- as.Date(paste(IDL$Db,IDL$Mb, IDL$Yb+105, sep = "-"), format = "%d-%m-%Y") 
# Calculate the beginning and end of observations for each country
# Beginning
b <- IDL[,min(DoE), by = list(Country)] 
names(b)[2] <- "b" 
# End
e <- IDL[,max(DoE), by = list(Country)]
names(e)[2] <- "e" 
# Merge them
IDL <- merge(IDL, b, by = "Country")
IDL <- merge(IDL, e, by = "Country")

# Here I calculate different ages and moments according to the observation schemes
# The age when they were observed (either dead or alive)
IDL$Age <-  time_length(difftime(IDL$DoE, IDL$DoB), "years")
# Remove those that didn't turn 105 years old
IDL$Age.at.b <- time_length(difftime(IDL$b, IDL$DoB), "years")

IDL <- subset(IDL, Age>= 105)

# Difference between the start of the observation window in each country
# minus the date when they turned 105 years
# There are some individuals that turned 105 years before b so they have negative "Lag"
IDL$Lag <-  time_length(difftime(IDL$DoH, IDL$b), "years")

# How much time do they stay in the observation window
# I took the absolute value because of the Lag is bigger than the Complete Stay
# They have been followed for very short times
# And they entered to the DB very late
IDL$Stay <- abs(time_length(difftime(IDL$DoE, IDL$b), "years") -pmax(IDL$Lag,0))


# 1.2 Put the data in the Surv format -------------------------------------
IDL$entry <- ifelse(IDL$Lag< 0, IDL$Age.at.b, 105)

# Time they spent in the observation window
IDL$exit <- IDL$entry + IDL$Stay

# If they experience the event (dead or alive)
IDL$delta <- ifelse( IDL$Status == "Dead", 1,0)
rIDL <- subset(IDL, Yb >= 1890)


# The actual database that we are going to use for the calculations
dati <- IDL
# Take just those that have some contribution to the likelihood
# This is to avoid warnings in the Surv Split
dati <- subset(dati,subset=dati$entry<dati$exit)
dati <- dati[order(dati$exit),]
dati$time <- dati$exit - dati$entry

# 2.1 Split the data into small intervals ---------------------------------

max(dati$exit)-min(dati$entry) # This is the total interval of ages
# Here the bins have to be split in relation to the ages
nx <-0.25 # Size of the interval
min <- nx
max <- round(max(dati$exit)-min(dati$entry),1)
tauj <- seq(min, max, nx)+105 # This is the age interval
M <- length(tauj)
range(tauj)


dati1 <- survSplit(Surv(entry, exit, delta)~.,data=dati, cut=tauj,
                    start="entry",  end="exit", 
                    event="delta", 
                    episode="interval")
 
dati1$y.new <- dati1$exit-dati1$entry
dati1$interval <- as.factor(dati1$interval)

dati1$Cohort5 <- ifelse(dati1$Yb>=1845 & dati1$Yb <=1849, "1845-49",dati1$Yb )

dati1$Cohort5 <- ifelse(dati1$Yb>=1850 & dati1$Yb <=1854, "1850-54",dati1$Cohort5 )
dati1$Cohort5 <- ifelse(dati1$Yb>=1855 & dati1$Yb <=1859, "1855-59",dati1$Cohort5 )


dati1$Cohort5 <- ifelse(dati1$Yb>=1860 & dati1$Yb <=1864, "1860-64",dati1$Cohort5 )
dati1$Cohort5 <- ifelse(dati1$Yb>=1865 & dati1$Yb <=1869, "1865-69",dati1$Cohort5 )

dati1$Cohort5 <- ifelse(dati1$Yb>=1870 & dati1$Yb <=1874, "1870-74",dati1$Cohort5 )
dati1$Cohort5 <- ifelse(dati1$Yb>=1875 & dati1$Yb <=1879, "1875-79",dati1$Cohort5 )

dati1$Cohort5 <- ifelse(dati1$Yb>=1880 & dati1$Yb <=1884, "1880-84",dati1$Cohort5 )
dati1$Cohort5 <- ifelse(dati1$Yb>=1885 & dati1$Yb <=1889, "1885-89",dati1$Cohort5 )

dati1$Cohort5 <- ifelse(dati1$Yb>=1890 & dati1$Yb <=1894, "1890-94",dati1$Cohort5 )
dati1$Cohort5 <- ifelse(dati1$Yb>=1895 & dati1$Yb <=1899, "1895-99",dati1$Cohort5 )

dati1$Cohort5 <- ifelse(dati1$Yb>=1900 & dati1$Yb <=1904, "1900-04",dati1$Cohort5 )
dati1$Cohort5 <- ifelse(dati1$Yb>=1905 & dati1$Yb <=1909, "1905-09",dati1$Cohort5 )

dati1$Cohort5 <- ifelse(dati1$Yb>=1910 & dati1$Yb <=1912, "1910-12",dati1$Cohort5 )


# Split country by country - All cohorts

ALL.t <- subset(dati1, Country != "USA")
EUR.t <- subset(dati1, Country %in% c("FIN","NOR", "SWE", "DNK")) # Nordic Coutries

EUF.t <- subset(dati1, Country %in% c("AUT", "ITA", "ESP", "CHE")) # Southern Europe
FRA.t <- subset(dati1, Country %in% c("FRA"))
GRB.t <- subset(dati1, Country %in% c("GRB"))
JPN.t <- subset(dati1, Country %in% c("JPN"))
AME.t <- subset(dati1, Country %in% c("CAN", "USA"))
USA.t <- subset(dati1, Country %in% c("USA"))
CAN.t <- subset(dati1, Country %in% c("CAN"))
BEL.t <- subset(dati1, Country %in% c("BEL"))
GER.t <- subset(dati1, Country %in% c("GER"))
AUT.t <- subset(dati1, Country %in% c("AUT"))
ITA.t <- subset(dati1, Country %in% c("ITA"))
FIN.t <- subset(dati1, Country %in% c("FIN"))
NOR.t <- subset(dati1, Country %in% c("NOR"))
ESP.t <- subset(dati1, Country %in% c("ESP"))
SWE.t <- subset(dati1, Country %in% c("SWE"))
CHE.t <- subset(dati1, Country %in% c("CHE"))
DNK.t <- subset(dati1, Country %in% c("DNK"))


ALL.f <- subset(ALL.t, Sex == "F")
EUR.f <- subset(EUR.t, Sex == "F")
EUF.f <- subset(EUF.t, Sex == "F")
FRA.f <- subset(FRA.t, Sex == "F")
GRB.f <- subset(GRB.t, Sex == "F")
JPN.f <- subset(JPN.t, Sex == "F")
AME.f <- subset(AME.t, Sex == "F")
USA.f <- subset(USA.t, Sex == "F")
CAN.f <- subset(CAN.t, Sex == "F")
BEL.f <- subset(BEL.t, Sex == "F")
GER.f <- subset(GER.t, Sex == "F")
AUT.f <- subset(AUT.t, Sex == "F")
ITA.f <- subset(ITA.t, Sex == "F")
FIN.f <- subset(FIN.t, Sex == "F")
NOR.f <- subset(NOR.t, Sex == "F")
ESP.f <- subset(ESP.t, Sex == "F")
SWE.f <- subset(SWE.t, Sex == "F")
CHE.f <- subset(CHE.t, Sex == "F")
DNK.f <- subset(DNK.t, Sex == "F")

ALL.m <- subset(ALL.t, Sex == "M")
EUR.m <- subset(EUR.t, Sex == "M")
EUF.m <- subset(EUF.t, Sex == "M")
FRA.m <- subset(FRA.t, Sex == "M")
GRB.m <- subset(GRB.t, Sex == "M")
JPN.m <- subset(JPN.t, Sex == "M")
AME.m <- subset(AME.t, Sex == "M")
USA.m <- subset(USA.t, Sex == "M")
CAN.m <- subset(CAN.t, Sex == "M")
BEL.m <- subset(BEL.t, Sex == "M")
GER.m <- subset(GER.t, Sex == "M")
AUT.m <- subset(AUT.t, Sex == "M")
ITA.m <- subset(ITA.t, Sex == "M")
FIN.m <- subset(FIN.t, Sex == "M")
NOR.m <- subset(NOR.t, Sex == "M")
ESP.m <- subset(ESP.t, Sex == "M")
SWE.m <- subset(SWE.t, Sex == "M")
CHE.m <- subset(CHE.t, Sex == "M")
DNK.m <- subset(DNK.t, Sex == "M")

ALL.r <- subset(ALL.t, Yb >= 1890)
EUR.r <- subset(EUR.t, Yb >= 1890)
EUF.r <- subset(EUF.t, Yb >= 1890)
FRA.r <- subset(FRA.t, Yb >= 1890)
GRB.r <- subset(GRB.t, Yb >= 1890)
JPN.r <- subset(JPN.t, Yb >= 1890)
AME.r <- subset(AME.t, Yb >= 1890)
USA.r <- subset(USA.t, Yb >= 1890)
CAN.r <- subset(CAN.t, Yb >= 1890)
BEL.r <- subset(BEL.t, Yb >= 1890)
GER.r <- subset(GER.t, Yb >= 1890)
AUT.r <- subset(AUT.t, Yb >= 1890)
ITA.r <- subset(ITA.t, Yb >= 1890)
FIN.r <- subset(FIN.t, Yb >= 1890)
NOR.r <- subset(NOR.t, Yb >= 1890)
ESP.r <- subset(ESP.t, Yb >= 1890)
SWE.r <- subset(SWE.t, Yb >= 1890)
CHE.r <- subset(CHE.t, Yb >= 1890)
DNK.r <- subset(DNK.t, Yb >= 1890)

FRA80.t <- subset(dati1, Country %in% c("FRA") & Cohort5== "1880-84")
FRA85.t <- subset(dati1, Country %in% c("FRA") & Cohort5== "1885-89")
FRA90.t <- subset(dati1, Country %in% c("FRA") & Cohort5== "1890-94")
FRA95.t <- subset(dati1, Country %in% c("FRA") & Cohort5== "1895-99")
FRA00.t <- subset(dati1, Country %in% c("FRA") & Cohort5== "1900-04")
FRA05.t <- subset(dati1, Country %in% c("FRA") & Cohort5== "1905-09")

# Estimate hazard


mx.ALL.t <- mx.calc(entry = ALL.t$entry, exit = ALL.t$exit, delta = ALL.t$delta, interval = ALL.t$interval,s.age = 105, nx = 0.25)
mx.EUR.t <- mx.calc(entry = EUR.t$entry, exit = EUR.t$exit, delta = EUR.t$delta, interval = EUR.t$interval,s.age = 105, nx = 0.25)
mx.EUF.t <- mx.calc(entry = EUF.t$entry, exit = EUF.t$exit, delta = EUF.t$delta, interval = EUF.t$interval,s.age = 105, nx = 0.25)
mx.FRA.t <- mx.calc(entry = FRA.t$entry, exit = FRA.t$exit, delta = FRA.t$delta, interval = FRA.t$interval,s.age = 105, nx = 0.25)
mx.GRB.t <- mx.calc(entry = GRB.t$entry, exit = GRB.t$exit, delta = GRB.t$delta, interval = GRB.t$interval,s.age = 105, nx = 0.25)
mx.JPN.t <- mx.calc(entry = JPN.t$entry, exit = JPN.t$exit, delta = JPN.t$delta, interval = JPN.t$interval,s.age = 105, nx = 0.25)
mx.AME.t <- mx.calc(entry = AME.t$entry, exit = AME.t$exit, delta = AME.t$delta, interval = AME.t$interval,s.age = 105, nx = 0.25)
mx.USA.t <- mx.calc(entry = USA.t$entry, exit = USA.t$exit, delta = USA.t$delta, interval = USA.t$interval,s.age = 105, nx = 0.25)
mx.CAN.t <- mx.calc(entry = CAN.t$entry, exit = CAN.t$exit, delta = CAN.t$delta, interval = CAN.t$interval,s.age = 105, nx = 0.25)
mx.BEL.t <- mx.calc(entry = BEL.t$entry, exit = BEL.t$exit, delta = BEL.t$delta, interval = BEL.t$interval,s.age = 105, nx = 0.25)
mx.GER.t <- mx.calc(entry = GER.t$entry, exit = GER.t$exit, delta = GER.t$delta, interval = GER.t$interval,s.age = 105, nx = 0.25)
mx.AUT.t <- mx.calc(entry = AUT.t$entry, exit = AUT.t$exit, delta = AUT.t$delta, interval = AUT.t$interval,s.age = 105, nx = 0.25)
mx.ITA.t <- mx.calc(entry = ITA.t$entry, exit = ITA.t$exit, delta = ITA.t$delta, interval = ITA.t$interval,s.age = 105, nx = 0.25)
mx.FIN.t <- mx.calc(entry = FIN.t$entry, exit = FIN.t$exit, delta = FIN.t$delta, interval = FIN.t$interval,s.age = 105, nx = 0.25)
mx.NOR.t <- mx.calc(entry = NOR.t$entry, exit = NOR.t$exit, delta = NOR.t$delta, interval = NOR.t$interval,s.age = 105, nx = 0.25)
mx.ESP.t <- mx.calc(entry = ESP.t$entry, exit = ESP.t$exit, delta = ESP.t$delta, interval = ESP.t$interval,s.age = 105, nx = 0.25)
mx.SWE.t <- mx.calc(entry = SWE.t$entry, exit = SWE.t$exit, delta = SWE.t$delta, interval = SWE.t$interval,s.age = 105, nx = 0.25)
mx.CHE.t <- mx.calc(entry = CHE.t$entry, exit = CHE.t$exit, delta = CHE.t$delta, interval = CHE.t$interval,s.age = 105, nx = 0.25)
mx.DNK.t <- mx.calc(entry = DNK.t$entry, exit = DNK.t$exit, delta = DNK.t$delta, interval = DNK.t$interval,s.age = 105, nx = 0.25)


mx.ALL.f <- mx.calc(entry = ALL.f$entry, exit = ALL.f$exit, delta = ALL.f$delta, interval = ALL.f$interval,s.age = 105, nx = 0.25)
mx.EUR.f <- mx.calc(entry = EUR.f$entry, exit = EUR.f$exit, delta = EUR.f$delta, interval = EUR.f$interval,s.age = 105, nx = 0.25)
mx.EUF.f <- mx.calc(entry = EUF.f$entry, exit = EUF.f$exit, delta = EUF.f$delta, interval = EUF.f$interval,s.age = 105, nx = 0.25)
mx.FRA.f <- mx.calc(entry = FRA.f$entry, exit = FRA.f$exit, delta = FRA.f$delta, interval = FRA.f$interval,s.age = 105, nx = 0.25)
mx.GRB.f <- mx.calc(entry = GRB.f$entry, exit = GRB.f$exit, delta = GRB.f$delta, interval = GRB.f$interval,s.age = 105, nx = 0.25)
mx.JPN.f <- mx.calc(entry = JPN.f$entry, exit = JPN.f$exit, delta = JPN.f$delta, interval = JPN.f$interval,s.age = 105, nx = 0.25)
mx.AME.f <- mx.calc(entry = AME.f$entry, exit = AME.f$exit, delta = AME.f$delta, interval = AME.f$interval,s.age = 105, nx = 0.25)
mx.USA.f <- mx.calc(entry = USA.f$entry, exit = USA.f$exit, delta = USA.f$delta, interval = USA.f$interval,s.age = 105, nx = 0.25)
mx.CAN.f <- mx.calc(entry = CAN.f$entry, exit = CAN.f$exit, delta = CAN.f$delta, interval = CAN.f$interval,s.age = 105, nx = 0.25)
mx.BEL.f <- mx.calc(entry = BEL.f$entry, exit = BEL.f$exit, delta = BEL.f$delta, interval = BEL.f$interval,s.age = 105, nx = 0.25)
mx.GER.f <- mx.calc(entry = GER.f$entry, exit = GER.f$exit, delta = GER.f$delta, interval = GER.f$interval,s.age = 105, nx = 0.25)
mx.AUT.f <- mx.calc(entry = AUT.f$entry, exit = AUT.f$exit, delta = AUT.f$delta, interval = AUT.f$interval,s.age = 105, nx = 0.25)
mx.ITA.f <- mx.calc(entry = ITA.f$entry, exit = ITA.f$exit, delta = ITA.f$delta, interval = ITA.f$interval,s.age = 105, nx = 0.25)
mx.FIN.f <- mx.calc(entry = FIN.f$entry, exit = FIN.f$exit, delta = FIN.f$delta, interval = FIN.f$interval,s.age = 105, nx = 0.25)
mx.NOR.f <- mx.calc(entry = NOR.f$entry, exit = NOR.f$exit, delta = NOR.f$delta, interval = NOR.f$interval,s.age = 105, nx = 0.25)
mx.ESP.f <- mx.calc(entry = ESP.f$entry, exit = ESP.f$exit, delta = ESP.f$delta, interval = ESP.f$interval,s.age = 105, nx = 0.25)
mx.SWE.f <- mx.calc(entry = SWE.f$entry, exit = SWE.f$exit, delta = SWE.f$delta, interval = SWE.f$interval,s.age = 105, nx = 0.25)
mx.CHE.f <- mx.calc(entry = CHE.f$entry, exit = CHE.f$exit, delta = CHE.f$delta, interval = CHE.f$interval,s.age = 105, nx = 0.25)
mx.DNK.f <- mx.calc(entry = DNK.f$entry, exit = DNK.f$exit, delta = DNK.f$delta, interval = DNK.f$interval,s.age = 105, nx = 0.25)



mx.ALL.m <- mx.calc(entry = ALL.m$entry, exit = ALL.m$exit, delta = ALL.m$delta, interval = ALL.m$interval,s.age = 105, nx = 0.25)
mx.EUR.m <- mx.calc(entry = EUR.m$entry, exit = EUR.m$exit, delta = EUR.m$delta, interval = EUR.m$interval,s.age = 105, nx = 0.25)
mx.EUF.m <- mx.calc(entry = EUF.m$entry, exit = EUF.m$exit, delta = EUF.m$delta, interval = EUF.m$interval,s.age = 105, nx = 0.25)
mx.FRA.m <- mx.calc(entry = FRA.m$entry, exit = FRA.m$exit, delta = FRA.m$delta, interval = FRA.m$interval,s.age = 105, nx = 0.25)
mx.GRB.m <- mx.calc(entry = GRB.m$entry, exit = GRB.m$exit, delta = GRB.m$delta, interval = GRB.m$interval,s.age = 105, nx = 0.25)
mx.JPN.m <- mx.calc(entry = JPN.m$entry, exit = JPN.m$exit, delta = JPN.m$delta, interval = JPN.m$interval,s.age = 105, nx = 0.25)
mx.AME.m <- mx.calc(entry = AME.m$entry, exit = AME.m$exit, delta = AME.m$delta, interval = AME.m$interval,s.age = 105, nx = 0.25)
mx.USA.m <- mx.calc(entry = USA.m$entry, exit = USA.m$exit, delta = USA.m$delta, interval = USA.m$interval,s.age = 105, nx = 0.25)
mx.CAN.m <- mx.calc(entry = CAN.m$entry, exit = CAN.m$exit, delta = CAN.m$delta, interval = CAN.m$interval,s.age = 105, nx = 0.25)
mx.BEL.m <- mx.calc(entry = BEL.m$entry, exit = BEL.m$exit, delta = BEL.m$delta, interval = BEL.m$interval,s.age = 105, nx = 0.25)
mx.GER.m <- mx.calc(entry = GER.m$entry, exit = GER.m$exit, delta = GER.m$delta, interval = GER.m$interval,s.age = 105, nx = 0.25)
mx.AUT.m <- mx.calc(entry = AUT.m$entry, exit = AUT.m$exit, delta = AUT.m$delta, interval = AUT.m$interval,s.age = 105, nx = 0.25)
mx.ITA.m <- mx.calc(entry = ITA.m$entry, exit = ITA.m$exit, delta = ITA.m$delta, interval = ITA.m$interval,s.age = 105, nx = 0.25)
#mx.FIN.m <- mx.calc(entry = FIN.m$entry, exit = FIN.m$exit, delta = FIN.m$delta, interval = FIN.m$interval,s.age = 105, nx = 0.25)
mx.NOR.m <- mx.calc(entry = NOR.m$entry, exit = NOR.m$exit, delta = NOR.m$delta, interval = NOR.m$interval,s.age = 105, nx = 0.25)
mx.ESP.m <- mx.calc(entry = ESP.m$entry, exit = ESP.m$exit, delta = ESP.m$delta, interval = ESP.m$interval,s.age = 105, nx = 0.25)
mx.SWE.m <- mx.calc(entry = SWE.m$entry, exit = SWE.m$exit, delta = SWE.m$delta, interval = SWE.m$interval,s.age = 105, nx = 0.25)
mx.CHE.m <- mx.calc(entry = CHE.m$entry, exit = CHE.m$exit, delta = CHE.m$delta, interval = CHE.m$interval,s.age = 105, nx = 0.25)
mx.DNK.m <- mx.calc(entry = DNK.m$entry, exit = DNK.m$exit, delta = DNK.m$delta, interval = DNK.m$interval,s.age = 105, nx = 0.25)


mx.ALL.r <- mx.calc(entry = ALL.r$entry, exit = ALL.r$exit, delta = ALL.r$delta, interval = ALL.r$interval,s.age = 105, nx = 0.25)
mx.EUR.r <- mx.calc(entry = EUR.r$entry, exit = EUR.r$exit, delta = EUR.r$delta, interval = EUR.r$interval,s.age = 105, nx = 0.25)
mx.EUF.r <- mx.calc(entry = EUF.r$entry, exit = EUF.r$exit, delta = EUF.r$delta, interval = EUF.r$interval,s.age = 105, nx = 0.25)
mx.FRA.r <- mx.calc(entry = FRA.r$entry, exit = FRA.r$exit, delta = FRA.r$delta, interval = FRA.r$interval,s.age = 105, nx = 0.25)
mx.GRB.r <- mx.calc(entry = GRB.r$entry, exit = GRB.r$exit, delta = GRB.r$delta, interval = GRB.r$interval,s.age = 105, nx = 0.25)
mx.JPN.r <- mx.calc(entry = JPN.r$entry, exit = JPN.r$exit, delta = JPN.r$delta, interval = JPN.r$interval,s.age = 105, nx = 0.25)
mx.AME.r <- mx.calc(entry = AME.r$entry, exit = AME.r$exit, delta = AME.r$delta, interval = AME.r$interval,s.age = 105, nx = 0.25)
mx.USA.r <- mx.calc(entry = USA.r$entry, exit = USA.r$exit, delta = USA.r$delta, interval = USA.r$interval,s.age = 105, nx = 0.25)
mx.CAN.r <- mx.calc(entry = CAN.r$entry, exit = CAN.r$exit, delta = CAN.r$delta, interval = CAN.r$interval,s.age = 105, nx = 0.25)
mx.BEL.r <- mx.calc(entry = BEL.r$entry, exit = BEL.r$exit, delta = BEL.r$delta, interval = BEL.r$interval,s.age = 105, nx = 0.25)
mx.GER.r <- mx.calc(entry = GER.r$entry, exit = GER.r$exit, delta = GER.r$delta, interval = GER.r$interval,s.age = 105, nx = 0.25)
mx.AUT.r <- mx.calc(entry = AUT.r$entry, exit = AUT.r$exit, delta = AUT.r$delta, interval = AUT.r$interval,s.age = 105, nx = 0.25)
mx.ITA.r <- mx.calc(entry = ITA.r$entry, exit = ITA.r$exit, delta = ITA.r$delta, interval = ITA.r$interval,s.age = 105, nx = 0.25)
mx.FIN.r <- mx.calc(entry = FIN.r$entry, exit = FIN.r$exit, delta = FIN.r$delta, interval = FIN.r$interval,s.age = 105, nx = 0.25)
mx.NOR.r <- mx.calc(entry = NOR.r$entry, exit = NOR.r$exit, delta = NOR.r$delta, interval = NOR.r$interval,s.age = 105, nx = 0.25)
mx.ESP.r <- mx.calc(entry = ESP.r$entry, exit = ESP.r$exit, delta = ESP.r$delta, interval = ESP.r$interval,s.age = 105, nx = 0.25)
mx.SWE.r <- mx.calc(entry = SWE.r$entry, exit = SWE.r$exit, delta = SWE.r$delta, interval = SWE.r$interval,s.age = 105, nx = 0.25)
mx.CHE.r <- mx.calc(entry = CHE.r$entry, exit = CHE.r$exit, delta = CHE.r$delta, interval = CHE.r$interval,s.age = 105, nx = 0.25)
mx.DNK.r <- mx.calc(entry = DNK.r$entry, exit = DNK.r$exit, delta = DNK.r$delta, interval = DNK.r$interval,s.age = 105, nx = 0.25)

mx.FRA80.t <- mx.calc(entry = FRA80.t$entry, exit = FRA80.t$exit, delta = FRA80.t$delta, interval = FRA80.t$interval,s.age = 105, nx = 0.25)
mx.FRA85.t <- mx.calc(entry = FRA85.t$entry, exit = FRA85.t$exit, delta = FRA85.t$delta, interval = FRA85.t$interval,s.age = 105, nx = 0.25)
mx.FRA90.t <- mx.calc(entry = FRA90.t$entry, exit = FRA90.t$exit, delta = FRA90.t$delta, interval = FRA90.t$interval,s.age = 105, nx = 0.25)
mx.FRA95.t <- mx.calc(entry = FRA95.t$entry, exit = FRA95.t$exit, delta = FRA95.t$delta, interval = FRA95.t$interval,s.age = 105, nx = 0.25)
mx.FRA00.t <- mx.calc(entry = FRA00.t$entry, exit = FRA00.t$exit, delta = FRA00.t$delta, interval = FRA00.t$interval,s.age = 105, nx = 0.25)
mx.FRA05.t <- mx.calc(entry = FRA05.t$entry, exit = FRA05.t$exit, delta = FRA05.t$delta, interval = FRA05.t$interval,s.age = 105, nx = 0.25)


mx.ALL.t$Sex.Country <- "Total.ALL"
mx.EUR.t$Sex.Country <- "Total.EUR"
mx.EUF.t$Sex.Country <- "Total.EUF"
mx.FRA.t$Sex.Country <- "Total.FRA"
mx.GRB.t$Sex.Country <- "Total.GRB"
mx.JPN.t$Sex.Country <- "Total.JPN"
mx.AME.t$Sex.Country <- "Total.AME"
mx.USA.t$Sex.Country <- "Total.USA"
mx.CAN.t$Sex.Country <- "Total.CAN"
mx.BEL.t$Sex.Country <- "Total.BEL"
mx.GER.t$Sex.Country <- "Total.GER"
mx.AUT.t$Sex.Country <- "Total.AUT"
mx.ITA.t$Sex.Country <- "Total.ITA"
mx.FIN.t$Sex.Country <- "Total.FIN"
mx.NOR.t$Sex.Country <- "Total.NOR"
mx.ESP.t$Sex.Country <- "Total.ESP"
mx.SWE.t$Sex.Country <- "Total.SWE"
mx.CHE.t$Sex.Country <- "Total.CHE"
mx.DNK.t$Sex.Country <- "Total.DNK"

mx.ALL.f$Sex.Country <- "Females.ALL"
mx.EUR.f$Sex.Country <- "Females.EUR"
mx.EUF.f$Sex.Country <- "Females.EUF"
mx.FRA.f$Sex.Country <- "Females.FRA"
mx.GRB.f$Sex.Country <- "Females.GRB"
mx.JPN.f$Sex.Country <- "Females.JPN"
mx.AME.f$Sex.Country <- "Females.AME"
mx.USA.f$Sex.Country <- "Females.USA"
mx.CAN.f$Sex.Country <- "Females.CAN"
mx.BEL.f$Sex.Country <- "Females.BEL"
mx.GER.f$Sex.Country <- "Females.GER"
mx.AUT.f$Sex.Country <- "Females.AUT"
mx.ITA.f$Sex.Country <- "Females.ITA"
mx.FIN.f$Sex.Country <- "Females.FIN"
mx.NOR.f$Sex.Country <- "Females.NOR"
mx.ESP.f$Sex.Country <- "Females.ESP"
mx.SWE.f$Sex.Country <- "Females.SWE"
mx.CHE.f$Sex.Country <- "Females.CHE"
mx.DNK.f$Sex.Country <- "Females.DNK"

mx.ALL.m$Sex.Country <- "Males.ALL"
mx.EUR.m$Sex.Country <- "Males.EUR"
mx.EUF.m$Sex.Country <- "Males.EUF"
mx.FRA.m$Sex.Country <- "Males.FRA"
mx.GRB.m$Sex.Country <- "Males.GRB"
mx.JPN.m$Sex.Country <- "Males.JPN"
mx.AME.m$Sex.Country <- "Males.AME"
mx.USA.m$Sex.Country <- "Males.USA"
mx.CAN.m$Sex.Country <- "Males.CAN"
mx.BEL.m$Sex.Country <- "Males.BEL"
mx.GER.m$Sex.Country <- "Males.GER"
mx.AUT.m$Sex.Country <- "Males.AUT"
mx.ITA.m$Sex.Country <- "Males.ITA"
#mx.FIN.m$Sex.Country <- "Males.FIN"
mx.NOR.m$Sex.Country <- "Males.NOR"
mx.ESP.m$Sex.Country <- "Males.ESP"
mx.SWE.m$Sex.Country <- "Males.SWE"
mx.CHE.m$Sex.Country <- "Males.CHE"
mx.DNK.m$Sex.Country <- "Males.DNK"

mx.ALL.r$Sex.Country <- "rTotal.ALL"
mx.EUR.r$Sex.Country <- "rTotal.EUR"
mx.EUF.r$Sex.Country <- "rTotal.EUF"
mx.FRA.r$Sex.Country <- "rTotal.FRA"
mx.GRB.r$Sex.Country <- "rTotal.GRB"
mx.JPN.r$Sex.Country <- "rTotal.JPN"
mx.AME.r$Sex.Country <- "rTotal.AME"
mx.USA.r$Sex.Country <- "rTotal.USA"
mx.CAN.r$Sex.Country <- "rTotal.CAN"
mx.BEL.r$Sex.Country <- "rTotal.BEL"
mx.GER.r$Sex.Country <- "rTotal.GER"
mx.AUT.r$Sex.Country <- "rTotal.AUT"
mx.ITA.r$Sex.Country <- "rTotal.ITA"
mx.FIN.r$Sex.Country <- "rTotal.FIN"
mx.NOR.r$Sex.Country <- "rTotal.NOR"
mx.ESP.r$Sex.Country <- "rTotal.ESP"
mx.SWE.r$Sex.Country <- "rTotal.SWE"
mx.CHE.r$Sex.Country <- "rTotal.CHE"
mx.DNK.r$Sex.Country <- "rTotal.DNK"

mx.FRA80.t$Sex.Country <- "FRA.1880-84"
mx.FRA85.t$Sex.Country <- "FRA.1885-89"
mx.FRA90.t$Sex.Country <- "FRA.1890-94"
mx.FRA95.t$Sex.Country <- "FRA.1895-99"
mx.FRA00.t$Sex.Country <- "FRA.1900-04"
mx.FRA05.t$Sex.Country <- "FRA.1905-09"


Mx <- rbind(mx.FRA.f,mx.CAN.f,mx.BEL.f,mx.GER.f,mx.AUT.f,mx.NOR.f,mx.DNK.f,
            mx.FRA.m,mx.CAN.m,mx.BEL.m,mx.GER.m,mx.AUT.m,mx.NOR.m,mx.DNK.m)

save(Mx, file="datiSplitFinal.RData")

