######################################################################################
# This code is to estimate Mx and split it in small intervals
# Updated by: Jesus-Adrian Alvarez
# Date: 8-09-2020
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

# Source code
source("Super100_FUN.R")

setwd("C:/Users/jmartinez/OneDrive - Syddansk Universitet/Super100/R")
#setwd("C:/Users/FVillavicencio/Dropbox/2020/Semi-superCent/Code")
# IDL <- data.table(read_csv("IDLOfficial.csv"))

# US DATA
IDL <- read.csv('idl_usa_dead_110.csv', sep = ';', header = T)
IDL <- data.table(IDL)

# Select and rename columns
IDL <- IDL[, c("IDNUMBER", "UPDATE", "VALIDATIONTYPE", "SEX", "DCOUNTRY",
               "BDATE", "DDATE", "AGEYEARS", "DAYSSINCEBD")]
names(IDL) <- c("id", "Update", "Validation", "Sex", "Country",
                'Yb', 'Ye', 'AgeYears', 'DaysSinceBd')
IDL$Status <- 'Dead'

# Select min age
minAge <- 110

# Select Sex
sex <- c('F', 'M')
IDL <- IDL[IDL$Sex %in% sex, ]


#--------------------#
# INTERVALS OF BIRTH #
#--------------------#

# Date of birth
IDL$DoB <- as.Date(NA)
# Date of entry to the DB
IDL$DoE <- as.Date(NA)
# Date when they turned minAge years old
IDL$DoH <- as.Date(NA)

# Individuals reaching AGE at DEATH in CURRENT YEAR
idCurrent <- which(IDL$AgeYears == IDL$Ye - IDL$Yb)
# Individuals reaching AGE at DEATH in PREVIOUS YEAR
idPrevious <- which(IDL$AgeYears != IDL$Ye - IDL$Yb)

# EARLIEST date of LAST BIRTHDAY
minBDay <- as.Date(paste(1, 1, IDL$Ye, sep = "-"), format = "%d-%m-%Y")
minBDay[idPrevious] <- minBDay[idPrevious] - IDL$DaysSinceBd[idPrevious]
  
# EARLIEST DATE OF BIRTH
IDL$minBirth <- minBDay
year(IDL$minBirth) <- IDL$Yb

# LATEST date of LAST BIRTHDAY
maxBDay <- as.Date(paste(31, 12, IDL$Ye, sep = "-"), format = "%d-%m-%Y")
maxBDay[idCurrent] <- maxBDay[idCurrent] - IDL$DaysSinceBd[idCurrent]
year(maxBDay[idPrevious]) <- year(maxBDay[idPrevious]) - 1

# LATEST DATE OF BIRTH
IDL$maxBirth <- maxBDay
year(IDL$maxBirth) <- IDL$Yb

# Checks
which(is.na(IDL$minBirth))
which(is.na(IDL$maxBirth))
which(IDL$Yb != year(IDL$minBirth))
which(IDL$Yb != year(IDL$maxBirth))
which(IDL$minBirth > IDL$maxBirth)

# Number of simulations (times hazard is estiamted)
nsim1 <- 1000

# Number of simulations (for confidence intervals)
nsim2 <- 100

# RANGE (in days) between possible BIRTH DATES
birth <- IDL$maxBirth - IDL$minBirth
range(birth)

#  Matrix with days to be added to all individuals
addBirth <- matrix(0, nrow = nsim1, ncol = length(birth))
for (i in 1:length(birth)) {
  addBirth[, i] <- sample(0:birth[i], size = nsim1, replace = T)
}
rm(birth)

# Checks
hist(addBirth[, 1])
hist(addBirth[, 2])


#-------------------------------------------------------#
# INITIATE SIMULATION WITH IMPUTATION OF TIMES OF BIRTH #
#-------------------------------------------------------#

# Time resolution
nx <- 0.25

# Age groups
nages <- seq(minAge, 120, nx)

# LISTS to store POINT ESTIMATE of the HAZARD (obtained after IMPUTING TIMES OF BIRTH)
mx.f <- list()
if ('M' %in% sex) {
  mx.m <- list()
  mx.t <- list()
}

# START SIMULATION
Start <- Sys.time()
for (i in 1:nsim1) {

  # IMPUT EXACT DATE OF BIRTH (specific to each individual)
  IDL$DoB <- IDL$minBirth + addBirth[i, ]
  
  # DATE they TURN minAge
  IDL$DoH <- IDL$DoB
  year(IDL$DoH) <- year(IDL$DoH) + minAge
  # Correct LEAP years
  upd <- which(is.na(IDL$DoH))
  if (length(upd) > 0) {
    IDL$DoH[upd] <- 
      as.Date(paste(1, 3, year(IDL$DoB[upd]) + minAge, sep = "-"),
              format = "%d-%m-%Y")
  }
  
  # DATE of ENTRY in the DATABASE
  IDL$DoE <- IDL$DoB
  year(IDL$DoE) <- year(IDL$DoE) + IDL$AgeYears
  # Correct LEAP years
  upd <- which(is.na(IDL$DoE))
  if (length(upd) > 0) {
    IDL$DoE[upd] <- 
      as.Date(paste(1, 3, year(IDL$DoB[upd]) + IDL$AgeYears, sep = "-"),
              format = "%d-%m-%Y")
  }
  IDL$DoE <- IDL$DoE + IDL$DaysSinceBd
  
  # Checks
  # which(year(IDL$DoB) != IDL$Yb)
  # which(year(IDL$DoE) != IDL$Ye)
  # which(IDL$maxBirth < IDL$DoB)
  # which(IDL$minBirth > IDL$DoB)
  
  # Calculate the beginning and end of observations for each country
  # Beginning
  b <- IDL[,min(DoE), by = list(Country)] 
  names(b)[2] <- "b" 
  # End
  e <- IDL[,max(DoE), by = list(Country)]
  names(e)[2] <- "e" 
  # Merge them
  IDL2 <- merge(IDL, b, by = "Country")
  IDL2 <- merge(IDL2, e, by = "Country")
  
  # Here I calculate different ages and moments according to the observation schemes
  # The age when they were observed (either dead or alive)
  
  # PANCHO 04.11.2020: Change the way compute Age to avoid Age < minAge
  # IDL2$Age <-  time_length(difftime(IDL2$DoE, IDL2$DoB), "years")
  IDL2$Age <- IDL2$AgeYears + IDL2$DaysSinceBd / 365.25
  IDL2$Age.at.b <- time_length(difftime(IDL2$b, IDL2$DoB), "years")
  # IDL2 <- subset(IDL2, Age >= minAge)
  
  # Difference between the start of the observation window in each country
  # minus the date when they turned minAge years
  # There are some individuals that turned minAge years before b so they have negative "Lag"
  IDL2$Lag <-  time_length(difftime(IDL2$DoH, IDL2$b), "years")
  
  # How much time do they stay in the observation window
  # I took the absolute value because of the Lag is bigger than the Complete Stay
  # They have been followed for very short times
  # And they entered to the DB very late
  IDL2$Stay <- abs(time_length(difftime(IDL2$DoE, IDL2$b), "years") -pmax(IDL2$Lag,0))
  
  
  # 1.2 Put the data in the Surv format -------------------------------------
  IDL2$entry <- ifelse(IDL2$Lag< 0, IDL2$Age.at.b, minAge)
  
  # Time they spent in the observation window
  IDL2$exit <- IDL2$entry + IDL2$Stay
  
  # If they experience the event (dead or alive)
  IDL2$delta <- ifelse( IDL2$Status == "Dead", 1,0)
  rIDL <- subset(IDL2, Yb >= 1890)
  
  # The actual database that we are going to use for the calculations
  dati <- IDL2[, c('id', 'Country', 'Sex', 'Yb', 'entry', 'exit', 'delta')]
  # dati <- IDL2[,c(2,1,5,8,22,23,24)]
  # Take just those that have some contribution to the likelihood
  # This is to avoid warnings in the Surv Split
  dati <- subset(dati,subset=dati$entry<dati$exit)
  dati <- dati[order(dati$exit),]
  dati$time <- dati$exit - dati$entry
  
  # 2.1 Split the data into small intervals ---------------------------------
  
  # Total interval of ages
  # max(dati$exit)-min(dati$entry)
  # Bins have to be split in relation to the ages
  maxInter <- round(max(dati$exit)-min(dati$entry),1)
  # This is the age interval
  tauj <- seq(nx, maxInter, nx) + minAge 
  M <- length(tauj)
  
  # SURVIVAL OBJECT
  dati1 <- survSplit(Surv(entry, exit, delta)~.,data=dati, cut=tauj,
                     start="entry",  end="exit", 
                     event="delta", 
                     episode="interval")
  
  dati1$y.new <- dati1$exit-dati1$entry
  dati1$interval <- as.factor(dati1$interval)
  
  # Split country by country - All cohorts
  USA.t <- subset(dati1, Country %in% c("USA"))
  USA.f <- subset(USA.t, Sex == "F")
  if ('M' %in% sex) {
    USA.m <- subset(USA.t, Sex == "M")
  }
  USA.r <- subset(USA.t, Yb >= 1890)
  
  # ESTIMATED HAZARD
  mx.USA.f <- mx.calc(entry = USA.f$entry, exit = USA.f$exit, delta = USA.f$delta, interval = USA.f$interval,s.age = minAge, nx = nx)
  mx.USA.r <- mx.calc(entry = USA.r$entry, exit = USA.r$exit, delta = USA.r$delta, interval = USA.r$interval,s.age = minAge, nx = nx)
  if ('M' %in% sex) {
    mx.USA.t <- mx.calc(entry = USA.t$entry, exit = USA.t$exit, delta = USA.t$delta, interval = USA.t$interval,s.age = minAge, nx = nx)
    mx.USA.m <- mx.calc(entry = USA.m$entry, exit = USA.m$exit, delta = USA.m$delta, interval = USA.m$interval,s.age = minAge, nx = nx)  
  }
  
  # Labels
  mx.USA.f$Sex.Country <- "Females.USA"
  mx.USA.r$Sex.Country <- "rTotal.USA"
  if ('M' %in% sex) {
    mx.USA.t$Sex.Country <- "Total.USA"
    mx.USA.m$Sex.Country <- "Males.USA"  
  }

  # Store results
  mx.f[[i]] <- mx.USA.f
  if ('M' %in% sex) {
    mx.t[[i]] <- mx.USA.t
    mx.m[[i]] <- mx.USA.m
  }
    
  if (i == 1) {
    # OUTPUT TABLES with Qx values from each simulation
    qxMat.f <- matrix(0, nrow = 0, ncol = length(nages),
                      dimnames = list(NULL, nages))
    if ('M' %in% IDL2$Sex) {
      qxMat.m <- matrix(0, nrow = 0, ncol = length(nages),
                        dimnames = list(NULL, nages))
      qxMat.t <- matrix(0, nrow = 0, ncol = length(nages),
                        dimnames = list(NULL, nages))
    }
  } 
  
  # FEMALES
  out <- simUSA(mx = mx.USA.f$mx, n = nrow(IDL2[IDL2$Sex == 'F']), nx = nx, 
                nage = nages, nsimCI = nsim2)
  qxMat.f <- rbind(qxMat.f, out)
  if ('M' %in% IDL2$Sex) {
    # MALES
    out <- simUSA(mx = mx.USA.m$mx, n = nrow(IDL2[IDL2$Sex == 'M']),
                  nx = nx, nage = nages, nsimCI = nsim2)
    qxMat.m <- rbind(qxMat.m, out)
    # TOTAL
    out <- simUSA(mx = mx.USA.t$mx, n = nrow(IDL2), 
                  nx = nx, nage = nages, nsimCI = nsim2)
    qxMat.t <- rbind(qxMat.t, out)
  }
    
  # Remove data
  rm(IDL2)
  
}
End <- Sys.time()
print(End - Start)


#------------------------------------------------#
# RATES, CONFIDENCE INTERVALS AND OTHER MEASURES #
#------------------------------------------------#

# FEMALES
ciFem <- ciUSA(qxMat = qxMat.f, Age = nages,nx = 0.25)
# MALES
ciMen <- ciUSA(qxMat = qxMat.m, Age = nages,nx = 0.25)
# TOTAL
ciTot <- ciUSA(qxMat = qxMat.t, Age = nages,nx = 0.25)

SAVE <- F
if (SAVE) {
  if ('M' %in% sex) {
    save(qxMat.f, qxMat.m, qxMat.t, mx.f, mx.m, mx.t, ciFem, ciMen, ciTot,
         file = paste0(format(Sys.Date(), "%Y%m%d"), '-TestUSA-CI.RData'))  
  } else {
    save(qxMat.f, mx.f, ciFem,
         file = paste0(format(Sys.Date(), "%Y%m%d"), '-TestUSA-CI.RData'))
  }
}


ggplot(subset(ciFem),
       aes(Age, mxMean))+
  geom_point(size = 0.2)+
  geom_line(size = 0.4)+
  scale_fill_manual(values = c("black"))+ ##5230fc #ff3729
  scale_color_manual(values = c("black"))+
  #facet_wrap(~Sex, ncol = 4, scales = "free_x")+
  geom_ribbon(aes(ymin = mx.low, ymax = mx.up), alpha = 0.1)+
  scale_y_continuous(expand = c(0,0),breaks = seq(0,1.8,by= 0.4))+
  scale_x_continuous(expand = c(0,0),breaks = seq(105,120,by = 1),
                     labels = c("105","'6","'7","'8","'9","110","'11","'12","'13","'14","115","","","","","120"))+
  coord_cartesian(ylim = c(0,1.6), xlim= c(105, 115))+
  theme_bw()+
  theme(strip.background = element_rect(fill="none"))+
  theme(strip.background = element_blank(),
        panel.border = element_blank(),
        axis.text.x  = element_text(hjust = 0.5),
        panel.spacing = unit(2.5, "lines"),
        panel.grid.major = element_line(size= .2),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 18,font_rc),
        axis.title.y = element_text(vjust = 2, size = 22),
        axis.title.x = element_text(vjust = -2, size = 22),
        axis.ticks = element_line(size = 0.3),
        aspect.ratio = 0.5,
        strip.text.x = element_text(size = 22, colour = "black"),
        legend.position = "none",
        plot.background = element_rect(fill = NA))+
  ylab("Risk of dying")
ggsave("USA_1fem.pdf", width = 12, height = 6, device = cairo_pdf)


save(ciFem,ciMen, file = "ciUSA.RData")
