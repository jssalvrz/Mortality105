######################################################################################
# This code is to estimate Confidence Intervals and plot the results
# Updated by: Jesus-Adrian Alvarez
# Date: 12-12-2020
# Interval: 0.25 years
######################################################################################

# Clear workspace
rm(list = ls())

# Libraries
library(lubridate)
library(tidyverse)
library(data.table)
library(survival)
library(changepoint)
library(lmvar)
library(hrbrthemes)
library(lemon)

# Load splitted data and functions
#load("datiSplit2.RData")
source("Super100_FUN.R")
setwd("C:/Users/jmartinez/OneDrive - Syddansk Universitet/Super100/R")
load("datiSplitFinal.RData") # Dataset including estimates for Mx

load("ciUSA.RData")


# Calculate CI for Mx using simulations -----------------------------------

Mx <- data.table(subset(Mx, expo>0))

# Create separate Mx dataframes and simulate CI with specific sample size
# We have to do it country by country :(


Mx.fraM <- data.table(subset(Mx, Sex.Country == "Males.FRA"))
mxSIM.fraM <- Mx.fraM[,CImx(Age = Age, mx = mx, nx=0.25, n = 863), by = list(Sex.Country)]
mxSIM.fraM$sex <- "Male"
mxSIM.fraM$country <- "France"
Mx.fraF <- data.table(subset(Mx, Sex.Country == "Females.FRA"))
mxSIM.fraF <- Mx.fraF[,CImx(Age = Age, mx = mx, nx=0.25, n = 8990), by = list(Sex.Country)]
mxSIM.fraF$sex <- "Female"
mxSIM.fraF$country <- "France"

Mx.gerM <- data.table(subset(Mx, Sex.Country == "Males.GER"))
mxSIM.gerM <- Mx.gerM[,CImx(Age = Age, mx = mx,nx=0.25, n = 111), by = list(Sex.Country)]
mxSIM.gerM$sex <- "Male"
mxSIM.gerM$country <- "Germany"
Mx.gerF <- data.table(subset(Mx, Sex.Country == "Females.GER"))
mxSIM.gerF <- Mx.gerF[,CImx(Age = Age, mx = mx, nx=0.25, n = 859), by = list(Sex.Country)]
mxSIM.gerF$sex <- "Female"
mxSIM.gerF$country <- "Germany"

# Mx.grbM <- data.table(subset(Mx, Sex.Country == "Males.GRB"))
# mxSIM.grbM <- Mx.grbM[,CImx(Age = Age, mx = mx,nx=0.5, n = 317), by = list(Sex.Country)]
# mxSIM.grbM$sex <- "Male"
# mxSIM.grbM$country <- "Great Britain"
# Mx.grbF <- data.table(subset(Mx, Sex.Country == "Females.GRB"))
# mxSIM.grbF <- Mx.grbF[,CImx(Age = Age, mx = mx, nx=0.5, n = 867), by = list(Sex.Country)]
# mxSIM.grbF$sex <- "Female"
# mxSIM.grbF$country <- "Great Britain"#867

Mx.belM <- data.table(subset(Mx, Sex.Country == "Males.BEL"))
mxSIM.belM <- Mx.belM[,CImx(Age = Age, mx = mx, nx=0.25, n = 82), by = list(Sex.Country)]
mxSIM.belM$sex <- "Male"
mxSIM.belM$country <- "Belgium"
Mx.belF <- data.table(subset(Mx, Sex.Country == "Females.BEL"))
mxSIM.belF <- Mx.belF[,CImx(Age = Age, mx = mx, nx=0.25, n = 777), by = list(Sex.Country)]
mxSIM.belF$sex <- "Female"
mxSIM.belF$country <- "Belgium"

Mx.dnkM <- data.table(subset(Mx, Sex.Country == "Males.DNK"))
mxSIM.dnkM <- Mx.dnkM[,CImx(Age = Age, mx = mx, nx=0.25, n = 67), by = list(Sex.Country)]
mxSIM.dnkM$sex <- "Male"
mxSIM.dnkM$country <- "Denmark"
Mx.dnkF <- data.table(subset(Mx, Sex.Country == "Females.DNK"))
mxSIM.dnkF <- Mx.dnkF[,CImx(Age = Age, mx = mx, nx=0.25, n = 417), by = list(Sex.Country)]
mxSIM.dnkF$sex <- "Female"
mxSIM.dnkF$country <- "Denmark"

Mx.autM <- data.table(subset(Mx, Sex.Country == "Males.AUT"))
mxSIM.autM <- Mx.autM[,CImx(Age = Age, mx = mx, nx=0.25, n = 44), by = list(Sex.Country)]
mxSIM.autM$sex <- "Male"
mxSIM.autM$country <- "Austria"
Mx.autF <- data.table(subset(Mx, Sex.Country == "Females.AUT"))
mxSIM.autF <- Mx.autF[,CImx(Age = Age, mx = mx, nx=0.25, n = 267), by = list(Sex.Country)]
mxSIM.autF$sex <- "Female"
mxSIM.autF$country <- "Austria"

# Mx.cheM <- data.table(subset(Mx, Sex.Country == "Males.CHE"))
# mxSIM.cheM <- Mx.cheM[,CImx(Age = Age, mx = mx, nx=0.5, n = 39), by = list(Sex.Country)]
# mxSIM.cheM$sex <- "Male"
# mxSIM.cheM$country <- "Switzerland"
# Mx.cheF <- data.table(subset(Mx, Sex.Country == "Females.CHE"))
# mxSIM.cheF <- Mx.cheF[,CImx(Age = Age, mx = mx, nx=0.5, n = 201), by = list(Sex.Country)]
# mxSIM.cheF$sex <- "Female"
# mxSIM.cheF$country <- "Switzerland"

Mx.norM <- data.table(subset(Mx, Sex.Country == "Males.NOR"))
mxSIM.norM <- Mx.norM[,CImx(Age = Age, mx = mx, nx=0.25, n = 41), by = list(Sex.Country)]
mxSIM.norM$sex <- "Male"
mxSIM.norM$country <- "Norway"
Mx.norF <- data.table(subset(Mx, Sex.Country == "Females.NOR"))
mxSIM.norF <- Mx.norF[,CImx(Age = Age, mx = mx, nx=0.25,n = 187), by = list(Sex.Country)]
mxSIM.norF$sex <- "Female"
mxSIM.norF$country <- "Norway"


Mx.canM <- data.table(subset(Mx, Sex.Country == "Males.CAN"))
mxSIM.canM <- Mx.canM[,CImx(Age = Age, mx = mx, nx=0.25, n = 46), by = list(Sex.Country)]
mxSIM.canM$sex <- "Male"
mxSIM.canM$country <- "Quebec"
Mx.canF <- data.table(subset(Mx, Sex.Country == "Females.CAN"))
mxSIM.canF <- Mx.canF[,CImx(Age = Age, mx = mx, nx=0.25, n = 287), by = list(Sex.Country)]
mxSIM.canF$sex <- "Female"
mxSIM.canF$country <- "Quebec"


# Mx.espM <- data.table(subset(Mx, Sex.Country == "Males.ESP"))
# mxSIM.espM <- Mx.espM[,CImx(Age = Age, mx = mx, n = 9), by = list(Sex.Country)]
# mxSIM.espM$sex <- "Male"
# mxSIM.espM$country <- "Spain"
# Mx.espF <- data.table(subset(Mx, Sex.Country == "Females.ESP"))
# mxSIM.espF <- Mx.espF[,CImx(Age = Age, mx = mx, n = 51), by = list(Sex.Country)]
# mxSIM.espF$sex <- "Female"
# mxSIM.espF$country <- "Spain"
# 
# Mx.itaM <- data.table(subset(Mx, Sex.Country == "Males.ITA"))
# mxSIM.itaM <- Mx.itaM[,CImx(Age = Age, mx = mx, n = 9), by = list(Sex.Country)]
# mxSIM.itaM$sex <- "Male"
# mxSIM.itaM$country <- "Italy"
# Mx.itaF <- data.table(subset(Mx, Sex.Country == "Females.ITA"))
# mxSIM.itaF <- Mx.itaF[,CImx(Age = Age, mx = mx, n = 51), by = list(Sex.Country)]
# mxSIM.itaF$sex <- "Female"
# mxSIM.itaF$country <- "Italy"
# 
# Mx.sweM <- data.table(subset(Mx, Sex.Country == "Males.SWE"))
# mxSIM.sweM <- Mx.sweM[,CImx(Age = Age, mx = mx, n = 1), by = list(Sex.Country)]
# mxSIM.sweM$sex <- "Male"
# mxSIM.sweM$country <- "Sweden"
# Mx.sweF <- data.table(subset(Mx, Sex.Country == "Females.SWE"))
# mxSIM.sweF <- Mx.sweF[,CImx(Age = Age, mx = mx, n = 11), by = list(Sex.Country)]
# mxSIM.sweF$sex <- "Female"
# mxSIM.sweF$country <- "Sweden"


# Mx.eurM <- data.table(subset(Mx, Sex.Country == "Males.EUR"))
# mxSIM.eurM <- Mx.eurM[,CImx(Age = Age, mx = mx, n = 110), by = list(Sex.Country)]
# mxSIM.eurM$sex <- "Male"
# mxSIM.eurM$country <- "Nordic Countries"
# Mx.eurF <- data.table(subset(Mx, Sex.Country == "Females.EUR"))
# mxSIM.eurF <- Mx.eurF[,CImx(Age = Age, mx = mx, n = 620), by = list(Sex.Country)]
# mxSIM.eurF$sex <- "Female"
# mxSIM.eurF$country <- "Nordic Countries"
# 
# Mx.eufM <- data.table(subset(Mx, Sex.Country == "Males.EUF"))
# mxSIM.eufM <- Mx.eufM[,CImx(Age = Age, mx = mx, n = 100), by = list(Sex.Country)]
# mxSIM.eufM$sex <- "Male"
# mxSIM.eufM$country <- "Southern Europe"
# Mx.eufF <- data.table(subset(Mx, Sex.Country == "Females.EUF"))
# mxSIM.eufF <- Mx.eufF[,CImx(Age = Age, mx = mx, n = 561), by = list(Sex.Country)]
# mxSIM.eufF$sex <- "Female"
# mxSIM.eufF$country <- "Southern Europe"


# Mx.finM <- data.table(subset(Mx, Sex.Country == "Males.FIN"))
# mxSIM.finM <- Mx.finM[,CImx(Age = Age, mx = mx, n = 82), by = list(Sex.Country)]
# mxSIM.finM$sex <- "Male"
# mxSIM.finM$country <- "Finland"
# Mx.finF <- data.table(subset(Mx, Sex.Country == "Females.FIN"))
# mxSIM.finF <- Mx.finF[,CImx(Age = Age, mx = mx, n = 777), by = list(Sex.Country)]
# mxSIM.finF$sex <- "Female"
# mxSIM.finF$country <- "Finland"



#Merge all the results to plot them together
mxSIM.red <- rbind(mxSIM.fraF, mxSIM.fraM, 
                   #mxSIM.grbF, mxSIM.grbM, 
                   mxSIM.gerF, mxSIM.gerM,
                   mxSIM.belF, mxSIM.belM,
                   mxSIM.dnkF, mxSIM.dnkM,
                   mxSIM.autF, mxSIM.autM,
                   #mxSIM.cheF, mxSIM.cheM,
                   mxSIM.canF, mxSIM.canM,
                   mxSIM.norF, mxSIM.norM)
                  # mxSIM.espF, mxSIM.espM,
                  # mxSIM.itaF, mxSIM.itaM,
                  # mxSIM.sweF, mxSIM.sweM,
                  # mxSIM.eurF, mxSIM.eurM,
                  # mxSIM.eufF, mxSIM.eufM)
                   #mxSIM.finF, mxSIM.finM)



mxest <- mxSIM.red[,c("country", "sex", "Age", "mx", "mx.low","mx.up")]

ciFem$country <- "United States"
ciFem$sex <- "Female"
mxUSA <- ciFem[,c("country", "sex", "Age", "mxMean", "mx.low", "mx.up")]
names(mxUSA)[4] <- "mx"

all <- subset(mxest, country == "France" & sex == "Female")
all$country <- "All"

mxest <- data.frame(rbind( mxest, mxUSA))





mxest$country<- factor(mxest$country, levels =(c("All","France", "Germany","Belgium", "United States",
                                                "Denmark", "Quebec","Austria",
                                                 "Norway","Great Britain" )))

# proportions according to survivorship function

prop <- read_csv("prop.csv")

mxest <- merge(mxest,prop)

ggplot(subset(mxest,  sex == "Female"),
       aes(Age+0.125, mx, alpha = p))+
  geom_hline(yintercept = 0.8, colour = "grey", linetype = "dashed" ,size = 0.3, alpha = 1)+
  geom_hline(yintercept = 0.6, colour = "grey",  linetype = "dashed",size = 0.3, alpha = 1)+
  geom_ribbon(aes(ymin = mx.low, ymax = mx.up), fill = "grey50", alpha = 0.1)+
  geom_line(size = 0.1, color = "black")+
  geom_point(size = 1.5, stroke = .00001, color = "black")+
  facet_rep_wrap(~country,ncol = 3, repeat.tick.labels = T)+

  scale_y_continuous(expand = c(0,0),breaks = seq(0,1.2,by= 0.2))+
   scale_x_continuous(expand = c(0,0),breaks = seq(105,120,by = 1),
                      labels = c("105","'6","'7","'8","'9","'10","'11","'12","'13","","115","","","","","120"))+
  coord_cartesian(ylim = c(0,1.2), xlim= c(105, 113))+
  theme_bw()+
  theme(strip.background = element_rect(fill="none"))+
  theme(strip.background = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(hjust = -2, size = 8),
        axis.text.x = element_text(vjust = -1, size = 8),
        panel.spacing = unit(2.5, "lines"),
        panel.grid.major = element_blank(),#element_line(size= .1, color = "grey90"),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 20,font_rc),
        axis.title.y = element_text(vjust = 2, size = 12),
        axis.title.x = element_text(vjust = -1, size = 12),
        axis.ticks = element_line(size = 0.3),
        aspect.ratio = 0.5,
        strip.text.x = element_text(size = 16, colour = "black"),
        legend.position = "none",
        plot.background = element_rect(fill = NA))+
  ylab("Risk of dying")+
  xlab("Age")

ggsave("Fig2__25.pdf", width = 8, height = 6, device = cairo_pdf)

