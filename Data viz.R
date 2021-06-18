library(tidyverse)
library(hrbrthemes)
library(skimr)
library(waffle)
library(lubridate)
library(lemon)
setwd("C:/Users/jmartinez/OneDrive - Syddansk Universitet/Super100/R")
IDL <-   data.frame(read_csv("IDLOfficialFinal.csv"))





# Date of birth
IDL$DoB <- as.Date(paste(IDL$Db,IDL$Mb, IDL$Yb, sep = "-"), format = "%d-%m-%Y") 
# Date of entry to the DB
IDL$DoE <- as.Date(paste(IDL$De,IDL$Me, IDL$Ye, sep = "-"), format = "%d-%m-%Y") 

IDL$Age <-  time_length(difftime(IDL$DoE, IDL$DoB), "years")


 IDL$Cohort5 <- ifelse(IDL$Yb>=1845 & IDL$Yb <=1849, "1845-49",IDL$Yb )
 
 IDL$Cohort5 <- ifelse(IDL$Yb>=1850 & IDL$Yb <=1854, "1850-54",IDL$Cohort5 )
 IDL$Cohort5 <- ifelse(IDL$Yb>=1855 & IDL$Yb <=1859, "1855-59",IDL$Cohort5 )
 
 
 IDL$Cohort5 <- ifelse(IDL$Yb>=1860 & IDL$Yb <=1864, "1860-64",IDL$Cohort5 )
 IDL$Cohort5 <- ifelse(IDL$Yb>=1865 & IDL$Yb <=1869, "1865-69",IDL$Cohort5 )
 
 IDL$Cohort5 <- ifelse(IDL$Yb>=1870 & IDL$Yb <=1874, "1870-74",IDL$Cohort5 )
 IDL$Cohort5 <- ifelse(IDL$Yb>=1875 & IDL$Yb <=1879, "1875-79",IDL$Cohort5 )
 
 IDL$Cohort5 <- ifelse(IDL$Yb>=1880 & IDL$Yb <=1884, "1880-84",IDL$Cohort5 )
 IDL$Cohort5 <- ifelse(IDL$Yb>=1885 & IDL$Yb <=1889, "1885-89",IDL$Cohort5 )
 
 IDL$Cohort5 <- ifelse(IDL$Yb>=1890 & IDL$Yb <=1894, "1890-94",IDL$Cohort5 )
 IDL$Cohort5 <- ifelse(IDL$Yb>=1895 & IDL$Yb <=1899, "1895-99",IDL$Cohort5 )
 
 IDL$Cohort5 <- ifelse(IDL$Yb>=1900 & IDL$Yb <=1904, "1900-04",IDL$Cohort5 )
 IDL$Cohort5 <- ifelse(IDL$Yb>=1905 & IDL$Yb <=1909, "1905-09",IDL$Cohort5 )
 
 IDL$Cohort5 <- ifelse(IDL$Yb>=1910 & IDL$Yb <=1912, "1910-12",IDL$Cohort5 )
 
 
 IDL$Country<- factor(IDL$Country, levels = rev(c("FRA", "DEU",
                                                  "BEL", "DNK", "CAN","AUT",
                                                  "CHE", "NOR", 
                                                  "ESP",  "SWE", "FIN",
                                                  "USA")))
 
  # IDL$Country<- factor(IDL$Country, levels = rev(c("AUT", "BEL", "CAN",
  #                                                  "CHE", "DEU", "DNK",
  #                                                  "ENW", "ESP", "FIN",
  #                                                  "FRA", "ITA", "JPN",
  #                                                  "NOR", "SWE", "USA")))

IDL$Sex[IDL$Sex == "F"] <- "Females"
IDL$Sex[IDL$Sex == "M"] <- "Males"

#IDL$Status[IDL$Country == "USA"]<- "Alive"

IDL$id <- paste(IDL$Sex,IDL$Status, sep="")

# 
 Cens <-data.frame(xtabs(~Country+Sex, data=IDL))
 Cens$id <- paste(Cens$Status,Cens$Sex, sep="")
 
 
 # Country <- Country %>% gather(Sex, Individuals,-Country, na.rm = TRUE)
 
 
                                                                                                                                    
# 
 # Cens <- Cens %>% spread(Status, -Country)
 # Cens$Total <- Cens$Alive+Cens$Dead
 # Cens$pAlive <- Cens$Alive / Cens$Total
 # Cens$pDead  <- Cens$Dead / Cens$Total
 
 
 
 
##################################################
 
 # Load data for the USA
 load("IDLUSA.RData")
 IDL2$Status <- "Alive"
 IDL2$Sex <- ifelse(IDL2$Sex == "M", "Males", "Females")
 IDL2$id <-  paste(IDL2$Sex,IDL2$Status, sep="")
 IDL2 <- IDL2[,c("id","Sex", "Country", "Status", "Age")]
 
 IDL3 <- IDL[,c("id","Sex", "Country", "Status", "Age")]
 
IDL3 <- rbind(IDL3, IDL2)
rename(IDL3$Country, replace=c( "NOR"= "Norway", "AUT"="Austria", "CAN"="Quebec",
                        "DNK"="Denmark",
                        "USA" = "United States", "BEL"="Belgium",
                        "DEU"= "Germany", "FRA"="France"))

 IDL3$Country<- factor(IDL3$Country, levels =(rev(c( "NOR", "AUT", "CAN", "DNK",
                                                    "USA", "BEL", "DEU", "FRA"))),
                   labels = rev(c( "Norway", "Austria","Quebec",
                          "Denmark","United States", "Belgium",
                           "Germany", "France")))



#ggsave("CohortDistr.pdf", width = 6, height = 4, device = cairo_pdf)



ggplot(IDL3)+
  geom_histogram(aes(Age-0.05, fill = id), position = "stack",binwidth = 0.1)+
  scale_fill_manual(values = c("#ff8f87","#ff3729","#b3a3ff", "#5230fc"))+
  scale_x_continuous(expand = c(0,0),breaks = seq(105,122,by = 1),
                     labels = c("105","","","","","110","","","","","115",
                                "","","","","120", "", "'22"))+
  scale_y_continuous(expand = c(0,0),
                     breaks = seq(0,1000,by= 200))+
                     #labels = c(0,"1","2","3","4","5","6","7","8","9","10"))+
  coord_cartesian(ylim = c(0,1000), xlim= c(105, 123))+
  theme_bw()+
  theme(strip.background = element_rect(fill="none"))+
  theme(strip.background = element_blank(),
        panel.border = element_blank(),
        axis.text.x  = element_text(hjust = 0.5),
        panel.spacing = unit(2, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 13,font_rc),
        axis.title.y = element_text(vjust = 2, size = 16),
        axis.title.x = element_text(vjust = -1, size = 16),
        aspect.ratio =1,
        strip.text.x = element_text(size = 17, colour = "black"),
        legend.position = "none",
        plot.background = element_rect(fill = NA))+
  ylab("Observed individuals")+
  xlab("Age")

ggsave("AgesDistr.pdf", width = 6, height = 4, device = cairo_pdf)



ggplot(subset(IDL3))+
  geom_histogram(aes(Age-0.05, fill = id), position = "stack",binwidth = 0.1)+
  scale_fill_manual(values = c("#ff8f87","#ff3729","#b3a3ff", "#5230fc"))+
  scale_x_continuous(expand = c(0,0),breaks = seq(105,120,by = 5),
                     labels = c("105","'10", "'15", "'20"))+
  scale_y_continuous(expand = c(0,0),
                     breaks = seq(0,1000,by= 200))+
  facet_rep_wrap(~Country, ncol = 3, scales = "fixed",  repeat.tick.labels = T)+
  #labels = c(0,"1","2","3","4","5","6","7","8","9","10"))+
  coord_cartesian(ylim = c(0,600), )+
  theme_bw()+
  theme(strip.background = element_rect(fill="none"))+
  theme(strip.background = element_blank(),
        panel.border = element_blank(),
        axis.text.x  = element_text(hjust = 0.5),
        axis.ticks = element_line(size = 0.1),
        panel.spacing = unit(1, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 10,font_rc),
        axis.title.y = element_text(vjust = 2, size = 16),
        axis.title.x = element_text(vjust = -1, size = 16),
        aspect.ratio =1,
        strip.text.x = element_text(size = 10, colour = "black"),
        legend.position = "none",
        plot.background = element_rect(fill = NA))+
  ylab("Observed individuals")+
  xlab("Age")

ggsave("CountryDistrallFra.pdf", width = 5, height = 5, device = cairo_pdf)


ggplot(subset(IDL3))+
  geom_histogram(aes(Age-0.05, fill = id), position = "stack",binwidth = 0.1)+
  scale_fill_manual(values = c("#ff8f87","#ff3729","#b3a3ff", "#5230fc"))+
  scale_x_continuous(expand = c(0,0),breaks = seq(105,120,by = 5),
                     labels = c("105","'10", "'15", "'20"))+
  scale_y_continuous(expand = c(0,0),
                     breaks = seq(0,1000,by= 200))+
  facet_rep_wrap(~Country, ncol = 2,  repeat.tick.labels = T)+
  #labels = c(0,"1","2","3","4","5","6","7","8","9","10"))+
  coord_cartesian(ylim = c(0,600) )+
  theme_bw()+
  theme(strip.background = element_rect(fill="none"))+
  theme(strip.background = element_blank(),
        panel.border = element_blank(),
        axis.text.x  = element_text(hjust = 0.5),
        axis.ticks = element_line(size = 0.1),
        panel.spacing = unit(1, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 10,font_rc),
        axis.title.y = element_text(vjust = 2, size = 16),
        axis.title.x = element_text(vjust = -1, size = 16),
        aspect.ratio =1,
        strip.text.x = element_text(size = 10, colour = "black"),
        legend.position = "none",
        plot.background = element_rect(fill = NA))+
  ylab("Observed individuals")+
  xlab("Age")

ggsave("CountryDistrfraAlt.pdf", width = 7, height = 7, device = cairo_pdf)




