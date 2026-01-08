#####Stinknet Analysis
#####8.9.2024

library(tidyverse)
library(car)
library(lmerTest)
library(emmeans)
library(readxl)
library(sjPlot)
library(jtools)
library(ggeffects)
library(ggplot2)
library(multcomp)
library(MuMIn)
library(gtsummary)
library(viridis)
install.packages("remotes")
remotes::install_github("ameztegui/Plasticity")
library(Plasticity)

height<-read_excel(file.choose())#read in greenhouse heights
height[height == 0] <- NA #Change 0's to NA

height$VPD<-height$population
height$VPD[height$VPD==1]<-1.22
height$VPD[height$VPD==2]<-2.06
height$VPD[height$VPD==3]<-3.62
height$VPD[height$VPD==4]<-3.81


height$population<-as.factor(height$population)
height$treatment<-as.factor(height$treatment)
height$RGR=log(height$height_cm_42/height$height_cm_0)
levels(height$treatment) <- c("Daily", "Pulsed")
height<-as.data.frame(height)
height<-within(height, population<-relevel(population, ref=4))

height.daily<-subset(height,treatment == "Daily")%>% droplevels()
height.daily<-within(height.daily, population<-relevel(population, ref=4))

height.pulse<-subset(height,treatment == "Pulsed")%>% droplevels()
height.pulse<-within(height.pulse, population<-relevel(population, ref=4))


bmass<-read_excel(file.choose())
bmass[bmass == 0] <- NA #Change 0's to NA
bmass<-spread(bmass, key=`plant part`, value=`total mass g`)#create new columsn for each tisuse type
bmass$root.shoot<-(bmass$root/bmass$shoot)#create new column with root:shoot ratio
bmass$total<-bmass$root+bmass$shoot

bmass$VPD<-bmass$population
bmass$VPD[bmass$VPD==1]<-1.22
bmass$VPD[bmass$VPD==2]<-2.06
bmass$VPD[bmass$VPD==3]<-3.62
bmass$VPD[bmass$VPD==4]<-3.81

bmass$population<-as.factor(bmass$population)
bmass$treatment<-as.factor(bmass$treatment)
levels(bmass$treatment) <- c("Daily", "Pulsed")
bmass<-as.data.frame(bmass)

#create a mortality column in the height file
height$mort<-height$height_cm_42
height$mort[is.na(height$mort)] <- 0
height$mort[height$mort > 0] <- 1 
#height<-within(height, population<-relevel(population, ref=4))

mort.by.treat<-aggregate(mort ~ treatment + population, data =height, FUN = sum)
mort.table<-flextable(data=mort.by.treat %>% rownames_to_column("Mortality"))
mort.table <- set_header_labels(mort.table, treatment = "Treatment", population="Population", mort="Surviving")
mort.table


bmass.daily<-filter(bmass,treatment == "Daily")%>% droplevels()
bmass.daily<-within(bmass.daily, population<-relevel(population, ref=4))

bmass.pulse<-subset(bmass,treatment == "Pulsed")%>% droplevels()
bmass.pulse<-within(bmass.pulse, population<-relevel(population, ref=4))


#inspect data
hist(height.daily$height_cm_42)#fairly normal, no adjustment needed
hist(bmass.daily$root.shoot)
hist(log(bmass.daily$root.shoot))#this is much less skewed
hist(bmass.daily$total)

############Mortality
mortality.out<-glm(mort~population*treatment, family="binomial"(link="logit"), data=height)
summary(mortality.out)###No significant effect of treatments..combine to look at pops

mortality.out.all<-glm(mort~population, family="binomial"(link="logit"), data=height)
summary(mortality.out.all)

tbl_regression(mortality.out)

#Need to create mort graph


###############Height

height.out<-aov(height_cm_42~population*treatment, data=height)
summary(height.out)#population and treatment main effects, no interaction
emmeans(height.out, pairwise~population)
emmeans(height.out, pairwise~treatment)


detach(package:plyr)

height.means<-height %>%
  group_by(population) %>%
  summarize(mean.height = mean(height_cm_42,na.rm=T), SD.height = sd(height_cm_42,na.rm=T),
            mean.RGR = mean(RGR,na.rm=T), SD.RGR = sd(RGR,na.rm=T))
control <- c("4", "1", "2", "3")

height.bar<-ggplot(height.means) +
  geom_bar( aes(x=population, y=mean.height, fill=population),stat="identity", alpha=0.7) +
  geom_errorbar( aes(x=population, ymin=mean.height-SD.height, ymax=mean.height+SD.height), width=0.4, colour="black", alpha=0.9)+
  labs(x = "Population", y = "Height (cm)")+
  scale_x_discrete(limits=control)+
  scale_fill_manual(values=c("chocolate", "chocolate", "chocolate","gray74"))+
  theme_classic()+
  theme(legend.position="none")
height.bar

height.box<-ggplot(height, aes(x=treatment, y=height_cm_42, color=treatment, fill=treatment))+
  geom_boxplot(alpha=0.3)+
  scale_color_manual(name = "Treatment",
                     labels=c("Daily", "Pulsed"),
                     values = c("darkblue","skyblue3")) +
  scale_fill_manual(name = "",
                    values = c("darkblue","skyblue3")) + 
  labs(x = "Treatment", y = "Height (cm)")+
  theme(axis.text=element_text(colour = 'black',size=10)) +
  theme(legend.text = element_text(size = 10)) +
  theme_classic()+
  theme(legend.position="none")
height.box




########Total Biomass

total.out<-aov(total~population*treatment, data=bmass)
summary(total.out)#treatment effect only
emmeans(total.out, pairwise~treatment)

total.box<-ggplot(bmass, aes(x=treatment, y=total, color=treatment, fill=treatment))+
  geom_boxplot(alpha=0.3)+
  scale_color_manual(name = "Treatment",
                     labels=c("Daily", "Pulsed"),
                     values = c("darkblue","skyblue3")) +
  scale_fill_manual(name = "",
                     values = c("darkblue","skyblue3")) + 
  labs(x = "Treatment", y = "Total Biomass (g)")+
  theme(axis.text=element_text(colour = 'black',size=10)) +
  theme(legend.text = element_text(size = 10)) +
  theme_classic()+
  theme(legend.position="none")
  
total.box

####Root:Shoot

bmass$l.root.shoot<-log(bmass$root.shoot)
bmass.daily$l.root.shoot<-log(bmass.daily$root.shoot)
bmass.pulse$l.root.shoot<-log(bmass.pulse$root.shoot)

rs.out<-aov(l.root.shoot~population*treatment, data=bmass)
summary(rs.out)#population differences only
emmeans(rs.out, pairwise~population)

rs.means<-bmass %>%
  group_by(population) %>%
  summarize(mean.rs = mean(root.shoot,na.rm=T), SD.rs = sd(root.shoot,na.rm=T))


rs.bar<-ggplot(rs.means) +
  geom_bar( aes(x=population, y=mean.rs, fill=population),stat="identity", alpha=0.7) +
  geom_errorbar( aes(x=population, ymin=mean.rs-SD.rs, ymax=mean.rs+SD.rs), width=0.4, colour="black", alpha=0.9)+
  labs(x = "Population", y = "Root to Shoot Ratio")+
  scale_x_discrete(limits=control)+
  scale_fill_manual(values=c("gray74", "chocolate", "chocolate","chocolate"))+
  theme_classic()+
  theme(legend.position="none")
rs.bar

#Growth Rate Analysis
growth.rate<-Anova(RGR~population*treatment, data=height)
summary(growth.rate)#treatment only
emmeans(growth.rate, pairwise~treatment)

RGR.box<-ggplot(height, aes(x=treatment, y=RGR, color=treatment, fill=treatment))+
  geom_boxplot(alpha=0.3)+
  scale_color_manual(name = "Treatment",
                     labels=c("Daily", "Pulsed"),
                     values = c("darkblue","skyblue3")) +
  scale_fill_manual(name = "",
                    values = c("darkblue","skyblue3")) + 
  labs(x = "Treatment", y = "Relative Growth Rate")+
  theme(axis.text=element_text(colour = 'black',size=10)) +
  theme(legend.text = element_text(size = 10)) +
  theme_classic()+
  theme(legend.position="none")
RGR.box

###########Make Results Tables
library(flextable)

height.out<-anova(height.out)
height.table<-flextable(data=height.out %>% rownames_to_column("Height (cm)"))
height.table

biomass.out<-anova(total.out)
biomass.table<-flextable(data=biomass.out %>% rownames_to_column("Biomass (g)"))
biomass.table

growth.out<-anova(total.rate)
rgr.table<-flextable(data=growth.out %>% rownames_to_column("Relative Growth Rate"))
rgr.table

rs.out<-anova(rs.out)
rs.table<-flextable(data=rs.out %>% rownames_to_column("Root:Shoot"))
rs.table

