#Figures for manuscript and supplemental information titled: 
#Light fuels while nitrogen suppresses symbiotic nitrogen fixation hotspots in Neotropical canopy gap seedlings
#published in New Phytologist 2021


#Packages####
library("tidyverse")
library("MASS")
library("DescTools")
library(ggplot2)
library(multcomp)
library(lsmeans)
library(grid)
library(ggpubr)
library(distreg.vis)


#Calculating N fix variables from Master datasheet ####
dat<-read.csv("LaSelva_MasterData_May2019.csv")

#Changing total N in isotope samples from micrograms to miligrams 
dat$Total.N.e <- dat$total.N.15n*.001
dat$Total.N.c <- dat$total.N.c*.001
#Calculating %N for each sample from the total N and sample weight of enriched samples ####
dat$XN.e<-((dat$Total.N.e)/dat$mg.tin.e)*100

#Calculating %N for each sample from the total N and sample weight of control samples ####
dat$XN.c<-((dat$total.N.c)/dat$mg.tin.c)*100

#Syringe Volume replaced (proportion), Enrichment of normal atmosphere, enrichment of added gas ####
svr<-.2
X15Nair<-.003663
X15Ngas<-.99

##svr.c 
svr.c <- 0

#Calculating the expected enrichment of syringe headspace during incubation.
#this will be very close to the volume replaced times the %15Ngas (how close depends on how close the gas enrichment is to 100%)
syr.enrich<-(X15Ngas*svr)+(X15Nair*(1-svr))
syr.enrich.c <-(X15Ngas*svr.c)+(X15Nair*(1-svr.c))

#Calculating the Percent of N atoms fixed in each sample ####
#calculated is by taking the difference between the enriched X15N and control X15N and (?dividing by expected enrichment)
dat$XNfixdX<-with(dat, (X15n.atX.15n-X15n.atX.c)/syr.enrich)

#Some samples only ended up having values for controls because enriched samples were below detection limit, 
#to avoid negatives I made them zeros
dat$XNfixdX[dat$XNfixdX < 0] = 0 

#Calculating the total amount of N (g) in the incubated sample amount. n15.mass is in grams
dat$tot.N<-with(dat, (n15.mass*XN.e))

#calculating the amount of N fixed in each sample (gN fixed per sample per incubation period) check units?
dat$Nfixd.sample<-with(dat, tot.N*XNfixdX)


#Calculating N fixation to units of (gN fixed per gram of nodule per hour) "SNF efficiency"####
#gN fixed per sample divided by the mass of nodules that were incubated (n15.mass)
dat$Nfixd<-with(dat, (Nfixd.sample/n15.mass)*2)
#Converting N fixation to units of (gN fixed per gram of nodule per year)
time.conv<-17520 #the number of "half hours" in a year
dat$Nfixd.yr<-with(dat, (Nfixd.sample/n15.mass)*time.conv)

#calculating total N fixed per plant ####
#Taking "SNF efficiency" and multiplying by the total nodule mass 
dat$c.mass[is.na(dat$c.mass)] <- 0
dat$n15.mass[is.na(dat$n15.mass)] <- 0
dat$nod.ex.mass[is.na(dat$nod.ex.mass)] <- 0
dat$totnodmass = with(dat, (nod.ex.mass+n15.mass+c.mass))
dat$tot.Nfixd <- with(dat, (Nfixd*totnodmass))

#calculating total N fixed per g of plant####
#Taking total N fixed per plant and dividing by the plants total mass
dat$c.mass[is.na(dat$l.mass)] <- 0
dat$n15.mass[is.na(dat$r.mass)] <- 0
dat$nod.ex.mass[is.na(dat$s.mass)] <- 0
dat$totmass = with(dat, (l.mass+r.mass+s.mass))
dat$tot.Nfixd.gplant <- with(dat, ((Nfixd*totnodmass)/totmass))

dat$tot.totnod = dat$totnodmass/dat$totmass
dat$nitro = as.factor(dat$nitro)
#Making light variable a percentage and adding x^2 for quadratic models  ####
dat$light.per2 = (dat$light.per*100)^2
dat$light.percent = dat$light.per*100

#Removing dead plants and making seperate dataframes for each species with and without zeros####
dat = dat %>%
  filter(totmass!=0)
dat$tot.Nfixd[is.na(dat$totnodmass)] <- 0
dat$tot.Nfixd[is.na(dat$tot.Nfixd)] <- 0
dat$tot.Nfixd.gplant[is.na(dat$tot.Nfixd.gplant)] <- 0
dat$nod<-ifelse(dat$totnodmass>0, yes=1, no=0)
dat$nodd<-ifelse(dat$tot.Nfixd>0, yes=1, no=0)
dat = dat %>%
  filter(!is.na(nitro))

#Making individual data frames from each legume species
dat1 = dat %>%
  filter(S.num == 1)
dat2 = dat %>%
  filter(S.num == 2)
dat3 = dat %>%
  filter(S.num == 3)
datSNF = dat %>%
  filter(!is.na(Nfixd))


dat$n = 0
dat$S.num[dat$S.num == 1] <- "P. macroloba"
dat$S.num[dat$S.num == "2"] <- "Z. longifolia"
dat$S.num[dat$S.num == "3"] <- "S. microstachyum"
dat$n[dat$nitro == 0] = "No Nitrogen Addition"
dat$n[dat$nitro == 1] = "Nitrogen Addition"
dat$light[dat$light == 0] = "Understory"
dat$light[dat$light == 1] = "Gap-Shade"
dat$light[dat$light == 2] = "Gap-Full Light"

dat$light = factor(dat$light, levels = c("Understory", "Gap-Shade", "Gap-Full Light"))

##Gap versus Intact forest graphs (figure 1)####
#Graphs of total biomass (g) ####

totbiomass.gi.all = ggplot(dat, aes(x=type, y=totmass))+
  geom_boxplot(fill = "grey56", notch = TRUE)+
  xlab ("") +
  ylab ("Total biomass (g)")+
  scale_y_continuous(expand= c(0,0),limits = c(0,35))+
  theme_pubr()+ 
  theme(text=element_text(size=15),axis.text.x=element_text(size=15, colour="black", angle=45,
                                                            vjust=1,hjust=.95))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  ggtitle("a")+
  theme(plot.title=element_text(hjust=-.1, face="bold"))+
  theme(plot.margin = unit(c(.1, .1, -.5, .25), "cm"))

totnodmass.gi.all = ggplot(dat, aes(x=type, y=totnodmass))+
  geom_boxplot(fill = "grey56", notch = TRUE)+
  xlab ("") +
  ylab(expression(atop("Nodule biomass",paste(~per~seedling~(g)))))+
  scale_y_continuous(expand= c(0,0),limits = c(0,1.5))+
  theme_pubr()+ 
  theme(text=element_text(size=15),axis.text.x=element_text(size=15, colour="black", angle=45,
                                                            vjust=1,hjust=.95))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  ggtitle("b")+
  theme(plot.title=element_text(hjust=-.1, face="bold"))+
  theme(plot.margin = unit(c(.1, .1, -.5, .25), "cm"))

totnodmass.gi.all.ins = ggplot(dat, aes(x=type, y=totnodmass))+
  geom_boxplot(fill = "grey56", notch = T)+
  xlab ("") +
  ylab ("")+
  scale_y_continuous(expand= c(0,0))+
  theme_pubr()+ 
  theme(text=element_text(size=15),axis.text.x=element_text(size=15, colour="black", angle=45,
                                                            vjust=1,hjust=.95))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  theme(axis.text.y=element_text(hjust=2.5,size=15, colour="black"))+
  coord_cartesian(xlim=c(1,2), ylim=c(0,.055), clip="on") +
  theme(plot.margin = unit(c(1, .1, -.5, -.5), "cm"))

totnodmass.gi.all. = totnodmass.gi.all + annotation_custom(ggplotGrob(totnodmass.gi.all.ins), xmin = 1.25, xmax = 2.5, 
                                                           ymin = .15, ymax = 1.5)



totnodtot.gi.all = ggplot(dat, aes(x=type, y=tot.totnod))+
  geom_boxplot(fill = "grey56", notch = T)+
  xlab ("") +
  ylab(expression(atop("Nodule allocation",paste(~(g~nodule~g^{-1}~seedling)))))+
  scale_y_continuous(expand= c(0,0),limits = c(0,.075))+
  theme_pubr()+ 
  theme(text=element_text(size=15),axis.text.x=element_text(size=15, colour="black", angle=45,
                                                            vjust=1,hjust=.95))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  ggtitle("c")+
  theme(plot.title=element_text(hjust=-.1, face="bold"))+
  theme(plot.margin = unit(c(.1, .1, .5, .5), "cm"))

totNfix.gi.all = ggplot(dat, aes(x=type, y=tot.Nfixd))+
  geom_boxplot(fill = "grey56", notch = T)+
  xlab("") +
  ylab(expression(atop("Total N fixed",paste(~(g~N~seedling^{-1}~hr^{-1})))))+
  scale_y_continuous(expand= c(0,0),limits = c(0,4))+
  theme_pubr()+ 
  theme(text=element_text(size=15),axis.text.x=element_text(size=15, colour="black", angle=45,
                                                            vjust=1,hjust=.95))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  theme(axis.text.y=element_text(hjust=2.5,size=15, colour="black"))+
  ggtitle("d")+
  theme(plot.title=element_text(hjust=-.1, face="bold"))+
  theme(plot.margin = unit(c(.1, .1, .25, .25), "cm"))

totNfix.gi.all.ins = ggplot(dat, aes(x=type, y=tot.Nfixd))+
  geom_boxplot(fill = "grey56", notch = TRUE)+
  xlab ("") +
  ylab ("")+
  scale_y_continuous(expand= c(0,0))+
  theme_pubr()+ 
  theme(text=element_text(size=15),axis.text.x=element_text(size=15, colour="black", angle=45,
                                                            vjust=1,hjust=.95))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  theme(axis.text.y=element_text(hjust=2.5,size=15, colour="black"))+
  coord_cartesian(xlim=c(1,2), ylim=c(0,.05), clip="on") +
  theme(plot.margin = unit(c(1, .1, -.5, -.5), "cm"))

totNfix.gi.all. = totNfix.gi.all + annotation_custom(ggplotGrob(totNfix.gi.all.ins), xmin = 1.25, xmax = 2.5, 
                                                     ymin = .5, ymax = 4)


ggarrange(totbiomass.gi.all,totnodmass.gi.all.,totnodtot.gi.all,totNfix.gi.all.,
          ncol=2,nrow=2,align="v",heights=c(1,1), widths=c(1,1))

quartz()
ggsave("Fig1.tiff",width = 9, height = 10, units = "in")


####Figure 2 : total biomass and nodule biomass for each legume species ####
#Calculating N fix variables from Master datasheet ####
dat<-read.csv("LaSelva_MasterData_May2019.csv")

tot.all<-ggplot(dat, aes(x=light.per, y=totmass, colour = nitro))+
  geom_point(size=3, shape=20)+
  xlab("% Light Transmittance")+
  theme_pubr()+
  theme(text=element_text(size=15),axis.text.x=element_text(size=15, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  theme(axis.text.y=element_text(size=15, colour="black"))+
  theme(axis.title.x=element_text(vjust=-4.5,size=15, colour="black"))+
  scale_color_manual(values = c("0" = "darkolivegreen3", "1" = "darkgreen"),labels = c("No N addition", "N addition"))+
  scale_x_continuous(breaks = c(.1, .2, .3, .4), labels = c("10", "20", "30", "40"))+
  theme(strip.text.x = element_text(face = "italic"))+
  ylab("Plant biomass (g)")+
  geom_smooth(data=dat,aes(x=light.per, y=totmass), method = "glm",identity = "log", se = FALSE)+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(.5,.1,1,.5),"cm"))+
  coord_cartesian(xlim=c(.05,.46), ylim=c(0,35), clip="off") +
  annotate("segment", x = 0.05, xend = 0.1, y = -4.5, yend = -4.5, colour = "grey20",size = 3)+
  annotate("segment", x = 0.09, xend = 0.2, y = -4.75, yend = -4.75, colour = "grey45",size = 3)+
  annotate("segment", x = 0.253, xend =  0.457, y =-4.75, yend = -4.75, colour = "grey70",size = 3)+
  facet_grid(.~S.num)

ggsave("Figure2_NP.tiff",width = 8, height = 6, units = "in")

tot.all.NFix<-ggplot(dat, aes(x=light.per, y=totnodmass, colour = nitro))+
  geom_point(size=3, shape=20)+
  xlab("% Light Transmittance")+
  theme_pubr()+
  theme(text=element_text(size=15),axis.text.x=element_text(size=15, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  theme(axis.text.y=element_text(size=15, colour="black"))+
  theme(axis.title.x=element_text(vjust=-4.5,size=15, colour="black"))+
  scale_color_manual(values = c("0" = "darkolivegreen3", "1" = "darkgreen"),labels = c("No N addition", "N addition"))+
  scale_x_continuous(breaks = c(.1, .2, .3, .4), labels = c("10", "20", "30", "40"))+
  theme(strip.text.x = element_text(face = "italic"))+
  ylab("Nodule biomass per seedling (g)")+
  geom_smooth(data=dat[dat$S.num=="P. macroloba",],aes(x=light.per, y=totnodmass), method = "gam", se = FALSE)+
  geom_smooth(data=dat[dat$S.num=="Z. longifolia",],aes(x=light.per, y=totnodmass, group = NA), color= "black", method = "gam", se = FALSE)+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(.5,.1,1,.5),"cm"))+
  coord_cartesian(xlim=c(.05,.46), ylim=c(0,1.3), clip="off") +
  annotate("segment", x = 0.05, xend = 0.1, y = -0.18, yend = -0.18, colour = "grey20",size = 3)+
  annotate("segment", x = 0.09, xend = 0.2, y = -.2, yend = -.2, colour = "grey45",size = 3)+
  annotate("segment", x = 0.25, xend =  0.47, y =-.2, yend = -.2, colour = "grey70",size = 3)+
  facet_grid(.~S.num)

quartz()

ggsave("Figure2_NP2.tiff",width = 8, height = 6, units = "in")



#### Figure 3: nodule biomass versus light and total N fixed versus light####

#Figure 3a: nodule biomass versus light

totnod.all<-ggplot(dat, aes(x=light.per, y=totnodmass, colour = nitro))+
  geom_point(size=3, shape=20)+
  xlab("")+
  ylab(expression(atop("Nodule biomass",paste(~per~seedling~(g)))))+
  theme_pubr()+
  theme(text=element_text(size=15),axis.text.x=element_text(size=15, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  theme(axis.text.y=element_text(hjust=2.5,size=15, colour="black"))+
  scale_color_manual(values = c("0" = "darkolivegreen3", "1" = "darkgreen"),labels = c("No N addition", "N addition"))+
  scale_x_continuous(breaks = c(.1, .2, .3, .4), labels = c("10", "20", "30", "40"))+
  geom_smooth(data=dat,aes(x=light.per, y=totnodmass), method = "gam", se = FALSE)+
  theme(legend.position = "none")+
  ggtitle("a")+
  theme(plot.title=element_text(hjust=-.1, face="bold"))+
  theme(plot.margin = unit(c(.5,.1,.1,.5),"cm"))+
  coord_cartesian(xlim=c(.05,.46), ylim=c(0,1.5), clip="off") +
  annotate("segment", x = 0.05, xend = 0.1, y = -0.18, yend = -0.18, colour = "grey20",size = 3)+
  annotate("segment", x = 0.09, xend = 0.2, y = -.2, yend = -.2, colour = "grey45",size = 3)+
  annotate("segment", x = 0.25, xend =  0.47, y =-.2, yend = -.2, colour = "grey70",size = 3)


totnod.all.ins<-ggplot(dat, aes(x=light.per, y=totnodmass, colour = nitro))+
  geom_point(size=3, shape=20)+
  xlab("")+
  theme_pubr()+
  theme(text=element_text(size=20),axis.text.x=element_text(size=20, colour="black"))+
  scale_color_manual(values = c("0" = "darkolivegreen3", "1" = "darkgreen"),labels = c("No N addition", "N addition"))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  theme(axis.text.y=element_text(hjust=2.5,size=15, colour="black"))+
  scale_x_continuous(breaks = c(.1, .2, .3, .4), labels = c("10", "20", "30", "40"))+
  scale_y_continuous(name="")+
  geom_smooth(data=dat,aes(x=light.per, y=totnodmass), method = "gam", se=FALSE)+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(.5,.1,.1,.5),"cm"))+
  coord_cartesian(xlim=c(.05,.46), ylim=c(0,.25), clip="on")+
  theme(panel.background = element_rect(fill = "grey95", colour = "grey95"))

totnodmass.all = totnod.all + annotation_custom(ggplotGrob(totnod.all.ins), xmin = .015, xmax = .35, 
                                                           ymin = .55, ymax = 1.6)

#Figure 3b: total N fixed and light
totNfixd.all<-ggplot(dat, aes(x=light.per, y=tot.Nfixd))+
  geom_point(size=3, shape=20)+
  xlab("% Light Transmittance")+
  theme_pubr()+
  theme(text=element_text(size=15),axis.text.x=element_text(size=15, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  theme(axis.text.y=element_text(hjust=2.5,size=15, colour="black"))+
  theme(axis.title.x=element_text(vjust=-1.5,size=15, colour="black"))+
  scale_x_continuous(breaks = c(.1, .2, .3, .4), labels = c("10", "20", "30", "40"))+
  ylab(expression(atop("Total N fixed",paste(~(g~N~seedling^{-1}~hr^{-1})))))+
  geom_smooth(data=dat,aes(x=light.per, y=tot.Nfixd), method = "gam",se = FALSE, colour = "black")+
  theme(legend.position = "none")+
  ggtitle("b")+
  theme(plot.title=element_text(hjust=-.1, face="bold"))+
  theme(plot.margin = unit(c(.5,.1,1,.5),"cm"))+
  coord_cartesian(xlim=c(.05,.46), ylim=c(0,4), clip="off") +
  annotate("segment", x = 0.05, xend = 0.1, y = -.5, yend = -.5, colour = "grey20",size = 3)+
  annotate("segment", x = 0.09, xend = 0.2, y = -.55, yend = -.55, colour = "grey45",size = 3)+
  annotate("segment", x = 0.253, xend =  0.457, y =-.55, yend = -.55, colour = "grey70",size = 3)


totNfixd.all.ins<-ggplot(dat, aes(x=light.per, y=tot.Nfixd))+
  geom_point(size=3, shape=20)+
  xlab("")+
  theme_pubr()+
  theme(text=element_text(size=15),axis.text.x=element_text(size=15, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  theme(axis.text.y=element_text(hjust=2.5,size=15, colour="black"))+
  scale_x_continuous(breaks = c(.1, .2, .3, .4), labels = c("10", "20", "30", "40"))+
  ylab("")+
  geom_smooth(data=dat,aes(x=light.per, y=tot.Nfixd), method = "gam", se = FALSE, colour = "black")+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(.5,.1,.1,.5),"cm"))+
  coord_cartesian(xlim=c(.05,.46), ylim=c(0,.5), clip="on")+
  theme(panel.background = element_rect(fill = "grey95", colour = "grey95"))

totNfixd.all = totNfixd.all + annotation_custom(ggplotGrob(totNfixd.all.ins), xmin = .015, xmax = .35, 
                                                ymin = 1.5, ymax = 4.25)

quartz()
ggarrange(totnodmass.all,totNfixd.all,
          ncol=1,nrow=2,align="v",heights=c(1,1), widths=c(1,1))

ggsave("Fig3.tiff",width = 6, height = 11, units = "in")


#Figure 4: nodule allocation and N fertilization

nodall.all = ggplot(dat, aes(x=nitro, y=tot.totnod))+
  geom_boxplot(fill = "grey56", notch = T)+
  xlab ("") +
  ylab(expression(atop("Nodule allocation",paste(~(g~nodule~g^{-1}~seedling)))))+
  scale_x_discrete(labels = c("Control", "N Fertilizer"))+
  scale_y_continuous(expand= c(0,0),limits = c(0,.075))+
  theme_pubr()+ 
  theme(text=element_text(size=15),axis.text.x=element_text(size=15, colour="black", angle=45,
                                                            vjust=1,hjust=.95))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  theme(axis.text.y=element_text(hjust=2.5,size=15, colour="black"))
ggsave("Figure4.tiff",width = 6, height = 7, units = "in")




####
####Supplemental Materials####


####Figure S2####
Phs<-ggplot(dat[dat$totnodmass>0,], aes(x=p.res.ppm, y=totnodmass, color = type))+
  geom_point(size=3, shape=20)+
  ylab("Nodule biomass per seedling (g)")+
  xlab(expression(atop("Phosphate availability",paste(~(mg~P~cm^{-2}~day^{-1})))))+
  theme_pubr()+
  scale_color_manual(values = c("Gap" = "darkblue", "Understory" = "lightblue"),labels = c("Gap", "Understory"))+
  theme(text=element_text(size=15),axis.text.x=element_text(size=15, colour="black"))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  theme(axis.text.y=element_text(hjust=2.5,size=15, colour="black"))+
  theme(axis.title.x=element_text(vjust=1,size=15, colour="black"))+
  theme(plot.margin = unit(c(1,1,1,1),"cm"))+
  geom_smooth(data=dat[dat$totnodmass>0,],aes(x=p.res.ppm, y=totnodmass, group = type), method = "glm", method.args = list(family = "Gamma"), se = FALSE)+
  theme(legend.position = c(.75,.75), legend.background = element_rect(fill="grey95", size =0.5, colour= "black"), legend.title = element_blank())
quartz()
Phs  

smoist<-ggplot(dat[dat$totnodmass>0,], aes(x=s.moist, y=totnodmass, color = type))+
  geom_point(size=3, shape=20)+
  ylab("Nodule biomass per seedling (g)")+
  xlab(expression(atop("Soil moisture",paste(~("%"~VWC)))))+
  theme_pubr()+
  theme(text=element_text(size=15),axis.text.x=element_text(size=15, colour="black"))+
  scale_color_manual(values = c("Gap" = "darkblue", "Understory" = "lightblue"),labels = c("Gap", "Understory"))+
  theme(panel.background=element_rect(fill="white", color="white"))+
  theme(axis.text.y=element_text(hjust=2.5,size=15, colour="black"))+
  theme(axis.title.x=element_text(vjust=1,size=15, colour="black"))+
  theme(plot.margin = unit(c(1,1,1,1),"cm"))+
  geom_smooth(data=dat[dat$totnodmass>0,],aes(x=s.moist, y=totnodmass, group = type), method = "glm", method.args = list(family = "Gamma"), se = FALSE)+
  theme(legend.position = c(.75,.75), legend.background = element_rect(fill="grey95", size =0.5, colour= "black"), legend.title = element_blank())

smoist


ggarrange(Phs, smoist, ncol = 2, nrow = 1, labels = c("a", "b"))
