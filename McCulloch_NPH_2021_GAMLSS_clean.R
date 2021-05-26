#Analyses using generalized additive models for location, scale and shape(GAMLSS)
#for manuscript titled: Light fuels while nitrogen suppresses symbiotic nitrogen fixation hotspots in Neotropical canopy gap seedlings
#published in New Phytologist 2021


#Packages####
library("tidyverse")
library("gamlss")
library("sjPlots")

#Total BIOMASS #############

#Calculating N fix variables from Master datasheet ####
sym<-read.csv("LaSelva_MasterData_May2019.csv")

#Calculating variables ####
sym$light <- as.factor(sym$light)
sym$tot <- sym$l.mass+sym$r.mass+sym$s.mass
#Removing dead plants
sym = sym %>%
  filter(!is.na(tot))
#changing resin strip units####
sym$n.res = (sym$n.res.ppm*30)/20/14
sym$p.res = (sym$p.res.ppm*30)/20/14
sym$nitro = as.factor(sym$nitro)
sym = data.frame(sym$P.num, sym$S.num, sym$tot, sym$n.res, sym$p.res, sym$s.moist,
                 sym$light.per, sym$nitro, sym$light, sym$type, sym$Herb)

names(sym)[1] = "P.num"
names(sym)[2] ="S.num"
names(sym)[3] = "tot"
names(sym)[4] = "n.res.ppm"
names(sym)[5] = "p.res.ppm"
names(sym)[6] = "s.moist"
names(sym)[7] = "light.per"
names(sym)[8] = "nitro"
names(sym)[9] = "light"
names(sym)[10] = "type"
names(sym)[11] = "Herb"


#legume species dataframes
sym1. = sym %>%
  filter(S.num == "1")
sym2. = sym %>%
  filter(S.num == "2")
sym3. = sym %>%
  filter(S.num == "3")


m1= gamlss(tot ~ light.per + nitro + n.res.ppm + p.res.ppm + s.moist+ S.num +type +Herb, 
            data =na.omit(sym), family = NO(mu.link="log"))

summary(m1)

tab_model(m1, show.se = TRUE, show.stat = TRUE , show.ci = FALSE)


#P.mac
m1 = gamlss(tot ~ light.per + nitro + n.res.ppm + p.res.ppm + s.moist +type + Herb, 
            data =na.omit(sym1.), family = NO(mu.link="log"))

summary(m1)
tab_model(m1, show.se = TRUE, show.stat = TRUE , show.ci = FALSE)

#Z.lon

m2 = gamlss(tot ~ light.per + nitro + n.res.ppm + p.res.ppm + s.moist +type + Herb, 
            data =na.omit(sym2.), family =NO(mu.link="log"))

summary(m2)
tab_model(m2, show.se = TRUE, show.stat = TRUE , show.ci = FALSE)

#S.mic

m3 =gamlss(tot ~ light.per + nitro + n.res.ppm + p.res.ppm + s.moist + type + Herb, 
           data =na.omit(sym3.), family =NO(mu.link="log"))
summary(m3)
tab_model(m3, show.se = TRUE, show.stat = TRUE , show.ci = FALSE)



###################################
#Total nodule biomass #############
###################################

#Calculating N fix variables from Master datasheet ####
sym<-read.csv("LaSelva_MasterData_May2019.csv")

#Calculating variables ####
sym$c.mass[is.na(sym$c.mass)] <- 0
sym$n15.mass[is.na(sym$n15.mass)] <- 0
sym$nod.ex.mass[is.na(sym$nod.ex.mass)] <- 0
sym$totnod <-sym$n15.mass+sym$c.mass+sym$nod.ex.mass
sym$light <- as.factor(sym$light)
sym$tot <- sym$l.mass+sym$r.mass+sym$s.mass
sym$tot.totnod <- sym$totnod/sym$tot
sym$nitro = as.factor(sym$nitro)
sym$ratio <- (sym$l.mass+sym$s.mass)/sym$r.mass
sym = sym %>%
  filter(!is.na(tot))
#changing resin strip units####
sym$n.res = (sym$n.res.ppm*30)/20/14
sym$p.res = (sym$p.res.ppm*30)/20/14
sym = data.frame(sym$P.num, sym$S.num, sym$totnod, sym$n.res, sym$p.res, sym$s.moist,
                 sym$light.per, sym$nitro, sym$light, sym$type)

names(sym)[1] = "P.num"
names(sym)[2] ="S.num"
names(sym)[3] = "totnod"
names(sym)[4] = "n.res.ppm"
names(sym)[5] = "p.res.ppm"
names(sym)[6] = "s.moist"
names(sym)[7] = "light.per"
names(sym)[8] = "nitro"
names(sym)[9] = "light"
names(sym)[10] = "type"

sym$S.num = as.factor(sym$S.num)

#legume species dataframes
sym1. = sym %>%
  filter(S.num == "1")
sym2. = sym %>%
  filter(S.num == "2")
sym3. = sym %>%
  filter(S.num == "3")


m1 = gamlss(totnod ~ light.per + nitro + n.res.ppm + p.res.ppm + s.moist+ S.num +type, 
            data =na.omit(sym), family = ZAGA)
summary(m1)

tab_model(m1, show.se = TRUE, show.stat = TRUE , show.ci = FALSE)

#P.mac
m1 = gamlss(totnod ~ light.per + nitro + n.res.ppm + p.res.ppm + s.moist +type, 
            data =na.omit(sym1.), family = ZAGA)

summary(m1)
tab_model(m1, show.se = TRUE, show.stat = TRUE , show.ci = FALSE)

#Z.lon

m2 = gamlss(totnod ~ light.per + nitro + n.res.ppm + p.res.ppm + s.moist +type, 
            data =na.omit(sym2.), family =ZAGA)
summary(m2)

tab_model(m2, show.se = TRUE, show.stat = TRUE , show.ci = FALSE)

#S.mic

m3 =gamlss(totnod ~ light.per + nitro + n.res.ppm + p.res.ppm + s.moist + type, 
           data =na.omit(sym3.), family =ZAGA)
summary(m3)
tab_model(m3, show.se = TRUE, show.stat = TRUE , show.ci = FALSE)



#######################################
#Nodule allocation #############
#######################################


#Calculating N fix variables from Master datasheet ####
sym<-read.csv("LaSelva_MasterData_May2019.csv")

#Calculating variables ####
sym$c.mass[is.na(sym$c.mass)] <- 0
sym$n15.mass[is.na(sym$n15.mass)] <- 0
sym$nod.ex.mass[is.na(sym$nod.ex.mass)] <- 0
sym$totnod <-sym$n15.mass+sym$c.mass+sym$nod.ex.mass
sym$light <- as.factor(sym$light)
sym$tot <- sym$l.mass+sym$r.mass+sym$s.mass
sym$tot.totnod <- sym$totnod/sym$tot
sym$ratio <- (sym$l.mass+sym$s.mass)/sym$r.mass
sym$nitro = as.factor(sym$nitro)
#changing resin strip units####
sym$n.res = (sym$n.res.ppm*30)/20/14
sym$p.res = (sym$p.res.ppm*30)/20/14

sym = data.frame(sym$P.num, sym$S.num, sym$tot.totnod, sym$n.res, sym$p.res, sym$s.moist,
                 sym$light.per, sym$nitro, sym$light, sym$type)
names(sym)[1] = "P.num"
names(sym)[2] ="S.num"
names(sym)[3] = "tot.totnod"
names(sym)[4] = "n.res.ppm"
names(sym)[5] = "p.res.ppm"
names(sym)[6] = "s.moist"
names(sym)[7] = "light.per"
names(sym)[8] = "nitro"
names(sym)[9] = "light"
names(sym)[10] = "type"

sym$S.num = as.factor(sym$S.num)

#legume species dataframes
sym1. = sym %>%
  filter(S.num == "1")
sym2. = sym %>%
  filter(S.num == "2")
sym3. = sym %>%
  filter(S.num == "3")


m1 = gamlss(tot.totnod ~ light.per + nitro + n.res.ppm + p.res.ppm + s.moist+ S.num +type, 
            data =na.omit(sym), family =ZAGA)

summary(m1)

tab_model(m1, show.se = TRUE, show.stat = TRUE , show.ci = FALSE)

#P.mac


m1 = gamlss(tot.totnod ~ light.per + nitro + n.res.ppm + p.res.ppm + s.moist +type, 
            data =na.omit(sym1.), family = ZAGA)

summary(m1)

#Z.lon

m2 = gamlss(tot.totnod ~ light.per + nitro + n.res.ppm + p.res.ppm + s.moist +type, 
            data =na.omit(sym2.), family =ZAGA)
summary(m2)

tab_model(m2, show.se = TRUE, show.stat = TRUE , show.ci = FALSE)

#S.mic

m3 =gamlss(tot.totnod ~ light.per + nitro + n.res.ppm + p.res.ppm + s.moist + type, 
           data =na.omit(sym3.), family =ZAGA)
summary(m3)

tab_model(m3, show.se = TRUE, show.stat = TRUE , show.ci = FALSE)

#######################################
#Total N fixed #############
#######################################


#Calculating N fix variables from Master datasheet ####
sym<-read.csv("LaSelva_MasterData_May2019.csv")

#Changing total N in isotope samples from micrograms to miligrams 
sym$Total.N.e <- sym$total.N.15n*.001
sym$Total.N.c <- sym$total.N.c*.001
#Calculating %N for each sample from the total N and sample weight of enriched samples ####
sym$XN.e<-((sym$Total.N.e)/sym$mg.tin.e)*100

#Calculating %N for each sample from the total N and sample weight of control samples ####
sym$XN.c<-((sym$total.N.c)/sym$mg.tin.c)*100

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
sym$XNfixdX<-with(sym, (X15n.atX.15n-X15n.atX.c)/syr.enrich)

#Some samples only ended up having values for controls because enriched samples were below detection limit, 
#to avoid negatives I made them zeros
sym$XNfixdX[sym$XNfixdX < 0] = 0 

#Calculating the total amount of N (g) in the incubated sample amount. n15.mass is in grams
sym$tot.N<-with(sym, (n15.mass*XN.e))

#calculating the amount of N fixed in each sample (gN fixed per sample per incubation period) check units?
sym$Nfixd.sample<-with(sym, tot.N*XNfixdX)


#Calculating N fixation to units of (gN fixed per gram of nodule per hour) "SNF efficiency"####
#gN fixed per sample divided by the mass of nodules that were incubated (n15.mass)
sym$Nfixd<-with(sym, (Nfixd.sample/n15.mass)*2)
#Converting N fixation to units of (gN fixed per gram of nodule per year)
time.conv<-17520 #the number of "half hours" in a year
sym$Nfixd.yr<-with(sym, (Nfixd.sample/n15.mass)*time.conv)

#calculating total N fixed per plant ####
#Taking "SNF efficiency" and multiplying by the total nodule mass 
sym$c.mass[is.na(sym$c.mass)] <- 0
sym$n15.mass[is.na(sym$n15.mass)] <- 0
sym$nod.ex.mass[is.na(sym$nod.ex.mass)] <- 0
sym$totnodmass = with(sym, (nod.ex.mass+n15.mass+c.mass))
sym$tot.Nfixd <- with(sym, (Nfixd*totnodmass))
sym$tot.Nfixd[(sym$totnod == 0)] <- 0

#calculating total N fixed per g of plant####
#Taking total N fixed per plant and dividing by the plants total mass
sym$c.mass[is.na(sym$l.mass)] <- 0
sym$n15.mass[is.na(sym$r.mass)] <- 0
sym$nod.ex.mass[is.na(sym$s.mass)] <- 0
sym$totmass = with(sym, (l.mass+r.mass+s.mass))
sym$tot.Nfixd.gplant <- with(sym, ((Nfixd*totnodmass)/totmass))

#changing resin strip units####
sym$n.res = (sym$n.res.ppm*30)/20/14
sym$p.res = (sym$p.res.ppm*30)/20/14

sym$nitro = as.factor(sym$nitro)

sym = data.frame(sym$P.num, sym$S.num, sym$tot.Nfixd, sym$n.res, sym$p.res, sym$s.moist,
                 sym$light.per, sym$nitro, sym$light, sym$type)
names(sym)[1] = "P.num"
names(sym)[2] ="S.num"
names(sym)[3] = "tot.Nfixd"
names(sym)[4] = "n.res.ppm"
names(sym)[5] = "p.res.ppm"
names(sym)[6] = "s.moist"
names(sym)[7] = "light.per"
names(sym)[8] = "nitro"
names(sym)[9] = "light"
names(sym)[10] = "type"

sym$S.num = as.factor(sym$S.num)

#legume species dataframes
sym1. = sym %>%
  filter(S.num == "1")
sym2. = sym %>%
  filter(S.num == "2")
sym3. = sym %>%
  filter(S.num == "3")


m1 = gamlss(tot.Nfixd ~ light.per + nitro + n.res.ppm + p.res.ppm + s.moist+ S.num +type, 
            data =na.omit(sym), family =ZAGA)

summary(m1)

tab_model(m1, show.se = TRUE, show.stat = TRUE , show.ci = FALSE)

#P.mac
m1 = gamlss(tot.Nfixd ~ light.per * nitro + n.res.ppm + p.res.ppm + s.moist +type, 
            data =na.omit(sym1.), family =ZAGA)

summary(m1)

#Z.lon

m2 = gamlss(tot.Nfixd ~ light.per + nitro + n.res.ppm + p.res.ppm + s.moist +type, 
            data =na.omit(sym2.), family =ZAGA)
summary(m2)

tab_model(m2, show.se = TRUE, show.stat = TRUE , show.ci = FALSE)

#S.mic

m3 =gamlss(tot.Nfixd ~ light.per + nitro + n.res.ppm + p.res.ppm + s.moist + type, 
           data =na.omit(sym3.), family =ZAGA)

summary(m3)
tab_model(m3, show.se = TRUE, show.stat = TRUE , show.ci = FALSE)


#######################################
#Total N fixed per g of plant#############
#######################################


#Calculating N fix variables from Master datasheet ####
sym<-read.csv("LaSelva_MasterData_May2019.csv")

#Changing total N in isotope samples from micrograms to miligrams 
sym$Total.N.e <- sym$total.N.15n*.001
sym$Total.N.c <- sym$total.N.c*.001
#Calculating %N for each sample from the total N and sample weight of enriched samples ####
sym$XN.e<-((sym$Total.N.e)/sym$mg.tin.e)*100

#Calculating %N for each sample from the total N and sample weight of control samples ####
sym$XN.c<-((sym$total.N.c)/sym$mg.tin.c)*100

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
sym$XNfixdX<-with(sym, (X15n.atX.15n-X15n.atX.c)/syr.enrich)

#Some samples only ended up having values for controls because enriched samples were below detection limit, 
#to avoid negatives I made them zeros
sym$XNfixdX[sym$XNfixdX < 0] = 0 

#Calculating the total amount of N (g) in the incubated sample amount. n15.mass is in grams
sym$tot.N<-with(sym, (n15.mass*XN.e))

#calculating the amount of N fixed in each sample (gN fixed per sample per incubation period) check units?
sym$Nfixd.sample<-with(sym, tot.N*XNfixdX)


#Calculating N fixation to units of (gN fixed per gram of nodule per hour) "SNF efficiency"####
#gN fixed per sample divided by the mass of nodules that were incubated (n15.mass)
sym$Nfixd<-with(sym, (Nfixd.sample/n15.mass)*2)
#Converting N fixation to units of (gN fixed per gram of nodule per year)
time.conv<-17520 #the number of "half hours" in a year
sym$Nfixd.yr<-with(sym, (Nfixd.sample/n15.mass)*time.conv)

#calculating total N fixed per plant ####
#Taking "SNF efficiency" and multiplying by the total nodule mass 
sym$c.mass[is.na(sym$c.mass)] <- 0
sym$n15.mass[is.na(sym$n15.mass)] <- 0
sym$nod.ex.mass[is.na(sym$nod.ex.mass)] <- 0
sym$totnodmass = with(sym, (nod.ex.mass+n15.mass+c.mass))
sym$tot.Nfixd <- with(sym, (Nfixd*totnodmass))
sym$tot.Nfixd[is.na(sym$tot.Nfixd)] <- 0

#calculating total N fixed per g of plant####
#Taking total N fixed per plant and dividing by the plants total mass
sym$c.mass[is.na(sym$l.mass)] <- 0
sym$n15.mass[is.na(sym$r.mass)] <- 0
sym$nod.ex.mass[is.na(sym$s.mass)] <- 0
sym$totmass = with(sym, (l.mass+r.mass+s.mass))
sym$tot.Nfixd.gplant <- with(sym, ((Nfixd*totnodmass)/totmass))
sym = sym %>%
  filter(!is.na(tot.Nfixd)) %>%
  filter(!is.na(tot))

#changing resin strip units####
sym$n.res = (sym$n.res.ppm*30)/20/14
sym$p.res = (sym$p.res.ppm*30)/20/14
sym$nitro = as.factor(sym$nitro)

sym = data.frame(sym$P.num, sym$S.num, sym$tot.Nfixd.gplant, sym$n.res, sym$p.res, sym$s.moist,
                 sym$light.per, sym$nitro, sym$light, sym$type)
names(sym)[1] = "P.num"
names(sym)[2] ="S.num"
names(sym)[3] = "tot.Nfixd.gplant"
names(sym)[4] = "n.res.ppm"
names(sym)[5] = "p.res.ppm"
names(sym)[6] = "s.moist"
names(sym)[7] = "light.per"
names(sym)[8] = "nitro"
names(sym)[9] = "light"
names(sym)[10] = "type"

sym$S.num = as.factor(sym$S.num)

#legume species dataframes
sym1. = sym %>%
  filter(S.num == "1")
sym2. = sym %>%
  filter(S.num == "2")
sym3. = sym %>%
  filter(S.num == "3")

m1 = gamlss(tot.Nfixd.gplant ~ light.per + nitro + n.res.ppm + p.res.ppm + s.moist+ S.num +type, 
            data =na.omit(sym), family =ZAGA)

summary(m1)

tab_model(m1, show.se = TRUE, show.stat = TRUE , show.ci = FALSE)

#P.mac
m1 = gamlss(tot.Nfixd.gplant ~ light.per + nitro + n.res.ppm + p.res.ppm + s.moist +type, 
            data =na.omit(sym1.), family =ZAGA)

summary(m1)


#Z.lon

m2 = gamlss(tot.Nfixd.gplant ~ light.per + nitro + n.res.ppm + p.res.ppm + s.moist +type, 
            data =na.omit(sym2.), family =ZAGA)
summary(m2)

tab_model(m2, show.se = TRUE, show.stat = TRUE , show.ci = FALSE)

#S.mic

m3 =gamlss(tot.Nfixd.gplant ~ light.per + nitro + n.res.ppm + p.res.ppm + s.moist + type, 
           data =na.omit(sym3.), family =ZAGA)
summary(m3)

tab_model(m3, show.se = TRUE, show.stat = TRUE , show.ci = FALSE)



#######################################
#SNF Efficiency#############
#######################################


#Calculating N fix variables from Master datasheet ####
sym<-read.csv("LaSelva_MasterData_May2019.csv")
sym$nitro = as.factor(sym$nitro)
#Changing total N in isotope samples from micrograms to miligrams 
sym$Total.N.e <- sym$total.N.15n*.001
sym$Total.N.c <- sym$total.N.c*.001
#Calculating %N for each sample from the total N and sample weight of enriched samples ####
sym$XN.e<-((sym$Total.N.e)/sym$mg.tin.e)*100

#Calculating %N for each sample from the total N and sample weight of control samples ####
sym$XN.c<-((sym$total.N.c)/sym$mg.tin.c)*100

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
sym$XNfixdX<-with(sym, (X15n.atX.15n-X15n.atX.c)/syr.enrich)

#Some samples only ended up having values for controls because enriched samples were below detection limit, 
#to avoid negatives I made them zeros
sym$XNfixdX[sym$XNfixdX < 0] = 0 

#Calculating the total amount of N (g) in the incubated sample amount. n15.mass is in grams
sym$tot.N<-with(sym, (n15.mass*XN.e))

#calculating the amount of N fixed in each sample (gN fixed per sample per incubation period) check units?
sym$Nfixd.sample<-with(sym, tot.N*XNfixdX)


#Calculating N fixation to units of (gN fixed per gram of nodule per hour) "SNF efficiency"####
#gN fixed per sample divided by the mass of nodules that were incubated (n15.mass)
sym$Nfixd<-with(sym, (Nfixd.sample/n15.mass)*2)
sym$totnodmass = with(sym, (nod.ex.mass+n15.mass+c.mass))

#changing resin strip units####
sym$n.res = (sym$n.res.ppm*30)/20/14
sym$p.res = (sym$p.res.ppm*30)/20/14

sym = data.frame(sym$P.num, sym$S.num, sym$Nfixd, sym$n.res, sym$p.res, sym$s.moist,
                 sym$light.per, sym$nitro, sym$light, sym$type)
names(sym)[1] = "P.num"
names(sym)[2] ="S.num"
names(sym)[3] = "Nfixd"
names(sym)[4] = "n.res.ppm"
names(sym)[5] = "p.res.ppm"
names(sym)[6] = "s.moist"
names(sym)[7] = "light.per"
names(sym)[8] = "nitro"
names(sym)[9] = "light"
names(sym)[10] = "type"


sym$S.num = as.factor(sym$S.num)

#legume species dataframes
sym1. = sym %>%
  filter(S.num == "1")
sym2. = sym %>%
  filter(S.num == "2")
sym3. = sym %>%
  filter(S.num == "3")

m1 = gamlss(Nfixd ~ light.per + nitro + n.res.ppm + p.res.ppm + s.moist+ S.num +type, 
            data =na.omit(sym), family =ZAGA)

summary(m1)

tab_model(m1, show.se = TRUE, show.stat = TRUE , show.ci = FALSE)

#P.mac
m1 = gamlss(Nfixd ~ light.per + nitro + n.res.ppm + p.res.ppm + s.moist +type, 
            data =na.omit(sym1.), family =ZAGA)

summary(m1)

#Z.lon

m2 = gamlss(Nfixd ~ light.per + nitro + n.res.ppm + p.res.ppm + s.moist +type, 
            data =na.omit(sym2.), family =ZAGA)
summary(m2)
tab_model(m2, show.se = TRUE, show.stat = TRUE , show.ci = FALSE)

#S.mic

m3 =gamlss(Nfixd ~ light.per + nitro + n.res.ppm + p.res.ppm + s.moist + type, 
           data =na.omit(sym3.), family =ZAGA)
summary(m3)
tab_model(m3, show.se = TRUE, show.stat = TRUE , show.ci = FALSE)

