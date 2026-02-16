#Metabolomics

rm(list=ls())

library(readxl)
library(vegan)

DF0=read_excel("Input/Experiment01_poop_quant.xlsx",sheet = 1)

DF1=DF0[,c(14:23)]

SAMPLES=gsub(pattern = "([A-Z0-9]*)_.*", replacement = "\\1", colnames(DF1))

DF2=t(DF1)
rownames(DF2) = SAMPLES
colnames(DF2)=DF0$`row ID`

DF3=DF2[-c(4,7),]

#####
#2)
#####

GROUPS=factor(c(1,1,2,1,1,2,2,2))

NMDS=metaMDS(DF3)

scores(NMDS,display = "sites")
plot(scores(NMDS,display = "sites"), col=GROUPS)

adonis2(DF3~GROUPS)

