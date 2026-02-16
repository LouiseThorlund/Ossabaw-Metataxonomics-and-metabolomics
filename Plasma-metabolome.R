#Metabolomics

rm(list=ls())

library(readxl)
library(vegan)

DP0=read_excel("Input/MZout_quant_plasma.xls",sheet = 1)

DP1=DP0[,c(14:18)]

PSAMPLES=gsub(pattern = "([A-Z0-9]*)_.*", replacement = "\\1", colnames(DP1))

DP2=t(DP1)
rownames(DP2) = PSAMPLES
colnames(DP2)=DP0$`row ID`

DP3=DP2[-c(3),]

#####
#2)
#####

GROUPS=factor(c(2,1,1,2))

NMDS=metaMDS(DP3)

scores(NMDS,display = "sites")
plot(scores(NMDS,display = "sites"), col=GROUPS)

adonis2(DP3~GROUPS)



