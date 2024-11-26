#Load in packages
library(dplyr)
library(ggplot2)
library(emmeans)
library(table1)
library(relaimpo)

adnimerge <- read.csv("ADNIMERGE.csv")
uwn <- read.csv("UWNPSYCHSUM_03_09_21.csv")
wmh <- read.csv("ADNI_UCD_WMH_04_30_19.csv", header=T)
nfl <- read.csv("ADNI_BLENNOWPLASMANFLLONG_10_03_18.csv", header=T)
gap43 <- read.csv("BLENNOW_LAB_CSF_GAP_43_06_08_21.csv",header=T)
strem2 <- read.csv("ADNI_HAASS_WASHU.csv", header=T)
upenn <- read.csv("UPENNBIOMK9_04_19_17.csv")
pet <- read.csv("UCBERKELEYFDG_8mm_02_17_23_22May2023.csv")
recmhist <- read.csv("RECMHIST.csv",header=T)


#-----------------------------------------------------------------------------#
#-----------------------CLEANING AND MERGING FILES----------------------------#
#-----------------------------------------------------------------------------#

adnimerge_filter <- adnimerge %>% filter(ORIGPROT=='ADNI2' | ORIGPROT=='ADNIGO', COLPROT=='ADNI2' | COLPROT=='ADNIGO')

uwn_filter <- uwn %>% filter(VISCODE2=='bl')
wmh_filter <- wmh %>% filter(VISCODE2=='scmri')
nfl_filter  <- nfl %>% filter(VISCODE2=='bl')
gap43_filter <- gap43 %>% filter(VISCODE2=='bl')
strem2_filter <- strem2 %>% filter(VISCODE2=='bl')
upenn_filter <- upenn %>% filter(VISCODE2=='bl')
pet_filter <- pet %>% filter(VISCODE2=='bl')

adnimerge_filter <- adnimerge_filter %>% dplyr::select(RID,VISCODE,DX_bl,AGE,PTGENDER,PTEDUCAT,MMSE,APOE4,ADAS13,Hippocampus_bl, ICV_bl)
wmh_filter <- wmh_filter %>% dplyr::select(RID, VISCODE2, CSF, GRAY, WHITE, WMH, TOTAL)
nfl_filter  <- nfl_filter %>% dplyr::select(RID, VISCODE2, PLASMA_NFL)
gap43_filter <- gap43_filter %>% dplyr::select(RID, VISCODE2, GAP_43)
uwn_filter <- uwn_filter %>% dplyr::select(RID, VISCODE2, ADNI_MEM, ADNI_EF)
strem2_filter <- strem2_filter %>% dplyr::select(RID, VISCODE2, WU_STREM2CORRECTED, MSD_STREM2CORRECTED)
upenn_filter <- upenn_filter %>% dplyr::select(RID,VISCODE2,ABETA,PTAU,TAU)
pet_filter <- pet_filter %>% dplyr::select(RID,ROINAME,MEAN,VISCODE2)

pet_filter <- aggregate(MEAN ~RID,pet_filter, function(x) x[1]/x[-1])
pet_filter$VISCODE2 <- "bl"
names(pet_filter)[names(pet_filter) == 'MEAN'] <- "metaROI"

names(adnimerge_filter)[names(adnimerge_filter) == 'VISCODE'] <- "VISCODE2"
wmh_filter$VISCODE2[wmh_filter$VISCODE2=='scmri'] <- 'bl'

all_data <- Reduce(function(x,y) merge(x,y, by=c("RID","VISCODE2")), list(adnimerge_filter,uwn_filter,nfl_filter,gap43_filter,wmh_filter,strem2_filter,upenn_filter, pet_filter))
all_data <- all_data %>% dplyr::select(RID,VISCODE2,AGE,PTEDUCAT,PTGENDER,ADNI_MEM,ADNI_EF,PLASMA_NFL,GAP_43, MSD_STREM2CORRECTED,MMSE,DX_bl,WMH,ABETA,TAU,PTAU,APOE4,TOTAL,Hippocampus_bl, ICV_bl, metaROI)

#Clean biomarker, WMH, APOE4, sex data
all_data$WMH <- all_data$WMH / all_data$ICV_bl

all_data$TAU[as.numeric(all_data$TAU) > 1300] <- 1300
all_data$TAU[as.numeric(all_data$TAU) < 80] <- 80

all_data$PTAU[all_data$PTAU=="<8"] <- 7
all_data$PTAU[as.numeric(all_data$PTAU) > 120] <- 120
all_data$PTAU[as.numeric(all_data$PTAU) < 8] <- 8

all_data$ABETA[all_data$ABETA==">1700"] <- 1800
all_data$ABETA[as.numeric(all_data$ABETA) > 1700] <- 1700
all_data$ABETA[as.numeric(all_data$ABETA) < 200] <- 200

all_data$APOE4[all_data$APOE4=="2"] <- 1

all_data$PTGENDER <- ifelse(all_data$PTGENDER=="Female", 1, 0)
all_data$PTGENDER <- as.factor(all_data$PTGENDER)

cor.test(all_data$GAP_43, all_data$PTAU) #collinearity; use FDG_PET as "N" biomarker


#Duplicates Check & Removal
all_data[all_data$RID %in% all_data$RID[duplicated(all_data$RID)],]

all_data <- all_data[!all_data$RID==4114 | all_data$PLASMA_NFL == 12.7,]
all_data <- all_data[!all_data$RID==4390 | all_data$PLASMA_NFL == 36.4,]
all_data <- all_data[!all_data$RID==5125 | all_data$PLASMA_NFL == 16.4,]
all_data <- all_data[!all_data$RID==4376,]


#Pull vascular pathology covariates
hx <- subset(recmhist,RID %in% all_data$RID)

#Hypertension data
h1 <- grepl('hypertension', hx$MHDESC) | grepl('HBP', hx$MHDESC)
hx$HYPT[h1] <- "1"
hx$HYPT[is.na(hx$HYPT)] <- 0

#Diabetes data
d1 <- grepl('diabetes', hx$MHDESC) | grepl('diabetic', hx$MHDESC) | grepl('insulin', hx$MHDESC)
hx$DIAB[d1] <- "1"
hx$DIAB[is.na(hx$DIAB)] <- 0

#Hyperlipidemia data
hl1 <- grepl('hyperlipidemia', hx$MHDESC)
hx$HYPL[hl1] <- "1"
hx$HYPL[is.na(hx$HYPL)] <- 0

#combine hypt, diabetes, hypl data
hx <- hx %>% dplyr::select(RID,HYPT,DIAB,HYPL)
hx <- hx %>% distinct(RID, .keep_all = TRUE)
length(hx$RID)

#Smoking data
medhist <- read.csv("MEDHIST.csv",header=T)
medhist <- medhist %>% dplyr::select(RID, MH16SMOK)
medhist <- subset(medhist,RID %in% all_data$RID)
medhist[medhist$RID %in% medhist$RID[duplicated(medhist$RID)],]
medhist <- medhist %>% distinct(RID, .keep_all = TRUE)

vasc_cov <- merge(medhist, hx , by="RID")
all_data <- merge(all_data,vasc_cov, by='RID')


#Missing data check
cbind(lapply(lapply(all_data, is.na), sum)) 
all_data$MSD_STREM2CORRECTED[is.na(all_data$MSD_STREM2CORRECTED)] <- median(all_data$MSD_STREM2CORRECTED, na.rm=TRUE)
all_data$PLASMA_NFL[is.na(all_data$PLASMA_NFL)] <- median(all_data$PLASMA_NFL, na.rm=TRUE)


#Outliers check
scattOut <- function(VOI) {
  upper3 = mean(VOI,na.rm=TRUE) + 3*sd(VOI,na.rm=TRUE)
  lower3 = mean(VOI,na.rm=TRUE) - 3*sd(VOI,na.rm=TRUE)
  upper5 = mean(VOI,na.rm=TRUE) + 5*sd(VOI,na.rm=TRUE)
  lower5 = mean(VOI,na.rm=TRUE) - 5*sd(VOI,na.rm=TRUE)
  plot(VOI)
  abline(h=upper3,col="blue")
  abline(h=lower3,col="blue")
  abline(h=upper5,col="red")
  abline(h=lower5,col="red")
}

scattOut(all_data$WMH)
scattOut(all_data$PLASMA_NFL)
scattOut(all_data$MSD_STREM2CORRECTED)
scattOut(all_data$GAP_43)
scattOut(all_data$ADNI_MEM)
scattOut(all_data$ADNI_EF)

#WMH outliers
uppermax_WMH <- mean(all_data$WMH) + 3*sd(all_data$WMH)
lowermax_WMH <- mean(all_data$WMH) - 3*sd(all_data$WMH)
WMH_out <- all_data$WMH[all_data$WMH > lowermax_WMH & all_data$WMH < uppermax_WMH]
all_data[which(all_data$WMH > uppermax_WMH | all_data$WMH < lowermax_WMH),] 
WMH_outliers <- c(2060, 2205, 4067, 4071, 4131, 4449, 4804, 5005, 5197)

#NFL outliers
uppermax_NFL <- mean(all_data$PLASMA_NFL) + 3*sd(all_data$PLASMA_NFL)
lowermax_NFL <- mean(all_data$PLASMA_NFL) - 3*sd(all_data$PLASMA_NFL)
all_data[which(all_data$PLASMA_NFL > uppermax_NFL | all_data$PLASMA_NFL < lowermax_NFL),] 
NFL_outliers <- c(2058, 2220, 4001, 4410, 4594, 4719, 4877, 5165)

#sTREM2 outliers
uppermax_strem <- mean(all_data$MSD_STREM2CORRECTED) + 3*sd(all_data$MSD_STREM2CORRECTED)
lowermax_strem <- mean(all_data$MSD_STREM2CORRECTED) - 3*sd(all_data$MSD_STREM2CORRECTED)
all_data[which(all_data$MSD_STREM2CORRECTED > uppermax_strem | all_data$MSD_STREM2CORRECTED < lowermax_strem),] 
strem_outliers <- c(2068, 2099, 4075, 4172, 4254, 4576, 4616, 4793, 4803, 5067)

#GAP43 outliers
uppermax_gap <- mean(all_data$GAP_43) + 3*sd(all_data$GAP_43)
lowermax_gap <- mean(all_data$GAP_43) - 3*sd(all_data$GAP_43)
all_data[which(all_data$GAP_43 > uppermax_gap | all_data$GAP_43 < lowermax_gap),]
gap_outliers <- c(2099, 2155, 4386, 4646, 476, 4853, 5018, 5251)

#ADNI_MEM outliers
uppermax_mem <- mean(all_data$ADNI_MEM) + 3*sd(all_data$ADNI_MEM)
lowermax_mem <- mean(all_data$ADNI_MEM) - 3*sd(all_data$ADNI_MEM)
all_data[which(all_data$ADNI_MEM > uppermax_mem | all_data$ADNI_MEM < lowermax_mem),] 
mem_outliers <- c(4599)

#ADNI_EF outliers
uppermax_ef <- mean(all_data$ADNI_EF) + 3*sd(all_data$ADNI_EF)
lowermax_ef <- mean(all_data$ADNI_EF) - 3*sd(all_data$ADNI_EF)
all_data[which(all_data$ADNI_EF > uppermax_ef | all_data$ADNI_EF < lowermax_ef),]
ef_outliers <- c(5005)

#Remove Outliers
outliers <- as.numeric(c(WMH_outliers,NFL_outliers,strem_outliers,gap_outliers,mem_outliers,ef_outliers))
outliers <- data.frame(outliers)
outliers <- outliers %>% distinct(outliers, .keep_all = TRUE)
all_data <- all_data[!(all_data$RID %in% outliers$outliers),]


#Checking data distributions
#continuous variables checked - normally distributed
shapiro.test(all_data$WMH)
    qqnorm(all_data$WMH)
shapiro.test(all_data$PLASMA_NFL)
    qqnorm(all_data$PLASMA_NFL)
shapiro.test(all_data$GAP_43)
    qqnorm(all_data$GAP_43)
shapiro.test(all_data$MSD_STREM2CORRECTED)
    qqnorm(all_data$MSD_STREM2CORRECTED)


#Data Log transformations
all_data$WMH_log <- log(all_data$WMH)
all_data$NFL_log <- log(all_data$PLASMA_NFL)
all_data$GAP43_log <- log(all_data$GAP_43)
all_data$MSD_STREM2_log <- log(all_data$MSD_STREM2CORRECTED)


#Standardize variables
all_data$PTEDUCAT_norm <- scale(all_data$PTEDUCAT)
all_data$AGE_norm <- scale(all_data$AGE)
all_data$MMSE_norm <- scale(all_data$MMSE)
all_data$NFL_log_norm <- scale(all_data$NFL_log)
all_data$GAP43_log_norm <- scale(all_data$GAP43_log)
all_data$MSD_STREM2_log_norm <- scale(all_data$MSD_STREM2_log)
all_data$WMH_log_norm <- scale(all_data$WMH_log)
all_data$ABETA_norm <- scale(as.numeric(all_data$ABETA))
all_data$PTAU_norm <- scale(as.numeric(all_data$PTAU))
all_data$metaROI_norm <- scale(as.numeric(all_data$metaROI))


#---------------------- AT/N Classification System ---------------------------#
#1. Binarize data
all_data$ABETA <- as.numeric(all_data$ABETA)
for (i in 1:length(all_data$ABETA)) {
  if (all_data$ABETA[i] < 977) {
    (all_data$ABETA_atn[i] <- 1)
  } else if (all_data$ABETA[i] >= 977) {
    (all_data$ABETA_atn[i] <- 0)
  }
}

all_data$PTAU <- as.numeric(all_data$PTAU)
for (i in 1:length(all_data$PTAU)) {
  if (all_data$PTAU[i] > 27) {
    (all_data$PTAU_atn[i] <- 1)
  } else if (all_data$PTAU[i] <= 27) {
    (all_data$PTAU_atn[i] <- 0)
  }
}

all_data$metaROI<- as.numeric(all_data$metaROI)
for (i in 1:length(all_data$metaROI)) {
  if (all_data$metaROI[i] <= 1.21) {
    (all_data$PET_atn[i] <- 1)
  } else if (all_data$metaROI[i] > 1.21) {
    (all_data$PET_atn[i] <- 0)
  }
}

#2. Create ATN classification strings - "ATN_class" [A(+/-), T(+/-), N(+/-)]
# [A+/-] for ABETA
for (i in 1:length(all_data$ABETA_atn)) {
  if (all_data$ABETA_atn[i] == 1) {
    (all_data$ATN_class[i] <- "A+")
  } else if (all_data$ABETA_atn[i] == 0) {
    (all_data$ATN_class[i] <- "A-")
  }
}

for (i in 1:length(all_data$PTAU_atn)) {
  if (all_data$PTAU_atn[i] == 1) {
    (all_data$ATN_class[i] <- paste(all_data$ATN_class[i],"T+"))
  } else if (all_data$PTAU_atn[i] == 0) {
    (all_data$ATN_class[i] <- paste(all_data$ATN_class[i],"T-"))
  }
}

for (i in 1:length(all_data$PET_atn)) {
  if (all_data$PET_atn[i] == 1) {
    (all_data$ATN_class[i] <- paste(all_data$ATN_class[i],"N+"))
  } else if (all_data$PET_atn[i] == 0) {
    (all_data$ATN_class[i] <- paste(all_data$ATN_class[i],"N-"))
  }
}

#3. Create ATN classification groups
for (i in 1:length(all_data$ATN_class)) {
  if (all_data$ATN_class[i] == "A- T- N-") {
    (all_data$ATN_group[i] <- "CN")
  } else if (all_data$ATN_class[i] == "A- T+ N+" | all_data$ATN_class[i] == "A- T- N+" | all_data$ATN_class[i] == "A- T+ N-") {
    (all_data$ATN_group[i] <- "SNAP")
  }
  else if (all_data$ATN_class[i] == "A+ T- N-" | all_data$ATN_class[i] =="A+ T+ N-" | all_data$ATN_class[i] == "A+ T+ N+" | all_data$ATN_class[i] == "A+ T- N+" )
    (all_data$ATN_group[i] <- "AD")
}


#------------------------ DESCRIPTIVE STATISTICS ------------------------------#
#Age
summary(aov(AGE ~ ATN_group, data=all_data))
TukeyHSD(aov(AGE ~ ATN_group, data=all_data))

#Education
all_data$PTEDUCAT <- as.numeric(all_data$PTEDUCAT)
summary(aov(PTEDUCAT ~ ATN_group, data=all_data))

#ABETA
summary(aov(ABETA ~ ATN_group, data=all_data))
TukeyHSD(aov(ABETA ~ ATN_group, data=all_data))

#T-TAU
summary(aov(TAU ~ ATN_group, data=all_data))
TukeyHSD(aov(TAU ~ ATN_group, data=all_data))

#P-TAU
summary(aov(PTAU ~ ATN_group, data=all_data))
TukeyHSD(aov(PTAU ~ ATN_group, data=all_data))
aggregate(all_data$PTAU, list(all_data$ATN_group), FUN=mean) 

#ADNI_EF
summary(aov(ADNI_EF ~ ATN_group, data=all_data))
TukeyHSD(aov(ADNI_EF ~ ATN_group, data=all_data))

#ADNI_MEM
summary(aov(ADNI_MEM ~ ATN_group, data=all_data))
TukeyHSD(aov(ADNI_MEM ~ ATN_group, data=all_data))

#NfL
summary(aov(PLASMA_NFL ~ ATN_group, data=all_data)) 
TukeyHSD(aov(PLASMA_NFL ~ ATN_group, data=all_data))

#GAP-43
summary(aov(GAP_43 ~ ATN_group, data=all_data))
TukeyHSD(aov(GAP_43 ~ ATN_group, data=all_data))

#sTREM2
summary(aov(MSD_STREM2CORRECTED ~ ATN_group, data=all_data)) 
TukeyHSD(aov(MSD_STREM2CORRECTED ~ ATN_group, data=all_data))

#WMH
summary(aov(WMH ~ ATN_group, data=all_data))
TukeyHSD(aov(WMH ~ ATN_group, data=all_data))

#MMSE
summary(aov(MMSE ~ ATN_group, data=all_data))
TukeyHSD(aov(MMSE ~ ATN_group, data=all_data))
aggregate(all_data$MMSE, list(all_data$ATN_group), FUN=mean) 

#FDG-PET
summary(aov(metaROI ~ ATN_group, data=all_data))
TukeyHSD(aov(metaROI ~ ATN_group, data=all_data))


#Categorical variables
#Sex
chisq.test(all_data$ATN_group, all_data$PTGENDER, correct=FALSE)
all_data %>% group_by(ATN_group,PTGENDER) %>% summarize((Freq=n()))

#Sex post-hoc
sex_1 <- all_data %>% filter(ATN_group=='CN' | ATN_group=='SNAP')
chisq.test(sex_1$ATN_group, sex_1$PTGENDER, correct=FALSE)

sex_2 <- all_data %>% filter(ATN_group=='CN' | ATN_group=='AD')
chisq.test(sex_2$ATN_group, sex_2$PTGENDER, correct=FALSE)

sex_3 <- all_data %>% filter(ATN_group=='SNAP' | ATN_group=='AD')
chisq.test(sex_3$ATN_group, sex_3$PTGENDER, correct=FALSE)

#APOE-4
all_data$APOE4 <- factor(all_data$APOE4, levels=c(1,0))
all_data$APOE4 <- as.numeric(all_data$APOE4)
all_data %>% group_by(ATN_group,APOE4) %>% summarize((Freq=n()))

chisq.test(all_data$ATN_group, all_data$APOE4, correct=FALSE)

#APOE4 post-hoc
apoe_1 <- all_data %>% filter(ATN_group=='CN' | ATN_group=='SNAP')
chisq.test(apoe_1$ATN_group, apoe_1$APOE4, correct=FALSE)

apoe_2 <- all_data %>% filter(ATN_group=='CN' | ATN_group=='AD')
chisq.test(apoe_2$ATN_group, apoe_2$APOE4, correct=FALSE)

apoe_3 <- all_data %>% filter(ATN_group=='SNAP' | ATN_group=='AD')
chisq.test(apoe_3$ATN_group, apoe_3$APOE4, correct=FALSE)



#----------------------------------------------------------------------------#
#------------------- Table 1. Participant Demographics ----------------------#
#----------------------------------------------------------------------------#
all_data$ATN_group <- factor(all_data$ATN_group, levels=c("CN","SNAP","AD"))    
levels(all_data$ATN_group) <-c("CN","SNAP","AD")
    
rndr <- function(x, name, ...) {
  if (!is.numeric(x)) return(render.categorical.default(x))
  what <- switch(name,
                 AGE = "Mean (SD)",
                 PTEDUCAT="Mean (SD)",
                 APOE4="Mean (SD)",
                 ABETA="Mean (SD)",
                 TAU="Mean (SD)",
                 metaROI = "Mean (SD)",
                 PTAU="Mean (SD)",
                 MMSE="Mean (SD)",
                 ADNI_EF="Mean (SD)",
                 ADNI_MEM="Mean (SD)",
                 WMH = "Mean (SD)",
                 PLASMA_NFL = "Mean (SD)",
                 GAP_43 = "Mean (SD)",
                 MSD_STREM2CORRECTED = "Mean (SD)")
  parse.abbrev.render.code(c("", what))(x)
}

labels <- list(
  variables=list(AGE="Age (in years)",
                 PTGENDER="Gender",
                 PTEDUCAT="Education (in years)",
                 APOE4="APOE4 allele status",
                 ABETA="Amyloid Beta",
                 TAU="Tau",
                 PTAU="P-Tau",
                 metaROI="FDG-PET",
                 MMSE="MMSE",
                 ADNI_EF="EF Composite Score",
                 ADNI_MEM="Memory Composite Score",
                 WMH = "Total WMH volume",
                 PLASMA_NFL = "Plasma NFL",
                 GAP_43 = "CSF GAP-43",
                 MSD_STREM2CORRECTED = "CSF sTREM2"),
  groups=list("CN","SNAP","AD"))

strata <- c(list(Total=all_data), split(all_data, all_data$ATN_group)) #list of data frames in order of display

my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=5), c("", "Mean (SD)"=sprintf("%s (&plusmn; %s)", MEAN, SD))) }
my.render.cat <- function(x) {
  c("", sapply(stats.default(x), function(y) with(y,sprintf("%d (%0.0f %%)", FREQ, PCT)))) }

table1(strata, labels, render.continuous=my.render.cont, render.categorical=my.render.cat, topclass= "Rtable1-times")


#==============================================================================#
#======================= Manuscript Aim 1 Analysis ============================#
#==============================================================================#

#SUPPLEMENTARY FIGURE 1.
#Supp Fig 1a. Unadjusted relaimpo for episodic memory
relimp.mem <- lm(all_data$ADNI_MEM ~ AGE_norm + PTGENDER + PTEDUCAT_norm + APOE4 + MH16SMOK + HYPT + DIAB + HYPL + ATN_group + WMH_log_norm + NFL_log_norm + GAP43_log_norm + MSD_STREM2_log_norm, data=all_data)
calc.relimp.mem <- calc.relimp(relimp.mem, type=c("lmg"))
calc.relimp.mem
plot(calc.relimp.mem)

#Supp Fig 1b. Unadjusted relaimpo for EF
relimp.ef <- lm(all_data$ADNI_EF ~ AGE_norm + PTGENDER + PTEDUCAT_norm + APOE4 + MH16SMOK + HYPT + DIAB + HYPL + ATN_group + WMH_log_norm + NFL_log_norm + GAP43_log_norm + MSD_STREM2_log_norm, data=all_data)
calc.relimp.ef <- calc.relimp(relimp.ef, type=c("lmg"))
calc.relimp.ef
plot(calc.relimp.ef)


#FIGURE 1.
#Fig 1a. Adjusted relaimpo for episodic memory
calc.relimp.memadj <- calc.relimp(relimp.mem, type=c("lmg"), always=c("AGE_norm","PTGENDER","PTEDUCAT_norm","APOE4","MH16SMOK","HYPT","DIAB","HYPL","ATN_group"))
calc.relimp.memadj
plot(calc.relimp.memadj)

#Fig 1b. Adjusted relaimpo for EF
calc.relimp.efadj <- calc.relimp(relimp.ef, type=c("lmg"), always=c("AGE_norm","PTGENDER","PTEDUCAT_norm","APOE4","MH16SMOK","HYPT","DIAB","HYPL","ATN_group"))
calc.relimp.efadj
plot(calc.relimp.efadj)



#==============================================================================#
#-====================== Manuscript Aim 2 Analysis ============================#
#==============================================================================#
#Simultaneous regression with relatively important biomarkers
all_data <- within(all_data, ATN_group <- relevel(factor(ATN_group), ref="CN"))
contrasts(all_data$ATN_group)

#Step 1.
mem_sim_step1 <- lm(all_data$ADNI_MEM ~ AGE_norm + PTGENDER + PTEDUCAT_norm + APOE4 + MH16SMOK + HYPT + DIAB + HYPL + ATN_group, data=all_data)        
summary(mem_sim_step1)

#Step 2.
mem_sim_step2 <- lm(all_data$ADNI_MEM ~ AGE_norm + PTGENDER + PTEDUCAT_norm + APOE4 + MH16SMOK + HYPT + DIAB + HYPL + ATN_group + NFL_log_norm + GAP43_log_norm, data=all_data)        
summary(mem_sim_step2)
vif(mem_sim_step2)

anova(mem_sim_step1, mem_sim_step2)
summary(mem_sim_step2)$r.squared - summary(mem_sim_step1)$r.squared
    


#==============================================================================#
#-====================== Manuscript Aim 3 Analysis ============================#
#==============================================================================#
#Biomarker x ATN group on episodic memory interaction analyses
all_data <- within(all_data, ATN_group <- relevel(factor(ATN_group), ref="CN"))
contrasts(all_data$ATN_group)

#Model 1. NfL
nfl_mem1 <- lm(all_data$ADNI_MEM ~ AGE_norm + PTGENDER + PTEDUCAT_norm + APOE4 + ATN_group, data=all_data)
nfl_mem2 <- lm(all_data$ADNI_MEM ~ AGE_norm + PTGENDER + PTEDUCAT_norm + APOE4 + ATN_group + NFL_log_norm, data=all_data)
nfl_mem_intx = lm(all_data$ADNI_MEM ~ AGE_norm + PTGENDER + PTEDUCAT_norm + APOE4 + NFL_log_norm*ATN_group, data=all_data)
summary(nfl_mem_intx)
anova(nfl_mem1, nfl_mem2, nfl_mem_intx)
vif(nfl_mem2)

#Figure 2a. Interaction Plot.
nfl_mem <- interact_plot(model=nfl_mem_intx, pred=NFL_log_norm, modx=ATN_group, plot.points=TRUE, x.label="Plasma NfL (pg/mL)", y.label="Composite episodic memory", legend.main="ATN Group",
                         line.thickness = 9, point.size = 13, point.alpha = 0.6)
nfl_mem <- nfl_mem + theme(plot.title = element_text(hjust = 0.5)) + theme(text = element_text(size=80)) + theme(axis.line = element_line(colour = 'black', size = 2)) + xlim(-4,3) + ylim(-3,3) + theme(legend.position = "none")


#Model 2. GAP-43
gap_mem1 <- lm(all_data$ADNI_MEM ~ AGE_norm + PTGENDER + PTEDUCAT_norm + APOE4 + ATN_group, data=all_data)
gap_mem2 <- lm(all_data$ADNI_MEM ~ AGE_norm + PTGENDER + PTEDUCAT_norm + APOE4 + ATN_group + GAP43_log_norm, data=all_data)
gap_mem_intx <- lm(all_data$ADNI_MEM ~ AGE_norm + PTGENDER + PTEDUCAT_norm + APOE4 + GAP43_log_norm*ATN_group, data=all_data)
summary(gap_mem_intx)
anova(gap_mem1, gap_mem2, gap_mem_intx)
vif(gap_mem2)

#Post-hoc
emtrends(gap_mem_intx, ~ATN_group, var="GAP43_log_norm")
emtrends(gap_mem_intx, pairwise~ATN_group, var="GAP43_log_norm")

#Figure 2b. Interaction Plot.
gap_mem <- interact_plot(model=gap_mem_intx, pred=GAP43_log_norm, modx=ATN_group, plot.points=TRUE, x.label="CSF GAP-43 (pg/mL)", y.label="Composite episodic memory", legend.main="ATN Group",
                         line.thickness = 9, point.size = 13, point.alpha = 0.6)
gap_mem <- gap_mem + theme(plot.title = element_text(hjust = 0.5)) + theme(text = element_text(size=80)) + theme(axis.line = element_line(colour = 'black', size = 2)) + xlim(-4,3) + ylim(-3,3) + theme(legend.position = "none")
    
