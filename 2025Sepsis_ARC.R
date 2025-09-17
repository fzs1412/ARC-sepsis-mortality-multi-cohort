library(ricu)
library(dplyr)
library(data.table)
library(stringr)
library(readxl)
library(openxlsx)
library(gtsummary)
setwd('d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis')
src <- c("mimic","eicu","miiv","aumc")
RICU_Cr <- load_concepts("crea", src)

dat <- read.csv("d:/OneDrive/R/drydata/RICU/drydata_crea.csv")
dat$charttime <- as.difftime(dat$charttime, units = "hours")
data <- rbind(dat_Cr, dat)
data <- data %>% filter(!is.na(icustay_id)) 
fwrite(data, "d:/OneDrive/R/0RICU/DataAnalysis/RICU_cbind_drydata_crea.csv")

dat <- fread("d:/OneDrive/R/0RICU/DataAnalysis/RICU_cbind_drydata_crea.csv")
## Sepsis with icustay >=24h
dat <- fread("d:/OneDrive/R/0RICU/DataAnalysis/ICU5data_2024x1.csv")
dat$death <- ifelse(is.na(dat$death), 0, dat$death)
dat <- dat[dat$age > 17, ]
dat <- dat[dat$sep3 == 1, ]
dat <- dat[dat$los_icu >= 1, ]
dat <- dat %>% select(-triglyceride) %>% select(-bmi)

## Remove all missing weight == dry_id_queshi
drydata <- dat[dat$source == "drydata", ]
dry_id_queshi <- drydata %>%
  filter(is.na(weight)) %>% select(icustay_id)

# RICU missing weight ids
id <- subset(dat, is.na(weight))
# Select non-drydata
id <- id[!id$source == "drydata", ]
id <- id %>% select(icustay_id) %>%
  mutate(id = 1)

dat0 <- dat[, 1:5]
library(mice)
aa <- mice(dat[, 6:52], seed = 123)
dat1 <- complete(aa, action = 3)
dat <- cbind(dat0, dat1)

## RICU missing weight ids
dat <- merge(dat, id, by = 'icustay_id', all.x = TRUE)
dat <- dat[!complete.cases(dat$id), ]
dat <- dat %>% select(-id)

fwrite(dat, "d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis/RICU_sepsis_mice.csv")
dat <- fread("d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis/RICU_sepsis_mice.csv")

# Outlier automation processing
fwrite(dat, "d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis/RICU_sepsis_mice1.csv")

dat <- fread("d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis/RICU_sepsis_mice1.csv")
# Exclude CRRT, AKI, cr 1.357 (120)
dat <- dat[dat$CKD == 0, ]
# Select charttime 0-72h with cr > 1.357 (120)
df <- fread("d:/OneDrive/R/0RICU/DataAnalysis/RICU_cbind_drydata_crea.csv")
df <- df %>%
  filter(charttime >= 0 & charttime <= 72) %>%
  filter(crea <= 1.357) %>%
  group_by(icustay_id) %>% 
  summarize(min_crea48h = min(crea))

df$icustay_id <- as.numeric(df$icustay_id)
dat$icustay_id <- as.numeric(dat$icustay_id)
# Merge
dat <- merge(dat, df, by = 'icustay_id', all.x = TRUE)

# Calculate creatinine clearance
# Define Cockcroft-Gault formula
calculate_crcl <- function(age, weight, min_crea48h, sex) {
  factor <- ifelse(sex == 0, 0.85, 1)
  crcl <- ((140 - age) * weight * factor) / (72 * min_crea48h)
  return(crcl)
}
# Calculate creatinine clearance (Cockcroft-Gault)
dat <- dat %>%
  filter(!is.na(min_crea48h)) %>%
  mutate(CrCl_CG = calculate_crcl(age, weight, min_crea48h, sex)) %>%
  mutate(ARC_CG = ifelse(CrCl_CG >= 130, 1, 0))



dat2411 <- fread("d:/OneDrive/R/0RICU/DataAnalysis/ICU5data_2024x2.csv") %>%
  select(icustay_id, adm, charlson, vaso_day1, cort_day1, rrt_day1, mv_day1, MV_time)

dat <- dat %>%
  mutate(icustay_id = as.numeric(icustay_id)) %>%
  left_join(dat2411, by = "icustay_id") %>%
  mutate(charlson = ifelse(is.na(charlson), 4, charlson))
dat <- dat %>%
  filter(rrt_day1 == 0) %>% select(-rrt_day1)

fwrite(dat, "d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis/RICU_sepsis_com.csv")
dat <- fread("d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis/RICU_sepsis_com.csv")

# Missing weight - remove dry_id_queshi
dry_id_queshi <- dry_id_queshi %>%
  mutate(dry_weight = 1)
dat <- fread("d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis/RICU_sepsis_com.csv") %>%
  mutate(icustay_id = as.numeric(icustay_id)) %>%
  left_join(dry_id_queshi, by = "icustay_id") %>%
  filter(is.na(dry_weight)) %>% select(-dry_weight)
fwrite(dat, "d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis/RICU_sepsis_com2.csv")


### CMAISE1.5v
data <- fread("d:/OneDrive/R/01sepsis/CMAISE1.5v_x0.csv")
data <- data %>%
  mutate(AKI = ifelse((sofa_cr + sofa_uo) > 0, 1, 0)) %>%
  filter(crrt == 0) %>%
  filter(renafailure == 0) %>%
  filter(Hospital_days >= 1) %>%
  filter(age > 17)
## AKI == 0
AKI0_id <- data %>%  
  group_by(PtID) %>% 
  summarize(sum_AKI = sum(AKI)) %>%
  filter(sum_AKI < 1) %>%
  select(PtID)
##
data <- data %>%
  inner_join(AKI0_id, by = "PtID")
## Exclude cr with NA's > 2
cr0_id <- data %>%
  mutate(cr0 = ifelse(is.na(cr), 1, 0)) %>%
  group_by(PtID) %>% 
  summarize(sum = sum(cr0)) %>%
  filter(sum < 2)
##
data <- data %>%
  inner_join(cr0_id, by = "PtID")
fwrite(data, "CMAISE.csv")
# Missing values
data <- fread("CMAISE.csv")
library(mice)
aa <- mice(data, seed = 123)
data <- complete(aa, action = 3)
# Outliers
data <- data %>%
  mutate(cr = cr / 88.4) %>%
  select(-sum, -crrt, -AKI)
fwrite(data, "CMAISEmice.csv")

dat <- fread("CMAISEmice.csv")
## 72h
# Define Cockcroft-Gault formula
calculate_crcl <- function(age, weight, cr, sex) {
  factor <- ifelse(sex == 0, 0.85, 1)
  crcl <- ((140 - age) * weight * factor) / (72 * cr)
  return(crcl)
}


ARC_id <- dat %>%
  mutate(CrCl = calculate_crcl(age, weight, cr, sex)) %>%
  mutate(ARC = ifelse(CrCl >= 130, 1, 0)) %>%
  filter(Days < 5) %>%
  group_by(PtID) %>% 
  summarize(sum_ARC = sum(ARC)) %>% 
  mutate(class_ARC = ifelse(sum_ARC > 0, 1, 0)) %>%
  select(-sum_ARC)


CMAISEday1 <- fread("CMAISEmice.csv") %>%
  filter(Days == 1) %>%
  mutate(CrCl = calculate_crcl(age, weight, cr, sex)) %>%
  left_join(ARC_id, by = "PtID") %>%
  rename(ARC = class_ARC)

fwrite(CMAISEday1, "CMAISEday1.csv")

CMAISEday135 <- fread("CMAISEmice.csv") %>%
  mutate(CrCl = calculate_crcl(age, weight, cr, sex)) %>%
  mutate(ARC = ifelse(CrCl >= 130, 1, 0))
fwrite(CMAISEday135, "CMAISEday135.csv")

dat <- fread("CMAISEday1.csv")


dat <- fread("CMAISEday1.csv") %>%
  select(ARC, mort)
tab1 <- twogrps(dat, gvar = "ARC")
tab1




# Remove drydata and merge CMAISEday1
# Merge vital_MaxMin
RICU_vital_MaxMin <- fread("d:/OneDrive/R/0RICU/DataAnalysis/ICU5data_2024x3.csv") %>%
  select(icustay_id, o2sat_min, o2sat_max, resp_min, hr_min, resp_max, hr_max, map_min, temp_min, map_max, temp_max)

dat <- fread("d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis/RICU_sepsis_com2.csv") %>%
  mutate(icustay_id = as.numeric(icustay_id)) %>%
  left_join(RICU_vital_MaxMin, by = "icustay_id")
# Merge CMAISE and select variables
dat <- dat %>%
  mutate(pf = (po2 / fio2) * 100) %>%
  select(source, CrCl, ARC, icustay_id, age, sex, weight, Hypertension,
         Diabete, CKD, MI, CHF, COPD,
         sofa, urine24h, vaso_day1, mv_day1,
         wbc, plt, hct,
         ph, po2, fio2, pco2, lact, ptt, bili,
         crea, bun, ca, k, na,
         death, los_hosp, MV_time,
         pf, hr_max, hr_min, map_max, map_min, o2sat_max, o2sat_min,
         resp_max, resp_min, temp_max, temp_min, charlson, los_icu) %>%
  filter(!source == "drydata") %>% 
  filter(los_hosp >= 1)
## mice vital_MaxMin
library(mice)
dat0 <- dat[, 1:5]  
library(mice)
aa <- mice(dat[, 6:48], seed = 123)
dat1 <- complete(aa, action = 3)
dat <- cbind(dat0, dat1)
## Outliers
dat0 <- dat[, 1:35]
# Outlier automation processing
data <- dat[, 36:48]
dat <- cbind(dat0, data)

# Process CMAISE variables
CMAISE <- fread("CMAISEday1.csv")
CMAISE <- CMAISE %>%
  select(CrCl, ARC, PtID, age, sex, height, weight, hyperten,
         diabete, myoinfarc, cardiofailure, copd, renafailure,
         SOFA, urine, sofa_vaso, mv,
         wbc, hct, plt,
         pha, pao, fio, paco, lac,
         aptt, bilirubin, bun, cr,
         k, na, ca,
         mort, MV_days, Hospital_days,
         pf, gcs, hrmax, hrmin, mapmax, mapmin, sapmax, sapmin, rrmax, rrmin, tmax, tmin) %>%
  mutate(bmi = 10000 * weight / (height^2), source = "CMAISE1.5v") %>%
  rename(icustay_id = PtID, Hypertension = hyperten, Diabete = diabete, CKD = renafailure,
         MI = myoinfarc, CHF = cardiofailure, COPD = copd, sofa = SOFA, urine24h = urine, vaso_day1 = sofa_vaso,
         mv_day1 = mv, ph = pha, po2 = pao, fio2 = fio, pco2 = paco, lact = lac, ptt = aptt, bili = bilirubin,
         crea = cr, death = mort, MV_time = MV_days, los_hosp = Hospital_days) %>%
  select(source, CrCl, ARC, icustay_id, age, sex, weight, Hypertension,
         Diabete, CKD, MI, CHF, COPD,
         sofa, urine24h, vaso_day1, mv_day1,
         wbc, plt, hct,
         ph, po2, fio2, pco2, lact, ptt, bili,
         crea, bun, ca, k, na,
         death, los_hosp, MV_time,
         pf, hrmax, hrmin, mapmax, mapmin, sapmax, sapmin, rrmax, rrmin, tmax, tmin) %>%
  rename(hr_max = hrmax, hr_min = hrmin, resp_max = rrmax,
         resp_min = rrmin, temp_max = tmax, temp_min = tmin,
         o2sat_max = sapmax, o2sat_min = sapmin, map_max = mapmax, map_min = mapmin) %>%
  mutate(vaso_day1 = ifelse(vaso_day1 > 0, 1, vaso_day1)) %>%
  mutate(bili = bili / 17.1, ca = ca * 4, charlson = NA, los_icu = NA)

# Merge
RICU_sepsis_com3 <- rbind(CMAISE, dat)
fwrite(RICU_sepsis_com3, "RICU_sepsis_com3.csv")

dat <- fread("d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis/RICU_sepsis_com3.csv") %>%
  select(source, age, sex, charlson,
         Hypertension, Diabete, CHF, MI, COPD, sofa, vaso_day1, 
         wbc, plt, ph, pf, lact, bili, crea, CrCl, bun,
         hr_max, hr_min, map_max, map_min,
         resp_max, resp_min, temp_max, temp_min,
         death, ARC, MV_time, los_icu, los_hosp)

# aumc, dat_CMAISE, eicu, miiv, mimic
dat_aumc <- dat[dat$source == "aumc", ]
dat_eicu <- dat[dat$source == "eicu", ]
dat_miiv <- dat[dat$source == "miiv", ]
dat_mimic <- dat[dat$source == "mimic", ]
dat_CMAISE <- dat[dat$source == "CMAISE1.5v", ] %>%
  select(source, age, sex,
         Hypertension, Diabete, CHF, MI, COPD, sofa, vaso_day1, 
         wbc, plt, ph, pf, lact, bili, crea, CrCl, bun,
         hr_max, hr_min, map_max, map_min,
         resp_max, resp_min, temp_max, temp_min,
         death, ARC, MV_time, los_hosp)
library(CBCgrps)
library(ggplot2)
dat <- dat %>% select(ARC, death)
tab1 <- twogrps(dat, gvar = "ARC")
tab1

skewvar2 <- c("bun", "los_icu", "los_hosp", "lact",
              "bili", "sofa", "crea", "CrCl")
tab1 <- multigrps(dat, gvar = "source", norm.rd = 1, cat.rd = 1, sk.rd = 1,
                  skewvar = skewvar2)

colnames(tab1) <- c("Col1", "Col2", "Col3", "Col4", "Col5", "Col6", "Col7", "Col8")
flextable::save_as_docx(flextable::flextable(as.data.frame(tab1)),
                        path = "../table/tab1_5data.docx")

# table2 ARC
## table for dat_mimic, dat_miiv, dat_eicu, dat_aumc
skewvar2 <- c("bun", "los_icu", "los_hosp", "lact",
              "bili", "sofa", "crea", "CrCl")
tab1 <- twogrps(dat_miiv, gvar = "ARC", norm.rd = 1, cat.rd = 1, sk.rd = 1,
                minfactorlevels = 5, skewvar = skewvar2)
colnames(tab1$Table) <- c("Col1", "Col2", "Col3", "Col4", "Col5")
flextable::save_as_docx(flextable::flextable(as.data.frame(tab1$Table)),
                        path = "../table/tab2_miiv_ARC.docx")

tab1 <- twogrps(dat_mimic, gvar = "ARC", norm.rd = 1, cat.rd = 1, sk.rd = 1,
                minfactorlevels = 5, skewvar = skewvar2)
colnames(tab1$Table) <- c("Col1", "Col2", "Col3", "Col4", "Col5")
flextable::save_as_docx(flextable::flextable(as.data.frame(tab1$Table)),
                        path = "../table/tabS_mimic_ARC.docx")

tab1 <- twogrps(dat_eicu, gvar = "ARC", norm.rd = 1, cat.rd = 1, sk.rd = 1,
                minfactorlevels = 5, skewvar = skewvar2)
colnames(tab1$Table) <- c("Col1", "Col2", "Col3", "Col4", "Col5")
flextable::save_as_docx(flextable::flextable(as.data.frame(tab1$Table)),
                        path = "../table/tabS_eicu_ARC.docx")
tab1 <- twogrps(dat_aumc, gvar = "ARC", norm.rd = 1, cat.rd = 1, sk.rd = 1,
                minfactorlevels = 5, skewvar = skewvar2)
colnames(tab1$Table) <- c("Col1", "Col2", "Col3", "Col4", "Col5")
flextable::save_as_docx(flextable::flextable(as.data.frame(tab1$Table)),
                        path = "../table/tabS_aumc_ARC.docx")
## table for dat_CMAISE (missing los_icu)
skewvar3 <- c("bun", "los_hosp", "lact",
              "bili", "sofa", "crea", "CrCl")
tab1 <- twogrps(dat_CMAISE, gvar = "ARC", norm.rd = 1, cat.rd = 1, sk.rd = 1,
                minfactorlevels = 5, skewvar = skewvar3)
colnames(tab1$Table) <- c("Col1", "Col2", "Col3", "Col4", "Col5")
flextable::save_as_docx(flextable::flextable(as.data.frame(tab1$Table)),
                        path = "../table/tabS_CMAISE_ARC.docx")

#
tab1 <- twogrps(dat, gvar = "death", norm.rd = 1, cat.rd = 1, sk.rd = 1,
                minfactorlevels = 5, skewvar = skewvar2)
colnames(tab1$Table) <- c("Col1", "Col2", "Col3", "Col4", "Col5")
flextable::save_as_docx(flextable::flextable(as.data.frame(tab1$Table)),
                        path = "../table/tabs_5data_mort.docx")





# Sensitivity analysis
dat <- fread("d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis/RICU_sepsis_com3.csv") %>%
  select(source, age, sex, charlson,
         Hypertension, Diabete, CHF, MI, COPD, sofa, vaso_day1, 
         wbc, plt, ph, pf, lact, bili, crea, CrCl, bun,
         hr_max, hr_min, map_max, map_min,
         resp_max, resp_min, temp_max, temp_min,
         death, ARC, MV_time, los_icu, los_hosp)

# aumc, dat_CMAISE, eicu, miiv, mimic
dat_aumc <- dat[dat$source == "aumc", ]
dat_eicu <- dat[dat$source == "eicu", ]
dat_miiv <- dat[dat$source == "miiv", ]
dat_mimic <- dat[dat$source == "mimic", ]
dat_CMAISE <- dat[dat$source == "CMAISE1.5v", ] %>%
  select(source, age, sex,
         Hypertension, Diabete, CHF, MI, COPD, sofa, vaso_day1, 
         wbc, plt, ph, pf, lact, bili, crea, CrCl, bun,
         hr_max, hr_min, map_max, map_min,
         resp_max, resp_min, temp_max, temp_min,
         death, ARC, MV_time, los_hosp)
library(CBCgrps)
library(ggplot2)



