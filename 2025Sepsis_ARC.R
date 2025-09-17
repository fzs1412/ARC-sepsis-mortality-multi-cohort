setwd('d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis')
library(ricu)
library(dplyr)
library(data.table)
library(stringr)
library(readxl)
library(openxlsx)
library(gtsummary)
#RICU
#排除RRT
dat=fread("d:/OneDrive/R/0RICU/DataAnalysis/data6_sepsis202504.csv") 


#sepsis  icustay >24h
dat=fread("d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis/RICU_sepsis_mice.csv")
# 异常值自动化后mice1
dat=fread("d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis/RICU_sepsis_mice1.csv")
##com
dat=fread("d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis/RICU_sepsis_com.csv")


#crea=crea / 88.4   ) # µmol/L to mg/dL

#合并RICU_crea  drydata_crea
setwd('d:/OneDrive/R/0RICU/DataAnalysis')
library(ricu)
library(dplyr)
library(data.table)
library(stringr)
src <- c("mimic","eicu","miiv","aumc")
RICU_Cr <- load_concepts("crea", src)

dat=read.csv("d:/OneDrive/R/drydata/RICU/drydata_crea.csv")
dat$charttime <- as.difftime(dat$charttime, units = "hours")
data=rbind(dat_Cr,dat)
data=data%>%filter(!is.na(icustay_id)) 
fwrite(data,"d:/OneDrive/R/0RICU/DataAnalysis/RICU_cbind_drydata_crea.csv")

dat=fread("d:/OneDrive/R/0RICU/DataAnalysis/RICU_cbind_drydata_crea.csv")
##sepsis  icustay 24h
dat=fread("d:/OneDrive/R/0RICU/DataAnalysis/ICU5data_2024x1.csv") #361922   完整的5个RICU数据  
dat$death=ifelse( is.na(dat$death), 0,dat$death )
dat=dat[dat$age>17,]   #354016
dat=dat[dat$sep3==1,]  #  80257
dat=dat[dat$los_icu>=1,] 
dat=dat%>%  select(-triglyceride)%>%  select(-bmi)

##删除所有缺失weight ==    dry_id_queshi
drydata=dat[dat$source=="drydata",]
dry_id_queshi=drydata%>%
  filter(is.na(weight))%>%  select(icustay_id)

#RICU  weight的缺失值id
id=subset(dat, is.na(weight)) 
#选择非drydata     4087
id=id[!id$source=="drydata",]   
id=id%>%  select( icustay_id  )%>%
  mutate( id = 1)

dat0<- dat[,1:5]  
library(mice)
aa <- mice(dat[,6:52], seed=123)  #??
dat1<-complete(aa, action=3)
dat=cbind(dat0,dat1)

##RICU  weight的缺失值id
dat=merge(dat,id,by = 'icustay_id',all.x = TRUE)
dat=dat[!complete.cases(dat$id),]
dat=dat%>%  select(-id)

fwrite(dat,"d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis/RICU_sepsis_mice.csv")
dat=fread("d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis/RICU_sepsis_mice.csv")

# 异常值自动化处理  689-727   714
fwrite(dat,"d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis/RICU_sepsis_mice1.csv")

dat=fread("d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis/RICU_sepsis_mice1.csv")
#排除  CRRT   AKI  cr  1.357 （120）
dat=dat[dat$CKD==0,]  #65037-56467===8570
#选择charttime 0-72h  cr > 1.357 （120）
df=fread("d:/OneDrive/R/0RICU/DataAnalysis/RICU_cbind_drydata_crea.csv")
df=df%>%
  filter( charttime>=0&charttime<=72 )%>%   #48
  filter(crea<=1.357)%>%
  group_by(icustay_id) %>% 
  summarize(min_crea48h = min(crea)) 

df$icustay_id <- as.numeric(df$icustay_id)
dat$icustay_id <- as.numeric(dat$icustay_id)
#合并
dat=merge(dat,df,by = 'icustay_id',all.x = TRUE)

#计算肌酐清除率
# 定义Cockcroft-Gault公式
calculate_crcl <- function(age, weight, min_crea48h, sex) {
  factor <- ifelse(sex == 0, 0.85, 1)                       #female    
  crcl <- ((140 - age) * weight * factor) / (72 * min_crea48h)
  return(crcl)
}
# 计算肌酐清除率（Cockcroft-Gault）
dat <- dat %>%
  filter(!is.na(min_crea48h)) %>%
  mutate(CrCl_CG = calculate_crcl(age, weight, min_crea48h, sex))%>%
  mutate(ARC_CG = ifelse(CrCl_CG >= 130, 1, 0))


# 定义计算CKD-EPI肌酐清除率的函数
CKD_EPI_eGFR <- function(min_crea48h, age, sex) {
  # 使用 ifelse 进行向量化判断，以避免 "length > 1" 错误
  result <- ifelse(sex == 1,
                   ifelse(min_crea48h <= 0.9,
                          141 * (min_crea48h / 0.9)^(-0.411) * 0.993^(age),
                          141 * (min_crea48h / 0.9)^(-1.209) * 0.993^(age)),
                   ifelse(min_crea48h <= 0.7,
                          144 * (min_crea48h / 0.7)^(-0.329) * 0.993^(age),
                          144 * (min_crea48h / 0.7)^(-1.209) * 0.993^(age)))
  return(result)
}

dat <- dat %>%
  mutate(CrCl_ckdepi = CKD_EPI_eGFR(min_crea48h, age, sex))%>%
  mutate(ARC_ckdepi = ifelse(CrCl_ckdepi >= 130, 1, 0))

# 定义计算MDRD简化公式的函数
# 输入: 血浆肌酐 (plasma_creatinine, 单位: mg/dL), 年龄 (age), 性别 (sex)
MDRD_175eGFR <- function(plasma_creatinine, age, sex) {
  # 根据数据中的性别编码（0为女性，1为男性）调整公式
  result <- 175 * (plasma_creatinine - 0.0113)^-1.154 * age^-0.203 *
    ifelse(sex == 0, 0.742, 1)
  return(result)
}
# 计算 eGFR 并判断是否为 ARC
dat <- dat %>%
  mutate(CrCl_MDRD = MDRD_175eGFR(min_crea48h, age, sex)) %>%
  mutate(ARC_MDRD = ifelse(CrCl_MDRD >= 130, 1, 0))


#处理异常值CrCl,min_crea48h,sofa
dat$sofa=ifelse( dat$sofa<2,2,dat$sofa)

q5<-quantile(dat$min_crea48h, 0.01)
dat[dat$ min_crea48h <q5,]$ min_crea48h <-q5

q99<-quantile(dat$CrCl, 0.95)
dat[dat$CrCl>q99,]$CrCl<-q99


dat2411=fread("d:/OneDrive/R/0RICU/DataAnalysis/ICU5data_2024x2.csv")%>%
  select(icustay_id,  adm,charlson,vaso_day1,cort_day1,rrt_day1,mv_day1, MV_time   )

dat=dat%>%
  mutate( icustay_id=as.numeric( icustay_id))%>%
  left_join(dat2411,by="icustay_id")%>%
  mutate(charlson= ifelse(is.na(charlson), 4, charlson)) 
dat=dat%>%
  filter( rrt_day1==0 )%>%select(-rrt_day1)

fwrite(dat,"d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis/RICU_sepsis_com.csv"    )
dat=fread("d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis/RICU_sepsis_com.csv")

#缺失weight   删除 dry_id_queshi
dry_id_queshi=dry_id_queshi%>%
  mutate(dry_weight=1)
dat=fread("d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis/RICU_sepsis_com.csv")%>%
  mutate(icustay_id = as.numeric(icustay_id))%>%
  left_join(dry_id_queshi,by = "icustay_id") %>%
  filter(is.na(dry_weight))%>%select(-dry_weight)
fwrite(dat,"d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis/RICU_sepsis_com2.csv"    )


###CMAISE1.5v
##AKI sofa
data=fread("d:/OneDrive/R/01sepsis/CMAISE1.5v_x0.csv")
data=data%>%
  mutate(AKI=ifelse(  (sofa_cr+sofa_uo)>0,1,0 ) )%>%
  filter( crrt==0 )%>%
  filter( renafailure==0 )%>%
  filter(Hospital_days>=1)%>%
  filter(age>17)
##AKI==0    369
AKI0_id=data%>%  
  group_by(PtID) %>% 
  summarize(sum_AKI = sum(AKI))%>%
  filter(sum_AKI<1)%>%
  select(  PtID )
##
data=data%>%
  inner_join(AKI0_id,by = "PtID")
##排除cr   NA's >2  :170 ??
cr0_id=data%>%
  mutate(cr0=ifelse(is.na(cr) ,1,0) )%>%
  group_by(PtID) %>% 
  summarize(sum = sum(cr0))%>%
  filter( sum<2)
##
data=data%>%
  inner_join(cr0_id,by = "PtID")
fwrite( data,"CMAISE.csv"         )
#  缺失值
data=fread(  "CMAISE.csv"         )
library(mice)
aa <- mice(data, seed=123)
data<-complete(aa, action=3)
#异常值  # 异常值自动化处理731
data=data%>%
  mutate(cr=cr/ 88.4)%>%
  select(-sum,-crrt,-AKI    )
fwrite( data,"CMAISEmice.csv"         )

dat=fread(  "CMAISEmice.csv"         )
##72h
# 定义Cockcroft-Gault公式
calculate_crcl <- function(age, weight, cr, sex) {
  factor <- ifelse(sex == 0, 0.85, 1)
  crcl <- ((140 - age) * weight * factor) / (72 * cr)
  return(crcl)
}
# 计算肌酐清除率（Cockcroft-Gault）

# 定义计算MDRD简化公式的函数
# 输入: 血浆肌酐 (plasma_creatinine, 单位: mg/dL), 年龄 (age), 性别 (sex)
MDRD_175eGFR <- function(plasma_creatinine, age, sex) {
  # 根据数据中的性别编码（0为女性，1为男性）调整公式
  result <- 175 * (plasma_creatinine - 0.0113)^-1.154 * age^-0.203 *
    ifelse(sex == 0, 0.742, 1)
  return(result)
}
# 计算 eGFR 并判断是否为 ARC
ARC_MDRD_id <- dat %>%
  mutate(CrCl_MDRD = MDRD_175eGFR(cr, age, sex)) %>%
  mutate(ARC_MDRD = ifelse(CrCl_MDRD >= 130, 1, 0))%>%
  filter( Days<5)%>%
  group_by(PtID) %>% 
  summarize(sum_ARC = sum(ARC_MDRD ))%>% 
  mutate(class_ARC_MDRD = ifelse(sum_ARC >0, 1, 0))%>%
  select(-sum_ARC)

ARC_id <- dat %>%
  mutate(CrCl = calculate_crcl(age, weight, cr, sex))%>%
  mutate(ARC = ifelse(CrCl >= 130, 1, 0))%>%
  filter( Days<5)%>%
  group_by(PtID) %>% 
  summarize(sum_ARC = sum(ARC ))%>% 
  mutate(class_ARC = ifelse(sum_ARC >0, 1, 0))%>%
  select(-sum_ARC)


CMAISEday1=fread(  "CMAISEmice.csv"         )%>%
  filter( Days==1    )%>%
   mutate(CrCl = calculate_crcl(age, weight, cr, sex))%>%
  left_join(ARC_id,by="PtID")%>%
  left_join(ARC_MDRD_id,by="PtID")%>%
  rename(ARC=class_ARC,ARC_MDRD=class_ARC_MDRD )

fwrite(CMAISEday1,"CMAISEday1.csv")

CMAISEday135=fread(  "CMAISEmice.csv"         )%>%
  mutate(CrCl = calculate_crcl(age, weight, cr, sex))%>%
  mutate(ARC = ifelse(CrCl >= 130, 1, 0))
fwrite(CMAISEday135,"CMAISEday135.csv")

dat=fread(  "CMAISEday1.csv"         )


dat=fread(  "CMAISEday1.csv"    ) %>%
  select(ARC,mort)
tab1<-twogrps(dat,gvar = "ARC")
tab1

dat=fread(  "CMAISEday1.csv")
dat=dat[,5:75]
mod1<-glm(mort~.,dat,
          family = "binomial")

logit.step<-step(mod1,direction = c("forward"))   
summary(logit.step)

mod1<-glm(mort~ARC,dat,
          family = "binomial")      
summary(mod1)


mod3<-glm(mort~ARC+age+wbc+SOFA+lac+pha  ,dat,
          family = "binomial")      
tbl_regression(mod3, exponentiate = TRUE)
summary(mod3)

cox_mode3 <- coxph(Surv(Hospital_days, mort) ~ ARC+age+wbc+SOFA+lac+pha    
                   ,dat)
tbl_regression(cox_mode3, exponentiate = TRUE)

mod3<-glm(mort~age+wbc+sofa+lymph+ARC+lact+ph+pt+temp_fir+Diabete#+urine24h
          ,dat_drydata,
          family = "binomial")      
tbl_regression(mod3, exponentiate = TRUE)


#删除drydata  合并CMAISEday1
#合并vital_MaxMin
RICU_vital_MaxMin=fread("d:/OneDrive/R/0RICU/DataAnalysis/ICU5data_2024x3.csv")%>%
  select(icustay_id,o2sat_min,o2sat_max,resp_min,hr_min,resp_max,hr_max,map_min,temp_min,map_max,temp_max)

dat=fread("d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis/RICU_sepsis_com2.csv")%>%
  mutate( icustay_id=as.numeric( icustay_id))%>%
  left_join(RICU_vital_MaxMin,by="icustay_id")
#合并CMAISE  选择变量
dat=dat%>%
  mutate(pf=(po2/fio2)*100          )%>%
  select( source, CrCl,ARC,icustay_id, age, sex,weight, Hypertension,
          Diabete, CKD, MI, CHF, COPD,
          sofa, urine24h, vaso_day1, mv_day1,
          wbc, plt, hct,
          ph, po2, fio2, pco2, lact, ptt, bili,
          crea, bun, ca, k, na,
          death, los_hosp, MV_time,
          pf,hr_max,hr_min,map_max,map_min,o2sat_max,o2sat_min,
          resp_max,resp_min,temp_max,temp_min,charlson,los_icu  )%>%#charlson
  filter(!source=="drydata"    )%>% 
  filter(los_hosp>=1)
##mice  vital_MaxMin
library(mice)
dat0<- dat[,1:5]  
library(mice)
aa <- mice(dat[,6:48], seed=123)  #??
dat1<-complete(aa, action=3)
dat=cbind(dat0,dat1)
##异常值
dat0<- dat[,1:35]
# 异常值自动化处理
data =dat[,36:48]
dat=cbind(dat0,data)

#处理CMAISE 变量
CMAISE=fread(  "CMAISEday1.csv"         )
CMAISE=CMAISE%>%
  select( CrCl,ARC,PtID,age,sex,height,weight,hyperten,              #CrCl,ARC,
          diabete,myoinfarc,cardiofailure,copd,renafailure,
          SOFA,urine,sofa_vaso,mv,
          wbc,hct,plt,
          pha,pao,fio,paco,lac,
          aptt,bilirubin,bun,cr,
          k,na,ca,
          mort,MV_days,Hospital_days,
          pf,gcs,hrmax,hrmin,mapmax,mapmin,sapmax,sapmin,rrmax,rrmin,tmax,tmin)%>%   #ricu缺
  mutate(bmi =10000* weight/(height^2), source="CMAISE1.5v")%>%  #
  rename(icustay_id=PtID, Hypertension=hyperten, Diabete=diabete, CKD=renafailure,
         MI=myoinfarc, CHF=cardiofailure, COPD=copd,sofa=SOFA,urine24h=urine,vaso_day1=sofa_vaso,
         mv_day1=mv,ph=pha,po2=pao,fio2=fio,pco2=paco,lact=lac,ptt=aptt,bili=bilirubin,
         crea=cr, death=mort, MV_time=MV_days,los_hosp=Hospital_days)%>%
  select( source, CrCl,ARC,icustay_id, age, sex,weight,  Hypertension,      #CrCl,ARC,
          Diabete, CKD, MI, CHF, COPD,
          sofa, urine24h, vaso_day1, mv_day1,
          wbc, plt, hct,
          ph, po2, fio2, pco2, lact, ptt, bili,
          crea, bun, ca, k, na,
          death, los_hosp, MV_time,
          pf,hrmax,hrmin,mapmax,mapmin,sapmax,sapmin,rrmax,rrmin,tmax,tmin)%>%
  rename(hr_max = hrmax, hr_min = hrmin, resp_max = rrmax,
         resp_min = rrmin,   temp_max = tmax,  temp_min = tmin,
         o2sat_max = sapmax,  o2sat_min = sapmin, map_max = mapmax, map_min = mapmin)%>%
  mutate(vaso_day1=ifelse(vaso_day1>0,1,vaso_day1))%>%
  mutate(bili=bili/17.1,ca=ca*4,charlson=NA,los_icu=NA)   #bili 17.1  crea 88.4 ca  4 钙离子 1 mmol/L=4 mg/dL crea=crea/88.42

#合并
RICU_sepsis_com3=rbind(CMAISE,dat)
fwrite(RICU_sepsis_com3,"RICU_sepsis_com3.csv"          )

dat=fread("d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis/RICU_sepsis_com3.csv")%>%
  select(source,age, sex, charlson,
         Hypertension,Diabete,CHF,MI,COPD,sofa,vaso_day1, 
         wbc,plt,ph,pf,lact,bili,crea,CrCl,bun,
         hr_max,hr_min,map_max,map_min,
         resp_max,resp_min,temp_max,temp_min,
         death, ARC,MV_time,los_icu,los_hosp)  #,ARC_MDRD

#aumc dat_CMAISE    eicu    miiv   mimic
dat_aumc=dat[dat$source=="aumc",]
dat_eicu=dat[dat$source=="eicu",]
dat_miiv=dat[dat$source=="miiv",]
dat_mimic=dat[dat$source=="mimic",]
dat_CMAISE=dat[dat$source=="CMAISE1.5v",]%>%   #缺失charlson los_icu
  select(source,age, sex,
         Hypertension,Diabete,CHF,MI,COPD,sofa,vaso_day1, 
         wbc,plt,ph,pf,lact,bili,crea,CrCl,bun,
         hr_max,hr_min,map_max,map_min,
         resp_max,resp_min,temp_max,temp_min,
         death, ARC,MV_time,los_hosp)
library(CBCgrps)
library(ggplot2)
dat <- dat %>%  select(ARC,death)
tab1<-twogrps(dat,gvar = "ARC")  #ARC
tab1
##table   dat_CMAISE  dat_mimic  dat_miiv  dat_eicu  dat_aumc

#library(DataExplorer)
#plot_qq(dat, sampled_rows = 1000L)  
skewvar2 = c("bun","los_icu","los_hosp","lact",
             "bili","sofa","crea","CrCl")
tab1<-multigrps(dat,gvar = "source",norm.rd = 1,cat.rd = 1, sk.rd=1,
                skewvar = skewvar2)

colnames(tab1) <- c("Col1", "Col2", "Col3", "Col4", "Col5", "Col6", "Col7", "Col8")
flextable::save_as_docx(flextable::flextable(as.data.frame(tab1)),
                        path = "../table/tab1_5data.docx")

#table2  ARC
##table     dat_mimic  dat_miiv  dat_eicu  dat_aumc
skewvar2 = c("bun","los_icu","los_hosp","lact",
             "bili","sofa","crea","CrCl")
tab1<-twogrps(dat_miiv,gvar = "ARC",norm.rd = 1,cat.rd = 1, sk.rd=1,
              minfactorlevels = 5,skewvar = skewvar2)
colnames(tab1$Table) <- c("Col1", "Col2", "Col3", "Col4", "Col5")
flextable::save_as_docx(flextable::flextable(as.data.frame(tab1$Table)),
                        path = "../table/tab2_miiv_ARC.docx")

tab1<-twogrps(dat_mimic,gvar = "ARC",norm.rd = 1,cat.rd = 1, sk.rd=1,
              minfactorlevels = 5,skewvar = skewvar2)
colnames(tab1$Table) <- c("Col1", "Col2", "Col3", "Col4", "Col5")
flextable::save_as_docx(flextable::flextable(as.data.frame(tab1$Table)),
                        path = "../table/tabS_mimic_ARC.docx")

tab1<-twogrps(dat_eicu,gvar = "ARC",norm.rd = 1,cat.rd = 1, sk.rd=1,
              minfactorlevels = 5,skewvar = skewvar2)
colnames(tab1$Table) <- c("Col1", "Col2", "Col3", "Col4", "Col5")
flextable::save_as_docx(flextable::flextable(as.data.frame(tab1$Table)),
                        path = "../table/tabS_eicu_ARC.docx")
tab1<-twogrps(dat_aumc,gvar = "ARC",norm.rd = 1,cat.rd = 1, sk.rd=1,
              minfactorlevels = 5,skewvar = skewvar2)
colnames(tab1$Table) <- c("Col1", "Col2", "Col3", "Col4", "Col5")
flextable::save_as_docx(flextable::flextable(as.data.frame(tab1$Table)),
                        path = "../table/tabS_aumc_ARC.docx")
##table   dat_CMAISE  缺失los_icu
skewvar3 = c("bun","los_hosp","lact",
             "bili","sofa","crea","CrCl")
tab1<-twogrps(dat_CMAISE,gvar = "ARC",norm.rd = 1,cat.rd = 1, sk.rd=1,
              minfactorlevels = 5,  skewvar = skewvar3)
colnames(tab1$Table) <- c("Col1", "Col2", "Col3", "Col4", "Col5")
flextable::save_as_docx(flextable::flextable(as.data.frame(tab1$Table)),
                        path = "../table/tabS_CMAISE_ARC.docx")

#
tab1<-twogrps(dat,gvar = "death",norm.rd = 1,cat.rd = 1, sk.rd=1,
              minfactorlevels = 5,skewvar = skewvar2)
colnames(tab1$Table) <- c("Col1", "Col2", "Col3", "Col4", "Col5")
flextable::save_as_docx(flextable::flextable(as.data.frame(tab1$Table)),
                        path = "../table/tabs_5data_mort.docx")


##logit
library(gtsummary)
dat=dat[,2:49]
mod1<-glm(death~.,dat,
          family = "binomial")
logit.step<-step(mod1,direction = c("forward"))   
summary(logit.step)
mod2<-glm(death~age+sex+charlson+vaso_day1+sofa+wbc+lact+ph+bili+crea+pf+
            hr_max+map_min+temp_min     ,dat,
          family = "binomial")
summary(mod2)

mod3<-glm(death~ARC+age+sex+charlson+vaso_day1+sofa+wbc+lact+ph+bili+crea+pf+
            hr_max+temp_min
          ,dat,
          family = "binomial")      
table1 <-tbl_regression(mod3, exponentiate = TRUE)
table1 


#dat_aumc
mod3<-glm(death~age+wbc+sofa+urine24h+lymph+pt+lact+ph+temp_fir+ARC
          ,dat_aumc,
          family = "binomial")      
tbl_regression(mod3, exponentiate = TRUE)
#dat_eicu
mod3<-glm(death~age+wbc+Diabete+sofa+urine24h+lymph+pt+lact+ph+temp_fir+ARC
          ,dat_eicu,
          family = "binomial")      
tbl_regression(mod3, exponentiate = TRUE)
#dat_miiv
mod3<-glm(death~age+wbc+Diabete+sofa+lymph+pt+lact+ph+temp_fir+ARC  #+urine24h
          ,dat_miiv,
          family = "binomial")      
tbl_regression(mod3, exponentiate = TRUE)
#dat_mimic
mod3<-glm(death~age+wbc+Diabete+sofa+lymph+pt+lact+ph+temp_fir+ARC  #+urine24h
          ,dat_mimic,
          family = "binomial")      
tbl_regression(mod3, exponentiate = TRUE)
#
 
#### 倾向性评分匹配


m.out=matchit(ARC~age+sex+vaso_day1+sofa+wbc+lact+ph+bili+crea+pf+
                hr_max+temp_min+charlson,              ,
              data=dat_mimic,method="nearest", ratio=1)
summary(m.out)
PSM_mimic<-match.data(m.out)

m.out=matchit(ARC~age+sex+vaso_day1+sofa+wbc+lact+ph+bili+crea+pf+
                hr_max+temp_min+charlson,              ,
              data=dat_miiv,method="nearest", ratio=1)
summary(m.out)
PSM_miiv<-match.data(m.out)
##后续meta_analysis_psm


# 加载必要的库
library(survival)
library(survminer)
# 绘制Kaplan-Meier生存曲线
dat=fread("d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis/RICU_sepsis_com3.csv")
dat_aumc=dat[dat$source=="aumc",]
dat_eicu=dat[dat$source=="eicu",]
dat_miiv=dat[dat$source=="miiv",]
dat_mimic=dat[dat$source=="mimic",]
dat_CMAISE <- dat[dat$source == "CMAISE1.5v", ]

##PSM_CMAISE
fit <- survfit(Surv(los_hosp, death) ~ ARC, data = PSM_CMAISE)   #data = dat
ggsurvplot(fit,
           xlab='days',    pval = T,   pval.size=3,       pval.coord = c(0, 0.6),     
           legend.labs=c('No_ARC','ARC'),
           surv.median.line = "hv",    ggtheme = theme_bw(), 
           palette = c("#FFA500", "#0000FF"),
           risk.table = TRUE, risk.table.col = "strata",   ylim = c(0.5, 1))

#mimic
fit <- survfit(Surv(los_hosp, death) ~ ARC, data = PSM_mimic)   #data = dat
ggsurvplot(fit,
           xlab='days',    pval = T,   pval.size=3,       pval.coord = c(0, 0.6),     
           legend.labs=c('No_ARC','ARC'),
           surv.median.line = "hv",    ggtheme = theme_bw(), 
           palette = c("#FFA500", "#0000FF"),
           risk.table = TRUE, risk.table.col = "strata",   ylim = c(0.5, 1))
#miiv
fit <- survfit(Surv(los_hosp, death) ~ ARC, data = PSM_miiv)   #data = dat
ggsurvplot(fit,
           xlab='days',    pval = T,   pval.size=3,       pval.coord = c(0, 0.6),     
           legend.labs=c('No_ARC','ARC'),
           surv.median.line = "hv",    ggtheme = theme_bw(), 
           palette = c("#FFA500", "#0000FF"),
           risk.table = TRUE, risk.table.col = "strata",   ylim = c(0.5, 1))
#eicu
fit <- survfit(Surv(los_hosp, death) ~ ARC, data = PSM_eicu)   #data = dat
ggsurvplot(fit,
           xlab='days',    pval = T,   pval.size=3,       pval.coord = c(0, 0.6),     
           legend.labs=c('No_ARC','ARC'),
           surv.median.line = "hv",    ggtheme = theme_bw(), 
           palette = c("#FFA500", "#0000FF"),
           risk.table = TRUE, risk.table.col = "strata",   ylim = c(0.5, 1))
#aumc
fit <- survfit(Surv(los_hosp, death) ~ ARC, data = PSM_aumc)   #data = dat
ggsurvplot(fit,
           xlab='days',    pval = T,   pval.size=3,       pval.coord = c(0, 0.6),     
           legend.labs=c('No_ARC','ARC'),
           surv.median.line = "hv",    ggtheme = theme_bw(), 
           palette = c("#FFA500", "#0000FF"),
           risk.table = TRUE, risk.table.col = "strata",   ylim = c(0.5, 1))



#预测ARC模型
dat=fread("d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis/RICU_sepsis_com.csv")
dat_aumc=dat[dat$source=="aumc",]
dat_eicu=dat[dat$source=="eicu",]
dat_miiv=dat[dat$source=="miiv",]
dat_mimic=dat[dat$source=="mimic",]
dat_drydata=dat[dat$source=="drydata",]

library(gtsummary)
dat=dat[,c(2:48,55)]

mod1<-glm(ARC~.,dat,
          family = "binomial")

logit.step<-step(mod1,direction = c("forward"))   
summary(logit.step)

mod2<-glm(ARC~age+sex+Hypertension+Diabete+MI +COPD+sofa+urine24h+plt+hct+
            lymph+pt+po2+ast+bicar+ca+crea+mg+na+glu+temp_fir            
          ,dat,
          family = "binomial")
summary(mod2)

mod3<-glm(ARC~age+sex+Hypertension+Diabete+MI +COPD+sofa+urine24h+plt+hct+
            lymph+pt+po2+ast+bicar+ca+mg+na+glu+temp_fir            
          ,dat,
          family = "binomial")      
tbl_regression(mod3, exponentiate = TRUE)

#dat_drydata
mod3<-glm(ARC~age+sex+sofa+urine24h+
            po2+ast+temp_fir            
          ,dat_drydata,
          family = "binomial")     
tbl_regression(mod3, exponentiate = TRUE)

#mimic
mod3<-glm(ARC~age+sex+Hypertension+Diabete +COPD+sofa+urine24h+hct+
            lymph+ast+bicar+ca+mg+na+glu            
          ,dat_mimic,
          family = "binomial")      
tbl_regression(mod3, exponentiate = TRUE)
#miiv
mod3<-glm(ARC~age+sex+Hypertension+Diabete+MI +COPD+sofa+urine24h+hct+
            lymph+ast+bicar+ca+mg+na+glu         
          ,dat_miiv,
          family = "binomial")      
tbl_regression(mod3, exponentiate = TRUE)
#eicu
mod3<-glm(ARC~age+sex+Diabete +COPD+sofa+urine24h+hct+
            lymph+ast+bicar+na+glu+temp_fir             
          ,dat_eicu,
          family = "binomial")      
tbl_regression(mod3, exponentiate = TRUE)
#aumc
mod3<-glm(ARC~age+sex+sofa+urine24h+plt+hct+
            bicar+mg+na+glu           
          ,dat_aumc,
          family = "binomial")      
tbl_regression(mod3, exponentiate = TRUE)

#敏感性分析
dat=fread("d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis/RICU_sepsis_com3.csv")%>%
  select(source,age, sex, charlson,
         Hypertension,Diabete,CHF,MI,COPD,sofa,vaso_day1, 
         wbc,plt,ph,pf,lact,bili,crea,CrCl,bun,
         hr_max,hr_min,map_max,map_min,
         resp_max,resp_min,temp_max,temp_min,
         death, ARC,MV_time,los_icu,los_hosp,ARC_MDRD)  #,ARC_MDRD

#aumc dat_CMAISE    eicu    miiv   mimic
dat_aumc=dat[dat$source=="aumc",]
dat_eicu=dat[dat$source=="eicu",]
dat_miiv=dat[dat$source=="miiv",]
dat_mimic=dat[dat$source=="mimic",]
dat_CMAISE=dat[dat$source=="CMAISE1.5v",]%>%   #缺失charlson los_icu
  select(source,age, sex,
         Hypertension,Diabete,CHF,MI,COPD,sofa,vaso_day1, 
         wbc,plt,ph,pf,lact,bili,crea,CrCl,bun,
         hr_max,hr_min,map_max,map_min,
         resp_max,resp_min,temp_max,temp_min,
         death, ARC,MV_time,los_hosp)
library(CBCgrps)
library(ggplot2)


mod3<-glm(death~ARC_MDRD+age+sex+charlson+vaso_day1+sofa+#+wbc+lact+ph+bili+crea+pf+
            hr_max+temp_min,dat_miiv,  family = "binomial")      
tbl_regression(mod3, exponentiate = TRUE)











##  Cockcroft-Gault）计算错误而重新写代码   加入MDRD简化公式的函数
id <- dat %>%
  select(  icustay_id,CrCl_CG,ARC_CG,ARC_MDRD   )
dat <- fread("RICU_sepsis_com3.csv") %>%
  left_join(mutate(id, icustay_id = as.character(icustay_id)), by = 'icustay_id')
dat=kk
CMAISE=fread(  "CMAISEday1.csv"    )%>%
  select(  PtID , ARC , CrCl,ARC_MDRD     )%>%
  rename( icustay_id=PtID, ARC1=ARC ,ARC_MDRD1=ARC_MDRD, CrCl1= CrCl        )
dat <- dat %>%
  left_join(CMAISE, by = 'icustay_id')
dat$ARC_CG=ifelse(is.na(dat$ARC_CG) ,dat$ARC1, dat$ARC_CG)
dat$ARC_MDRD=ifelse(is.na(dat$ARC_MDRD) ,dat$ARC_MDRD1, dat$ARC_MDRD)
dat$ARC_MDRD=ifelse(is.na(dat$ARC_MDRD) ,dat$ARC_CG, dat$ARC_MDRD)
dat$CrCl_CG=ifelse(is.na(dat$CrCl_CG) ,dat$CrCl1, dat$CrCl_CG)

q99<-quantile(dat$CrCl_CG, 0.95)
dat[dat$CrCl_CG>q99,]$CrCl_CG<-q99
dat <- dat %>%
  select(-ARC,-ARC1,-CrCl,-CrCl1,-ARC_MDRD1  ) %>%
  rename( CrCl=CrCl_CG,     ARC=ARC_CG  )
fwrite( dat, "RICU_sepsis_com3.csv"       )
