#This script relates to the manuscript "Juvenile idiopathic arthritis polygenic risk scores are associated with cardiovascular phenotypes in early adulthood". 
#This script details the ALSPAC variables used in this analysis, details of data cleaning and variable recoding.
#This script excludes the variable extraction from the ALSPAC data catalogue as this differs between internal/external data users. 
#This script was written by Dr SLN Clarke (University of Bristol) with assistance from Dr GC Sharp (University of Exeter) in February 2022.
#For queries related to this script please contact sarah.clarke@bristol.ac.uk.

####Load packages####
library(devtools)
library(tidyverse)
library(alspac)
library(haven)
library(childsds)
library(R.utils)
library(data.table)
library(ggplot2)
library(forestplot)
library(meta)

####Set up workspace####
#remove all existing variables
rm(list=ls(all=TRUE)) 

#Set environment
setwd() #path to local working directory

####load functions####
#outlier removal by IQR
remove.iqr3.outliers <- function(x){x[(x<(quantile(x,0.25,na.rm=T)-(IQR(x,na.rm=T)*3)))|(x>(quantile(x,0.75,na.rm=T)+(IQR(x,na.rm=T)*3)))]<-NA;x}

#summarise results
mysummary <- function(x){
  funs <- c(mean,sd,median,IQR)
  lapply(funs,function(f) f(x,na.rm = T))
}

#regression analysis
run_analysis<-function(exposure,covariates,outcomes,df){
  
  # run logistic regression for all binary outcomes:
  binary_key <- key[key$var.type=="binary",]
  if(nrow(binary_key)==0){
    binary_res <- NA
  }else{
    binary_res <- apply(binary_key,1,function(x){
      outcome <- as.character(x[1])
      fit <-try(glm(formula=as.formula(paste( outcome, "~", paste(c(exposure,covariates),collapse="+") )),
                    data=df,family="binomial"),silent=F)
      # extract results
      res<-try(data.frame(summary(fit)$coef,1-(fit$deviance/fit$null.deviance),nrow(model.frame(fit))),silent=F)
      if(class(res)=="try-error"){res<-as.data.frame(t(data.frame(error=rep(NA,6))))
      }else{
        if(is.na(fit$coefficients[2])){res<-as.data.frame(t(data.frame(error=rep(NA,6))))
        }else{
          res <- res
        }}
      res$exposure <-rownames(res)
      try(names(res)<-c("est","se","t.z","p","r2","n","exposure"))
      res$regression_type <-"glm"
      if(is.na(res$est[1])){
        res$ci.l <- NA
        res$ci.u <- NA
        res$intercept_est <- NA
        res$intercept_se <- NA
        res$intercept_p <- NA
      }else{
        res$ci.l <- res$est-(res$se*1.96)
        res$ci.u <- res$est+(res$se*1.96)
        res$intercept_est <- res[1,1]
        res$intercept_se <- res[1,2]
        res$intercept_p <- res[1,4]
      }
      # remove unneeded results
      if(any(c(is.na(covariates),class(fit)=="try-error")==TRUE)){
        res<-res
      }else{
        if(is.na(fit$coefficients[2])){
          res<-res
        }else{
          res<-res[-grep(res$exposure,pattern=paste0(c("(Intercept)",covariates),collapse="|")),]
        }
      }
      res
    })
    
    binary_res <-Map(cbind, binary_res, outcome = binary_key$outcome)
    binary_res <- as.data.frame(do.call(rbind,binary_res))
  }
  
  # run linear regression for all continuous/integer outcomes:
  cont_key <- key[key$var.type=="continuous",]
  if(nrow(cont_key)==0){cont_res <-NA
  }else{
    cont_res <- apply(cont_key,1,function(x){
      outcome <- as.character(x[1])
      fit <-try(lm(formula=as.formula(paste( outcome, "~", paste(c(exposure,covariates),collapse="+") )),data=df),silent=F)
      # extract results
      res<-try(data.frame(summary(fit)$coef,summary(fit)$r.squared,nrow(model.frame(fit))),silent=F)
      if(class(res)=="try-error"){res<-as.data.frame(t(data.frame(error=rep(NA,6))))
      }else{
        if(is.na(fit$coefficients[2])){res<-as.data.frame(t(data.frame(error=rep(NA,6))))
        }else{
          res <- res
        }}
      res$exposure <-rownames(res)
      try(names(res)<-c("est","se","t.z","p","r2","n","exposure"))
      res$regression_type <-"lm"
      if(is.na(res$est[1])){
        res$ci.l <- NA
        res$ci.u <- NA
        res$intercept_est <- NA
        res$intercept_se <- NA
        res$intercept_p <- NA
      }else{
        res$ci.l <- res$est-(res$se*1.96)
        res$ci.u <- res$est+(res$se*1.96)
        res$intercept_est <- res[1,1]
        res$intercept_se <- res[1,2]
        res$intercept_p <- res[1,4]
      }
      # remove unneeded results
      if(any(c(is.na(covariates),class(fit)=="try-error")==TRUE)){
        res<-res
      }else{
        if(is.na(fit$coefficients[2])){
          res<-res
        }else{
          res<-res[-grep(res$exposure,pattern=paste0(c("(Intercept)",covariates),collapse="|")),]
        }
      }
      res
    })
    
    cont_res <-Map(cbind, cont_res, outcome = cont_key$outcome)
    cont_res <- as.data.frame(do.call(rbind,cont_res))
  }
  
  # Combine logistic and linear regression results
  all_res <- rbind.data.frame(binary_res,cont_res)
  row.names(all_res)<-NULL
  all_res[,c("outcome","exposure","n","est","se","ci.l","ci.u","p","r2","t.z","intercept_est","intercept_se","intercept_p")]
}


#####Generate dataframe of ALSPAC offspring from cohort profile and maternal data####
#Get cohort profile (file cp_2b.dta)
cp<-read_dta() #path to ALSPAC data catalogue cohort profile

#the ALSPAC sample is described in the cohort profile file.
cp$covs_sex_binary <- NA
cp$covs_sex_binary[cp$kz021==2]<-"female" #Female
cp$covs_sex_binary[cp$kz021==1]<-"male" #Male

cp$survive_one_year_binary <-NA
cp$survive_one_year_binary[cp$kz011b==1] <-1
cp$survive_one_year_binary[cp$kz011b==2] <-0

nrow(cp[cp$in_alsp==1,])
nrow(cp[cp$in_alsp==1&cp$survive_one_year_binary==1,])

#extract ALSPAC offspring data (file kz_5c.dta)
kz<-read_dta() #path to ALSPAC data catalogue offspring data file

#merge cohort profile with ALSPAC offspring from kz_5c.dta
dat<-NULL
dat<-left_join(cp,kz,by="aln_qlet")

#Get maternal pregnancy data (file mz_5a.dta)
mz<-read_dta() #path to ALSPAC data catalogue maternal data file
mz <- mz[mz$mz001==1,]
mz <- mz[-which(mz$mz010<1),]

#link siblings to mothers in ALSPAC if more than one sibling in core cohort (file mult_pregs in alspac_CORE.dta):
mult_pregs_core<-read_dta() #path to ALSPAC data catalogue multiple pregnancies data file

#idenitfy pregnancy prder for multiple sibings in core cohort (pregnancy_order_no):
mult_pregs_core<-left_join(mult_pregs_core,mz,by="aln")
mult_pregs_core<-mult_pregs_core[order(mult_pregs_core$mz024b,mult_pregs_core$mz024a),]#order by year and month of delivery
mult_pregs_core$pregnancy_order_no_discrete <-1 #first child
mult_pregs_core$pregnancy_order_no_discrete[duplicated(mult_pregs_core$mult_no)]<-2 #second child
mult_pregs_core[-which(is.na(mult_pregs_core$aln3)),][3,"pregnancy_order_no_discrete"]<-3 #third child
mult_pregs_core<-mult_pregs_core[,c("aln","pregnancy_order_no_discrete","mult_no","aln2","aln3")]
mz<-full_join(mz,mult_pregs_core,by="aln")
mz$pregnancy_order_no_discrete[is.na(mz$pregnancy_order_no_discrete)]<-0 #0 for women with only one pregnancy in core sample

dat<-full_join(kz,mz,by="aln")
dat<-left_join(cp,dat,by=c("aln","qlet","aln_qlet"))

#generate unique identifier for child
dat$mother_aln_qlet <- paste0(dat$pregnancy_order_no_discrete,"_",dat$aln_qlet)

####Extract variables for analysis####
####Time of questionnaire completion####
time_vars <- c("a902","b924","c991","d990","e699","pa900","pa901","pb900","pb901")
time_meta_data <- findVars(time_vars) #create a table of meta data for each variable
time_meta_data <- subset(time_meta_data,subset=tolower(name) %in% time_vars)

#extract variables
time_data <- extractVars(time_meta_data) #extracting data and perform withdrawal of consent
time_data<-time_data[,c("aln","qlet",time_vars)]
dat <- left_join(dat,time_data,by=c("aln","qlet")) #merge with main dataframe 

#derive clean variable
dat$time_qA_mother_gestation_discrete<-dat$a902
dat$time_qA_mother_gestation_discrete[dat$time_qA_mother_gestation_discrete<0] <- NA

dat$time_qB_mother_gestation_discrete<-dat$b924
dat$time_qB_mother_gestation_discrete[dat$time_qB_mother_gestation_discrete<0] <- NA

dat$time_qC_mother_gestation_discrete<-dat$c991
dat$time_qC_mother_gestation_discrete[dat$time_qC_mother_gestation_discrete<0] <- NA

dat$time_qD_mother_gestation_discrete<-dat$d990
dat$time_qD_mother_gestation_discrete[dat$time_qD_mother_gestation_discrete<0] <- NA

dat$time_qPA_partner_gestation_discrete<-dat$pa900
dat$time_qPA_partner_gestation_discrete[dat$time_qPA_partner_gestation_discrete<0] <- NA

dat$time_qPB_partner_gestation<-dat$pb900
dat$time_qPB_partner_gestation[dat$time_qPB_partner_gestation<0] <- NA

dat$time_qE_mother_sincebirth_discrete<-dat$e699
dat$time_qE_mother_sincebirth_discrete[dat$time_qE_mother_sincebirth_discrete<0] <- NA

dat$time_qPA_partner_sincebirth_discrete<-dat$pa901
dat$time_qPA_partner_sincebirth_discrete[dat$time_qPA_partner_sincebirth_discrete<0] <- NA

dat$time_qPB_partner_sincebirth_discrete<-dat$pa901
dat$time_qPB_partner_sincebirth_discrete[dat$time_qPB_partner_sincebirth_discrete<0] <- NA

####Anthropometry variables age 7 years to age 24 years####
anthropometry_meta_data <- findVars("sga", "lga", "gestational", "height", "weight", "waist", "circumference", "BMI", "obese", "underweight", "fat", "fat mass percentage")
anthropometry_meta_data_child <-anthropometry_meta_data[anthropometry_meta_data$cat3 %in% c("Child","Child Based","Obstetric"),]
anthropometry_meta_data_child <-anthropometry_meta_data_child[grep(lapply(anthropometry_meta_data_child$obj,substring,1,3),pattern= "f07|f08|f09|f10|F11|tf1|tf2|tf3|tf4|F24|cif|ka|kb|kc|kd|ke|kf|kg|kh|ki|kj|kl|km|kn|ko|kp|kr|ks|kt|ku|kv|kq|kr|ks|kt|ku|kv|kw|ccs|ccb|ccc|ccd|cce|ccf|ccg|cch|cci|ccj|cck|OA|OB|OC"),]

#set up meta-data
child_anthropometry_vars <-c("bestgest","centiles")

#Height - clinic measures only
child_anthropometry_vars <-c(child_anthropometry_vars,c("cf054", "cf055", "cf056", "cf057", "cf058", "f7ms010","f8lf020", "f9ms010","fdms010","fems010",
                                                        "fg3100","fh3000","FJMR020","FKMS1000"))
#WC- clinic measures only
child_anthropometry_vars <-c(child_anthropometry_vars,c("cf075","cf076","cf077","cf078","cf079","f7ms018","f9ms018","fems018","fg3120","fh4020","FKMS1052"))

#Weight- clinic measures only
child_anthropometry_vars <-c(child_anthropometry_vars,c("cf040","cf041","cf042","cf043","cf044","cf045","cf046","cf047","cf048","cf049","f7ms026","f9ms026",
                                                        "fems026","fg3130","fh3010","FJMR022","FKMS1030"))
#BMI
child_anthropometry_vars <-c(child_anthropometry_vars,c("cf060","cf061","cf062","cf063","cf064","cf065","cf066","cf067","cf068","cf069","f7ms026a","f9ms026a","fems026a",
                                                        "fg3139","fh3019","FJMR022a","FKMS1040"))
#Fat mass 
child_anthropometry_vars <-c(child_anthropometry_vars,c("f9dx135","fedx135","fg3254","fh2254","FJDX135","FKDX1001"))

anthropometry_meta_data <- findVars(unique(child_anthropometry_vars)) #create a table ofmeta data for each variable

anthropometry_meta_data <- subset(anthropometry_meta_data,subset=name %in% child_anthropometry_vars) #to ensure only variables of interest are extracted and not variables with longer names that include the same strings

#Some of the variables (bestgest, cf046, cf056, cf066, cf076) have the same name in multiple files, so to remove duplicates:
anthropometry_meta_data <- anthropometry_meta_data[which(anthropometry_meta_data$cat2 %in% c("Clinic","Quest","bestgest","birthweight_centile")),]

#extract variables
child_anthropometry_data <- extractVars(anthropometry_meta_data[anthropometry_meta_data$name %in% unique(child_anthropometry_vars),]) #extract data and perform withdrawal of consent

child_anthropometry_data<-child_anthropometry_data[,c("aln","qlet",child_anthropometry_vars)]
dat <- left_join(dat,child_anthropometry_data,by=c("aln","qlet")) 

#Clean variables
#Variables are coded according to age - stage 3 (7yrs), stage 4 (9yrs), stage 5 (11yrs), stage 6 (13yrs), stage 7 (15yrs), stage 8 (17yrs) and stage 9 (24yrs)
#Height
#stage 3
dat$anthro_height_stage3_continuous<-dat$f7ms010
dat$anthro_height_stage3_continuous[dat$anthro_height_stage3_continuous<0]<-NA
dat$anthro_height_stage3_zscore_continuous <- scale(dat$anthro_height_stage3_continuous)

#stage 4
dat$anthro_height_stage4_continuous<-dat$fdms010
dat$anthro_height_stage4_continuous[dat$anthro_height_stage4_continuous<0]<-NA
dat$anthro_height_stage4_zscore_continuous <- scale(dat$anthro_height_stage4_continuous)

#stage 5
dat$anthro_height_stage5_continuous<-dat$fems010
dat$anthro_height_stage5_continuous[dat$anthro_height_stage5_continuous<0]<-NA
dat$anthro_height_stage5_zscore_continuous <- scale(dat$anthro_height_stage5_continuous)

#stage 6
dat$anthro_height_stage6_continuous<-dat$fg3100
dat$anthro_height_stage6_continuous[dat$anthro_height_stage6_continuous<0]<-NA
dat$anthro_height_stage6_zscore_continuous <- scale(dat$anthro_height_stage6_continuous)

#stage 7
dat$anthro_height_stage7_continuous<-dat$fh3000
dat$anthro_height_stage7_continuous[dat$anthro_height_stage7_continuous<0]<-NA
dat$anthro_height_stage7_zscore_continuous <- scale(dat$anthro_height_stage7_continuous)

#stage 8
dat$anthro_height_stage8_continuous<-dat$FJMR020
dat$anthro_height_stage8_continuous[dat$anthro_height_stage8_continuous<0]<-NA
dat$anthro_height_stage8_zscore_continuous <- scale(dat$anthro_height_stage8_continuous)

#stage 9
dat$anthro_height_stage9_continuous<-dat$FKMS1000/10
dat$anthro_height_stage9_continuous[dat$anthro_height_stage9_continuous<0]<-NA
dat$anthro_height_stage9_zscore_continuous <- scale(dat$anthro_height_stage9_continuous)

#Waist circumference
#stage 3
dat$anthro_waist_stage3_continuous<-dat$f7ms018
dat$anthro_waist_stage3_continuous[dat$anthro_waist_stage3_continuous<0]<-NA
dat$anthro_waist_stage3_zscore_continuous <- scale(dat$anthro_waist_stage3_continuous)

#stage 4
dat$anthro_waist_stage4_continuous<-dat$f9ms018
dat$anthro_waist_stage4_continuous[dat$anthro_waist_stage4_continuous<0]<-NA
dat$anthro_waist_stage4_zscore_continuous <- scale(dat$anthro_waist_stage4_continuous)

#stage 5
dat$anthro_waist_stage5_continuous<-dat$fems018
dat$anthro_waist_stage5_continuous[dat$anthro_waist_stage5_continuous<0]<-NA
dat$anthro_waist_stage5_zscore_continuous <- scale(dat$anthro_waist_stage5_continuous)

#stage 6
dat$anthro_waist_stage6_continuous<-dat$fg3120
dat$anthro_waist_stage6_continuous[dat$anthro_waist_stage6_continuous<0]<-NA
dat$anthro_waist_stage6_zscore_continuous <- scale(dat$anthro_waist_stage6_continuous)

#stage 7
dat$anthro_waist_stage7_continuous<-dat$fh4020
dat$anthro_waist_stage7_continuous[dat$anthro_waist_stage7_continuous<0]<-NA
dat$anthro_waist_stage7_zscore_continuous <- scale(dat$anthro_waist_stage7_continuous)

#stage 8 - no waist circumference measurement
dat$anthro_waist_stage8_continuous<-NA
dat$anthro_waist_stage8_continuous[dat$anthro_waist_stage8_continuous<0]<-NA
dat$anthro_waist_stage8_zscore_continuous <- scale(dat$anthro_waist_stage8_continuous)

#stage 9
dat$anthro_waist_stage9_continuous<-dat$FKMS1052
dat$anthro_waist_stage9_continuous[dat$anthro_waist_stage9_continuous<0]<-NA
dat$anthro_waist_stage9_zscore_continuous <- scale(dat$anthro_waist_stage9_continuous)

#Weight
#stage 3
dat$anthro_weight_stage3_continuous<-dat$f7ms026
dat$anthro_weight_stage3_continuous[dat$anthro_weight_stage3_continuous<0]<-NA
dat$anthro_weight_stage3_zscore_continuous <- scale(dat$anthro_weight_stage3_continuous)

#stage 4
dat$anthro_weight_stage4_continuous<-dat$f9ms026
dat$anthro_weight_stage4_continuous[dat$anthro_weight_stage4_continuous<0]<-NA
dat$anthro_weight_stage4_zscore_continuous <- scale(dat$anthro_weight_stage4_continuous)

#stage 5
dat$anthro_weight_stage5_continuous<-dat$fems026
dat$anthro_weight_stage5_continuous[dat$anthro_weight_stage5_continuous<0]<-NA
dat$anthro_weight_stage5_zscore_continuous <- scale(dat$anthro_weight_stage5_continuous)

#stage 6
dat$anthro_weight_stage6_continuous<-dat$fg3130
dat$anthro_weight_stage6_continuous[dat$anthro_weight_stage6_continuous<0]<-NA
dat$anthro_weight_stage6_zscore_continuous <- scale(dat$anthro_weight_stage6_continuous)

#stage 7
dat$anthro_weight_stage7_continuous<-dat$fh3010
dat$anthro_weight_stage7_continuous[dat$anthro_weight_stage7_continuous<0]<-NA
dat$anthro_weight_stage7_zscore_continuous <- scale(dat$anthro_weight_stage7_continuous)

#stage 8
dat$anthro_weight_stage8_continuous<-dat$FJMR022
dat$anthro_weight_stage8_continuous[dat$anthro_weight_stage8_continuous<0]<-NA
dat$anthro_weight_stage8_zscore_continuous <- scale(dat$anthro_weight_stage8_continuous)

#stage 9
dat$anthro_weight_stage9_continuous<-dat$FKMS1030
dat$anthro_weight_stage9_continuous[dat$anthro_weight_stage9_continuous<0]<-NA
dat$anthro_weight_stage9_zscore_continuous <- scale(dat$anthro_weight_stage9_continuous)

#Age at clinic visit
#stage 3
f7 <-read_dta() #path to ALSPAC file f07_5a.dta
#stage 4
f9 <-read_dta() #path to ALSPAC file f09_4c.dta
#stage 5
f11 <-read_dta() #path to ALSPAC file F11_5d.dta
#stage 6
tf2 <-read_dta() #path to ALSPAC file tf2_5a.dta
#stage 7
tf3 <-read_dta() #path to ALSPAC file tf3_4c.dta
#stage 8
tf4 <-read_dta() #path to ALSPAC file tf4_6a.dta
#stage 9
f24 <-read_dta() #path to ALSPAC file F24_5a.dta

ages <- merge(f7[,c("aln","qlet","f7003c")],f9[,c("aln","qlet","f9003c")],by=c("aln","qlet"),all=T)
ages <- merge(ages,f11[,c("aln","qlet","fe003c")],by=c("aln","qlet"),all=T)
ages <- merge(ages,tf2[,c("aln","qlet","fg0011a")],by=c("aln","qlet"),all=T)
ages <- merge(ages,tf3[,c("aln","qlet","fh0011a")],by=c("aln","qlet"),all=T)
ages <- merge(ages,tf4[,c("aln","qlet","FJ003a")],by=c("aln","qlet"),all=T)
ages <- merge(ages,f24[,c("aln","qlet","FKAR0010")],by=c("aln","qlet"),all=T)

ages$f7003c[ages$f7003c<0]<-NA
ages$f9003c[ages$f9003c<0]<-NA
ages$fe003c[ages$fe003c<0]<-NA
ages$fg0011a[ages$fg0011a<0]<-NA
ages$fh0011a[ages$fh0011a<0]<-NA
ages$FJ003a[ages$FJ003a<0]<-NA
ages$FKAR0010[ages$FKAR0010<0]<-NA

dat <- left_join(dat,ages,by=c("aln","qlet")) 

dat$covs_age_child_stage3_f7_continuous <- dat$f7003c
dat$covs_age_child_stage4_f9_continuous <- dat$f9003c
dat$covs_age_child_stage5_f11_continuous <- dat$fe003c
dat$covs_age_child_stage6_tf2_continuous <- dat$fg0011a
dat$covs_age_child_stage7_tf3_continuous <- dat$fh0011a
dat$covs_age_child_stage8_tf4_continuous <- dat$FJ003a
dat$covs_age_child_stage9_f24_continuous <- dat$FKAR0010


#BMI
#stage 3
dat$anthro_bmi_stage3_continuous <- dat$anthro_weight_stage3_continuous / ((dat$anthro_height_stage3_continuous/100)^2)

#stage 4
dat$anthro_bmi_stage4_continuous <- dat$anthro_weight_stage4_continuous / ((dat$anthro_height_stage4_continuous/100)^2)

#stage 5
dat$anthro_bmi_stage5_continuous <- dat$anthro_weight_stage5_continuous / ((dat$anthro_height_stage5_continuous/100)^2)

#stage 6
dat$anthro_bmi_stage6_continuous <- dat$anthro_weight_stage6_continuous / ((dat$anthro_height_stage6_continuous/100)^2)

#stage 7
dat$anthro_bmi_stage7_continuous <- dat$anthro_weight_stage7_continuous / ((dat$anthro_height_stage7_continuous/100)^2)

#stage 8
dat$anthro_bmi_stage8_continuous <- dat$anthro_weight_stage8_continuous / ((dat$anthro_height_stage8_continuous/100)^2)

#stage 9
dat$anthro_bmi_stage9_continuous <- dat$anthro_weight_stage9_continuous / ((dat$anthro_height_stage9_continuous/100)^2)

#WHO obese and overweight - child values take into account age and sex and use SDS scores, stage 8 and 9 use adult BMI cut offs
#stage 3
dat$anthro_bmi_stage3_sds_continuous <- sds(value = dat$anthro_bmi_stage3_continuous,
                                            age = dat$covs_age_child_stage3_f7_continuous/12,
                                            sex = dat$covs_sex_binary, male = "male", female = "female",
                                            ref = ukwho.ref, item = "bmi")
#stage 4
dat$anthro_bmi_stage4_sds_continuous <- sds(value = dat$anthro_bmi_stage4_continuous,
                                            age = dat$covs_age_child_stage4_f9_continuous/12,
                                            sex = dat$covs_sex_binary, male = "male", female = "female",
                                            ref = ukwho.ref, item = "bmi")

#stage 5
dat$anthro_bmi_stage5_sds_continuous <- sds(value = dat$anthro_bmi_stage5_continuous,
                                            age = dat$covs_age_child_stage5_f11_continuous/12,
                                            sex = dat$covs_sex_binary, male = "male", female = "female",
                                            ref = ukwho.ref, item = "bmi")
#stage 6
dat$anthro_bmi_stage6_sds_continuous <- sds(value = dat$anthro_bmi_stage6_continuous,
                                            age = dat$covs_age_child_stage6_tf2_continuous/12,
                                            sex = dat$covs_sex_binary, male = "male", female = "female",
                                            ref = ukwho.ref, item = "bmi")

#stage 7
dat$anthro_bmi_stage7_sds_continuous <- sds(value = dat$anthro_bmi_stage7_continuous,
                                            age = dat$covs_age_child_stage7_tf3_continuous/12,
                                            sex = dat$covs_sex_binary, male = "male", female = "female",
                                            ref = ukwho.ref, item = "bmi")


#stage 8 SDS not used as some participants >18yrs

#stage 9 SDS not used as all participants >18yrs

#define BMI category (for <18yrs - underweight <-2SD, overweight >+1SD, obese >+2SD, for stage 8 and 9 adult BMI cut offs used)
#stage 3
dat$anthro_overweightobese_stage3_binary <- NA
dat$anthro_overweightobese_stage3_binary[dat$anthro_bmi_stage3_sds_continuous>1] <-1
dat$anthro_overweightobese_stage3_binary[which(dat$anthro_bmi_stage3_sds_continuous>=(-2) &dat$anthro_bmi_stage3_sds_continuous<=1)] <-0

dat$anthro_obese_stage3_binary <- NA
dat$anthro_obese_stage3_binary[which(dat$anthro_bmi_stage3_sds_continuous>2)] <-1
dat$anthro_obese_stage3_binary[which(dat$anthro_bmi_stage3_sds_continuous<=1 &dat$anthro_bmi_stage3_sds_continuous>=(-2))] <-0

dat$anthro_underweight_stage3_binary <- NA
dat$anthro_underweight_stage3_binary[which(dat$anthro_bmi_stage3_sds_continuous<(-2))] <-1
dat$anthro_underweight_stage3_binary[which((dat$anthro_bmi_stage3_sds_continuous>=(-2) &dat$anthro_bmi_stage3_sds_continuous<=1))] <-0

#stage 4
dat$anthro_overweightobese_stage4_binary <- NA
dat$anthro_overweightobese_stage4_binary[dat$anthro_bmi_stage4_sds_continuous>1] <-1
dat$anthro_overweightobese_stage4_binary[which(dat$anthro_bmi_stage4_sds_continuous>=(-2) &dat$anthro_bmi_stage4_sds_continuous<=1)] <-0

dat$anthro_obese_stage4_binary <- NA
dat$anthro_obese_stage4_binary[which(dat$anthro_bmi_stage4_sds_continuous>2)] <-1
dat$anthro_obese_stage4_binary[which(dat$anthro_bmi_stage4_sds_continuous<=1 &dat$anthro_bmi_stage4_sds_continuous>=(-2))] <-0

dat$anthro_underweight_stage4_binary <- NA
dat$anthro_underweight_stage4_binary[which(dat$anthro_bmi_stage4_sds_continuous<(-2))] <-1
dat$anthro_underweight_stage4_binary[which((dat$anthro_bmi_stage4_sds_continuous>=(-2) &dat$anthro_bmi_stage4_sds_continuous<=1))] <-0

#stage 5
dat$anthro_overweightobese_stage5_binary <- NA
dat$anthro_overweightobese_stage5_binary[dat$anthro_bmi_stage5_sds_continuous>1] <-1
dat$anthro_overweightobese_stage5_binary[which(dat$anthro_bmi_stage5_sds_continuous>=(-2) &dat$anthro_bmi_stage5_sds_continuous<=1)] <-0

dat$anthro_obese_stage5_binary <- NA
dat$anthro_obese_stage5_binary[which(dat$anthro_bmi_stage5_sds_continuous>2)] <-1
dat$anthro_obese_stage5_binary[which(dat$anthro_bmi_stage5_sds_continuous>=(-2) &dat$anthro_bmi_stage5_sds_continuous<=1)] <-0

dat$anthro_underweight_stage5_binary <- NA
dat$anthro_underweight_stage5_binary[which(dat$anthro_bmi_stage5_sds_continuous<(-2))] <-1
dat$anthro_underweight_stage5_binary[which(dat$anthro_bmi_stage5_sds_continuous>=(-2) &dat$anthro_bmi_stage5_sds_continuous<=1)] <-0

#stage 6
dat$anthro_overweightobese_stage6_binary <- NA
dat$anthro_overweightobese_stage6_binary[dat$anthro_bmi_stage6_sds_continuous>1] <-1
dat$anthro_overweightobese_stage6_binary[which(dat$anthro_bmi_stage6_sds_continuous>=(-2) &dat$anthro_bmi_stage6_sds_continuous<=1)] <-0

dat$anthro_obese_stage6_binary <- NA
dat$anthro_obese_stage6_binary[which(dat$anthro_bmi_stage6_sds_continuous>2)] <-1
dat$anthro_obese_stage6_binary[which((dat$anthro_bmi_stage6_sds_continuous>=(-2) &dat$anthro_bmi_stage6_sds_continuous<=1))] <-0

dat$anthro_underweight_stage6_binary <- NA
dat$anthro_underweight_stage6_binary[which(dat$anthro_bmi_stage6_sds_continuous<(-2))] <-1
dat$anthro_underweight_stage6_binary[which((dat$anthro_bmi_stage6_sds_continuous>=(-2) &dat$anthro_bmi_stage6_sds_continuous<=1))] <-0

#stage 7
dat$anthro_overweightobese_stage7_binary <- NA
dat$anthro_overweightobese_stage7_binary[dat$anthro_bmi_stage7_sds_continuous>1] <-1
dat$anthro_overweightobese_stage7_binary[which(dat$anthro_bmi_stage7_sds_continuous>=(-2) &dat$anthro_bmi_stage7_sds_continuous<=1)] <-0

dat$anthro_obese_stage7_binary <- NA
dat$anthro_obese_stage7_binary[which(dat$anthro_bmi_stage7_sds_continuous>2)] <-1
dat$anthro_obese_stage7_binary[which((dat$anthro_bmi_stage7_sds_continuous>=(-2) &dat$anthro_bmi_stage6_sds_continuous<=1))] <-0

dat$anthro_underweight_stage7_binary <- NA
dat$anthro_underweight_stage7_binary[which(dat$anthro_bmi_stage7_sds_continuous<(-2))] <-1
dat$anthro_underweight_stage7_binary[which((dat$anthro_bmi_stage7_sds_continuous>=(-2) &dat$anthro_bmi_stage7_sds_continuous<=1))] <-0

#stage 8
dat$anthro_overweightobese_stage8_binary <- NA
dat$anthro_overweightobese_stage8_binary[which(dat$anthro_bmi_stage8_continuous>=25)] <-1
dat$anthro_overweightobese_stage8_binary[which(dat$anthro_bmi_stage8_continuous<25& dat$anthro_bmi_stage8_continuous>=18.5)] <-0

dat$anthro_obese_stage8_binary <- NA
dat$anthro_obese_stage8_binary[which(dat$anthro_bmi_stage8_continuous>=30)] <-1
dat$anthro_obese_stage8_binary[which(dat$anthro_bmi_stage8_continuous<25& dat$anthro_bmi_stage8_continuous>=18.5)] <-0

dat$anthro_underweight_stage8_binary<-NA
dat$anthro_underweight_stage8_binary[which(dat$anthro_bmi_stage8_continuous<18.5)] <-1
dat$anthro_underweight_stage8_binary[which(dat$anthro_bmi_stage8_continuous<25& dat$anthro_bmi_stage8_continuous>=18.5)] <-0

#stage 9
dat$anthro_overweightobese_stage9_binary <- NA
dat$anthro_overweightobese_stage9_binary[which(dat$anthro_bmi_stage9_continuous>=25)] <-1
dat$anthro_overweightobese_stage9_binary[which(dat$anthro_bmi_stage9_continuous<25& dat$anthro_bmi_stage9_continuous>=18.5)] <-0

dat$anthro_obese_stage9_binary <- NA
dat$anthro_obese_stage9_binary[which(dat$anthro_bmi_stage9_continuous>=30)] <-1
dat$anthro_obese_stage9_binary[which(dat$anthro_bmi_stage9_continuous<25& dat$anthro_bmi_stage9_continuous>=18.5)] <-0

dat$anthro_underweight_stage9_binary<-NA
dat$anthro_underweight_stage9_binary[which(dat$anthro_bmi_stage9_continuous<18.5)] <-1
dat$anthro_underweight_stage9_binary[which(dat$anthro_bmi_stage9_continuous<25& dat$anthro_bmi_stage9_continuous>=18.5)] <-0

#Fat mass index
#stage 4
dat$anthro_fmi_stage4_continuous <- dat$f9dx135
dat$anthro_fmi_stage4_continuous[dat$f9dx135<0]<-NA
dat$anthro_fmi_stage4_continuous <- dat$anthro_fmi_stage4_continuous / ((dat$anthro_height_stage4_continuous/100)^2)
dat$anthro_fmi_stage4_zscore_continuous <- scale(dat$anthro_fmi_stage4_continuous)

#stage 5
dat$anthro_fmi_stage5_continuous <- dat$fedx135
dat$anthro_fmi_stage5_continuous[dat$fedx135<0]<-NA
dat$anthro_fmi_stage5_continuous <- dat$anthro_fmi_stage5_continuous / ((dat$anthro_height_stage5_continuous/100)^2)
dat$anthro_fmi_stage5_zscore_continuous <- scale(dat$anthro_fmi_stage5_continuous)

#stage 6
dat$anthro_fmi_stage6_continuous <- dat$fg3254
dat$anthro_fmi_stage6_continuous[dat$fg3254<0]<-NA
dat$anthro_fmi_stage6_continuous <- dat$anthro_fmi_stage6_continuous / ((dat$anthro_height_stage6_continuous/100)^2)
dat$anthro_fmi_stage6_zscore_continuous <- scale(dat$anthro_fmi_stage6_continuous)

#stage 8
dat$anthro_fmi_stage8_continuous <- dat$FJDX135
dat$anthro_fmi_stage8_continuous[dat$FJDX135<0]<-NA
dat$anthro_fmi_stage8_continuous <- dat$anthro_fmi_stage8_continuous / ((dat$anthro_height_stage8_continuous/100)^2)
dat$anthro_fmi_stage8_zscore_continuous <- scale(dat$anthro_fmi_stage8_continuous)

#stage 9
dat$anthro_fmi_stage9_continuous <- dat$FKDX1001
dat$anthro_fmi_stage9_continuous[dat$FKDX1001<0]<-NA
dat$anthro_fmi_stage9_continuous <- dat$anthro_fmi_stage9_continuous / ((dat$anthro_height_stage9_continuous/100)^2)
dat$anthro_fmi_stage9_zscore_continuous <- scale(dat$anthro_fmi_stage9_continuous)

#####Immunological variables age 7 years to age 24 years####
#set up meta-data
#get a list of immunological measures asked for offspring
immunological_meta_data <- findVars("asthma", "eczema", "allergies", "allergic", "pollen","df","dp","dust","cat","dog","grass","sesame","almond","nut","peanut","milk","latex", "wheeze", "whistle", "whilstling", "wheezing", "wheezed","reaction","bee","mite","atopy","rash","egg","allergy","wasp")
immunological_meta_data <-immunological_meta_data[immunological_meta_data$cat3 %in% c("Child","Child Based"),]

#select measures up to age 24
immunological_meta_data <-immunological_meta_data[grep(lapply(immunological_meta_data$obj,substring,1,3),pattern="f07|f08|f09|f10|F11|tf1|tf2|tf3|tf4|F24|cif|ka|kb|kc|kd|ke|kf|kg|kh|ki|kj|kk|kl|km|kn|ko|kp|kq|kr|ks|kt|ku|kv|kw|ccs|ccb|ccc|ccd|cce|ccf|ccg|cch|cci|ccj|cck|ta|tb|tc|ccs|YPB"),]

#list variables to extract
immuno_vars<-unique(c("YPB1212","YPB1218","YPB1219","YPB1220","YPB1221","YPB1222","YPB1223","YPB1224","YPB1225","YPB1226","YPB1227","YPB1228","YPB1229","CRP_cord","CRP_f9","crp_TF3","CRP_TF4","CRP_F24","YPB9992"))

#set up final meta-data
immuno_meta_data <- findVars(unique(immuno_vars)) 
immuno_meta_data <- subset(immuno_meta_data,subset=name %in% immuno_vars) 

#extract variables
immuno_data <- extractVars(immuno_meta_data[immuno_meta_data$name %in% unique(immuno_vars),]) 
immuno_data<-immuno_data[,c("aln","qlet",immuno_vars)]

dat <- left_join(dat,immuno_data,by=c("aln","qlet")) 

#autoimmunity
#psoriasis
dat$immuno_psoriasis_stage9_binary_ever<-NA
dat$immuno_psoriasis_stage9_binary_ever[which(dat$YPB1212==3)] <-0
dat$immuno_psoriasis_stage9_binary_ever[which(dat$YPB1212%in%c(1,2))] <-1

#Crohn's
dat$immuno_cd_stage9_binary_ever<-NA
dat$immuno_cd_stage9_binary_ever[which(dat$YPB1218==3)] <-0
dat$immuno_cd_stage9_binary_ever[which(dat$YPB1218%in%c(1,2))] <-1

#ulcerative colitis
dat$immuno_uc_stage9_binary_ever<-NA
dat$immuno_uc_stage9_binary_ever[which(dat$YPB1219==3)] <-0
dat$immuno_uc_stage9_binary_ever[which(dat$YPB1219%in%c(1,2))] <-1

#ankylosing spondylitis
dat$immuno_as_stage9_binary_ever<-NA
dat$immuno_as_stage9_binary_ever[which(dat$YPB1220==3)] <-0
dat$immuno_as_stage9_binary_ever[which(dat$YPB1220%in%c(1,2))] <-1

#psoriatic arthritis
dat$immuno_psa_stage9_binary_ever<-NA
dat$immuno_psa_stage9_binary_ever[which(dat$YPB1221==3)] <-0
dat$immuno_psa_stage9_binary_ever[which(dat$YPB1221%in%c(1,2))] <-1

#spondyloarthropathy
dat$immuno_sarth_stage9_binary_ever<-NA
dat$immuno_sarth_stage9_binary_ever[which(dat$YPB1222==3)] <-0
dat$immuno_sarth_stage9_binary_ever[which(dat$YPB1222%in%c(1,2))] <-1

#rheumatoid arthritis
dat$immuno_RA_stage9_binary_ever<-NA
dat$immuno_RA_stage9_binary_ever[which(dat$YPB1223==3)] <-0
dat$immuno_RA_stage9_binary_ever[which(dat$YPB1223%in%c(1,2))] <-1

#Sjogrens
dat$immuno_sjs_stage9_binary_ever<-NA
dat$immuno_sjs_stage9_binary_ever[which(dat$YPB1224==3)] <-0
dat$immuno_sjs_stage9_binary_ever[which(dat$YPB1224%in%c(1,2))] <-1

#systemic lupus erythematosus
dat$immuno_sle_stage9_binary_ever<-NA
dat$immuno_sle_stage9_binary_ever[which(dat$YPB1225==3)] <-0
dat$immuno_sle_stage9_binary_ever[which(dat$YPB1225%in%c(1,2))] <-1

#Grave's
dat$immuno_graves_stage9_binary_ever<-NA
dat$immuno_graves_stage9_binary_ever[which(dat$YPB1226==3)] <-0
dat$immuno_graves_stage9_binary_ever[which(dat$YPB1226%in%c(1,2))] <-1

#Multiple sclerosis
dat$immuno_ms_stage9_binary_ever<-NA
dat$immuno_ms_stage9_binary_ever[which(dat$YPB1227==3)] <-0
dat$immuno_ms_stage9_binary_ever[which(dat$YPB1227%in%c(1,2))] <-1

#Hashimoto's thyroiditis
dat$immuno_hashi_stage9_binary_ever<-NA
dat$immuno_hashi_stage9_binary_ever[which(dat$YPB1228==3)] <-0
dat$immuno_hashi_stage9_binary_ever[which(dat$YPB1228%in%c(1,2))] <-1

#Type 1 diabetes
dat$immuno_t1d_stage9_binary_ever<-NA
dat$immuno_t1d_stage9_binary_ever[which(dat$YPB1229==3)] <-0
dat$immuno_t1d_stage9_binary_ever[which(dat$YPB1229%in%c(1,2))] <-1

#Any autoimmune disease
dat$immuno_autoimmune_any_stage9_binary_ever<-NA

dat$immuno_autoimmune_any_stage9_binary_ever<- rowSums(dat[,c("immuno_psoriasis_stage9_binary_ever","immuno_cd_stage9_binary_ever","immuno_uc_stage9_binary_ever","immuno_as_stage9_binary_ever",
                                                              "immuno_psa_stage9_binary_ever","immuno_sarth_stage9_binary_ever","immuno_RA_stage9_binary_ever","immuno_sjs_stage9_binary_ever",
                                                              "immuno_sle_stage9_binary_ever","immuno_graves_stage9_binary_ever","immuno_ms_stage9_binary_ever","immuno_hashi_stage9_binary_ever",
                                                              "immuno_t1d_stage9_binary_ever")],na.rm=T)
dat$immuno_autoimmune_any_stage9_binary_ever[dat$immuno_autoimmune_any_stage9_binary_ever>0]<-1
dat$immuno_autoimmune_any_stage9_binary_ever[which(is.na(dat$immuno_psoriasis_stage9_binary_ever) & is.na(dat$immuno_cd_stage9_binary_ever) & is.na(dat$immuno_uc_stage9_binary_ever) & is.na(dat$immuno_as_stage9_binary_ever)
                                                   & is.na(dat$immuno_psa_stage9_binary_ever) & is.na(dat$immuno_sarth_stage9_binary_ever) & is.na(dat$immuno_RA_stage9_binary_ever) & is.na(dat$immuno_sjs_stage9_binary_ever)
                                                   & is.na(dat$immuno_sle_stage9_binary_ever) & is.na(dat$immuno_graves_stage9_binary_ever) & is.na(dat$immuno_ms_stage9_binary_ever) & is.na(dat$immuno_hashi_stage9_binary_ever)
                                                   & is.na(dat$immuno_t1d_stage9_binary_ever))] <- NA

#CRP
dat$immuno_CRP_stage4_continuous <- dat$CRP_f9
dat$immuno_CRP_stage4_continuous[dat$CRP_f9<0]<-NA
dat$immuno_CRP_stage4_zscore_continuous <- scale(dat$immuno_CRP_stage4_continuous)
dat$immuno_CRP_stage4_no_outliers_continuous<- remove.iqr3.outliers(dat$immuno_CRP_stage4_continuous)
dat$immuno_CRP_stage4_no_outliers_zscore_continuous <- scale(dat$immuno_CRP_stage4_no_outliers_continuous)

dat$immuno_CRP_stage7_continuous <- dat$crp_TF3
dat$immuno_CRP_stage7_continuous[dat$crp_TF3<0]<-NA
dat$immuno_CRP_stage7_zscore_continuous <- scale(dat$immuno_CRP_stage7_continuous)
dat$immuno_CRP_stage7_no_outliers_continuous<- remove.iqr3.outliers(dat$immuno_CRP_stage7_continuous)
dat$immuno_CRP_stage7_no_outliers_zscore_continuous <- scale(dat$immuno_CRP_stage7_no_outliers_continuous)

dat$immuno_CRP_stage8_continuous <- dat$CRP_TF4
dat$immuno_CRP_stage8_continuous[dat$CRP_TF4<0]<-NA
dat$immuno_CRP_stage8_zscore_continuous <- scale(dat$immuno_CRP_stage8_continuous)
dat$immuno_CRP_stage8_no_outliers_continuous<- remove.iqr3.outliers(dat$immuno_CRP_stage8_continuous)
dat$immuno_CRP_stage8_no_outliers_zscore_continuous <- scale(dat$immuno_CRP_stage8_no_outliers_continuous)

dat$immuno_CRP_stage9_continuous <- dat$CRP_F24
dat$immuno_CRP_stage9_continuous[dat$CRP_F24<0]<-NA
dat$immuno_CRP_stage9_zscore_continuous <- scale(dat$immuno_CRP_stage9_continuous)
dat$immuno_CRP_stage9_no_outliers_continuous<- remove.iqr3.outliers(dat$immuno_CRP_stage9_continuous)
dat$immuno_CRP_stage9_no_outliers_zscore_continuous <- scale(dat$immuno_CRP_stage9_no_outliers_continuous)

#######Cardiovascular variables age 7 years to age 24 years####
#set up meta-data
cardiomet_vars<-c("f7sa022", "f7sa021", "TRIG_F7", "HDL_F7", "LDL_F7", "CHOL_F7", "Glc_F7", "ApoA1_F7", "ApoB_F7","ApoBApoA1_F7", "f9sa022", "f9sa021", "trig_f9", "HDL_f9", 
                  "LDL_f9", "CHOL_F9", "insulin_F9", "ApoAI_f9", "ApoB_f9", "fesa022", "fesa021", "fg1022", "fg1021","fg6120","fg6125","fg6121","fg6126", "fh2031", "fh2036", "fh2030", "fh2035", "trig_TF3", "hdl_TF3", 
                  "ldl_TF3", "chol_TF3", "Glc_TF3", "insulin_TF3", "ApoA1_TF3", "ApoB_TF3","ApoBApoA1_TF3", "FKBP1031","FKBP1030", "Trig_F24", "HDL_F24", "LDL_F24", "Chol_F24", "Glc_F24", "Insulin_F24", 
                  "ApoA1_F24", "ApoB_F24","ApoBApoA1_F24","FJAR019b","FJAR019a", "TRIG_TF4","HDL_TF4","LDL_TF4","CHOL_TF4","Glc_TF4","insulin_TF4","ApoA1_TF4","ApoB_TF4" ,"ApoBApoA1_TF4","Adiponectin_TF3","Adiponectin_f9","LEPTIN_f9",
                  "Gp_F7", "Gp_TF3","Gp_TF4","Gp_F24","FKCV1131","FKCV2131","FJAR083d","FKCV4200","FKEC5330","FKEC5350","FKEC5250","FKEC5050","FKEC5180","FKEC5080","FKEC5090","FKEC5250","FKEC5310",
                  "FJGR114",	"FJGR113",	"FJGR089",	"FJGR101",	"FJGR031",	"FJGR053",	"FJGR043",	"FJGR047","FKEC5070","FJGR039")

cardiomet_meta_data <- findVars(cardiomet_vars)

cardiomet_meta_data <- subset(cardiomet_meta_data,subset=name %in% cardiomet_vars) # to ensure only extract variables of interest and not variables with longer names that include the same strings
cardiomet_meta_data <- cardiomet_meta_data[cardiomet_meta_data$cat1=="Current",]

#Extract variables
cardiomet_data <- extractVars(cardiomet_meta_data[cardiomet_meta_data$name %in% unique(cardiomet_vars),]) #extracting data and perform withdrawal of consent

cardiomet_data<-cardiomet_data[,c("aln","qlet",cardiomet_vars)]
dat <- left_join(dat,cardiomet_data,by=c("aln","qlet"))

#Clean variables
#systolic BP
#stage 3
dat$cardio_sbp_stage3_continuous <- dat$f7sa021
dat$cardio_sbp_stage3_continuous[dat$f7sa021<0]<-NA
dat$cardio_sbp_stage3_continuous<- remove.iqr3.outliers(dat$cardio_sbp_stage3_continuous)
dat$cardio_sbp_stage3_zscore_continuous <- scale(dat$cardio_sbp_stage3_continuous)

#stage4
dat$cardio_sbp_stage4_continuous <- dat$f9sa021
dat$cardio_sbp_stage4_continuous[dat$f9sa021<0]<-NA
dat$cardio_sbp_stage4_continuous<- remove.iqr3.outliers(dat$cardio_sbp_stage4_continuous)
dat$cardio_sbp_stage4_zscore_continuous <- scale(dat$cardio_sbp_stage4_continuous)

#stage5
dat$cardio_sbp_stage5_continuous <- dat$fesa021
dat$cardio_sbp_stage5_continuous[dat$fesa021<0]<-NA
dat$cardio_sbp_stage5_continuous<- remove.iqr3.outliers(dat$cardio_sbp_stage5_continuous)
dat$cardio_sbp_stage5_zscore_continuous <- scale(dat$cardio_sbp_stage5_continuous)

#stage 6
dat$cardio_sbp_1_stage6_continuous<-NA
dat$cardio_sbp_1_stage6_continuous[which(dat$fg6120>=0)]<-dat$fg6120[which(dat$fg6120>=0)]
dat$cardio_sbp_2_stage6_continuous<-NA
dat$cardio_sbp_2_stage6_continuous[which(dat$fg6125>=0)]<-dat$fg6125[which(dat$fg6125>=0)]
dat$cardio_sbp_stage6_continuous<-NA
dat$cardio_sbp_stage6_continuous[is.na(dat$cardio_sbp_1_stage6_continuous)]<-dat$cardio_sbp_2_stage6_continuous[is.na(dat$cardio_sbp_1_stage6_continuous)]
dat$cardio_sbp_stage6_continuous[is.na(dat$cardio_sbp_2_stage6_continuous)]<-dat$cardio_sbp_1_stage6_continuous[is.na(dat$cardio_sbp_2_stage6_continuous)]
dat$cardio_sbp_stage6_continuous[which(!is.na(dat$cardio_sbp_2_stage6_continuous)&!is.na(dat$cardio_sbp_1_stage6_continuous))]<-(dat$cardio_sbp_1_stage6_continuous[which(!is.na(dat$cardio_sbp_2_stage6_continuous)&!is.na(dat$cardio_sbp_1_stage6_continuous))]+dat$cardio_sbp_2_stage6_continuous[which(!is.na(dat$cardio_sbp_2_stage6_continuous)&!is.na(dat$cardio_sbp_1_stage6_continuous))])/2
dat$cardio_sbp_stage6_continuous<- remove.iqr3.outliers(dat$cardio_sbp_stage6_continuous)
dat$cardio_sbp_stage6_zscore_continuous <- scale(dat$cardio_sbp_stage6_continuous)

#stage 7
dat$cardio_sbp_1_stage7_continuous<-NA
dat$cardio_sbp_1_stage7_continuous[which(dat$fh2030>=0)]<-dat$fh2030[which(dat$fh2030>=0)]
dat$cardio_sbp_2_stage7_continuous<-NA
dat$cardio_sbp_2_stage7_continuous[which(dat$fh2035>=0)]<-dat$fh2035[which(dat$fh2035>=0)]
dat$cardio_sbp_stage7_continuous<-NA
dat$cardio_sbp_stage7_continuous[is.na(dat$cardio_sbp_1_stage7_continuous)]<-dat$cardio_sbp_2_stage7_continuous[is.na(dat$cardio_sbp_1_stage7_continuous)]
dat$cardio_sbp_stage7_continuous[is.na(dat$cardio_sbp_2_stage7_continuous)]<-dat$cardio_sbp_1_stage7_continuous[is.na(dat$cardio_sbp_2_stage7_continuous)]
dat$cardio_sbp_stage7_continuous[which(!is.na(dat$cardio_sbp_2_stage7_continuous)&!is.na(dat$cardio_sbp_1_stage7_continuous))]<-(dat$cardio_sbp_1_stage7_continuous[which(!is.na(dat$cardio_sbp_2_stage7_continuous)&!is.na(dat$cardio_sbp_1_stage7_continuous))]+dat$cardio_sbp_2_stage7_continuous[which(!is.na(dat$cardio_sbp_2_stage7_continuous)&!is.na(dat$cardio_sbp_1_stage7_continuous))])/2 #
dat$cardio_sbp_stage7_continuous<- remove.iqr3.outliers(dat$cardio_sbp_stage7_continuous)
dat$cardio_sbp_stage7_zscore_continuous <- scale(dat$cardio_sbp_stage7_continuous)

#stage 8
dat$cardio_sbp_stage8_continuous <- dat$FJAR019a
dat$cardio_sbp_stage8_continuous[dat$FJAR019a<0]<-NA
dat$cardio_sbp_stage8_continuous<- remove.iqr3.outliers(dat$cardio_sbp_stage8_continuous)
dat$cardio_sbp_stage8_zscore_continuous <- scale(dat$cardio_sbp_stage8_continuous)

#stage 9
dat$cardio_sbp_stage9_continuous <- dat$FKBP1030
dat$cardio_sbp_stage9_continuous[dat$FKBP1030<0]<-NA
dat$cardio_sbp_stage9_continuous<- remove.iqr3.outliers(dat$cardio_sbp_stage9_continuous)
dat$cardio_sbp_stage9_zscore_continuous <- scale(dat$cardio_sbp_stage9_continuous)

#Diastolic BP
#stage 3
dat$cardio_dbp_stage3_continuous <- dat$f7sa022
dat$cardio_dbp_stage3_continuous[dat$f7sa022<0]<-NA
dat$cardio_dbp_stage3_continuous<- remove.iqr3.outliers(dat$cardio_dbp_stage3_continuous)
dat$cardio_dbp_stage3_zscore_continuous <- scale(dat$cardio_dbp_stage3_continuous)

#stage 4
dat$cardio_dbp_stage4_continuous <- dat$f9sa022
dat$cardio_dbp_stage4_continuous[dat$f9sa022<0]<-NA
dat$cardio_dbp_stage4_continuous<- remove.iqr3.outliers(dat$cardio_dbp_stage4_continuous)
dat$cardio_dbp_stage4_zscore_continuous <- scale(dat$cardio_dbp_stage4_continuous)

#stage 5
dat$cardio_dbp_stage5_continuous <- dat$fesa022
dat$cardio_dbp_stage5_continuous[dat$fesa022<0]<-NA
dat$cardio_dbp_stage5_continuous<- remove.iqr3.outliers(dat$cardio_dbp_stage5_continuous)
dat$cardio_dbp_stage5_zscore_continuous <- scale(dat$cardio_dbp_stage5_continuous)

#stage 6
dat$cardio_dbp_1_stage6_continuous<-NA
dat$cardio_dbp_1_stage6_continuous[which(dat$fg6121>=0)]<-dat$fg6121[which(dat$fg6121>=0)]
dat$cardio_dbp_2_stage6_continuous<-NA
dat$cardio_dbp_2_stage6_continuous[which(dat$fg6126>=0)]<-dat$fg6126[which(dat$fg6126>=0)]
dat$cardio_dbp_stage6_continuous<-NA
dat$cardio_dbp_stage6_continuous[is.na(dat$cardio_dbp_1_stage6_continuous)]<-dat$cardio_dbp_2_stage6_continuous[is.na(dat$cardio_dbp_1_stage6_continuous)]
dat$cardio_dbp_stage6_continuous[is.na(dat$cardio_dbp_2_stage6_continuous)]<-dat$cardio_dbp_1_stage6_continuous[is.na(dat$cardio_dbp_2_stage6_continuous)]
dat$cardio_dbp_stage6_continuous[which(!is.na(dat$cardio_dbp_2_stage6_continuous)&!is.na(dat$cardio_dbp_1_stage6_continuous))]<-(dat$cardio_dbp_1_stage6_continuous[which(!is.na(dat$cardio_dbp_2_stage6_continuous)&!is.na(dat$cardio_dbp_1_stage6_continuous))]+dat$cardio_dbp_2_stage6_continuous[which(!is.na(dat$cardio_dbp_2_stage6_continuous)&!is.na(dat$cardio_dbp_1_stage6_continuous))])/2
dat$cardio_dbp_stage6_continuous<- remove.iqr3.outliers(dat$cardio_dbp_stage6_continuous)
dat$cardio_dbp_stage6_zscore_continuous <- scale(dat$cardio_dbp_stage6_continuous)

#stage 7
dat$cardio_dbp_1_stage7_continuous<-NA
dat$cardio_dbp_1_stage7_continuous[which(dat$fh2031>=0)]<-dat$fh2031[which(dat$fh2031>=0)]
dat$cardio_dbp_2_stage7_continuous<-NA
dat$cardio_dbp_2_stage7_continuous[which(dat$fh2036>=0)]<-dat$fh2036[which(dat$fh2036>=0)]
dat$cardio_dbp_stage7_continuous<-NA
dat$cardio_dbp_stage7_continuous[is.na(dat$cardio_dbp_1_stage7_continuous)]<-dat$cardio_dbp_2_stage7_continuous[is.na(dat$cardio_dbp_1_stage7_continuous)]
dat$cardio_dbp_stage7_continuous[is.na(dat$cardio_dbp_2_stage7_continuous)]<-dat$cardio_dbp_1_stage7_continuous[is.na(dat$cardio_dbp_2_stage7_continuous)]
dat$cardio_dbp_stage7_continuous[which(!is.na(dat$cardio_dbp_2_stage7_continuous)&!is.na(dat$cardio_dbp_1_stage7_continuous))]<-(dat$cardio_dbp_1_stage7_continuous[which(!is.na(dat$cardio_dbp_2_stage7_continuous)&!is.na(dat$cardio_dbp_1_stage7_continuous))]+dat$cardio_dbp_2_stage7_continuous[which(!is.na(dat$cardio_dbp_2_stage7_continuous)&!is.na(dat$cardio_dbp_1_stage7_continuous))])/2
dat$cardio_dbp_stage7_continuous<- remove.iqr3.outliers(dat$cardio_dbp_stage7_continuous)
dat$cardio_dbp_stage7_zscore_continuous <- scale(dat$cardio_dbp_stage7_continuous)

#stage 8
dat$cardio_dbp_stage8_continuous <- dat$FJAR019b
dat$cardio_dbp_stage8_continuous[dat$FJAR019b<0]<-NA
dat$cardio_dbp_stage8_continuous<- remove.iqr3.outliers(dat$cardio_dbp_stage8_continuous)
dat$cardio_dbp_stage8_zscore_continuous <- scale(dat$cardio_dbp_stage8_continuous)

#stage 9
dat$cardio_dbp_stage9_continuous <- dat$FKBP1031
dat$cardio_dbp_stage9_continuous[dat$FKBP1031<0]<-NA
dat$cardio_dbp_stage9_continuous<- remove.iqr3.outliers(dat$cardio_dbp_stage9_continuous)
dat$cardio_dbp_stage9_zscore_continuous <- scale(dat$cardio_dbp_stage9_continuous)

#hypertension
dat$cardio_htn_stage9_binary<-NA
dat$cardio_htn_stage9_binary[dat$cardio_sbp_stage9_continuous>=140&dat$cardio_dbp_stage9_continuous>=90]<-1
dat$cardio_htn_stage9_binary[dat$cardio_sbp_stage9_continuous<140|dat$cardio_dbp_stage9_continuous<90]<-0

#isolated diastolic hypertension
dat$cardio_idh_stage9_binary<-NA
dat$cardio_idh_stage9_binary[dat$cardio_dbp_stage9_continuous>=90&dat$cardio_sbp_stage9_continuous<140]<-1
dat$cardio_idh_stage9_binary[dat$cardio_dbp_stage9_continuous<90]<-0

#pre-hypertensive
dat$cardio_prehtn_stage9_binary<-NA
dat$cardio_prehtn_stage9_binary[dat$cardio_sbp_stage9_continuous>=120&dat$cardio_dbp_stage9_continuous>=80]<-1
dat$cardio_prehtn_stage9_binary[dat$cardio_sbp_stage9_continuous<120|dat$cardio_dbp_stage9_continuous<80]<-0

#isolated systolic hypertension
dat$cardio_ish_stage9_binary<-NA
dat$cardio_ish_stage9_binary[dat$cardio_sbp_stage9_continuous>=160&dat$cardio_dbp_stage9_continuous<=90]<-1
dat$cardio_ish_stage9_binary[dat$cardio_sbp_stage9_continuous<160&dat$cardio_dbp_stage9_continuous<=90]<-0

#Triglycerides
#stage 3
dat$cardio_trig_stage3_continuous <- dat$TRIG_F7
dat$cardio_trig_stage3_continuous[dat$TRIG_F7<0]<-NA
dat$cardio_trig_stage3_continuous<- remove.iqr3.outliers(dat$cardio_trig_stage3_continuous)
dat$cardio_trig_stage3_zscore_continuous <- scale(dat$cardio_trig_stage3_continuous)

#stage 4
dat$cardio_trig_stage4_continuous <- dat$trig_f9
dat$cardio_trig_stage4_continuous[dat$trig_f9<0]<-NA
dat$cardio_trig_stage4_continuous<- remove.iqr3.outliers(dat$cardio_trig_stage4_continuous)
dat$cardio_trig_stage4_zscore_continuous <- scale(dat$cardio_trig_stage4_continuous)

#stage 7
dat$cardio_trig_stage7_continuous <- dat$trig_TF3
dat$cardio_trig_stage7_continuous[dat$trig_TF3<0]<-NA
dat$cardio_trig_stage7_continuous<- remove.iqr3.outliers(dat$cardio_trig_stage7_continuous)
dat$cardio_trig_stage7_zscore_continuous <- scale(dat$cardio_trig_stage7_continuous)

#stage 8
dat$cardio_trig_stage8_continuous <- dat$TRIG_TF4
dat$cardio_trig_stage8_continuous[dat$TRIG_TF4<0]<-NA
dat$cardio_trig_stage8_continuous<- remove.iqr3.outliers(dat$cardio_trig_stage8_continuous)
dat$cardio_trig_stage8_zscore_continuous <- scale(dat$cardio_trig_stage8_continuous)

#stage 9
dat$cardio_trig_stage9_continuous <- dat$Trig_F24
dat$cardio_trig_stage9_continuous[dat$Trig_F24<0]<-NA
dat$cardio_trig_stage9_continuous<- remove.iqr3.outliers(dat$cardio_trig_stage9_continuous)
dat$cardio_trig_stage9_zscore_continuous <- scale(dat$cardio_trig_stage9_continuous)

#HDL
#stage 3
dat$cardio_HDL_stage3_continuous <- dat$HDL_F7
dat$cardio_HDL_stage3_continuous[dat$HDL_F7<0]<-NA
dat$cardio_HDL_stage3_continuous<- remove.iqr3.outliers(dat$cardio_HDL_stage3_continuous)
dat$cardio_HDL_stage3_zscore_continuous <- scale(dat$cardio_HDL_stage3_continuous)

#stage 4
dat$cardio_HDL_stage4_continuous <- dat$HDL_f9
dat$cardio_HDL_stage4_continuous[dat$HDL_f9<0]<-NA
dat$cardio_HDL_stage4_continuous<- remove.iqr3.outliers(dat$cardio_HDL_stage4_continuous)
dat$cardio_HDL_stage4_zscore_continuous <- scale(dat$cardio_HDL_stage4_continuous)

#stage 7
dat$cardio_HDL_stage7_continuous <- dat$hdl_TF3
dat$cardio_HDL_stage7_continuous[dat$hdl_TF3<0]<-NA
dat$cardio_HDL_stage7_continuous<- remove.iqr3.outliers(dat$cardio_HDL_stage7_continuous)
dat$cardio_HDL_stage7_zscore_continuous <- scale(dat$cardio_HDL_stage7_continuous)

#stage 8
dat$cardio_HDL_stage8_continuous <- dat$HDL_TF4
dat$cardio_HDL_stage8_continuous[dat$HDL_TF4<0]<-NA
dat$cardio_HDL_stage8_continuous<- remove.iqr3.outliers(dat$cardio_HDL_stage8_continuous)
dat$cardio_HDL_stage8_zscore_continuous <- scale(dat$cardio_HDL_stage8_continuous)

#stage 9
dat$cardio_HDL_stage9_continuous <- dat$HDL_F24
dat$cardio_HDL_stage9_continuous[dat$HDL_F24<0]<-NA
dat$cardio_HDL_stage9_continuous<- remove.iqr3.outliers(dat$cardio_HDL_stage9_continuous)
dat$cardio_HDL_stage9_zscore_continuous <- scale(dat$cardio_HDL_stage9_continuous)

#LDL
#stage 3
dat$cardio_LDL_stage3_continuous <- dat$LDL_F7
dat$cardio_LDL_stage3_continuous[dat$LDL_F7<0]<-NA
dat$cardio_LDL_stage3_continuous<- remove.iqr3.outliers(dat$cardio_LDL_stage3_continuous)
dat$cardio_LDL_stage3_zscore_continuous <- scale(dat$cardio_LDL_stage3_continuous)

#stage 4
dat$cardio_LDL_stage4_continuous <- dat$LDL_f9
dat$cardio_LDL_stage4_continuous[dat$LDL_f9<0]<-NA
dat$cardio_LDL_stage4_continuous<- remove.iqr3.outliers(dat$cardio_LDL_stage4_continuous)
dat$cardio_LDL_stage4_zscore_continuous <- scale(dat$cardio_LDL_stage4_continuous)

#stage 7
dat$cardio_LDL_stage7_continuous <- dat$ldl_TF3
dat$cardio_LDL_stage7_continuous[dat$ldl_TF3<0]<-NA
dat$cardio_LDL_stage7_continuous<- remove.iqr3.outliers(dat$cardio_LDL_stage7_continuous)
dat$cardio_LDL_stage7_zscore_continuous <- scale(dat$cardio_LDL_stage7_continuous)

#stage 8
dat$cardio_LDL_stage8_continuous <- dat$LDL_TF4
dat$cardio_LDL_stage8_continuous[dat$LDL_TF4<0]<-NA
dat$cardio_LDL_stage8_continuous<- remove.iqr3.outliers(dat$cardio_LDL_stage8_continuous)
dat$cardio_LDL_stage8_zscore_continuous <- scale(dat$cardio_LDL_stage8_continuous)

#stage 9
dat$cardio_LDL_stage9_continuous <- dat$LDL_F24
dat$cardio_LDL_stage9_continuous[dat$LDL_F24<0]<-NA
dat$cardio_LDL_stage9_continuous<- remove.iqr3.outliers(dat$cardio_LDL_stage9_continuous)
dat$cardio_LDL_stage9_zscore_continuous <- scale(dat$cardio_LDL_stage9_continuous)

#Total cholesterol
#stage 3
dat$cardio_chol_stage3_continuous <- dat$CHOL_F7
dat$cardio_chol_stage3_continuous[dat$CHOL_F7<0]<-NA
dat$cardio_chol_stage3_continuous<- remove.iqr3.outliers(dat$cardio_chol_stage3_continuous)
dat$cardio_chol_stage3_zscore_continuous <- scale(dat$cardio_chol_stage3_continuous)

#stage 4
dat$cardio_chol_stage4_continuous <- dat$CHOL_F9
dat$cardio_chol_stage4_continuous[dat$CHOL_F9<0]<-NA
dat$cardio_chol_stage4_continuous<- remove.iqr3.outliers(dat$cardio_chol_stage4_continuous)
dat$cardio_chol_stage4_zscore_continuous <- scale(dat$cardio_chol_stage4_continuous)

#stage 7
dat$cardio_chol_stage7_continuous <- dat$chol_TF3
dat$cardio_chol_stage7_continuous[dat$chol_TF3<0]<-NA
dat$cardio_chol_stage7_continuous<- remove.iqr3.outliers(dat$cardio_chol_stage7_continuous)
dat$cardio_chol_stage7_zscore_continuous <- scale(dat$cardio_chol_stage7_continuous)

#stage8
dat$cardio_chol_stage8_continuous <- dat$CHOL_TF4
dat$cardio_chol_stage8_continuous[dat$CHOL_TF4<0]<-NA
dat$cardio_chol_stage8_continuous<- remove.iqr3.outliers(dat$cardio_chol_stage8_continuous)
dat$cardio_chol_stage8_zscore_continuous <- scale(dat$cardio_chol_stage8_continuous)

#stage 9
dat$cardio_chol_stage9_continuous <- dat$Chol_F24
dat$cardio_chol_stage9_continuous[dat$Chol_F24<0]<-NA
dat$cardio_chol_stage9_continuous<- remove.iqr3.outliers(dat$cardio_chol_stage9_continuous)
dat$cardio_chol_stage9_zscore_continuous <- scale(dat$cardio_chol_stage9_continuous)


#Glucose
#stage 3
dat$cardio_glucose_stage3_continuous <- dat$Glc_F7
dat$cardio_glucose_stage3_continuous[dat$Glc_F7<0]<-NA
dat$cardio_glucose_stage3_continuous<- remove.iqr3.outliers(dat$cardio_glucose_stage3_continuous)
dat$cardio_glucose_stage3_zscore_continuous <- scale(dat$cardio_glucose_stage3_continuous)

#stage 7
dat$cardio_glucose_stage7_continuous <- dat$Glc_TF3
dat$cardio_glucose_stage7_continuous[dat$Glc_TF3<0]<-NA
dat$cardio_glucose_stage7_continuous<- remove.iqr3.outliers(dat$cardio_glucose_stage7_continuous)
dat$cardio_glucose_stage7_zscore_continuous <- scale(dat$cardio_glucose_stage7_continuous)

#stage 8
dat$cardio_glucose_stage8_continuous <- dat$Glc_TF4
dat$cardio_glucose_stage8_continuous[dat$Glc_TF4<0]<-NA
dat$cardio_glucose_stage8_continuous<- remove.iqr3.outliers(dat$cardio_glucose_stage8_continuous)
dat$cardio_glucose_stage8_zscore_continuous <- scale(dat$cardio_glucose_stage8_continuous)

#stage 9
dat$cardio_glucose_stage9_continuous <- dat$Glc_F24
dat$cardio_glucose_stage9_continuous[dat$Glc_F24<0]<-NA
dat$cardio_glucose_stage9_continuous<- remove.iqr3.outliers(dat$cardio_glucose_stage9_continuous)
dat$cardio_glucose_stage9_zscore_continuous <- scale(dat$cardio_glucose_stage9_continuous)

#Insulin
#stage 4
dat$cardio_insulin_stage4_continuous <- dat$insulin_F9
dat$cardio_insulin_stage4_continuous[dat$insulin_F9<0]<-NA
dat$cardio_insulin_stage4_continuous<- remove.iqr3.outliers(dat$cardio_insulin_stage4_continuous)
dat$cardio_insulin_stage4_zscore_continuous <- scale(dat$cardio_insulin_stage4_continuous)

#stage 7
dat$cardio_insulin_stage7_continuous <- dat$insulin_TF3
dat$cardio_insulin_stage7_continuous[dat$insulin_TF3<0]<-NA
dat$cardio_insulin_stage7_continuous<- remove.iqr3.outliers(dat$cardio_insulin_stage7_continuous)
dat$cardio_insulin_stage7_zscore_continuous <- scale(dat$cardio_insulin_stage7_continuous)

#stage 8
dat$cardio_insulin_stage8_continuous <- dat$insulin_TF4
dat$cardio_insulin_stage8_continuous[dat$insulin_TF4<0]<-NA
dat$cardio_insulin_stage8_continuous<- remove.iqr3.outliers(dat$cardio_insulin_stage8_continuous)
dat$cardio_insulin_stage8_zscore_continuous <- scale(dat$cardio_insulin_stage8_continuous)

#stage 9
dat$cardio_insulin_stage9_continuous <- dat$Insulin_F24
dat$cardio_insulin_stage9_continuous[dat$Insulin_F24<0]<-NA
dat$cardio_insulin_stage9_continuous<- remove.iqr3.outliers(dat$cardio_insulin_stage9_continuous)
dat$cardio_insulin_stage9_zscore_continuous <- scale(dat$cardio_insulin_stage9_continuous)

#APO-A
#stage 3
dat$cardio_apoA_stage3_continuous <- dat$ApoA1_F7
dat$cardio_apoA_stage3_continuous[dat$ApoA1_F7<0]<-NA
dat$cardio_apoA_stage3_continuous<- remove.iqr3.outliers(dat$cardio_apoA_stage3_continuous)
dat$cardio_apoA_stage3_zscore_continuous <- scale(dat$cardio_apoA_stage3_continuous)

#stage 4
dat$cardio_apoA_stage4_continuous <- dat$ApoAI_f9
dat$cardio_apoA_stage4_continuous[dat$ApoAI_f9<0]<-NA
dat$cardio_apoA_stage4_continuous<- remove.iqr3.outliers(dat$cardio_apoA_stage4_continuous)
dat$cardio_apoA_stage4_zscore_continuous <- scale(dat$cardio_apoA_stage4_continuous)

#stage 7
dat$cardio_apoA_stage7_continuous <- dat$ApoA1_TF3
dat$cardio_apoA_stage7_continuous[dat$ApoA1_TF3<0]<-NA
dat$cardio_apoA_stage7_continuous<- remove.iqr3.outliers(dat$cardio_apoA_stage7_continuous)
dat$cardio_apoA_stage7_zscore_continuous <- scale(dat$cardio_apoA_stage7_continuous)

#stage 8
dat$cardio_apoA_stage8_continuous <- dat$ApoA1_TF4
dat$cardio_apoA_stage8_continuous[dat$ApoA1_TF4<0]<-NA
dat$cardio_apoA_stage8_continuous<- remove.iqr3.outliers(dat$cardio_apoA_stage8_continuous)
dat$cardio_apoA_stage8_zscore_continuous <- scale(dat$cardio_apoA_stage8_continuous)

#stage 9
dat$cardio_apoA_stage9_continuous <- dat$ApoA1_F24
dat$cardio_apoA_stage9_continuous[dat$ApoA1_F24<0]<-NA
dat$cardio_apoA_stage9_continuous<- remove.iqr3.outliers(dat$cardio_apoA_stage9_continuous)
dat$cardio_apoA_stage9_zscore_continuous <- scale(dat$cardio_apoA_stage9_continuous)

#APO-B
#stage 3
dat$cardio_apoB_stage3_continuous <- dat$ApoB_F7
dat$cardio_apoB_stage3_continuous[dat$ApoB_F7<0]<-NA
dat$cardio_apoB_stage3_continuous<- remove.iqr3.outliers(dat$cardio_apoB_stage3_continuous)
dat$cardio_apoB_stage3_zscore_continuous <- scale(dat$cardio_apoB_stage3_continuous)

#stage 4
dat$cardio_apoB_stage4_continuous <- dat$ApoB_f9
dat$cardio_apoB_stage4_continuous[dat$ApoB_f9<0]<-NA
dat$cardio_apoB_stage4_continuous<- remove.iqr3.outliers(dat$cardio_apoB_stage4_continuous)
dat$cardio_apoB_stage4_zscore_continuous <- scale(dat$cardio_apoB_stage4_continuous)

#stage 7
dat$cardio_apoB_stage7_continuous <- dat$ApoB_TF3
dat$cardio_apoB_stage7_continuous[dat$ApoB_TF3<0]<-NA
dat$cardio_apoB_stage7_continuous<- remove.iqr3.outliers(dat$cardio_apoB_stage7_continuous)
dat$cardio_apoB_stage7_zscore_continuous <- scale(dat$cardio_apoB_stage7_continuous)

#stage 8
dat$cardio_apoB_stage8_continuous <- dat$ApoB_TF4
dat$cardio_apoB_stage8_continuous[dat$ApoB_TF4<0]<-NA
dat$cardio_apoB_stage8_continuous<- remove.iqr3.outliers(dat$cardio_apoB_stage8_continuous)
dat$cardio_apoB_stage8_zscore_continuous <- scale(dat$cardio_apoB_stage8_continuous)

#stage 9
dat$cardio_apoB_stage9_continuous <- dat$ApoB_F24
dat$cardio_apoB_stage9_continuous[dat$ApoB_F24<0]<-NA
dat$cardio_apoB_stage9_continuous<- remove.iqr3.outliers(dat$cardio_apoB_stage9_continuous)
dat$cardio_apoB_stage9_zscore_continuous <- scale(dat$cardio_apoB_stage9_continuous)

#ApoB:A
#stage 3
dat$cardio_apoba_stage3_continuous<-dat$ApoBApoA1_F7
dat$cardio_apoba_stage3_continuous[dat$ApoBApoA1_F7<0]<-NA
dat$cardio_apoba_stage3_continuous<-remove.iqr3.outliers(dat$cardio_apoba_stage3_continuous)
dat$cardio_apoba_stage3_zscore_continuous<-scale(dat$cardio_apoba_stage3_continuous)

#stage 4
dat$cardio_apoba_stage4_continuous<-dat$cardio_apoB_stage4_continuous/dat$cardio_apoA_stage4_continuous
dat$cardio_apoba_stage4_continuous[dat$cardio_apoba_stage4_continuous<0]<-NA
dat$cardio_apoba_stage4_continuous<-remove.iqr3.outliers(dat$cardio_apoba_stage4_continuous)
dat$cardio_apoba_stage4_zscore_continuous<-scale(dat$cardio_apoba_stage4_continuous)

#stage 7
dat$cardio_apoba_stage7_continuous<-dat$ApoBApoA1_TF3
dat$cardio_apoba_stage7_continuous[dat$ApoBApoA1_TF3<0]<-NA
dat$cardio_apoba_stage7_continuous<-remove.iqr3.outliers(dat$cardio_apoba_stage7_continuous)
dat$cardio_apoba_stage7_zscore_continuous<-scale(dat$cardio_apoba_stage7_continuous)

#stage 8
dat$cardio_apoba_stage8_continuous<-dat$ApoBApoA1_TF4
dat$cardio_apoba_stage8_continuous[dat$ApoBApoA1_TF4<0]<-NA
dat$cardio_apoba_stage8_continuous<-remove.iqr3.outliers(dat$cardio_apoba_stage8_continuous)
dat$cardio_apoba_stage8_zscore_continuous<-scale(dat$cardio_apoba_stage8_continuous)

#stage 9
dat$cardio_apoba_stage9_continuous<-dat$ApoBApoA1_F24
dat$cardio_apoba_stage9_continuous[dat$ApoBApoA1_F24<0]<-NA
dat$cardio_apoba_stage9_continuous<-remove.iqr3.outliers(dat$cardio_apoba_stage9_continuous)
dat$cardio_apoba_stage9_zscore_continuous<-scale(dat$cardio_apoba_stage9_continuous)

#Glyoprotein acetyls
#stage 4
dat$cardio_GlyAce_stage4_continuous<-dat$Gp_F7
dat$cardio_GlyAce_stage4_continuous[which(dat$Gp_F7<0)]<-NA
dat$cardio_GlyAce_stage4_zscore_continuous<-scale(dat$cardio_GlyAce_stage4_continuous)
dat$cardio_GlyAce_stage4_no_outliers_continuous<-remove.iqr3.outliers(dat$cardio_GlyAce_stage4_continuous)
dat$cardio_GlyAce_stage4_no_outliers_zscore_continuous<-scale(dat$cardio_GlyAce_stage4_no_outliers_continuous)

#stage 7
dat$cardio_GlyAce_stage7_continuous<-dat$Gp_TF3
dat$cardio_GlyAce_stage7_continuous[which(dat$Gp_TF3<0)]<-NA
dat$cardio_GlyAce_stage7_zscore_continuous<-scale(dat$cardio_GlyAce_stage7_continuous)
dat$cardio_GlyAce_stage7_no_outliers_continuous<-remove.iqr3.outliers(dat$cardio_GlyAce_stage7_continuous)
dat$cardio_GlyAce_stage7_no_outliers_zscore_continuous<-scale(dat$cardio_GlyAce_stage7_no_outliers_continuous)

#stage 8
dat$cardio_GlyAce_stage8_continuous<-dat$Gp_TF4
dat$cardio_GlyAce_stage8_continuous[which(dat$Gp_TF4<0)]<-NA
dat$cardio_GlyAce_stage8_zscore_continuous<-scale(dat$cardio_GlyAce_stage8_continuous)
dat$cardio_GlyAce_stage8_no_outliers_continuous<-remove.iqr3.outliers(dat$cardio_GlyAce_stage8_continuous)
dat$cardio_GlyAce_stage8_no_outliers_zscore_continuous<-scale(dat$cardio_GlyAce_stage8_no_outliers_continuous)

#stage 9
dat$cardio_GlyAce_stage9_continuous<-dat$Gp_F24
dat$cardio_GlyAce_stage9_continuous[which(dat$Gp_F24<0)]<-NA
dat$cardio_GlyAce_stage9_zscore_continuous<-scale(dat$cardio_GlyAce_stage9_continuous)
dat$cardio_GlyAce_stage9_no_outliers_continuous<-remove.iqr3.outliers(dat$cardio_GlyAce_stage9_continuous)
dat$cardio_GlyAce_stage9_no_outliers_zscore_continuous<-scale(dat$cardio_GlyAce_stage9_no_outliers_continuous)

#cIMT
dat$cardio_rCIMT_stage9_continuous<-dat$FKCV1131
dat$cardio_rCIMT_stage9_continuous[which(dat$FKCV1131<0)]<-NA
dat$cardio_rCIMT_stage9_continuous<-remove.iqr3.outliers(dat$cardio_rCIMT_stage9_continuous)
dat$cardio_rCIMT_stage9_zscore_continuous<-scale(dat$cardio_rCIMT_stage9_continuous)

dat$cardio_lCIMT_stage9_continuous<-dat$FKCV2131
dat$cardio_lCIMT_stage9_continuous[which(dat$FKCV2131<0)]<-NA
dat$cardio_lCIMT_stage9_continuous<-remove.iqr3.outliers(dat$cardio_lCIMT_stage9_continuous)
dat$cardio_lCIMT_stage9_zscore_continuous<-scale(dat$cardio_lCIMT_stage9_continuous)

dat$cardio_CIMT_stage9_continuous[which(is.na(dat$cardio_rCIMT_stage9_continuous))]<-dat$cardio_lCIMT_stage9_continuous[which(is.na(dat$cardio_rCIMT_stage9_continuous))]
dat$cardio_CIMT_stage9_continuous[which(is.na(dat$cardio_lCIMT_stage9_continuous))]<-dat$cardio_rCIMT_stage9_continuous[which(is.na(dat$cardio_lCIMT_stage9_continuous))]
dat$cardio_CIMT_stage9_continuous[which(!is.na(dat$cardio_lCIMT_stage9_continuous)&!is.na(dat$cardio_rCIMT_stage9_continuous))]<-
  (dat$cardio_lCIMT_stage9_continuous[which(!is.na(dat$cardio_rCIMT_stage9_continuous)&!is.na(dat$cardio_lCIMT_stage9_continuous))]+dat$cardio_rCIMT_stage9_continuous[which(!is.na(dat$cardio_rCIMT_stage9_continuous)&!is.na(dat$cardio_lCIMT_stage9_continuous))])/2
dat$cardio_CIMT_stage9_zscore_continuous<-scale(dat$cardio_CIMT_stage9_continuous)

#pulse wave velocity
#stage 8
dat$cardio_PWV_stage8_continuous<-dat$FJAR083d
dat$cardio_PWV_stage8_continuous<-remove.iqr3.outliers(dat$cardio_PWV_stage8_continuous)
dat$cardio_PWV_stage8_zscore_continuous<-dat$cardio_PWV_stage8_continuous

#stage 9
dat$cardio_PWV_stage9_continuous<-dat$FKCV4200
dat$cardio_PWV_stage9_continuous[which(dat$FKCV4200<0)]<-NA
dat$cardio_PWV_stage9_continuous<-remove.iqr3.outliers(dat$cardio_PWV_stage9_continuous)
dat$cardio_PWV_stage9_zscore_continuous<-dat$cardio_PWV_stage9_continuous

#Mitral A wave peak
#stage 9
dat$cardio_Awave_stage9_continuous<-dat$FKEC5330
dat$cardio_Awave_stage9_continuous[which(dat$FKEC5330<0)]<-NA
dat$cardio_Awave_stage9_continuous<-remove.iqr3.outliers(dat$cardio_Awave_stage9_continuous)
dat$cardio_Awave_stage9_zscore_continuous<-scale(dat$cardio_Awave_stage9_continuous)

#Mitral E wave peak
#stage 8
dat$cardio_Ewave_stage8_continuous<-dat$FJGR114
dat$cardio_Ewave_stage8_continuous[which(dat$FJGR114<0)]<-NA
dat$cardio_Ewave_stage8_continuous<-remove.iqr3.outliers(dat$cardio_Ewave_stage8_continuous)
dat$cardio_Ewave_stage8_zscore_continuous<-scale(dat$cardio_Ewave_stage8_continuous)

#stage 9
dat$cardio_Ewave_stage9_continuous<-dat$FKEC5350
dat$cardio_Ewave_stage9_continuous[which(dat$FKEC5350<0)]<-NA
dat$cardio_Ewave_stage9_continuous<-remove.iqr3.outliers(dat$cardio_Ewave_stage9_continuous)
dat$cardio_Ewave_stage9_zscore_continuous<-scale(dat$cardio_Ewave_stage9_continuous)

#E/A
#stage 8
dat$cardio_EAratio_stage8_continuous<-dat$FJGR113
dat$cardio_EAratio_stage8_continuous[which(dat$FJGR113<0)]<-NA
dat$cardio_EAratio_stage8_continuous<-remove.iqr3.outliers(dat$cardio_EAratio_stage8_continuous)
dat$cardio_EAratio_stage8_zscore_continuous<-scale(dat$cardio_EAratio_stage8_continuous)

#stage 9
dat$cardio_EAratio_stage9_continuous<-dat$cardio_Ewave_stage9_continuous/dat$cardio_Awave_stage9_continuous
dat$cardio_EAratio_stage9_continuous[is.na(dat$cardio_Ewave_stage9_continuous)]<-NA
dat$cardio_EAratio_stage9_continuous[is.na(dat$cardio_Awave_stage9_continuous)]<-NA
dat$cardio_EAratio_stage9_continuous<-remove.iqr3.outliers(dat$cardio_EAratio_stage9_continuous)
dat$cardio_EAratio_stage9_zscore_continuous<-scale(dat$cardio_EAratio_stage9_continuous)

#Derive e' =(e'lateral+e' septal)/2
#stage 8
dat$cardio_elat_stage8_continuous<-dat$FJGR089
dat$cardio_elat_stage8_continuous[which(dat$FJGR089<0)]<-NA
dat$cardio_elat_stage8_zscore_continuous<-scale(dat$cardio_elat_stage8_continuous)
dat$cardio_esep_stage8_continuous<-dat$FJGR101
dat$cardio_esep_stage8_continuous[which(dat$FJGR101<0)]<-NA
dat$cardio_esep_stage8_zscore_continuous<-scale(dat$cardio_esep_stage8_continuous)
dat$cardio_e_stage8_continuous<-(dat$cardio_elat_stage8_continuous+dat$cardio_esep_stage8_continuous)/2
dat$cardio_e_stage8_continuous[is.na(dat$cardio_elat_stage8_continuous)]<-NA
dat$cardio_e_stage8_continuous[is.na(dat$cardio_esep_stage8_continuous)]<-NA
dat$cardio_e_stage8_continuous<-remove.iqr3.outliers(dat$cardio_e_stage8_continuous)
dat$cardio_e_stage8_zscore_continuous<-scale(dat$cardio_e_stage8_continuous)

#stage 9
dat$cardio_elat_stage9_continuous<-dat$FKEC5310
dat$cardio_elat_stage9_continuous[which(dat$FKEC5310<0)]<-NA
dat$cardio_elat_stage9_zscore_continuous<-scale(dat$cardio_elat_stage9_continuous)
dat$cardio_esep_stage9_continuous<-dat$FKEC5250
dat$cardio_esep_stage9_continuous[which(dat$FKEC5250<0)]<-NA
dat$cardio_esep_stage9_zscore_continuous<-scale(dat$cardio_esep_stage9_continuous)
dat$cardio_e_stage9_continuous<-(dat$cardio_elat_stage9_continuous+dat$cardio_esep_stage9_continuous)/2
dat$cardio_e_stage9_continuous[is.na(dat$cardio_elat_stage9_continuous)]<-NA
dat$cardio_e_stage9_continuous[is.na(dat$cardio_esep_stage9_continuous)]<-NA
dat$cardio_e_stage9_continuous<-remove.iqr3.outliers(dat$cardio_e_stage9_continuous)
dat$cardio_e_stage9_zscore_continuous<-scale(dat$cardio_e_stage9_continuous)

#E/e'
#stage 8
dat$cardio_Ee_stage8_continuous<-dat$cardio_Ewave_stage8_continuous/dat$cardio_e_stage8_continuous
dat$cardio_Ee_stage8_continuous[is.na(dat$cardio_Ewave_stage8_continuous)]<-NA
dat$cardio_Ee_stage8_continuous[which(dat$cardio_Ewave_stage8_continuous<0)]<-NA
dat$cardio_Ee_stage8_continuous[is.na(dat$cardio_e_stage8_continuous)]<-NA
dat$cardio_Ee_stage8_continuous[which(dat$cardio_e_stage8_continuous<0)]<-NA
dat$cardio_Ee_stage8_continuous<-remove.iqr3.outliers(dat$cardio_Ee_stage8_continuous)
dat$cardio_Ee_stage8_zscore_continuous<-scale(dat$cardio_Ee_stage8_continuous)

#stage 9
dat$cardio_Ee_stage9_continuous<-dat$cardio_Ewave_stage9_continuous/dat$cardio_e_stage9_continuous
dat$cardio_Ee_stage9_continuous[is.na(dat$cardio_Ewave_stage9_continuous)]<-NA
dat$cardio_Ee_stage9_continuous[which(dat$cardio_Ewave_stage9_continuous<0)]<-NA
dat$cardio_Ee_stage9_continuous[is.na(dat$cardio_e_stage9_continuous)]<-NA
dat$cardio_Ee_stage9_continuous[which(dat$cardio_e_stage9_continuous<0)]<-NA
dat$cardio_Ee_stage9_continuous<-remove.iqr3.outliers(dat$cardio_Ee_stage9_continuous)
dat$cardio_Ee_stage9_zscore_continuous<-scale(dat$cardio_Ee_stage9_continuous)

#LAD
#stage 8
dat$cardio_LAD_stage8_continuous<-dat$FJGR039
dat$cardio_LAD_stage8_continuous[which(dat$FJGR039<0)]<-NA
dat$cardio_LAD_stage8_continuous<-remove.iqr3.outliers(dat$cardio_LAD_stage8_continuous)
dat$cardio_LAD_stage8_zscore_continuous<-scale(dat$cardio_LAD_stage8_continuous)

#stage 9
dat$cardio_LAD_stage9_continuous<-dat$FKEC5070
dat$cardio_LAD_stage9_continuous[which(dat$FKEC5070<0)]<-NA
dat$cardio_LAD_stage9_continuous<-remove.iqr3.outliers(dat$cardio_LAD_stage9_continuous)
dat$cardio_LAD_stage9_zscore_continuous<-scale(dat$cardio_LAD_stage9_continuous)

#IVSd
#stage 8
dat$cardio_IVSd_stage8_continuous<-dat$FJGR031
dat$cardio_IVSd_stage8_continuous[which(dat$FJGR031<0)]<-NA
dat$cardio_IVSd_stage8_continuous<-remove.iqr3.outliers(dat$cardio_IVSd_stage8_continuous)
dat$cardio_IVSd_stage8_zscore_continuous<-scale(dat$cardio_IVSd_stage8_continuous)

#stage 9
dat$cardio_IVSd_stage9_continuous<-dat$FKEC5050
dat$cardio_IVSd_stage9_continuous[which(dat$FKEC5050<0)]<-NA
dat$cardio_IVSd_stage9_continuous<-remove.iqr3.outliers(dat$cardio_IVSd_stage9_continuous)
dat$cardio_IVSd_stage9_zscore_continuous<-scale(dat$cardio_IVSd_stage9_continuous)

#PWd
#stage 8
dat$cardio_PWd_stage8_continuous<-dat$FJGR053
dat$cardio_PWd_stage8_continuous[which(dat$FJGR053<0)]<-NA
dat$cardio_PWd_stage8_continuous<-remove.iqr3.outliers(dat$cardio_PWd_stage8_continuous)
dat$cardio_PWd_stage8_zscore_continuous<-scale(dat$cardio_PWd_stage8_continuous)

#stage 9
dat$cardio_PWd_stage9_continuous<-dat$FKEC5180
dat$cardio_PWd_stage9_continuous[which(dat$FKEC5180<0)]<-NA
dat$cardio_PWd_stage9_continuous<-remove.iqr3.outliers(dat$cardio_PWd_stage9_continuous)
dat$cardio_PWd_stage9_zscore_continuous<-scale(dat$cardio_PWd_stage9_continuous)

#LVIDd
#stage 8
dat$cardio_LVIDd_stage8_continuous<-dat$FJGR043
dat$cardio_LVIDd_stage8_continuous[which(dat$FJGR043<0)]<-NA
dat$cardio_LVIDd_stage8_continuous<-remove.iqr3.outliers(dat$cardio_LVIDd_stage8_continuous)
dat$cardio_LVIDd_stage8_zscore_continuous<-scale(dat$cardio_LVIDd_stage8_continuous)

#stage 9
dat$cardio_LVIDd_stage9_continuous<-dat$FKEC5080
dat$cardio_LVIDd_stage9_continuous[which(dat$FKEC5080<0)]<-NA
dat$cardio_LVIDd_stage9_continuous<-remove.iqr3.outliers(dat$cardio_LVIDd_stage9_continuous)
dat$cardio_LVIDd_stage9_zscore_continuous<-scale(dat$cardio_LVIDd_stage9_continuous)

#LVIDs
#stage 8
dat$cardio_LVIDs_stage8_continuous<-dat$FJGR047
dat$cardio_LVIDs_stage8_continuous[which(dat$FJGR047<0)]<-NA
dat$cardio_LVIDs_stage8_continuous<-remove.iqr3.outliers(dat$cardio_LVIDs_stage8_continuous)
dat$cardio_LVIDs_stage8_zscore_continuous<-scale(dat$cardio_LVIDs_stage8_continuous)

#stage 9
dat$cardio_LVIDs_stage9_continuous<-dat$FKEC5090
dat$cardio_LVIDs_stage9_continuous[which(dat$FKEC5090<0)]<-NA
dat$cardio_LVIDs_stage9_continuous<-remove.iqr3.outliers(dat$cardio_LVIDs_stage9_continuous)
dat$cardio_LVIDs_stage9_zscore_continuous<-scale(dat$cardio_LVIDs_stage9_continuous)

#Ejection fraction - Teicholz
#Derive Vs and Vd =[7/(2.4 + LVIDs)] * LVIDs^3 and =[7/(2.4 + LVIDd)] * LVIDd^3
#stage 8
dat$cardio_Vs_stage8_continuous<-NA
dat$cardio_Vs_stage8_continuous<- (7/(2.4 + dat$cardio_LVIDs_stage8_continuous))* (dat$cardio_LVIDs_stage8_continuous^3)
dat$cardio_Vs_stage8_continuous[is.na(dat$cardio_LVIDs_stage8_continuous)]<-NA
dat$cardio_Vs_stage8_continuous<-remove.iqr3.outliers(dat$cardio_Vs_stage8_continuous)
dat$cardio_Vs_stage8_zscore_continuous<-scale(dat$cardio_Vs_stage8_continuous)

dat$cardio_Vd_stage8_continuous<-NA
dat$cardio_Vd_stage8_continuous<- (7/(2.4 + dat$cardio_LVIDd_stage8_continuous)) * (dat$cardio_LVIDd_stage8_continuous^3)
dat$cardio_Vd_stage8_continuous[is.na(dat$cardio_LVIDd_stage8_continuous)]<-NA
dat$cardio_Vd_stage8_continuous<-remove.iqr3.outliers(dat$cardio_Vd_stage8_continuous)
dat$cardio_Vd_stage8_zscore_continuous<-scale(dat$cardio_Vd_stage8_continuous)

#stage 9
dat$cardio_Vs_stage9_continuous<-NA
dat$cardio_Vs_stage9_continuous<- (7/(2.4 + dat$cardio_LVIDs_stage9_continuous))* (dat$cardio_LVIDs_stage9_continuous^3)
dat$cardio_Vs_stage9_continuous[is.na(dat$cardio_LVIDs_stage9_continuous)]<-NA
dat$cardio_Vs_stage9_continuous<-remove.iqr3.outliers(dat$cardio_Vs_stage9_continuous)
dat$cardio_Vs_stage9_zscore_continuous<-scale(dat$cardio_Vs_stage9_continuous)

dat$cardio_Vd_stage9_continuous<-NA
dat$cardio_Vd_stage9_continuous<- (7/(2.4 + dat$cardio_LVIDd_stage9_continuous)) * (dat$cardio_LVIDd_stage9_continuous^3)
dat$cardio_Vd_stage9_continuous<-remove.iqr3.outliers(dat$cardio_Vd_stage9_continuous)
dat$cardio_Vd_stage9_continuous[is.na(dat$cardio_LVIDd_stage9_continuous)]<-NA
dat$cardio_Vd_stage9_zscore_continuous<-scale(dat$cardio_Vd_stage9_continuous)

#Derive EF =(Vd-Vs)/Vd*100
#stage 8
dat$cardio_ef_stage8_continuous<-(dat$cardio_Vd_stage8_continuous-dat$cardio_Vs_stage8_continuous)/dat$cardio_Vd_stage8_continuous*100
dat$cardio_ef_stage8_continuous[is.na(dat$cardio_Vd_stage8_continuous)]<-NA
dat$cardio_ef_stage8_continuous[is.na(dat$cardio_Vs_stage8_continuous)]<-NA
dat$cardio_ef_stage8_continuous<-remove.iqr3.outliers(dat$cardio_ef_stage8_continuous)
dat$cardio_ef_stage8_zscore_continuous<-scale(dat$cardio_ef_stage8_continuous)

#stage 9
dat$cardio_ef_stage9_continuous<-(dat$cardio_Vd_stage9_continuous-dat$cardio_Vs_stage9_continuous)/dat$cardio_Vd_stage9_continuous*100
dat$cardio_ef_stage9_continuous[is.na(dat$cardio_Vd_stage9_continuous)]<-NA
dat$cardio_ef_stage9_continuous[is.na(dat$cardio_Vs_stage9_continuous)]<-NA
dat$cardio_ef_stage9_continuous<-remove.iqr3.outliers(dat$cardio_ef_stage9_continuous)
dat$cardio_ef_stage9_zscore_continuous<-scale(dat$cardio_ef_stage9_continuous)

#Fractional shortening =(LVIDd-LVIDs)/LVIDd*100
#stage 8
dat$cardio_FS_stage8_continuous<-(dat$cardio_LVIDd_stage8_continuous-dat$cardio_LVIDs_stage8_continuous)/dat$cardio_LVIDd_stage8_continuous*100
dat$cardio_FS_stage8_continuous[is.na(dat$cardio_LVIDd_stage8_continuous)]<-NA
dat$cardio_FS_stage8_continuous[is.na(dat$cardio_LVIDs_stage8_continuous)]<-NA
dat$cardio_FS_stage8_continuous<-remove.iqr3.outliers(dat$cardio_FS_stage8_continuous)
dat$cardio_FS_stage8_zscore_continuous<-scale(dat$cardio_FS_stage8_continuous)

#stage 9
dat$cardio_FS_stage9_continuous<-(dat$cardio_LVIDd_stage9_continuous-dat$cardio_LVIDs_stage9_continuous)/dat$cardio_LVIDd_stage9_continuous*100
dat$cardio_FS_stage9_continuous[is.na(dat$cardio_LVIDd_stage9_continuous)]<-NA
dat$cardio_FS_stage9_continuous[is.na(dat$cardio_LVIDs_stage9_continuous)]<-NA
dat$cardio_FS_stage9_continuous<-remove.iqr3.outliers(dat$cardio_FS_stage9_continuous)
dat$cardio_FS_stage9_zscore_continuous<-scale(dat$cardio_FS_stage9_continuous)

#LAP =(E/e´ x 1.25) + 1.9
#stage 8
dat$cardio_lap_stage8_continuous<-(dat$cardio_Ee_stage8_continuous*1.25) + 1.9
dat$cardio_lap_stage8_continuous<-remove.iqr3.outliers(dat$cardio_lap_stage8_continuous)
dat$cardio_lap_stage8_zscore_continuous<-scale(dat$cardio_lap_stage8_continuous)

#stage 9
dat$cardio_lap_stage9_continuous<-(dat$cardio_Ee_stage9_continuous*1.25) + 1.9
dat$cardio_lap_stage9_continuous<-remove.iqr3.outliers(dat$cardio_lap_stage9_continuous)
dat$cardio_lap_stage9_zscore_continuous<-scale(dat$cardio_lap_stage9_continuous)

#LVMI- Deveruex and Reichek = 0.8*104*([IVSd+LVIDd+PWd]^3 - LVIDd^3)+0.6
#stage 8
dat$cardio_lvm_stage8_continuous<- 0.8*1.04*((dat$cardio_IVSd_stage8_continuous+dat$cardio_LVIDd_stage8_continuous+dat$cardio_PWd_stage8_continuous)^3 - dat$cardio_LVIDd_stage8_continuous^3)+0.6
dat$cardio_lvm_stage8_continuous[is.na(dat$cardio_IVSd_stage8_continuous)]<-NA
dat$cardio_lvm_stage8_continuous[is.na(dat$cardio_LVIDd_stage8_continuous)]<-NA
dat$cardio_lvm_stage8_continuous[is.na(dat$cardio_PWd_stage8_continuous)]<-NA
dat$cardio_lvm_stage8_zscore_continuous<-scale(dat$cardio_lvm_stage8_continuous)
dat$cardio_lvmi_stage8_continuous<-dat$cardio_lvm_stage8_continuous/(dat$anthro_height_stage8_continuous/100)^2.7
dat$cardio_lvmi_stage8_continuous[is.na(dat$anthro_height_stage8_continuous)]<-NA
dat$cardio_lvmi_stage8_continuous[is.na(dat$cardio_lvm_stage8_continuous)]<-NA
dat$cardio_lvmi_stage8_continuous<-remove.iqr3.outliers(dat$cardio_lvmi_stage8_continuous)
dat$cardio_lvmi_stage8_zscore_continuous<-scale(dat$cardio_lvmi_stage8_continuous)

#stage 9
dat$cardio_lvm_stage9_continuous<- 0.8*1.04*((dat$cardio_IVSd_stage9_continuous+dat$cardio_LVIDd_stage9_continuous+dat$cardio_PWd_stage9_continuous)^3 - dat$cardio_LVIDd_stage9_continuous^3)+0.6
dat$cardio_lvm_stage9_continuous[is.na(dat$cardio_IVSd_stage9_continuous)]<-NA
dat$cardio_lvm_stage9_continuous[is.na(dat$cardio_LVIDd_stage9_continuous)]<-NA
dat$cardio_lvm_stage9_continuous[is.na(dat$cardio_PWd_stage9_continuous)]<-NA
dat$cardio_lvm_stage9_zscore_continuous<-scale(dat$cardio_lvm_stage9_continuous)
dat$cardio_lvmi_stage9_continuous<-dat$cardio_lvm_stage9_continuous/(dat$anthro_height_stage9_continuous/100)^2.7
dat$cardio_lvmi_stage9_continuous[is.na(dat$anthro_height_stage9_continuous)]<-NA
dat$cardio_lvmi_stage9_continuous[is.na(dat$cardio_lvm_stage9_continuous)]<-NA
dat$cardio_lvmi_stage9_continuous<-remove.iqr3.outliers(dat$cardio_lvmi_stage9_continuous)
dat$cardio_lvmi_stage9_zscore_continuous<-scale(dat$cardio_lvmi_stage9_continuous)

#####Extract negative controls####
#set-up meta-data
negcon_vars <- c("g572","g572a","h462","h462a","g571","g571a","h461","h461a","kk140a")

negcon_meta_data <- findVars(unique(negcon_vars)) #create a table wof meta data for each variable

negcon_meta_data<- subset(negcon_meta_data,subset=tolower(name) %in% negcon_vars) #to ensure only variables of interest are extracted and not variables with longer names that include the same strings

#extractvariables
negcon_data <- extractVars(negcon_meta_data[negcon_meta_data$name %in% unique(negcon_vars),]) #extracting data and perform withdrawal of consent

negcon_data<-negcon_data[,c("aln","qlet",negcon_vars)]

dat <- left_join(dat,negcon_data,by=c("aln","qlet")) 

#pigeon infestation
dat$negcon_pigeons_binary <-NA
dat$negcon_pigeons_binary[which(dat$g572a==2|dat$h462a==2)]<-0
dat$negcon_pigeons_binary[which(dat$g572a==1|dat$h462a==1)]<-1

#mouse infestation
dat$negcon_mice_binary <-NA
dat$negcon_mice_binary[which(dat$g571a==2|dat$h461a==2)]<-0
dat$negcon_mice_binary[which(dat$g571a==1|dat$h461a==1)]<-1

#sting
dat$negcon_sting_ever_binary<-NA
dat$negcon_sting_ever_binary[which(dat$kk140a==2)]<-0
dat$negcon_sting_ever_binary[which(dat$kk140a==1)]<-1

####Polygenic risk scores####
#Read in ancestry informed principal components
pcs<-read.table() #file path to ALSPAC principal components file
pcs$aln <-as.numeric(substr(pcs$V1, 1, 5))
pcs$qlet<-substr(pcs$V1, 6, 6)
pcs[,c(1:2)]<-NULL
colnames(pcs)<-c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9",'PC10',
                 "PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20",
                 "aln","qlet")
dat<-left_join(dat,pcs, by =c("aln", "qlet"))


#Read in PRS generated in PRSice
#Read in PRSs with no MHC SNPs
prsMHC1<-fread() #file path to PRS scores
prsMHC1<-prsMHC1[,c(1,2,3,4,6,8,10,12)]
prsnoMHC2<-fread() #file path to PRS scores
names(prsnoMHC2)
prsnoMHC2<-prsnoMHC2[,c(1,2,3,12,102)]
prsnoMHC3<-fread() #file path to PRS scores
names(prsnoMHC3)
prsnoMHC3<-prsnoMHC3[,c(1,2,3,12)]
prsnoMHC4<-fread() #file path to PRS scores
names(prsnoMHC4)
prsnoMHC4<-prsnoMHC4[,c(1,2,3,4)]

prsnoMHC<-left_join(prsnoMHC2,prsnoMHC3,by=c("FID","IID"))
prsnoMHC<-left_join(prsnoMHC,prsnoMHC4,by=c("FID","IID"))
prsnoMHC<-prsnoMHC[,c(1,2,5,4,3,7,6,9,8)]
prsnoMHC$aln <-as.numeric(substr(prsnoMHC$FID, 1, 5))
prsnoMHC$qlet<-substr(prsnoMHC$FID, 6, 6)
prsnoMHC<-data.frame(prsnoMHC)
prsnoMHC<-prsnoMHC[,-which(names(prsnoMHC) %in% c("FID","IID"))]

#Read in PRSs with 1 MHC SNP scores
prsMHC1<-fread() #file path to PRS scores
prsMHC1<-prsMHC1[,c(1,2,3,4,6,8,10,12)]
prsMHC2<-fread() #file path to PRS scores
names(prsMHC2)
prsMHC2<-prsMHC2[,c(1,2,3,12,102)]
prsMHC3<-fread() #file path to PRS scores
names(prsMHC3)
prsMHC3<-prsMHC3[,c(1,2,3,12)]
prsMHC4<-fread() #file path to PRS scores
names(prsMHC4)
prsMHC4<-prsMHC4[,c(1,2,3,4)]

prsMHC<-left_join(prsMHC2,prsMHC3,by=c("FID","IID"))
prsMHC<-left_join(prsMHC,prsMHC4,by=c("FID","IID"))
prsMHC<-prsMHC[,c(1,2,5,4,3,7,6,9,8)]

prsMHC$aln <-as.numeric(substr(prsMHC$FID, 1, 5))
prsMHC$qlet<-substr(prsMHC$FID, 6, 6)
prsMHC<-data.frame(prsMHC)
prsMHC<-prsMHC[,-which(names(prsMHC) %in% c("FID","IID"))]

names(prsMHC)<-c("Pt_0.01MHC","Pt_0.001MHC","Pt_0.0001MHC","Pt_1e.05MHC","Pt_1e.06MHC","Pt_1e.07MHC","Pt_5e.08MHC","aln","qlet")

dat<-left_join(dat,prsnoMHC, by =c("aln","qlet"))
dat<-left_join(dat,prsMHC, by =c("aln","qlet"))

####Remove original variables from dataframe####
dat<-dat[,-which(names(dat) %in% time_vars)]
dat<-dat[,-which(names(dat)%in%c("f9003c","f7003c","fe003c","fg0011a","fh0011a","FKAR0010"))]
dat<- dat[,-which(names(dat) %in% names(ages)[3:nrow(ages)])] 
dat<-dat[,-which(colnames(dat)%in%child_anthropometry_vars)]
dat<-dat[,-which(colnames(dat)%in%immuno_vars)]
dat<-dat[,-which(colnames(dat)%in%negcon_vars)]
dat<-dat[,-which(colnames(dat)%in%cardiomet_vars)]

####Save data file####
saveRDS(dat,()) #add file path


####Data cleaning#####
#recode covariate sex
dat$covs_sex_binary[which(dat$covs_sex_binary=="female")]<-1
dat$covs_sex_binary[which(dat$covs_sex_binary=="male")]<-0
dat$covs_sex_binary<-as.numeric(dat$covs_sex_binary)

#Make covariate ages numeric
dat$covs_age_child_stage9_f24_continuous<-as.numeric(dat$covs_age_child_stage9_f24_continuous)
dat$covs_age_child_stage8_tf4_continuous<-as.numeric(dat$covs_age_child_stage8_tf4_continuous)
dat$covs_age_child_stage7_tf3_continuous<-as.numeric(dat$covs_age_child_stage7_tf3_continuous)
dat$covs_age_child_stage6_tf2_continuous<-as.numeric(dat$covs_age_child_stage6_tf2_continuous)
dat$covs_age_child_stage5_f11_continuous<-as.numeric(dat$covs_age_child_stage5_f11_continuous)
dat$covs_age_child_stage4_f9_continuous<-as.numeric(dat$covs_age_child_stage4_f9_continuous)
dat$covs_age_child_stage3_f7_continuous<-as.numeric(dat$covs_age_child_stage3_f7_continuous)
dat$covs_age_child_stage2_cif_continuous<-as.numeric(dat$covs_age_child_stage2_cif_continuous)
dat$covs_age_child_stage1_cif_continuous<-as.numeric(dat$covs_age_child_stage1_cif_continuous)

#Create PRS Z-scores
#PRS with no MHC SNPs
dat$Pt_0.0001_zscore_continuous<-scale(dat$Pt_0.0001)
dat$Pt_0.001_zscore_continuous<-scale(dat$Pt_0.001)
dat$Pt_0.01_zscore_continuous<-scale(dat$Pt_0.01)
dat$Pt_1e.06_zscore_continuous<-scale(dat$Pt_1e.06)
dat$Pt_1e.07_zscore_continuous<-scale(dat$Pt_1e.07)
dat$Pt_5e.08_zscore_continuous<-scale(dat$Pt_5e.08)
dat$Pt_1e.05_zscore_continuous<-scale(dat$Pt_1e.05)

#PRS with MHC SNPs
dat$Pt_0.0001MHC_zscore_continuous<-scale(dat$Pt_0.0001MHC)
dat$Pt_0.001MHC_zscore_continuous<-scale(dat$Pt_0.001MHC)
dat$Pt_0.01MHC_zscore_continuous<-scale(dat$Pt_0.01MHC)
dat$Pt_1e.06MHC_zscore_continuous<-scale(dat$Pt_1e.06MHC)
dat$Pt_1e.07MHC_zscore_continuous<-scale(dat$Pt_1e.07MHC)
dat$Pt_5e.08MHC_zscore_continuous<-scale(dat$Pt_5e.08MHC)
dat$Pt_1e.05MHC_zscore_continuous<-scale(dat$Pt_1e.05MHC)

#Create Z scores for BMI
dat$anthro_bmi_stage1_zscore_continuous<-scale(dat$anthro_bmi_stage1_continuous)
dat$anthro_bmi_stage2_zscore_continuous<-scale(dat$anthro_bmi_stage2_continuous)
dat$anthro_bmi_stage3_zscore_continuous<-scale(dat$anthro_bmi_stage3_continuous)
dat$anthro_bmi_stage4_zscore_continuous<-scale(dat$anthro_bmi_stage4_continuous)
dat$anthro_bmi_stage5_zscore_continuous<-scale(dat$anthro_bmi_stage5_continuous)
dat$anthro_bmi_stage6_zscore_continuous<-scale(dat$anthro_bmi_stage6_continuous)
dat$anthro_bmi_stage7_zscore_continuous<-scale(dat$anthro_bmi_stage7_continuous)
dat$anthro_bmi_stage8_zscore_continuous<-scale(dat$anthro_bmi_stage8_continuous)
dat$anthro_bmi_stage9_zscore_continuous<-scale(dat$anthro_bmi_stage9_continuous)

#convert fmi to kg/m2 from g/m2
#stage 4
dat$anthro_fmi_stage4_continuous<-dat$anthro_fmi_stage4_continuous/1000
dat$anthro_fmi_stage4_zscore_continuous<-scale(dat$anthro_fmi_stage4_continuous)

#stage 5
dat$anthro_fmi_stage5_continuous<-dat$anthro_fmi_stage5_continuous/1000
dat$anthro_fmi_stage5_zscore_continuous<-scale(dat$anthro_fmi_stage5_continuous)

#stage 6
dat$anthro_fmi_stage6_continuous<-dat$anthro_fmi_stage6_continuous/1000
dat$anthro_fmi_stage6_zscore_continuous<-scale(dat$anthro_fmi_stage6_continuous)

#stage 8
dat$anthro_fmi_stage8_continuous<-dat$anthro_fmi_stage8_continuous/1000
dat$anthro_fmi_stage8_zscore_continuous<-scale(dat$anthro_fmi_stage8_continuous)

#stage 9
dat$anthro_fmi_stage9_continuous<-dat$anthro_fmi_stage9_continuous/1000
dat$anthro_fmi_stage9_zscore_continuous<-scale(dat$anthro_fmi_stage9_continuous)

#Convert stage 9 waist circumference to cm
dat$anthro_waist_stage9_continuous<-dat$anthro_waist_stage9_continuous/10
dat$anthro_waist_stage9_zscore_continuous<-scale(dat$anthro_waist_stage9_continuous)

#Convert stage 4 Apo A and B to g/l from mg/dl
dat$cardio_apoA_stage4_continuous<-dat$cardio_apoA_stage4_continuous/100
dat$cardio_apoA_stage4_zscore_continuous<-scale(dat$cardio_apoA_stage4_continuous)
dat$cardio_apoB_stage4_continuous<-dat$cardio_apoB_stage4_continuous/100
dat$cardio_apoB_stage4_zscore_continuous<-scale(dat$cardio_apoB_stage4_continuous)

#generate combined cIMT
dat$cardio_CIMT_stage9_continuous[which(is.na(dat$cardio_rCIMT_stage9_continuous))]<-dat$cardio_lCIMT_stage9_continuous[which(is.na(dat$cardio_rCIMT_stage9_continuous))]
dat$cardio_CIMT_stage9_continuous[which(is.na(dat$cardio_lCIMT_stage9_continuous))]<-dat$cardio_rCIMT_stage9_continuous[which(is.na(dat$cardio_lCIMT_stage9_continuous))]
dat$cardio_CIMT_stage9_continuous[which(!is.na(dat$cardio_lCIMT_stage9_continuous)&!is.na(dat$cardio_rCIMT_stage9_continuous))]<-
  (dat$cardio_lCIMT_stage9_continuous[which(!is.na(dat$cardio_rCIMT_stage9_continuous)&!is.na(dat$cardio_lCIMT_stage9_continuous))]+dat$cardio_rCIMT_stage9_continuous[which(!is.na(dat$cardio_rCIMT_stage9_continuous)&!is.na(dat$cardio_lCIMT_stage9_continuous))])/2
dat$cardio_CIMT_stage9_zscore_continuous<-scale(dat$cardio_CIMT_stage9_continuous)

#Remove participants with no PRS scores and who didn't attend Focus24 clinic
dat<-dat[!is.na(dat$Pt_1e.05),] 
dat<-dat[!is.na(dat$covs_age_child_stage9_f24_continuous),]

#Change insulin from U/L to pmol/mL
dat$cardio_insulin_pmolL_stage9_continuous<-dat$cardio_insulin_stage9_continuous*6
dat$cardio_insulin_pmolL_stage8_continuous<-dat$cardio_insulin_stage8_continuous*6
dat$cardio_insulin_pmolL_stage7_continuous<-dat$cardio_insulin_stage7_continuous*6
dat$cardio_insulin_pmolL_stage4_continuous<-dat$cardio_insulin_stage4_continuous*6

#Export offspring insulin data to HOMA2_IR calculator
HOMA<-dat[grep("glucose_|insulin_pmol|aln_qlet",names(dat))]
HOMA<-HOMA[-grep("mother",names(HOMA))]
HOMA<-HOMA[-grep("zscore",names(HOMA))]
write.table(HOMA,"",quote=FALSE,row.names=FALSE) #file path to add

#import HOMA2_IR variables
HOMA2_IR<-read_csv("",
                   col_types=list(col_character(),col_number(),col_number(),col_number(),col_number(),col_number(),col_number(),col_number(),col_number(),col_number(),col_number(),col_number(),col_number(),col_number(),col_number(),col_number())) #add file path
names(HOMA2_IR)<-c("aln_qlet","cardio_glucose_stage9_continuous","cardio_insulin_pmolL_stage9_continuous","cardio_HOMAB_stage9_continuous","cardio_HOMAS_stage9_continuous","cardio_IR_stage9_continuous",
                   "cardio_glucose_stage7_continuous","cardio_insulin_pmolL_stage7_continuous","cardio_HOMAB_stage7_continuous","cardio_HOMAS_stage7_continuous","cardio_IR_stage7_continuous",
                   "cardio_glucose_stage8_continuous","cardio_insulin_pmolL_stage8_continuous","cardio_HOMAB_stage8_continuous","cardio_HOMAS_stage8_continuous","cardio_IR_stage8_continuous")
#stage 9
HOMA2_IR$cardio_HOMAB_stage9_zscore_continuous<-scale(HOMA2_IR$cardio_HOMAB_stage9_continuous)
HOMA2_IR$cardio_HOMAS_stage9_zscore_continuous<-scale(HOMA2_IR$cardio_HOMAS_stage9_continuous)
HOMA2_IR$cardio_IR_stage9_zscore_continuous<-scale(HOMA2_IR$cardio_IR_stage9_continuous)

#stage 8
HOMA2_IR$cardio_HOMAB_stage8_zscore_continuous<-scale(HOMA2_IR$cardio_HOMAB_stage8_continuous)
HOMA2_IR$cardio_HOMAS_stage8_zscore_continuous<-scale(HOMA2_IR$cardio_HOMAS_stage8_continuous)
HOMA2_IR$cardio_IR_stage8_zscore_continuous<-scale(HOMA2_IR$cardio_IR_stage8_continuous)

#stage 7
HOMA2_IR$cardio_HOMAB_stage7_zscore_continuous<-scale(HOMA2_IR$cardio_HOMAB_stage7_continuous)
HOMA2_IR$cardio_HOMAS_stage7_zscore_continuous<-scale(HOMA2_IR$cardio_HOMAS_stage7_continuous)
HOMA2_IR$cardio_IR_stage7_zscore_continuous<-scale(HOMA2_IR$cardio_IR_stage7_continuous)

#add HOMA2_IR data to main data file
HOMA2_final<-HOMA2_IR[,c(1,4,5,6,9,10,11,14:25)]
HOMA2_final<-data.frame(HOMA2_final)
dat <- left_join(dat,HOMA2_final,by="aln_qlet")


#Remove variables not needed for downstream analysis and convert to dataframe
card<-dat[,grep("waist|fmi|bmi|sting|mice|pigeons|left|allstages|PC|sex|ish_|idh_|htn|cardio|autoimmune|covs_age_child_stage|prs|zscore|obese|sds|CRP", names(dat))]
card<-card[,-grep("psor|_as_|_uc_|_cd_|psa|_ms_|sjs|sarth|_sle_|_RA_|hashi|graves|t1d|_month_continuous|height|_weight|_lvm_|_Vd_|_Vs_|_LVIDs_|_LVIDd_|PWd_|_IVSd_|_e_|_esep_|_elat_|_Ewave_|_Awave_", names(card))]
card<-data.frame(card)

####Descriptive stats for cohort####
#Histogram of PRSs
plot<-ggplot(data=card, aes(Pt_1e.05MHC_zscore_continuous)) + 
  geom_histogram(binwidth=0.25,fill="gray",color="black") + 
  scale_x_continuous(limits = c(-4,5), expand = c(0, 0),breaks=c(-4-3,-2,-1,0,1,2,3,4,5))+ 
  scale_y_continuous(limits = c(0,350), expand = c(0, 0),breaks=c(0,50,100,150,200,250,300,350)) + 
  labs(x="JIA PRS Z score", y = "Density") +
  theme_classic()
plot

#Descriptive stats for PRS
x<-data.frame(card$Pt_1e.05MHC_zscore_continuous)
a<-sapply(x, sd, na.rm=TRUE)
b<-sapply(x, mean, na.rm=TRUE)
c<-sapply(x, min, na.rm=TRUE)
d<-sapply(x, max, na.rm=TRUE)
e<- "JIA PRS"
df<-data.frame(e,2815,a,b,c,d,row.names=NULL)
names(df)<-c("Cohort","N","Mean","SD","Min","Max")

#Descriptive stats for cohort
x<-data.frame(card$covs_age_child_stage9_f24_continuous)
a<-as.numeric(sapply(x, sd, na.rm=TRUE)/12)
b<-as.numeric(sapply(x, mean, na.rm=TRUE)/12)
c<-as.numeric(sapply(x, min, na.rm=TRUE)/12)
d<-as.numeric(sapply(x, max, na.rm=TRUE)/12)
e<-"Focus 24+"
f<-(sum(card$covs_sex_binary)/nrow(card))*100
df2<-data.frame(e,2815,b,a,c,d,f)
names(df2)<-c("Cohort","N","Mean","SD","Min","Max","Female (%)")

#####Descriptive stats for cardiovascular phenotypes#####
#Separate stage 9 phenotypes into binary and continuous categories 
cardstage9<-card[grep("stage9|sex|neg",names(card))]
cardstage9<-cardstage9[-grep("zscore|lefth",names(cardstage9))]
cardstage9b<-cardstage9[grep("binary|covs_sex",names(cardstage9))] #binary variables
cardstage9c<-cardstage9[grep("continuous|covs_age",names(cardstage9))] #continuous variables

#Descriptive stats continuous variables
a<-lapply(cardstage9c,mysummary)

b<-c(mean(cardstage9c[,1],na.rm=T),mean(cardstage9c[,2],na.rm=T),mean(cardstage9c[,3],na.rm=T),mean(cardstage9c[,4],na.rm=T),mean(cardstage9c[,5],na.rm=T),
     mean(cardstage9c[,6],na.rm=T),mean(cardstage9c[,7],na.rm=T),mean(cardstage9c[,8],na.rm=T),mean(cardstage9c[,9],na.rm=T),mean(cardstage9c[,10],na.rm=T),
     mean(cardstage9c[,11],na.rm=T),mean(cardstage9c[,12],na.rm=T),mean(cardstage9c[,13],na.rm=T),mean(cardstage9c[,14],na.rm=T),mean(cardstage9c[,15],na.rm=T),
     mean(cardstage9c[,16],na.rm=T),mean(cardstage9c[,17],na.rm=T),mean(cardstage9c[,18],na.rm=T),mean(cardstage9c[,19],na.rm=T),mean(cardstage9c[,20],na.rm=T),
     mean(cardstage9c[,21],na.rm=T),mean(cardstage9c[,22],na.rm=T),mean(cardstage9c[,23],na.rm=T),mean(cardstage9c[,24],na.rm=T),mean(cardstage9c[,25],na.rm=T),
     mean(cardstage9c[,26],na.rm=T),mean(cardstage9c[,27],na.rm=T),mean(cardstage9c[,28],na.rm=T),mean(cardstage9c[,29],na.rm=T),mean(cardstage9c[,30],na.rm=T),
     mean(cardstage9c[,31],na.rm=T),mean(cardstage9c[,32],na.rm=T),mean(cardstage9c[,33],na.rm=T),mean(cardstage9c[,34],na.rm=T))

c<-c(sd(cardstage9c[,1],na.rm=T),sd(cardstage9c[,2],na.rm=T),sd(cardstage9c[,3],na.rm=T),sd(cardstage9c[,4],na.rm=T),sd(cardstage9c[,5],na.rm=T),
     sd(cardstage9c[,6],na.rm=T),sd(cardstage9c[,7],na.rm=T),sd(cardstage9c[,8],na.rm=T),sd(cardstage9c[,9],na.rm=T),sd(cardstage9c[,10],na.rm=T),
     sd(cardstage9c[,11],na.rm=T),sd(cardstage9c[,12],na.rm=T),sd(cardstage9c[,13],na.rm=T),sd(cardstage9c[,14],na.rm=T),sd(cardstage9c[,15],na.rm=T),
     sd(cardstage9c[,16],na.rm=T),sd(cardstage9c[,17],na.rm=T),sd(cardstage9c[,18],na.rm=T),sd(cardstage9c[,19],na.rm=T),sd(cardstage9c[,20],na.rm=T),
     sd(cardstage9c[,21],na.rm=T),sd(cardstage9c[,22],na.rm=T),sd(cardstage9c[,23],na.rm=T),sd(cardstage9c[,24],na.rm=T),sd(cardstage9c[,25],na.rm=T),
     sd(cardstage9c[,26],na.rm=T),sd(cardstage9c[,27],na.rm=T),sd(cardstage9c[,28],na.rm=T),sd(cardstage9c[,29],na.rm=T),sd(cardstage9c[,30],na.rm=T),
     sd(cardstage9c[,31],na.rm=T),sd(cardstage9c[,32],na.rm=T),sd(cardstage9c[,33],na.rm=T),sd(cardstage9c[,34],na.rm=T))

summ<-data.frame(b,c)
row.names(summ)<-names(cardstage9c)

summ$b<-as.numeric(summ$b)
summ$c<-as.numeric(summ$c)
summ$b<-formatC(summ$b,digits=2,format="f")
summ$c<-formatC(summ$c,digits=2,format="f")
summ$d<-paste0(summ$b," (",summ$c,")")
summ$b<-NULL
summ$c<-NULL

summ$labels<-(row.names(summ))
summ$newname[summ$labels=="anthro_waist_stage9_continuous"]<-"Waist circumference, mean (SD), mm"
summ$newname[summ$labels=="covs_age_child_stage9_f24_continuous"]<-"Age,  mean (SD), months"
summ$newname[summ$labels=="anthro_bmi_stage9_continuous"]<-"Body mass index,  mean (SD), kg/m2"
summ$newname[summ$labels=="anthro_fmi_stage9_continuous"]<-"Fat mass index,  mean (SD), kg/m2"
summ$newname[summ$labels=="immuno_CRP_stage9_no_outliers_continuous"]<-"hsCRP,  mean (SD), mg/L"
summ$newname[summ$labels=="cardio_sbp_stage9_continuous"]<-"Systolic BP,  mean (SD), mmHg"
summ$newname[summ$labels=="cardio_dbp_stage9_continuous"]<-"Diastolic BP,  mean (SD), mmHg"
summ$newname[summ$labels=="cardio_trig_stage9_continuous"]<-"Triglycerides,  mean (SD), mmol/L"
summ$newname[summ$labels=="cardio_HDL_stage9_continuous"]<-"HDL,  mean (SD), mmol/L"
summ$newname[summ$labels=="cardio_LDL_stage9_continuous"]<-"LDL,  mean (SD), mmol/L"
summ$newname[summ$labels=="cardio_chol_stage9_continuous"]<-"Total cholesterol,  mean (SD), mmol/L"
summ$newname[summ$labels=="cardio_glucose_stage9_continuous"]<-"Glucose,  mean (SD), mmol/L"
summ$newname[summ$labels=="cardio_insulin_stage9_continuous"]<-"Insulin,  mean (SD), mU/L"
summ$newname[summ$labels=="cardio_apoA_stage9_continuous"]<-"Apo-AI,  mean (SD), g/L"
summ$newname[summ$labels=="cardio_apoB_stage9_continuous"]<-"Apo-B,  mean (SD), g/L"
summ$newname[summ$labels=="cardio_apoba_stage9_continuous"]<-"ApoB:AI,  mean (SD)"
summ$newname[summ$labels=="cardio_GlyAce_stage9_no_outliers_continuous"]<-"Glycoprotein acetylation,  mean (SD), mmol/L"
summ$newname[summ$labels=="cardio_rCIMT_stage9_continuous"]<-"Right cIMT,  mean (SD), mm"
summ$newname[summ$labels=="cardio_lCIMT_stage9_continuous"]<-"Left cIMT,  mean (SD), mm"
summ$newname[summ$labels=="cardio_PWV_stage9_continuous"]<-"Pulse wave velocity,  mean (SD), m/s"
summ$newname[summ$labels=="cardio_EAratio_stage9_continuous"]<-"Mitral E/A,  mean (SD)"
summ$newname[summ$labels=="cardio_Ee_stage9_continuous"]<-"E/e',  mean (SD)"
summ$newname[summ$labels=="cardio_FS_stage9_continuous"]<-"Fractional shortening,  mean (SD), %"
summ$newname[summ$labels=="cardio_ef_stage9_continuous"]<-"Ejection fraction,  mean (SD), %"
summ$newname[summ$labels=="cardio_lap_stage9_continuous"]<-"Left atrial pressure,  mean (SD), mmHg"
summ$newname[summ$labels=="cardio_lvmi_stage9_continuous"]<-"LVMI,  mean (SD), g/m2.7"
summ$newname[summ$labels=="cardio_IR_stage9_continuous"]<-"HOMA2_IR,  mean (SD)"
summ$newname[summ$labels=="cardio_CIMT_stage9_continuous"]<-"cIMT,  mean (SD), mm"
summ$newname[summ$labels=="cardio_LAD_stage9_continuous"]<-"Left atrial diameter,  mean (SD), cm"

#Descriptive stats binary variables
x<-names(cardstage9b)

summ2<-missing%>%
  filter(row.names(missing)%in% x)
summ2$per<-NA
summ2$per<-summ2$cases/summ2$total

summ2$labels<-(row.names(summ2))

summ2$newname[summ2$labels=="covs_sex_binary"]<-"Female, n (%)"
summ2$newname[summ2$labels=="anthro_overweightobese_stage9_binary"]<-"Overweight or obese, n (%)"
summ2$newname[summ2$labels=="anthro_obese_stage9_binary"]<-"Obese, n (%)"
summ2$newname[summ2$labels=="immuno_autoimmune_any_stage9_binary_ever"]<-"Autoimmune disease, n (%)"
summ2$newname[summ2$labels=="cardio_htn_stage9_binary"]<-"Hypertension (>140/90mmHg), n (%)"
summ2$newname[summ2$labels=="cardio_idh_stage9_binary"]<-"Isolated diastolic hyerptension (<140/>90mmHg), n (%)"
summ2$newname[summ2$labels=="cardio_prehtn_stage9_binary"]<-"High normal blood pressure (>120/80mmHg), n (%)"
summ2$newname[summ2$labels=="negcon_pigeons_binary"]<-"Pigeon infestation, n (%)"
summ2$newname[summ2$labels=="negcon_mice_binary"]<-"Mice infestation, n (%)"
summ2$newname[summ2$labels=="negcon_sting_ever_binary"]<-"Wasp sting, n (%)"

summ2$per<-formatC(summ2$per,digits=3,format="f")
summ2$d<-paste0(summ2$cases, " (",summ2$per,")")
summ2<-summ2[,c(7:9)]

#complete descriptive stats
summary<-rbind(summ,summ2)
summary<-na.omit(summary)
row.names(summary)<-summary$newname
summary9<-summary[1]

write.table(summary9,"", quote=FALSE) #dd file path

#####Histograms for continuous phenotpyes####
datconf24<-card[,-grep("binary|zscore", names(card))]
datconf24<-datconf24[,grep("stage9", names(datconf24))]
for(i in 1:ncol(datconf24)) {
  print(i)
  x<-names(datconf24)[i]
  print(x)
  d<-density(datconf24[, i], na.rm=TRUE)
  plot(d,main=x)
}

#####log transformation of CRP####
#stage 4
card$immuno_CRP_log_stage4_continuous<-log(card$immuno_CRP_stage4_continuous)
card$immuno_CRP_log_stage4_zscore_continuous<-scale(card$immuno_CRP_log_stage4_continuous)
card$immuno_CRP_log_stage4_no_outliers_continuous<-log(card$immuno_CRP_stage4_no_outliers_continuous)
card$immuno_CRP_log_stage4_no_outliers_zscore_continuous<-scale(card$immuno_CRP_log_stage4_no_outliers_continuous)

#stage 7
card$immuno_CRP_log_stage7_continuous<-log(card$immuno_CRP_stage7_continuous)
card$immuno_CRP_log_stage7_zscore_continuous<-scale(card$immuno_CRP_log_stage7_continuous)
card$immuno_CRP_log_stage7_no_outliers_continuous<-log(card$immuno_CRP_stage7_no_outliers_continuous)
card$immuno_CRP_log_stage7_no_outliers_zscore_continuous<-scale(card$immuno_CRP_log_stage7_no_outliers_continuous)

#stage 8
card$immuno_CRP_log_stage8_continuous<-log(card$immuno_CRP_stage8_continuous)
card$immuno_CRP_log_stage8_zscore_continuous<-scale(card$immuno_CRP_log_stage8_continuous)
card$immuno_CRP_log_stage8_no_outliers_continuous<-log(card$immuno_CRP_stage8_no_outliers_continuous)
card$immuno_CRP_log_stage8_no_outliers_zscore_continuous<-scale(card$immuno_CRP_log_stage8_no_outliers_continuous)

#stage 9
card$immuno_CRP_log_stage9_continuous<-log(card$immuno_CRP_stage9_continuous)
card$immuno_CRP_log_stage9_zscore_continuous<-scale(card$immuno_CRP_log_stage9_continuous)
card$immuno_CRP_log_stage9_no_outliers_continuous<-log(card$immuno_CRP_stage9_no_outliers_continuous)
card$immuno_CRP_log_stage9_no_outliers_zscore_continuous<-scale(card$immuno_CRP_log_stage9_no_outliers_continuous)

#boxplots for binary variables
datbin<-card[,grep("binary|Pt_1e.05MHC_zscore_continuous", names(card))]
for(i in 3:ncol(datbin)) {
  print(i)
  x<-names(datbin)[i]
  print(x)
  p<-ggplot(datbin, aes(as.factor(datbin[,i]), Pt_1e.05MHC_zscore_continuous)) + 
    geom_boxplot() +
    stat_summary(fun=mean, geom="point", shape=23, size=4) +
    #geom_jitter(shape=16, position=position_jitter(0.2)) +
    scale_x_discrete(limits=c("0", "1")) +
    xlab(x)
  plot(p)
}


####Regression analysis####
##Stage 9 positive and negative controls####
outcomes<-names(card[grepl("neg|autoimm",names(card))])
outcomes<-outcomes[!grepl("left", outcomes)]
covariates<-names(card[,c(316:325)]) #10 pcs
covariates1<-names(card[,c(316:325,1)]) #10 pcs, sex
covariates2<-names(card[,c(316:325,1,29)]) #10 pcs, sex, age
key <- data.frame(outcome=outcomes,
                  var.type=ifelse(!grepl("binary",outcomes),"continuous","binary")
)

#Run PheWAS
phewas_res_unadj <- run_analysis(exposure="Pt_1e.05MHC_zscore_continuous",
                                 outcomes=outcomes,
                                 covariates = NULL,
                                 df=card)
phewas_res_adj <- run_analysis(exposure="Pt_1e.05MHC_zscore_continuous",
                               outcomes=outcomes,
                               covariates = covariates,
                               df=card)
phewas_res_adj1 <- run_analysis(exposure="Pt_1e.05MHC_zscore_continuous",
                                outcomes=outcomes,
                                covariates = covariates1,
                                df=card)
phewas_res_adj2 <- run_analysis(exposure="Pt_1e.05MHC_zscore_continuous",
                                outcomes=outcomes,
                                covariates = covariates2,
                                df=card)

phewas_all<-NULL
phewas_all<-left_join(phewas_res_unadj,phewas_res_adj,by=c("outcome","exposure"))
phewas_all<-left_join(phewas_all,phewas_res_adj1,by=c("outcome","exposure"))
phewas_all<-left_join(phewas_all,phewas_res_adj2,by=c("outcome","exposure"))

phewas_all$label[grep("immuno_autoimmune_",phewas_all$outcome)]<- "Autoimmune disease"
phewas_all$label[grep("negcon_sting_ever",phewas_all$outcome)]<- "Bee/wasp sting"
phewas_all$label[grep("negcon_pigeons_binary",phewas_all$outcome)]<- "Pigeon infestation"
phewas_all$label[grep("negcon_mice_binary",phewas_all$outcome)]<- "Mouse infestation"
phewas_all$label[grep("negcon_lefthand_binary",phewas_all$outcome)]<- "Left handedness"

phewas_bin<-phewas_all[grep("binary",phewas_all$outcome),]

phewas_bin$or.x<-exp(phewas_bin$est.x)
phewas_bin$or.y<-exp(phewas_bin$est.y)
phewas_bin$cil.x<-exp((phewas_bin$est.x-(1.96*phewas_bin$se.x)))
phewas_bin$ciu.x<-exp((phewas_bin$est.x+(1.96*phewas_bin$se.x)))
phewas_bin$cil.y<-exp((phewas_bin$est.y-(1.96*phewas_bin$se.y)))
phewas_bin$ciu.y<-exp((phewas_bin$est.y+(1.96*phewas_bin$se.y)))
phewas_bin$or.x.x<-exp(phewas_bin$est.x.x)
phewas_bin$or.y.y<-exp(phewas_bin$est.y.y)
phewas_bin$cil.x.x<-exp((phewas_bin$est.x.x-(1.96*phewas_bin$se.x.x)))
phewas_bin$ciu.x.x<-exp((phewas_bin$est.x.x+(1.96*phewas_bin$se.x.x)))
phewas_bin$cil.y.y<-exp((phewas_bin$est.y.y-(1.96*phewas_bin$se.y.y)))
phewas_bin$ciu.y.y<-exp((phewas_bin$est.y.y+(1.96*phewas_bin$se.y.y)))

#####Forest plot results
tabletext <- cbind(c("Outcome", phewas_bin$label),c("Sample size",phewas_bin$n.x))
forestplot(tabletext, 
           mean = cbind(c(NA,phewas_bin$or.x), c(NA, phewas_bin$or.y),c(NA, phewas_bin$or.x.x),c(NA, phewas_bin$or.y.y)),
           lower = cbind(c(NA,phewas_bin$cil.x), c(NA, phewas_bin$cil.y), c(NA, phewas_bin$cil.x.x), c(NA, phewas_bin$cil.y.y)),
           upper = cbind(c(NA,phewas_bin$ciu.x), c(NA,phewas_bin$ciu.y), c(NA,phewas_bin$ciu.x.x), c(NA,phewas_bin$ciu.y.y)),
           col = fpColors(box = c("gray60","gray60","black","black")),
           fn.ci_norm =c(fpDrawCircleCI,fpDrawDiamondCI,fpDrawCircleCI,fpDrawDiamondCI),
           boxsize = 0.1,
           graphwidth=unit(4,"cm"),
           xlog=TRUE,
           xticks=c(0.4,1.0,2.5),
           clip=c(0.4,1.0,2.5),
           txt_gp = fpTxtGp(xlab=gpar(cex=1.0),ticks=gpar(cex=0.75)),
           xlab="OR per SD increase in JIA PRS",
           hrzl_lines = list("2" = gpar(lty = 1)),
           lineheight = unit(1.5,"cm"),
           legend = c("Unadjusted", "10PCs","10 PCs, age","10 PCs, age, sex"),
           legend_args = fpLegend(pos = list(x=1.5, y=0.5)))

write.table(phewas_bin,"",row.names=FALSE, quote=FALSE) #add file path


###Stage9 all cardiovascular phenotypes####
outcomes<-names(card[grepl("neg|cardio|anthro|autoimm|bmi|CRP",names(card))])
outcomes<-outcomes[grepl("stage9|neg",outcomes)]
outcomes<-outcomes[grepl("neg|zscore|binary",outcomes)]
outcomes<-outcomes[!grepl("RTA|covs|autoimmune|ish|HOMA|rCIMT|lCIMT|lap", outcomes)]
print(outcomes)
covariates<-names(card[,c(316:325)]) #10 pcs
covariates1<-names(card[,c(316:325,1)]) #10 pcs, sex
covariates2<-names(card[,c(316:325,1,29)]) #10 pcs, sex, age
key <- data.frame(outcome=outcomes,
                  var.type=ifelse(!grepl("binary",outcomes),"continuous","binary")
)

#Run PheWAS
phewas_res_unadj <- run_analysis(exposure="Pt_1e.05MHC_zscore_continuous",
                                 outcomes=outcomes,
                                 covariates = NULL,
                                 df=card)
phewas_res_adj <- run_analysis(exposure="Pt_1e.05MHC_zscore_continuous",
                               outcomes=outcomes,
                               covariates = covariates,
                               df=card)
phewas_res_adj1 <- run_analysis(exposure="Pt_1e.05MHC_zscore_continuous",
                                outcomes=outcomes,
                                covariates = covariates1,
                                df=card) 
phewas_res_adj2 <- run_analysis(exposure="Pt_1e.05MHC_zscore_continuous",
                                outcomes=outcomes,
                                covariates = covariates2,
                                df=card) 

phewas_all<-NULL
phewas_all<-left_join(phewas_res_unadj,phewas_res_adj,by=c("outcome","exposure"))
phewas_all<-left_join(phewas_all,phewas_res_adj1,by=c("outcome","exposure"))
phewas_all<-left_join(phewas_all,phewas_res_adj2,by=c("outcome","exposure"))
phewas_all$label[grep("anthro_overweightobese",phewas_all$outcome)]<- "Overweight or obese"
phewas_all$label[grep("anthro_obese_",phewas_all$outcome)]<- "Obese"
phewas_all$label[grep("immuno_autoimmune_",phewas_all$outcome)]<- "Ever had autoimmune disease"
phewas_all$label[grep("negcon_pigeons_binary",phewas_all$outcome)]<- "Pigeon infestation"
phewas_all$label[grep("negcon_mice_binary",phewas_all$outcome)]<- "Mice infestation"
phewas_all$label[grep("negcon_lefthand_binary",phewas_all$outcome)]<- "Left handedness by age 3yrs"
phewas_all$label[grep("anthro_waist_",phewas_all$outcome)]<- "Waist circumference"
phewas_all$label[grep("anthro_bmi_stage",phewas_all$outcome)]<- "BMI"
phewas_all$label[grep("anthro_fmi_stage",phewas_all$outcome)]<- "FMI"
phewas_all$label[grep("immuno_CRP_stage",phewas_all$outcome)]<- "CRP"
phewas_all$label[grep("cardio_sbp_stage",phewas_all$outcome)]<- "Systolic BP"
phewas_all$label[grep("cardio_dbp_stage",phewas_all$outcome)]<- "Diastolic BP"
phewas_all$label[grep("cardio_trig_stage",phewas_all$outcome)]<- "Triglycerides"
phewas_all$label[grep("cardio_HDL_stage",phewas_all$outcome)]<- "HDL"
phewas_all$label[grep("cardio_LDL_stage",phewas_all$outcome)]<- "LDL"
phewas_all$label[grep("cardio_chol_stage",phewas_all$outcome)]<- "Total cholesterol"
phewas_all$label[grep("cardio_glucose_stage",phewas_all$outcome)]<- "Blood glucose"
phewas_all$label[grep("cardio_insulin_stage",phewas_all$outcome)]<- "Blood insulin"
phewas_all$label[grep("cardio_apoA_stage",phewas_all$outcome)]<- "Apo-AI"
phewas_all$label[grep("cardio_apoB_stage",phewas_all$outcome)]<- "Apo-B"
phewas_all$label[grep("cardio_apoba_stage",phewas_all$outcome)]<- "Apo-B:AI"
phewas_all$label[grep("cardio_prehtn_stage",phewas_all$outcome)]<- "High normal BP (>120/80mmHg)"
phewas_all$label[grep("negcon_sting_ever",phewas_all$outcome)]<- "Ever been stung"
phewas_all$label[grep("cardio_htn_stage",phewas_all$outcome)]<- "Hypertensive (>140/90mmHg)"
phewas_all$label[grep("cardio_idh_stage",phewas_all$outcome)]<- "IDH (<140/>90mmHg)"
phewas_all$label[grep("cardio_ish_stage",phewas_all$outcome)]<- "Isolated systolic hypertension (>140/<80mmHg)"
phewas_all$label[grep("cardio_GlyAce_stage9_no_outliers_zscore_continuous",phewas_all$outcome)]<- "   Glycoprotein acetylation"
phewas_all$label[grep("cardio_rCIMT_stage",phewas_all$outcome)]<- "Right cIMT"
phewas_all$label[grep("cardio_lCIMT_stage",phewas_all$outcome)]<- "Left cIMT"
phewas_all$label[grep("cardio_PWV_stage",phewas_all$outcome)]<- "PWV"
phewas_all$label[grep("cardio_EAratio_stage",phewas_all$outcome)]<- "Mitral E/A ratio"
phewas_all$label[grep("cardio_Ee_stage",phewas_all$outcome)]<- "E/e'"
phewas_all$label[grep("cardio_ef_stage",phewas_all$outcome)]<- "Ejection fraction"
phewas_all$label[grep("cardio_FS_stage",phewas_all$outcome)]<- "Fractional shortening"
phewas_all$label[grep("cardio_lap_stage",phewas_all$outcome)]<- "Left atrial pressure"
phewas_all$label[grep("cardio_lvmi_stage",phewas_all$outcome)]<- "LVMI"
phewas_all$label[grep("immuno_CRP_stage9_zscore_continuous",phewas_all$outcome)]<- "CRP"
phewas_all$label[grep("immuno_CRP_stage9_no_outliers_zscore_continuous",phewas_all$outcome)]<- "CRP (no outliers)"
phewas_all$label[grep("immuno_CRP_log_stage9_zscore_continuous",phewas_all$outcome)]<- "CRP (log)"
phewas_all$label[grep("immuno_CRP_log_stage9_no_outliers_zscore_continuous",phewas_all$outcome)]<- "log hsCRP"
phewas_all$label[grep("anthro_waist_log_stage9_",phewas_all$outcome)]<- "Waist circumference (log)"
phewas_all$label[grep("anthro_bmi_log_stage9_",phewas_all$outcome)]<- "BMI (log)"
phewas_all$label[grep("anthro_fmi_log_stage9_",phewas_all$outcome)]<- "FMI (log)"
phewas_all$label[grep("cardio_trig_log_stage9_",phewas_all$outcome)]<- "Triglycerides (log)"
phewas_all$label[grep("cardio_insulin_log_stage9_",phewas_all$outcome)]<- "Insulin (log)"
phewas_all$label[grep("cardio_IR_stage9_",phewas_all$outcome)]<- "HOMA2_IR"
phewas_all$label[grep("cardio_CIMT_stage",phewas_all$outcome)]<- "cIMT"
phewas_all$label[grep("cardio_LAD_stage9",phewas_all$outcome)]<- "Left atrial diameter"

phewas_bin<-phewas_all[grep("binary",phewas_all$outcome),]

phewas_con<-phewas_all[-grep("binary",phewas_all$outcome),]

phewas_bin$or.x<-exp(phewas_bin$est.x)
phewas_bin$or.y<-exp(phewas_bin$est.y)
phewas_bin$cil.x<-exp((phewas_bin$est.x-(1.96*phewas_bin$se.x)))
phewas_bin$ciu.x<-exp((phewas_bin$est.x+(1.96*phewas_bin$se.x)))
phewas_bin$cil.y<-exp((phewas_bin$est.y-(1.96*phewas_bin$se.y)))
phewas_bin$ciu.y<-exp((phewas_bin$est.y+(1.96*phewas_bin$se.y)))
phewas_bin$or.x.x<-exp(phewas_bin$est.x.x)
phewas_bin$or.y.y<-exp(phewas_bin$est.y.y)
phewas_bin$cil.x.x<-exp((phewas_bin$est.x.x-(1.96*phewas_bin$se.x.x)))
phewas_bin$ciu.x.x<-exp((phewas_bin$est.x.x+(1.96*phewas_bin$se.x.x)))
phewas_bin$cil.y.y<-exp((phewas_bin$est.y.y-(1.96*phewas_bin$se.y.y)))
phewas_bin$ciu.y.y<-exp((phewas_bin$est.y.y+(1.96*phewas_bin$se.y.y)))

x<-c("Ever had autoimmune disease", "Ever been stung","Pigeon infestation","Mice infestation","Left handedness by age 3yrs",
     "Overweight or obese","Obese","High normal BP (>120/80mmHg)","Hypertensive (>140/90mmHg)","IDH (<140/>90mmHg)")
rownames(phewas_bin)<-phewas_bin$label
phewas_bin<-phewas_bin[x,]

x<-c("Systolic BP", "Diastolic BP","Triglycerides","HDL","LDL","Total cholesterol","Apo-AI","Apo-B",
     "Apo-B:AI","Blood glucose","Blood insulin","HOMA2_IR","   Glycoprotein acetylation","log hsCRP",
     "Waist circumference","FMI","BMI","cIMT","PWV","Mitral E/A ratio","E/e'","Left atrial diameter","Ejection fraction","Fractional shortening",
     "LVMI")
rownames(phewas_con)<-phewas_con$label
phewas_con<-phewas_con[x,]

#Forest plot results
#unadj and PC/sex adjusted binary
tabletext <- cbind(c("Category","Anthropometry      ",NA,"Blood pressure",NA,NA),c("Outcome (24yrs)", phewas_bin$label[6:10]),c("Sample size",phewas_bin$n.x[6:10]))
forestplot(tabletext, 
           mean = cbind(c(NA,phewas_bin$or.x[6:10]),c(NA, phewas_bin$or.x.x[6:10])),
           lower = cbind(c(NA,phewas_bin$cil.x[6:10]),c(NA, phewas_bin$cil.x.x[6:10])),
           upper = cbind(c(NA,phewas_bin$ciu.x[6:10]),c(NA,phewas_bin$ciu.x.x[6:10])),
           col = fpColors(box = c("grey", "black")),
           boxsize = 0.12,
           graphwidth=unit(6,"cm"),
           xlog=TRUE,
           xticks=c(0.25,1.0,4.0),
           clip=c(0.25,1.0,4.0),
           txt_gp = fpTxtGp(xlab=gpar(cex=1.0),ticks=gpar(cex=0.75)),
           xlab="OR per SD increase in JIA PRS",
           hrzl_lines = list("2" = gpar(lty = 1),
                             "4" = gpar(lty = 1)),
           lineheight = unit(0.8,"cm"))

#unadj and PCS/sex adj continuous
tabletext <- cbind(c("Category", "Blood pressure", NA, "Blood biomarkers",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,"Anthropometry",NA,NA,"Early atherosclerosis/","arteriosclerosis","Cardiac structure/function",NA,NA,NA,NA,NA),
                   c("Outcome (24yrs)", phewas_con$label),c("Sample size",phewas_con$n.x))
forestplot(tabletext, 
           mean = cbind(c(NA,phewas_con$est.x),c(NA,phewas_con$est.x.x)),
           lower = cbind(c(NA,phewas_con$ci.l.x),c(NA, phewas_con$ci.l.x.x)),
           upper = cbind(c(NA,phewas_con$ci.u.x),c(NA,phewas_con$ci.u.x.x)),
           col = fpColors(box = c("grey", "black")),
           boxsize = 0.12,
           graphwidth=unit(6,"cm"),
           xlog=FALSE,
           txt_gp = fpTxtGp(xlab=gpar(cex=1.0),ticks=gpar(cex=0.75)),
           xlab="SD increase per SD increase in JIA PRS",
           hrzl_lines = list("2" = gpar(lty = 1),
                             "4" = gpar(lty = 1),
                             "16" = gpar(lty = 1),
                             "19" = gpar(lty = 1),
                             "21" = gpar(lty = 1)),
           xticks=c(-0.10,-0.05,0,0.05,0.10),
           clip=c(-0.10,0.10),
           lineheight = unit(0.6,"cm"),
           graph.pos = 4)

write.table(phewas_con, "",row.names=FALSE, quote=FALSE) #add file path
write.table(phewas_bin, "",row.names=FALSE, quote=FALSE) #add file path

#####Sensitivity analysis 1 - remove MHC####
#Stage 9 no MHC
outcomesnoMHC<-names(card[grepl("neg|cardio|anthro|autoimm|bmi|CRP",names(card))])
outcomesnoMHC<-outcomesnoMHC[grepl("stage9|neg",outcomesnoMHC)]
outcomesnoMHC<-outcomesnoMHC[grepl("neg|zscore|binary",outcomesnoMHC)]
outcomesnoMHC<-outcomesnoMHC[!grepl("RTA|covs|left|ish|HOMA|lap", outcomesnoMHC)]
print(outcomesnoMHC)

exposures <- c("Pt_1e.05MHC_zscore_continuous","Pt_1e.05_zscore_continuous")
key <- data.frame(outcome=outcomesnoMHC,
                  var.type=ifelse(!grepl("binary",outcomesnoMHC),"continuous","binary")
)


phewas_res_list <- lapply(exposures,
                          run_analysis,
                          outcomes=outcomesnoMHC,covariates=covariates1,df=card)

phewas_res_all <- dplyr::bind_rows(phewas_res_list)

phewas_all<-phewas_res_all
phewas_all$label<-NA
phewas_all$label[grep("anthro_overweightobese",phewas_all$outcome)]<- "Overweight or obese"
phewas_all$label[grep("anthro_obese_",phewas_all$outcome)]<- "Obese"
phewas_all$label[grep("immuno_autoimmune_",phewas_all$outcome)]<- "Ever had autoimmune disease"
phewas_all$label[grep("negcon_pigeons_binary",phewas_all$outcome)]<- "Pigeon infestation"
phewas_all$label[grep("negcon_mice_binary",phewas_all$outcome)]<- "Mice infestation"
phewas_all$label[grep("negcon_lefthand_binary",phewas_all$outcome)]<- "Left handedness by age 3yrs"
phewas_all$label[grep("anthro_waist_",phewas_all$outcome)]<- "Waist circumference"
phewas_all$label[grep("anthro_bmi_stage",phewas_all$outcome)]<- "BMI"
phewas_all$label[grep("anthro_fmi_stage",phewas_all$outcome)]<- "Fat mass index"
phewas_all$label[grep("immuno_CRP_stage",phewas_all$outcome)]<- "CRP"
phewas_all$label[grep("cardio_sbp_stage",phewas_all$outcome)]<- "Systolic BP"
phewas_all$label[grep("cardio_dbp_stage",phewas_all$outcome)]<- "Diastolic BP"
phewas_all$label[grep("cardio_trig_stage",phewas_all$outcome)]<- "Triglycerides"
phewas_all$label[grep("cardio_HDL_stage",phewas_all$outcome)]<- "HDL cholesterol"
phewas_all$label[grep("cardio_LDL_stage",phewas_all$outcome)]<- "LDL cholesterol"
phewas_all$label[grep("cardio_chol_stage",phewas_all$outcome)]<- "Total cholesterol"
phewas_all$label[grep("cardio_glucose_stage",phewas_all$outcome)]<- "Blood glucose"
phewas_all$label[grep("cardio_insulin_stage",phewas_all$outcome)]<- "Blood insulin"
phewas_all$label[grep("cardio_apoA_stage",phewas_all$outcome)]<- "Apo-lipoprotein A"
phewas_all$label[grep("cardio_apoB_stage",phewas_all$outcome)]<- "Apo-lipoprotein B"
phewas_all$label[grep("cardio_apoba_stage",phewas_all$outcome)]<- "Apo-lipopotein B:A"
phewas_all$label[grep("cardio_prehtn_stage",phewas_all$outcome)]<- "High normal BP (>120/80mmHg)"
phewas_all$label[grep("negcon_sting_ever",phewas_all$outcome)]<- "Ever been stung"
phewas_all$label[grep("cardio_htn_stage",phewas_all$outcome)]<- "Hypertensive (>140/90mmHg)"
phewas_all$label[grep("cardio_idh_stage",phewas_all$outcome)]<- "Isolated diastolic hypertension (<140/>90mmHg)"
phewas_all$label[grep("cardio_ish_stage",phewas_all$outcome)]<- "Isolated systolic hypertension (>140/<80mmHg)"
phewas_all$label[grep("cardio_GlyAce_stage9_no_outliers_zscore_continuous",phewas_all$outcome)]<- "   Glycoprotein acetylation"
phewas_all$label[grep("cardio_rCIMT_stage",phewas_all$outcome)]<- "Right CIMT"
phewas_all$label[grep("cardio_lCIMT_stage",phewas_all$outcome)]<- "Left CIMT"
phewas_all$label[grep("cardio_PWV_stage",phewas_all$outcome)]<- "Pulse wave velocity"
phewas_all$label[grep("cardio_EAratio_stage",phewas_all$outcome)]<- "Mitral E/A ratio"
phewas_all$label[grep("cardio_Ee_stage",phewas_all$outcome)]<- "E/e'"
phewas_all$label[grep("cardio_ef_stage",phewas_all$outcome)]<- "Ejection fraction"
phewas_all$label[grep("cardio_FS_stage",phewas_all$outcome)]<- "Fractional shortening"
phewas_all$label[grep("cardio_lap_stage",phewas_all$outcome)]<- "Left atrial pressure"
phewas_all$label[grep("cardio_lvmi_stage",phewas_all$outcome)]<- "Left ventricular mass index"
phewas_all$label[grep("immuno_CRP_stage9_zscore_continuous",phewas_all$outcome)]<- "CRP"
phewas_all$label[grep("immuno_CRP_stage9_no_outliers_zscore_continuous",phewas_all$outcome)]<- "CRP (no outliers)"
phewas_all$label[grep("immuno_CRP_log_stage9_zscore_continuous",phewas_all$outcome)]<- "CRP (log)"
phewas_all$label[grep("immuno_CRP_log_stage9_no_outliers_zscore_continuous",phewas_all$outcome)]<- "logCRP"
phewas_all$label[grep("anthro_waist_log_stage9_",phewas_all$outcome)]<- "Waist circumference (log)"
phewas_all$label[grep("anthro_bmi_log_stage9_",phewas_all$outcome)]<- "BMI (log)"
phewas_all$label[grep("anthro_fmi_log_stage9_",phewas_all$outcome)]<- "FMI (log)"
phewas_all$label[grep("cardio_trig_log_stage9_",phewas_all$outcome)]<- "Triglycerides (log)"
phewas_all$label[grep("cardio_insulin_log_stage9_",phewas_all$outcome)]<- "Insulin (log)"
phewas_all$label[grep("cardio_IR_stage9_",phewas_all$outcome)]<- "Insulin resistance"
phewas_all$label[grep("cardio_CIMT_stage",phewas_all$outcome)]<- "cIMT"
phewas_all$label[grep("cardio_LAD_stage9",phewas_all$outcome)]<- "Left atrial diameter"

phewas_all$sort<-NA
phewas_all$sort[phewas_all$exposure=="Pt_1e.05MHC_zscore_continuous"]<-paste0(phewas_all$label,"1")
phewas_all$sort[phewas_all$exposure=="Pt_1e.05_zscore_continuous"]<-paste0(phewas_all$label,"2")

phewas_bin<-phewas_all[grep("binary",phewas_all$outcome),]

phewas_con<-phewas_all[-grep("binary",phewas_all$outcome),]

phewas_bin$or<-exp(phewas_bin$est)
phewas_bin$or<-exp(phewas_bin$est)
phewas_bin$cil<-exp((phewas_bin$est-(1.96*phewas_bin$se)))
phewas_bin$ciu<-exp((phewas_bin$est+(1.96*phewas_bin$se)))
phewas_bin$cil<-exp((phewas_bin$est-(1.96*phewas_bin$se)))
phewas_bin$ciu<-exp((phewas_bin$est+(1.96*phewas_bin$se)))
phewas_bin$or<-exp(phewas_bin$est)
phewas_bin$or<-exp(phewas_bin$est)
phewas_bin$cil<-exp((phewas_bin$est-(1.96*phewas_bin$se)))
phewas_bin$ciu<-exp((phewas_bin$est+(1.96*phewas_bin$se)))
phewas_bin$cil<-exp((phewas_bin$est-(1.96*phewas_bin$se)))
phewas_bin$ciu<-exp((phewas_bin$est+(1.96*phewas_bin$se)))

x<-c("Ever had autoimmune disease1", "Ever had autoimmune disease2", 
     "Ever been stung1","Ever been stung2",
     "Pigeon infestation1","Pigeon infestation2",
     "Mice infestation1","Mice infestation2",
     "Overweight or obese1","Overweight or obese2",
     "Obese1","Obese2",
     "High normal BP (>120/80mmHg)1","High normal BP (>120/80mmHg)2",
     "Hypertensive (>140/90mmHg)1","Hypertensive (>140/90mmHg)2",
     "Isolated diastolic hypertension (<140/>90mmHg)1","Isolated diastolic hypertension (<140/>90mmHg)2")
rownames(phewas_bin)<-phewas_bin$sort
phewas_bin<-phewas_bin[x,]

#Binary forest plot
##unadj and PC/sex adjusted

tabletext <- cbind(c("Category","Controls",NA,NA,NA,"Anthropometry",NA,"Blood pressure",NA,NA),
                   c("Outcome (24yrs)", phewas_bin$label[1],phewas_bin$label[3],phewas_bin$label[5],phewas_bin$label[7],
                     phewas_bin$label[9],phewas_bin$label[11],phewas_bin$label[13],phewas_bin$label[15],phewas_bin$label[17]),
                   c("Sample size", phewas_bin$n[1],phewas_bin$n[3],phewas_bin$n[5],phewas_bin$n[7],
                     phewas_bin$n[9],phewas_bin$n[11],phewas_bin$n[13],phewas_bin$n[15],phewas_bin$n[17]))

forestplot(tabletext, 
           mean = cbind(c(NA,phewas_bin$or[c(1,3,5,7,9,11,13,15,17)]),c(NA,phewas_bin$or[c(2,4,6,8,10,12,14,16,18)])),
           lower = cbind(c(NA,phewas_bin$cil[c(1,3,5,7,9,11,13,15,17)]),c(NA,phewas_bin$cil[c(2,4,6,8,10,12,14,16,18)])),
           upper = cbind(c(NA,phewas_bin$ciu[c(1,3,5,7,9,11,13,15,17)]),c(NA,phewas_bin$ciu[c(2,4,6,8,10,12,14,16,18)])),
           col = fpColors(box = c("grey", "black")),
           boxsize = 0.1,
           graphwidth=unit(6,"cm"),
           xlog=TRUE,
           xticks=c(0.25,1.0,4.0),
           clip=c(0.25,1.0,4.0),
           txt_gp = fpTxtGp(xlab=gpar(cex=1.0),ticks=gpar(cex=0.75)),
           xlab="OR per SD increase in JIA PRS",
           hrzl_lines = list("2" = gpar(lty = 1),
                             "6" = gpar(lty = 1),
                             "8" = gpar(lty = 1)),
           lineheight = unit(0.75,"cm"))

#Continuous forest plot 
x<-c("Systolic BP1","Systolic BP2", 
     "Diastolic BP1","Diastolic BP2",
     "Triglycerides1","Triglycerides2",
     "HDL cholesterol1","HDL cholesterol2",
     "LDL cholesterol1","LDL cholesterol2",
     "Total cholesterol1","Total cholesterol2",
     "Apo-lipoprotein A1","Apo-lipoprotein A2",
     "Apo-lipoprotein B1","Apo-lipoprotein B2",
     "Apo-lipopotein B:A1","Apo-lipopotein B:A2",
     "Blood glucose1","Blood glucose2",
     "Blood insulin1","Blood insulin2",
     "Insulin resistance1","Insulin resistance2",
     "   Glycoprotein acetylation1","   Glycoprotein acetylation2",
     "logCRP1","logCRP2",
     "Waist circumference1","Waist circumference2",
     "Fat mass index1","Fat mass index2",
     "BMI1","BMI2",
     "cIMT1","cIMT2",
     "Pulse wave velocity1","Pulse wave velocity2",
     "Mitral E/A ratio1","Mitral E/A ratio2",
     "E/e'1","E/e'2",
     "Left atrial diameter1","Left atrial diameter2",
     "Ejection fraction1","Ejection fraction2",
     "Fractional shortening1","Fractional shortening2",
     "Left ventricular mass index1","Left ventricular mass index2")
rownames(phewas_con)<-phewas_con$sort
phewas_con<-phewas_con[x,]

#unadj and PCS/sex adj only
tabletext <- cbind(c("Category", "Blood pressure", NA,
                     "Blood biomarkers",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                     "Anthropometry",NA,NA,"Early atherosclerosis",NA,"Cardiac structure and function",NA,NA,NA,NA,NA),
                   c("Outcome (24yrs)",phewas_con$label[1],phewas_con$label[3],phewas_con$label[5],phewas_con$label[7],phewas_con$label[9],
                     phewas_con$label[11],phewas_con$label[13],phewas_con$label[15],phewas_con$label[17],phewas_con$label[19],phewas_con$label[21],
                     phewas_con$label[23],phewas_con$label[25],phewas_con$label[27],phewas_con$label[29],phewas_con$label[31],phewas_con$label[33],
                     phewas_con$label[35],phewas_con$label[37],phewas_con$label[39],phewas_con$label[41],phewas_con$label[43],phewas_con$label[45],
                     phewas_con$label[47],phewas_con$label[49]),
                   c("Sample size",phewas_con$n[1],phewas_con$n[3],phewas_con$n[5],phewas_con$n[7],phewas_con$n[9],
                     phewas_con$n[11],phewas_con$n[13],phewas_con$n[15],phewas_con$n[17],phewas_con$n[19],phewas_con$n[21],
                     phewas_con$n[23],phewas_con$n[25],phewas_con$n[27],phewas_con$n[29],phewas_con$n[31],phewas_con$n[33],
                     phewas_con$n[35],phewas_con$n[37],phewas_con$n[39],phewas_con$n[41],phewas_con$n[43],phewas_con$n[45],
                     phewas_con$n[47],phewas_con$n[49]))


forestplot(tabletext, 
           mean = cbind(c(NA,phewas_con$est[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49)]),
                        c(NA,phewas_con$est[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50)])),
           lower = cbind(c(NA,phewas_con$ci.l[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49)]),
                         c(NA,phewas_con$ci.l[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50)])),
           upper = cbind(c(NA,phewas_con$ci.u[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49)]),
                         c(NA,phewas_con$ci.u[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50)])),
           col = fpColors(box = c("grey", "black")),
           boxsize = 0.1,
           graphwidth=unit(6,"cm"),
           xlog=FALSE,
           txt_gp = fpTxtGp(xlab=gpar(cex=1.0),ticks=gpar(cex=0.75)),
           xlab="SD increase per SD increase in JIA PRS",
           hrzl_lines = list("2" = gpar(lty = 1),
                             "4" = gpar(lty = 1),
                             "16" = gpar(lty = 1),
                             "19" = gpar(lty = 1),
                             "21" = gpar(lty = 1)),
           xticks=c(-0.10,-0.05,0,0.05,0.10),
           clip=c(-0.10,0.10),
           lineheight = unit(0.75,"cm"),
           graph.pos = 4)

write.table(phewas_con,"",row.names=FALSE, quote=FALSE)
write.table(phewas_bin,"",row.names=FALSE, quote=FALSE)

#####Sensitivity analysis 2 - remove AI participants#### 
noAI<-card%>%
  filter(immuno_autoimmune_any_stage9_binary_ever!=1|is.na(immuno_autoimmune_any_stage9_binary_ever))
outcomesnoAI<-names(noAI[grepl("neg|cardio|anthro|bmi|CRP",names(noAI))])
outcomesnoAI<-outcomesnoAI[grepl("stage9|neg",outcomesnoAI)]
outcomesnoAI<-outcomesnoAI[grepl("neg|zscore|binary",outcomesnoAI)]
outcomesnoAI<-outcomesnoAI[!grepl("RTA|covs|ish|HOMA|rCIMT|lCIMT|lap", outcomesnoAI)]
print(outcomesnoAI)
covariates1<-names(card[,c(316:325,1)]) #10 pcs, sex

exposures <- c("Pt_1e.05MHC_zscore_continuous")
key <- data.frame(outcome=outcomesnoAI,
                  var.type=ifelse(!grepl("binary",outcomesnoAI),"continuous","binary")
)

phewas_res_list <- lapply(exposures,
                          run_analysis,
                          outcomes=outcomesnoAI,covariates=covariates1,df=noAI)

phewas_res_all <- dplyr::bind_rows(phewas_res_list)


outcomes<-names(card[grepl("neg|cardio|anthro|autoimm|bmi|CRP",names(card))])
outcomes<-outcomes[grepl("stage9|neg",outcomes)]
outcomes<-outcomes[grepl("neg|zscore|binary",outcomes)]
outcomes<-outcomes[!grepl("RTA|covs|ish|HOMA|rCIMT|lCIMT|lap", outcomes)]
print(outcomes)
covariates<-names(card[,c(316:325)]) #10 pcs
covariates1<-names(card[,c(316:325,1)]) #10 pcs, sex
covariates2<-names(card[,c(316:325,1,29)]) #10 pcs, sex, age
key <- data.frame(outcome=outcomes,
                  var.type=ifelse(!grepl("binary",outcomes),"continuous","binary")
)

phewas_res_adj1 <- run_analysis(exposure="Pt_1e.05MHC_zscore_continuous",
                                outcomes=outcomes,
                                covariates = covariates1,
                                df=card)

phewas_all<-NULL
phewas_all<-left_join(phewas_res_adj1,phewas_res_all,by=c("outcome","exposure"))

phewas_all$label<-NA
phewas_all$label[grep("anthro_overweightobese",phewas_all$outcome)]<- "Overweight or obese"
phewas_all$label[grep("anthro_obese_",phewas_all$outcome)]<- "Obese"
phewas_all$label[grep("immuno_autoimmune_",phewas_all$outcome)]<- "Ever had autoimmune disease"
phewas_all$label[grep("negcon_pigeons_binary",phewas_all$outcome)]<- "Pigeon infestation"
phewas_all$label[grep("negcon_mice_binary",phewas_all$outcome)]<- "Mice infestation"
phewas_all$label[grep("negcon_lefthand_binary",phewas_all$outcome)]<- "Left handedness by age 3yrs"
phewas_all$label[grep("anthro_waist_",phewas_all$outcome)]<- "Waist circumference"
phewas_all$label[grep("anthro_bmi_stage",phewas_all$outcome)]<- "BMI"
phewas_all$label[grep("anthro_fmi_stage",phewas_all$outcome)]<- "Fat mass index"
phewas_all$label[grep("immuno_CRP_stage",phewas_all$outcome)]<- "CRP"
phewas_all$label[grep("cardio_sbp_stage",phewas_all$outcome)]<- "Systolic BP"
phewas_all$label[grep("cardio_dbp_stage",phewas_all$outcome)]<- "Diastolic BP"
phewas_all$label[grep("cardio_trig_stage",phewas_all$outcome)]<- "Triglycerides"
phewas_all$label[grep("cardio_HDL_stage",phewas_all$outcome)]<- "HDL cholesterol"
phewas_all$label[grep("cardio_LDL_stage",phewas_all$outcome)]<- "LDL cholesterol"
phewas_all$label[grep("cardio_chol_stage",phewas_all$outcome)]<- "Total cholesterol"
phewas_all$label[grep("cardio_glucose_stage",phewas_all$outcome)]<- "Blood glucose"
phewas_all$label[grep("cardio_insulin_stage",phewas_all$outcome)]<- "Blood insulin"
phewas_all$label[grep("cardio_apoA_stage",phewas_all$outcome)]<- "Apo-lipoprotein A"
phewas_all$label[grep("cardio_apoB_stage",phewas_all$outcome)]<- "Apo-lipoprotein B"
phewas_all$label[grep("cardio_apoba_stage",phewas_all$outcome)]<- "Apo-lipopotein B:A"
phewas_all$label[grep("cardio_prehtn_stage",phewas_all$outcome)]<- "High normal BP (>120/80mmHg)"
phewas_all$label[grep("negcon_sting_ever",phewas_all$outcome)]<- "Ever been stung"
phewas_all$label[grep("cardio_htn_stage",phewas_all$outcome)]<- "Hypertensive (>140/90mmHg)"
phewas_all$label[grep("cardio_idh_stage",phewas_all$outcome)]<- "Isolated diastolic hypertension (<140/>90mmHg)"
phewas_all$label[grep("cardio_ish_stage",phewas_all$outcome)]<- "Isolated systolic hypertension (>140/<80mmHg)"
phewas_all$label[grep("cardio_GlyAce_stage9_no_outliers_zscore_continuous",phewas_all$outcome)]<- "Glycoprotein acetylation"
phewas_all$label[grep("cardio_GlyAce_stage9_zscore_continuous",phewas_all$outcome)]<-"Glycoprotein acetylation outliers"
phewas_all$label[grep("cardio_rCIMT_stage",phewas_all$outcome)]<- "Right CIMT"
phewas_all$label[grep("cardio_lCIMT_stage",phewas_all$outcome)]<- "Left CIMT"
phewas_all$label[grep("cardio_PWV_stage",phewas_all$outcome)]<- "Pulse wave velocity"
phewas_all$label[grep("cardio_EAratio_stage",phewas_all$outcome)]<- "Mitral E/A ratio"
phewas_all$label[grep("cardio_Ee_stage",phewas_all$outcome)]<- "E/e'"
phewas_all$label[grep("cardio_ef_stage",phewas_all$outcome)]<- "Ejection fraction"
phewas_all$label[grep("cardio_FS_stage",phewas_all$outcome)]<- "Fractional shortening"
phewas_all$label[grep("cardio_lap_stage",phewas_all$outcome)]<- "Left atrial pressure"
phewas_all$label[grep("cardio_lvmi_stage",phewas_all$outcome)]<- "Left ventricular mass index"
phewas_all$label[grep("immuno_CRP_stage9_zscore_continuous",phewas_all$outcome)]<- "CRP"
phewas_all$label[grep("immuno_CRP_stage9_no_outliers_zscore_continuous",phewas_all$outcome)]<- "CRP (no outliers)"
phewas_all$label[grep("immuno_CRP_log_stage9_zscore_continuous",phewas_all$outcome)]<- "CRP (log)"
phewas_all$label[grep("immuno_CRP_log_stage9_no_outliers_zscore_continuous",phewas_all$outcome)]<- "logCRP"
phewas_all$label[grep("anthro_waist_log_stage9_",phewas_all$outcome)]<- "Waist circumference (log)"
phewas_all$label[grep("anthro_bmi_log_stage9_",phewas_all$outcome)]<- "BMI (log)"
phewas_all$label[grep("anthro_fmi_log_stage9_",phewas_all$outcome)]<- "FMI (log)"
phewas_all$label[grep("cardio_trig_log_stage9_",phewas_all$outcome)]<- "Triglycerides (log)"
phewas_all$label[grep("cardio_insulin_log_stage9_",phewas_all$outcome)]<- "Insulin (log)"
phewas_all$label[grep("cardio_IR_stage9_",phewas_all$outcome)]<- "Insulin resistance"
phewas_all$label[grep("cardio_CIMT_stage",phewas_all$outcome)]<- "cIMT"
phewas_all$label[grep("cardio_LAD_stage9",phewas_all$outcome)]<- "Left atrial diameter"

phewas_bin<-phewas_all[grep("binary",phewas_all$outcome),]

phewas_con<-phewas_all[-grep("binary",phewas_all$outcome),]

phewas_bin$or.x<-exp(phewas_bin$est.x)
phewas_bin$cil.x<-exp((phewas_bin$est.x-(1.96*phewas_bin$se.x)))
phewas_bin$ciu.x<-exp((phewas_bin$est.x+(1.96*phewas_bin$se.x)))
phewas_bin$or.y<-exp(phewas_bin$est.y)
phewas_bin$cil.y<-exp((phewas_bin$est.y-(1.96*phewas_bin$se.y)))
phewas_bin$ciu.y<-exp((phewas_bin$est.y+(1.96*phewas_bin$se.y)))


x<-c("Ever had autoimmune disease", "Ever been stung","Pigeon infestation","Mice infestation","Left handedness by age 3yrs",
     "Overweight or obese","Obese","High normal BP (>120/80mmHg)","Hypertensive (>140/90mmHg)","Isolated diastolic hypertension (<140/>90mmHg)")
rownames(phewas_bin)<-phewas_bin$label
phewas_bin<-phewas_bin[x,]

#Binary forest plot #Autoimmune variable removed from plot as meaningless with AI participants removed
##unadj and PC/sex adjusted
tabletext <- cbind(c("Category","Controls",NA,NA,NA,NA,NA,"Anthropometry",NA,NA,NA,"Blood pressure",NA,NA,NA,NA,NA),
                   c("Outcome (24yrs)", phewas_bin$label[2],NA,phewas_bin$label[3],NA,phewas_bin$label[4],NA,
                     phewas_bin$label[6],NA,phewas_bin$label[7],NA,phewas_bin$label[8],NA,phewas_bin$label[9],NA,phewas_bin$label[10],NA),
                   c("Sample size", phewas_bin$n.x[2],phewas_bin$n.y[2],phewas_bin$n.x[3],phewas_bin$n.y[3],
                     phewas_bin$n.x[4],phewas_bin$n.y[4],phewas_bin$n.x[6],phewas_bin$n.y[6],phewas_bin$n.x[7],
                     phewas_bin$n.y[7],phewas_bin$n.x[8],phewas_bin$n.y[8],phewas_bin$n.x[9],phewas_bin$n.y[9],phewas_bin$n.x[10],phewas_bin$n.y[10]))
styles <- fpShapesGp(
  lines = list(
    gpar(col = "black"),
    gpar(col = "grey"),
    gpar(col = "black"),
    gpar(col = "grey"),
    gpar(col = "black"),
    gpar(col = "grey"),
    gpar(col = "black"),
    gpar(col = "grey"),
    gpar(col = "black"),
    gpar(col = "grey"),
    gpar(col = "black"),
    gpar(col = "grey"),
    gpar(col = "black"),
    gpar(col = "grey"),
    gpar(col = "black"),
    gpar(col = "grey"),
    gpar(col = "black")
  ),
  box = list(
    gpar(fill = "black"),
    gpar(fill = "grey"),
    gpar(fill = "black"),
    gpar(fill = "grey"),
    gpar(fill = "black"),
    gpar(fill = "grey"),
    gpar(fill = "black"),
    gpar(fill = "grey"),
    gpar(fill = "black"),
    gpar(fill = "grey"),
    gpar(fill = "black"),
    gpar(fill = "grey"),
    gpar(fill = "black"),
    gpar(fill = "grey"),
    gpar(fill = "black"),
    gpar(fill = "grey"),
    gpar(fill = "black")
  ) 
)
forestplot(tabletext, 
           mean = cbind(c(NA,phewas_bin$or.x[2],phewas_bin$or.y[2],phewas_bin$or.x[3],phewas_bin$or.y[3],
                          phewas_bin$or.x[4],phewas_bin$or.y[4],phewas_bin$or.x[6],phewas_bin$or.y[6],phewas_bin$or.x[7],phewas_bin$or.y[7],
                          phewas_bin$or.x[8],phewas_bin$or.y[8],phewas_bin$or.x[9],phewas_bin$or.y[9],phewas_bin$or.x[10],phewas_bin$or.y[10])),
           lower = cbind(c(NA,phewas_bin$cil.x[2],phewas_bin$cil.y[2],phewas_bin$cil.x[3],phewas_bin$cil.y[3],
                           phewas_bin$cil.x[4],phewas_bin$cil.y[4],phewas_bin$cil.x[6],phewas_bin$cil.y[6],phewas_bin$cil.x[7],phewas_bin$cil.y[7],
                           phewas_bin$cil.x[8],phewas_bin$cil.y[8],phewas_bin$cil.x[9],phewas_bin$cil.y[9],phewas_bin$cil.x[10],phewas_bin$cil.y[10])),
           upper = cbind(c(NA,phewas_bin$ciu.x[2],phewas_bin$ciu.y[2],phewas_bin$ciu.x[3],phewas_bin$ciu.y[3],
                           phewas_bin$ciu.x[4],phewas_bin$ciu.y[4],phewas_bin$ciu.x[6],phewas_bin$ciu.y[6],phewas_bin$ciu.x[7],phewas_bin$ciu.y[7],
                           phewas_bin$ciu.x[8],phewas_bin$ciu.y[8],phewas_bin$ciu.x[9],phewas_bin$ciu.y[9],phewas_bin$ciu.x[10],phewas_bin$ciu.y[10])),
           #col = fpColors(box = c("grey", "black")),
           boxsize = 0.1,
           graphwidth=unit(6,"cm"),
           xlog=TRUE,
           xticks=c(0.25,1.0,4.0),
           clip=c(0.25,1.0,4.0),
           txt_gp = fpTxtGp(xlab=gpar(cex=1.0),ticks=gpar(cex=0.75)),
           xlab="OR per SD increase in JIA PRS",
           hrzl_lines = list("2" = gpar(lty = 1),
                             "8" = gpar(lty = 1),
                             "12" = gpar(lty = 1)),
           lineheight = unit(0.75,"cm"),
           shapes_gp = styles)

#Continuous forest plot 
x<-c("Systolic BP", "Diastolic BP","Triglycerides","HDL cholesterol","LDL cholesterol","Total cholesterol","Apo-lipoprotein A","Apo-lipoprotein B",
     "Apo-lipopotein B:A","Blood glucose","Blood insulin","Insulin resistance","Glycoprotein acetylation","logCRP",
     "Waist circumference","Fat mass index","BMI","cIMT","Pulse wave velocity","Mitral E/A ratio","E/e'","Left atrial diameter","Ejection fraction","Fractional shortening",
     "Left ventricular mass index")
rownames(phewas_con)<-phewas_con$label
phewas_con<-phewas_con[x,]

#unadj and PCS/sex adj only
tabletext <- cbind(c("Category", "Blood pressure", NA,NA,NA, 
                     "Blood biomarkers",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                     "Anthropometry",NA,NA,NA,NA,NA,"Early atherosclerosis",NA,NA,NA,"Cardiac structure and function",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
                   c("Outcome (24yrs)",phewas_con$label[1],NA,phewas_con$label[2],NA,phewas_con$label[3],NA,phewas_con$label[4],NA,phewas_con$label[5],NA,
                     phewas_con$label[6],NA,phewas_con$label[7],NA,phewas_con$label[8],NA,phewas_con$label[9],NA,phewas_con$label[10],NA,
                     phewas_con$label[11],NA,phewas_con$label[12],NA,phewas_con$label[13],NA,phewas_con$label[14],NA,phewas_con$label[15],NA,
                     phewas_con$label[16],NA,phewas_con$label[17],NA,phewas_con$label[18],NA,phewas_con$label[19],NA,phewas_con$label[20],NA,
                     phewas_con$label[21],NA,phewas_con$label[22],NA,phewas_con$label[23],NA,phewas_con$label[24],NA,phewas_con$label[25],NA),
                   c("Sample size",phewas_con$n.x[1],phewas_con$n.y[1],phewas_con$n.x[2],phewas_con$n.y[2],phewas_con$n.x[3],phewas_con$n.y[3],
                     phewas_con$n.x[4],phewas_con$n.y[4],phewas_con$n.x[5],phewas_con$n.y[5],phewas_con$n.x[6],phewas_con$n.y[6],phewas_con$n.x[7],
                     phewas_con$n.y[7],phewas_con$n.x[8],phewas_con$n.y[8],phewas_con$n.x[9],phewas_con$n.y[9],phewas_con$n.x[10],phewas_con$n.y[10],
                     phewas_con$n.x[11],phewas_con$n.y[11],phewas_con$n.x[12],phewas_con$n.y[12],phewas_con$n.x[13],phewas_con$n.y[13],phewas_con$n.x[14],phewas_con$n.y[14],
                     phewas_con$n.x[15],phewas_con$n.y[15],phewas_con$n.x[16],phewas_con$n.y[16],phewas_con$n.x[17],phewas_con$n.y[17],phewas_con$n.x[18],phewas_con$n.y[18],
                     phewas_con$n.x[19],phewas_con$n.y[19],phewas_con$n.x[20],phewas_con$n.y[20],phewas_con$n.x[21],phewas_con$n.y[21],phewas_con$n.x[22],phewas_con$n.y[22],
                     phewas_con$n.x[23],phewas_con$n.y[23],phewas_con$n.x[24],phewas_con$n.y[24],phewas_con$n.x[25],phewas_con$n.y[25]))

styles <- fpShapesGp(
  lines = list(
    gpar(col = "black"),
    gpar(col = "grey"),
    gpar(col = "black"),
    gpar(col = "grey"),
    gpar(col = "black"),
    gpar(col = "grey"),
    gpar(col = "black"),
    gpar(col = "grey"),
    gpar(col = "black"),
    gpar(col = "grey"),
    gpar(col = "black"),
    gpar(col = "grey"),
    gpar(col = "black"),
    gpar(col = "grey"),
    gpar(col = "black"),
    gpar(col = "grey"),
    gpar(col = "black"),
    gpar(col = "grey"),
    gpar(col = "black"),
    gpar(col = "grey"),
    gpar(col = "black"),
    gpar(col = "grey"),
    gpar(col = "black"),
    gpar(col = "grey"),
    gpar(col = "black"),
    gpar(col = "grey"),
    gpar(col = "black"),
    gpar(col = "grey"),
    gpar(col = "black"),
    gpar(col = "grey"),
    gpar(col = "black"),
    gpar(col = "grey"),
    gpar(col = "black"),
    gpar(col = "grey"),
    gpar(col = "black"),
    gpar(col = "grey"),
    gpar(col = "black"),
    gpar(col = "grey"),
    gpar(col = "black"),
    gpar(col = "grey"),
    gpar(col = "black"),
    gpar(col = "grey"),
    gpar(col = "black"),
    gpar(col = "grey"),
    gpar(col = "black"),
    gpar(col = "grey"),
    gpar(col = "black"),
    gpar(col = "grey"),
    gpar(col = "black"),
    gpar(col = "grey"),
    gpar(col = "black")
  ),
  box = list(
    gpar(fill = "black"),
    gpar(fill = "grey"),
    gpar(fill = "black"),
    gpar(fill = "grey"),
    gpar(fill = "black"),
    gpar(fill = "grey"),
    gpar(fill = "black"),
    gpar(fill = "grey"),
    gpar(fill = "black"),
    gpar(fill = "grey"),
    gpar(fill = "black"),
    gpar(fill = "grey"),
    gpar(fill = "black"),
    gpar(fill = "grey"),
    gpar(fill = "black"),
    gpar(fill = "grey"),
    gpar(fill = "black"),
    gpar(fill = "grey"),
    gpar(fill = "black"),
    gpar(fill = "grey"),
    gpar(fill = "black"),
    gpar(fill = "grey"),
    gpar(fill = "black"),
    gpar(fill = "grey"),
    gpar(fill = "black"),
    gpar(fill = "grey"),
    gpar(fill = "black"),
    gpar(fill = "grey"),
    gpar(fill = "black"),
    gpar(fill = "grey"),
    gpar(fill = "black"),
    gpar(fill = "grey"),
    gpar(fill = "black"),
    gpar(fill = "grey"),
    gpar(fill = "black"),
    gpar(fill = "grey"),
    gpar(fill = "black"),
    gpar(fill = "grey"),
    gpar(fill = "black"),
    gpar(fill = "grey"),
    gpar(fill = "black"),
    gpar(fill = "grey"),
    gpar(fill = "black"),
    gpar(fill = "grey"),
    gpar(fill = "black"),
    gpar(fill = "grey"),
    gpar(fill = "black"),
    gpar(fill = "grey"),
    gpar(fill = "black"),
    gpar(fill = "grey"),
    gpar(fill = "black")
  ) 
)

forestplot(tabletext, 
           mean = cbind(c(NA,phewas_con$est.x[1],phewas_con$est.y[1],phewas_con$est.x[2],phewas_con$est.y[2],phewas_con$est.x[3],phewas_con$est.y[3],
                          phewas_con$est.x[4],phewas_con$est.y[4],phewas_con$est.x[5],phewas_con$est.y[5],phewas_con$est.x[6],phewas_con$est.y[6],
                          phewas_con$est.x[7],phewas_con$est.y[7],phewas_con$est.x[8],phewas_con$est.y[8],phewas_con$est.x[9],phewas_con$est.y[9],
                          phewas_con$est.x[10],phewas_con$est.y[10],phewas_con$est.x[11],phewas_con$est.y[11],phewas_con$est.x[12],phewas_con$est.y[12],
                          phewas_con$est.x[13],phewas_con$est.y[13],phewas_con$est.x[14],phewas_con$est.y[14],phewas_con$est.x[15],phewas_con$est.y[15],
                          phewas_con$est.x[16],phewas_con$est.y[16],phewas_con$est.x[17],phewas_con$est.y[17],phewas_con$est.x[18],phewas_con$est.y[18],
                          phewas_con$est.x[19],phewas_con$est.y[19],phewas_con$est.x[20],phewas_con$est.y[20],phewas_con$est.x[21],phewas_con$est.y[21],
                          phewas_con$est.x[22],phewas_con$est.y[22],phewas_con$est.x[23],phewas_con$est.y[23],phewas_con$est.x[24],phewas_con$est.y[24],
                          phewas_con$est.x[25],phewas_con$est.y[25])),
           lower = cbind(c(NA,phewas_con$ci.l.x[1],phewas_con$ci.l.y[1],phewas_con$ci.l.x[2],phewas_con$ci.l.y[2],phewas_con$ci.l.x[3],phewas_con$ci.l.y[3],
                           phewas_con$ci.l.x[4],phewas_con$ci.l.y[4],phewas_con$ci.l.x[5],phewas_con$ci.l.y[5],phewas_con$ci.l.x[6],phewas_con$ci.l.y[6],
                           phewas_con$ci.l.x[7],phewas_con$ci.l.y[7],phewas_con$ci.l.x[8],phewas_con$ci.l.y[8],phewas_con$ci.l.x[9],phewas_con$ci.l.y[9],
                           phewas_con$ci.l.x[10],phewas_con$ci.l.y[10],phewas_con$ci.l.x[11],phewas_con$ci.l.y[11],phewas_con$ci.l.x[12],phewas_con$ci.l.y[12],
                           phewas_con$ci.l.x[13],phewas_con$ci.l.y[13],phewas_con$ci.l.x[14],phewas_con$ci.l.y[14],phewas_con$ci.l.x[15],phewas_con$ci.l.y[15],
                           phewas_con$ci.l.x[16],phewas_con$ci.l.y[16],phewas_con$ci.l.x[17],phewas_con$ci.l.y[17],phewas_con$ci.l.x[18],phewas_con$ci.l.y[18],
                           phewas_con$ci.l.x[19],phewas_con$ci.l.y[19],phewas_con$ci.l.x[20],phewas_con$ci.l.y[20],phewas_con$ci.l.x[21],phewas_con$ci.l.y[21],
                           phewas_con$ci.l.x[22],phewas_con$ci.l.y[22],phewas_con$ci.l.x[23],phewas_con$ci.l.y[23],phewas_con$ci.l.x[24],phewas_con$ci.l.y[24],
                           phewas_con$ci.l.x[25],phewas_con$ci.l.y[25])),
           upper = cbind(c(NA,phewas_con$ci.u.x[1],phewas_con$ci.u.y[1],phewas_con$ci.u.x[2],phewas_con$ci.u.y[2],phewas_con$ci.u.x[3],phewas_con$ci.u.y[3],
                           phewas_con$ci.u.x[4],phewas_con$ci.u.y[4],phewas_con$ci.u.x[5],phewas_con$ci.u.y[5],phewas_con$ci.u.x[6],phewas_con$ci.u.y[6],
                           phewas_con$ci.u.x[7],phewas_con$ci.u.y[7],phewas_con$ci.u.x[8],phewas_con$ci.u.y[8],phewas_con$ci.u.x[9],phewas_con$ci.u.y[9],
                           phewas_con$ci.u.x[10],phewas_con$ci.u.y[10],phewas_con$ci.u.x[11],phewas_con$ci.u.y[11],phewas_con$ci.u.x[12],phewas_con$ci.u.y[12],
                           phewas_con$ci.u.x[13],phewas_con$ci.u.y[13],phewas_con$ci.u.x[14],phewas_con$ci.u.y[14],phewas_con$ci.u.x[15],phewas_con$ci.u.y[15],
                           phewas_con$ci.u.x[16],phewas_con$ci.u.y[16],phewas_con$ci.u.x[17],phewas_con$ci.u.y[17],phewas_con$ci.u.x[18],phewas_con$ci.u.y[18],
                           phewas_con$ci.u.x[19],phewas_con$ci.u.y[19],phewas_con$ci.u.x[20],phewas_con$ci.u.y[20],phewas_con$ci.u.x[21],phewas_con$ci.u.y[21],
                           phewas_con$ci.u.x[22],phewas_con$ci.u.y[22],phewas_con$ci.u.x[23],phewas_con$ci.u.y[23],phewas_con$ci.u.x[24],phewas_con$ci.u.y[24],
                           phewas_con$ci.u.x[25],phewas_con$ci.u.y[25])),
           #col = fpColors(box = c("grey", "black")),
           boxsize = 0.1,
           graphwidth=unit(6,"cm"),
           xlog=FALSE,
           txt_gp = fpTxtGp(xlab=gpar(cex=1.0),ticks=gpar(cex=0.75)),
           xlab="SD increase per SD increase in JIA PRS",
           hrzl_lines = list("2" = gpar(lty = 1),
                             "6" = gpar(lty = 1),
                             "30" = gpar(lty = 1),
                             "36" = gpar(lty = 1),
                             "40" = gpar(lty = 1)),
           xticks=c(-0.10,-0.05,0,0.05,0.10),
           clip=c(-0.10,0.10),
           lineheight = unit(0.5,"cm"),
           graph.pos = 4,
           shapes_gp = styles)

write.table(phewas_con,"",row.names=FALSE, quote=FALSE) #add file path
write.table(phewas_bin,"",row.names=FALSE, quote=FALSE) #add file path

#####Sensitivity analysis 3 - different PRS P value thresholds####
outcomessens3<-names(card[grepl("neg|cardio|anthro|autoimm|bmi|CRP",names(card))])
outcomessens3<-outcomessens3[grepl("stage9|neg",outcomessens3)]
outcomessens3<-outcomessens3[grepl("neg|zscore|binary",outcomessens3)]
outcomessens3<-outcomessens3[!grepl("RTA|covs|left|ish|HOMA|rCIMT|lCIMT|lap", outcomessens3)]
print(outcomessens3)
covariates1<-names(card[,c(316:325,1)]) #10 pcs, sex, age

exposures <- c("Pt_0.0001MHC_zscore_continuous","Pt_0.001MHC_zscore_continuous","Pt_0.01MHC_zscore_continuous",
               "Pt_1e.06MHC_zscore_continuous","Pt_1e.07MHC_zscore_continuous","Pt_5e.08MHC_zscore_continuous",
               "Pt_1e.05MHC_zscore_continuous")
key <- data.frame(outcome=outcomessens3,
                  var.type=ifelse(!grepl("binary",outcomessens3),"continuous","binary")
)

phewas_res_list <- lapply(exposures,
                          run_analysis,
                          outcomes=outcomessens3,covariates=covariates1,df=card)

phewas_res_allsens3 <- dplyr::bind_rows(phewas_res_list)

phewas_allsens3<-phewas_res_allsens3
phewas_allsens3$label<-NA
phewas_allsens3$label[grep("anthro_overweightobese",phewas_allsens3$outcome)]<- "Overweight or obese"
phewas_allsens3$label[grep("anthro_obese_",phewas_allsens3$outcome)]<- "Obese"
phewas_allsens3$label[grep("immuno_autoimmune_",phewas_allsens3$outcome)]<- "Ever had autoimmune disease"
phewas_allsens3$label[grep("negcon_pigeons_binary",phewas_allsens3$outcome)]<- "Pigeon infestation"
phewas_allsens3$label[grep("negcon_mice_binary",phewas_allsens3$outcome)]<- "Mice infestation"
phewas_allsens3$label[grep("negcon_lefthand_binary",phewas_allsens3$outcome)]<- "Left handedness by age 3yrs"
phewas_allsens3$label[grep("anthro_waist_",phewas_allsens3$outcome)]<- "Waist circumference"
phewas_allsens3$label[grep("anthro_bmi_stage",phewas_allsens3$outcome)]<- "BMI"
phewas_allsens3$label[grep("anthro_fmi_stage",phewas_allsens3$outcome)]<- "FMI"
phewas_allsens3$label[grep("immuno_CRP_stage",phewas_allsens3$outcome)]<- "CRP"
phewas_allsens3$label[grep("cardio_sbp_stage",phewas_allsens3$outcome)]<- "Systolic BP"
phewas_allsens3$label[grep("cardio_dbp_stage",phewas_allsens3$outcome)]<- "Diastolic BP"
phewas_allsens3$label[grep("cardio_trig_stage",phewas_allsens3$outcome)]<- "Triglycerides"
phewas_allsens3$label[grep("cardio_HDL_stage",phewas_allsens3$outcome)]<- "HDL"
phewas_allsens3$label[grep("cardio_LDL_stage",phewas_allsens3$outcome)]<- "LDL"
phewas_allsens3$label[grep("cardio_chol_stage",phewas_allsens3$outcome)]<- "Total cholesterol"
phewas_allsens3$label[grep("cardio_glucose_stage",phewas_allsens3$outcome)]<- "Blood glucose"
phewas_allsens3$label[grep("cardio_insulin_stage",phewas_allsens3$outcome)]<- "Blood insulin"
phewas_allsens3$label[grep("cardio_apoA_stage",phewas_allsens3$outcome)]<- "Apo-AI"
phewas_allsens3$label[grep("cardio_apoB_stage",phewas_allsens3$outcome)]<- "Apo-B"
phewas_allsens3$label[grep("cardio_apoba_stage",phewas_allsens3$outcome)]<- "Apo-B:AI"
phewas_allsens3$label[grep("cardio_prehtn_stage",phewas_allsens3$outcome)]<- "High normal BP (>120/80mmHg)"
phewas_allsens3$label[grep("negcon_sting_ever",phewas_allsens3$outcome)]<- "Ever been stung"
phewas_allsens3$label[grep("cardio_htn_stage",phewas_allsens3$outcome)]<- "Hypertensive (>140/90mmHg)"
phewas_allsens3$label[grep("cardio_idh_stage",phewas_allsens3$outcome)]<- "IDH (<140/>90mmHg)"
phewas_allsens3$label[grep("cardio_ish_stage",phewas_allsens3$outcome)]<- "Isolated systolic hypertension (>140/<80mmHg)"
phewas_allsens3$label[grep("cardio_GlyAce_stage9_no_outliers_zscore_continuous",phewas_allsens3$outcome)]<- "   Glycoprotein acetylation"
phewas_allsens3$label[grep("cardio_rCIMT_stage",phewas_allsens3$outcome)]<- "Right CIMT"
phewas_allsens3$label[grep("cardio_lCIMT_stage",phewas_allsens3$outcome)]<- "Left CIMT"
phewas_allsens3$label[grep("cardio_PWV_stage",phewas_allsens3$outcome)]<- "Pulse wave velocity"
phewas_allsens3$label[grep("cardio_EAratio_stage",phewas_allsens3$outcome)]<- "Mitral E/A ratio"
phewas_allsens3$label[grep("cardio_Ee_stage",phewas_allsens3$outcome)]<- "E/e'"
phewas_allsens3$label[grep("cardio_ef_stage",phewas_allsens3$outcome)]<- "Ejection fraction"
phewas_allsens3$label[grep("cardio_FS_stage",phewas_allsens3$outcome)]<- "Fractional shortening"
phewas_allsens3$label[grep("cardio_lap_stage",phewas_allsens3$outcome)]<- "Left atrial pressure"
phewas_allsens3$label[grep("cardio_lvmi_stage",phewas_allsens3$outcome)]<- "LVMI"
phewas_allsens3$label[grep("immuno_CRP_stage9_zscore_continuous",phewas_allsens3$outcome)]<- "CRP"
phewas_allsens3$label[grep("immuno_CRP_stage9_no_outliers_zscore_continuous",phewas_allsens3$outcome)]<- "CRP (no outliers)"
#phewas_allsens3$label[grep("cardio_GlyAce_stage9_no_outliers_zscore_continuous",phewas_allsens3$outcome)]<- "Glycoprotein acetylation (no outliers)"
phewas_allsens3$label[grep("immuno_CRP_log_stage9_zscore_continuous",phewas_allsens3$outcome)]<- "CRP (log)"
phewas_allsens3$label[grep("immuno_CRP_log_stage9_no_outliers_zscore_continuous",phewas_allsens3$outcome)]<- "log hsCRP"
phewas_allsens3$label[grep("anthro_waist_log_stage9_",phewas_allsens3$outcome)]<- "Waist circumference (log)"
phewas_allsens3$label[grep("anthro_bmi_log_stage9_",phewas_allsens3$outcome)]<- "BMI (log)"
phewas_allsens3$label[grep("anthro_fmi_log_stage9_",phewas_allsens3$outcome)]<- "FMI (log)"
phewas_allsens3$label[grep("cardio_trig_log_stage9_",phewas_allsens3$outcome)]<- "Triglycerides (log)"
phewas_allsens3$label[grep("cardio_insulin_log_stage9_",phewas_allsens3$outcome)]<- "Insulin (log)"
phewas_allsens3$label[grep("cardio_IR_stage9_",phewas_allsens3$outcome)]<- "HOMA2_IR"
phewas_allsens3$label[grep("cardio_CIMT_stage",phewas_allsens3$outcome)]<- "cIMT"
phewas_allsens3$label[grep("cardio_LAD_stage9",phewas_allsens3$outcome)]<- "Left atrial diameter"

phewas_allsens3$sort<-NA
phewas_allsens3$sort[phewas_allsens3$exposure=="Pt_0.0001MHC_zscore_continuous"]<-paste0(phewas_allsens3$label[phewas_allsens3$exposure=="Pt_0.0001MHC_zscore_continuous"],"3")
phewas_allsens3$sort[phewas_allsens3$exposure=="Pt_0.001MHC_zscore_continuous"]<-paste0(phewas_allsens3$label[phewas_allsens3$exposure=="Pt_0.001MHC_zscore_continuous"],"2")
phewas_allsens3$sort[phewas_allsens3$exposure=="Pt_0.01MHC_zscore_continuous"]<-paste0(phewas_allsens3$label[phewas_allsens3$exposure=="Pt_0.01MHC_zscore_continuous"],"1")
phewas_allsens3$sort[phewas_allsens3$exposure=="Pt_1e.05MHC_zscore_continuous"]<-paste0(phewas_allsens3$label[phewas_allsens3$exposure=="Pt_1e.05MHC_zscore_continuous"],"4")
phewas_allsens3$sort[phewas_allsens3$exposure=="Pt_1e.06MHC_zscore_continuous"]<-paste0(phewas_allsens3$label[phewas_allsens3$exposure=="Pt_1e.06MHC_zscore_continuous"],"5")
phewas_allsens3$sort[phewas_allsens3$exposure=="Pt_1e.07MHC_zscore_continuous"]<-paste0(phewas_allsens3$label[phewas_allsens3$exposure=="Pt_1e.07MHC_zscore_continuous"],"6")
phewas_allsens3$sort[phewas_allsens3$exposure=="Pt_5e.08MHC_zscore_continuous"]<-paste0(phewas_allsens3$label[phewas_allsens3$exposure=="Pt_5e.08MHC_zscore_continuous"],"7")

phewas_bin<-phewas_allsens3[grep("binary",phewas_allsens3$outcome),]

phewas_con<-phewas_allsens3[-grep("binary",phewas_allsens3$outcome),]

phewas_bin$or<-exp(phewas_bin$est)
phewas_bin$or<-exp(phewas_bin$est)
phewas_bin$cil<-exp((phewas_bin$est-(1.96*phewas_bin$se)))
phewas_bin$ciu<-exp((phewas_bin$est+(1.96*phewas_bin$se)))
phewas_bin$cil<-exp((phewas_bin$est-(1.96*phewas_bin$se)))
phewas_bin$ciu<-exp((phewas_bin$est+(1.96*phewas_bin$se)))
phewas_bin$or<-exp(phewas_bin$est)
phewas_bin$or<-exp(phewas_bin$est)
phewas_bin$cil<-exp((phewas_bin$est-(1.96*phewas_bin$se)))
phewas_bin$ciu<-exp((phewas_bin$est+(1.96*phewas_bin$se)))
phewas_bin$cil<-exp((phewas_bin$est-(1.96*phewas_bin$se)))
phewas_bin$ciu<-exp((phewas_bin$est+(1.96*phewas_bin$se)))

x<-c("Ever had autoimmune disease1", "Ever had autoimmune disease2", "Ever had autoimmune disease3", "Ever had autoimmune disease4", "Ever had autoimmune disease5", "Ever had autoimmune disease6", "Ever had autoimmune disease7", 
     "Ever been stung1","Ever been stung2","Ever been stung3","Ever been stung4","Ever been stung5","Ever been stung6","Ever been stung7",
     "Pigeon infestation1","Pigeon infestation2","Pigeon infestation3","Pigeon infestation4","Pigeon infestation5","Pigeon infestation6","Pigeon infestation7",
     "Mice infestation1","Mice infestation2","Mice infestation3","Mice infestation4","Mice infestation5","Mice infestation6","Mice infestation7",
     "Overweight or obese1","Overweight or obese2","Overweight or obese3","Overweight or obese4","Overweight or obese5","Overweight or obese6","Overweight or obese7",
     "Obese1","Obese2","Obese3","Obese4","Obese5","Obese6","Obese7",
     "High normal BP (>120/80mmHg)1","High normal BP (>120/80mmHg)2","High normal BP (>120/80mmHg)3","High normal BP (>120/80mmHg)4","High normal BP (>120/80mmHg)5","High normal BP (>120/80mmHg)6","High normal BP (>120/80mmHg)7",
     "Hypertensive (>140/90mmHg)1","Hypertensive (>140/90mmHg)2","Hypertensive (>140/90mmHg)3","Hypertensive (>140/90mmHg)4","Hypertensive (>140/90mmHg)5","Hypertensive (>140/90mmHg)6","Hypertensive (>140/90mmHg)7",
     "IDH (<140/>90mmHg)1","IDH (<140/>90mmHg)2","IDH (<140/>90mmHg)3","IDH (<140/>90mmHg)4","IDH (<140/>90mmHg)5","IDH (<140/>90mmHg)6","IDH (<140/>90mmHg)7")
rownames(phewas_bin)<-phewas_bin$sort
phewas_bin<-phewas_bin[x,]

#Binary forest plot
##unadj and PC/sex adjusted

tabletext <- cbind(c("Category","Controls",NA,NA,NA,"Anthropometry   ",NA,"Blood pressure",NA,NA),
                   c("Outcome (24yrs)", phewas_bin$label[1],phewas_bin$label[8],phewas_bin$label[15],phewas_bin$label[22],
                     phewas_bin$label[29],phewas_bin$label[36],phewas_bin$label[43],phewas_bin$label[50],phewas_bin$label[57]),
                   c("Sample size", phewas_bin$n[1],phewas_bin$n[8],phewas_bin$n[15],phewas_bin$n[22],
                     phewas_bin$n[29],phewas_bin$n[36],phewas_bin$n[43],phewas_bin$n[50],phewas_bin$n[57]))

forestplot(tabletext, 
           mean = cbind(c(NA,phewas_bin$or[c(1,8,15,22,29,36,43,50,57)]),c(NA,phewas_bin$or[c(2,9,16,23,30,37,44,51,58)]),c(NA,phewas_bin$or[c(3,10,17,4,31,38,45,52,59)]),c(NA,phewas_bin$or[c(4,11,18,25,32,39,46,53,60)]),c(NA,phewas_bin$or[c(6,13,20,27,34,41,48,55,62)]),c(NA,phewas_bin$or[c(6,13,20,27,34,41,48,55,62)]),c(NA,phewas_bin$or[c(7,14,21,28,35,42,49,56,63)])),
           lower = cbind(c(NA,phewas_bin$cil[c(1,8,15,22,29,36,43,50,57)]),c(NA,phewas_bin$cil[c(2,9,16,23,30,37,44,51,58)]),c(NA,phewas_bin$cil[c(3,10,17,4,31,38,45,52,59)]),c(NA,phewas_bin$cil[c(4,11,18,25,32,39,46,53,60)]),c(NA,phewas_bin$cil[c(6,13,20,27,34,41,48,55,62)]),c(NA,phewas_bin$cil[c(6,13,20,27,34,41,48,55,62)]),c(NA,phewas_bin$cil[c(7,14,21,28,35,42,49,56,63)])),
           upper = cbind(c(NA,phewas_bin$ciu[c(1,8,15,22,29,36,43,50,57)]),c(NA,phewas_bin$ciu[c(2,9,16,23,30,37,44,51,58)]),c(NA,phewas_bin$ciu[c(3,10,17,4,31,38,45,52,59)]),c(NA,phewas_bin$ciu[c(4,11,18,25,32,39,46,53,60)]),c(NA,phewas_bin$ciu[c(6,13,20,27,34,41,48,55,62)]),c(NA,phewas_bin$ciu[c(6,13,20,27,34,41,48,55,62)]),c(NA,phewas_bin$ciu[c(7,14,21,28,35,42,49,56,63)])),
           col = fpColors(box = c("grey90","grey75","grey60","grey45","grey30","grey15","black")),
           boxsize = 0.05,
           graphwidth=unit(6,"cm"),
           xlog=TRUE,
           xticks=c(0.25,1.0,4.0),
           clip=c(0.25,1.0,4.0),
           txt_gp = fpTxtGp(xlab=gpar(cex=1.0),ticks=gpar(cex=0.75)),
           xlab="OR per SD increase in JIA PRS",
           hrzl_lines = list("2" = gpar(lty = 1),
                             "6" = gpar(lty = 1),
                             "8" = gpar(lty = 1)),
           lineheight = unit(1.1,"cm"),
           legend = c("1e-2", "1e-3","1e-4","1e-5","1e-6","1e-7","5e-8"),
           legend_args = fpLegend(title="P value threshold",pos=list(x=1.25,y=0.75)))

#Continuous forest plot 
x<-c("Systolic BP1","Systolic BP2", "Systolic BP3","Systolic BP4","Systolic BP5","Systolic BP6","Systolic BP7",
     "Diastolic BP1","Diastolic BP2","Diastolic BP3","Diastolic BP4","Diastolic BP5","Diastolic BP6","Diastolic BP7",
     "Triglycerides1","Triglycerides2","Triglycerides3","Triglycerides4","Triglycerides5","Triglycerides6","Triglycerides7",
     "HDL1","HDL2","HDL3","HDL4","HDL5","HDL6","HDL7",
     "LDL1","LDL2","LDL3","LDL4","LDL5","LDL6","LDL7",
     "Total cholesterol1","Total cholesterol2","Total cholesterol3","Total cholesterol4","Total cholesterol5","Total cholesterol6","Total cholesterol7",
     "Apo-AI1","Apo-AI2","Apo-AI3","Apo-AI4","Apo-AI5","Apo-AI6","Apo-AI7",
     "Apo-B1","Apo-B2","Apo-B3","Apo-B4","Apo-B5","Apo-B6","Apo-B7",
     "Apo-B:AI1","Apo-B:AI2","Apo-B:AI3","Apo-B:AI4","Apo-B:AI5","Apo-B:AI6","Apo-B:AI7",
     "Blood glucose1","Blood glucose2","Blood glucose3","Blood glucose4","Blood glucose5","Blood glucose6","Blood glucose7",
     "Blood insulin1","Blood insulin2","Blood insulin3","Blood insulin4","Blood insulin5","Blood insulin6","Blood insulin7",
     "HOMA2_IR1","HOMA2_IR2","HOMA2_IR3","HOMA2_IR4","HOMA2_IR5","HOMA2_IR6","HOMA2_IR7",
     "   Glycoprotein acetylation1","   Glycoprotein acetylation2","   Glycoprotein acetylation3","   Glycoprotein acetylation4","   Glycoprotein acetylation5","   Glycoprotein acetylation6","   Glycoprotein acetylation7",
     "log hsCRP1","log hsCRP2","log hsCRP3","log hsCRP4","log hsCRP5","log hsCRP6","log hsCRP7",
     "Waist circumference1","Waist circumference2","Waist circumference3","Waist circumference4","Waist circumference5","Waist circumference6","Waist circumference7",
     "FMI1","FMI2","FMI3","FMI4","FMI5","FMI6","FMI7",
     "BMI1","BMI2","BMI3","BMI4","BMI5","BMI6","BMI7",
     "cIMT1","cIMT2","cIMT3","cIMT4","cIMT5","cIMT6","cIMT7",
     "Pulse wave velocity1","Pulse wave velocity2","Pulse wave velocity3","Pulse wave velocity4","Pulse wave velocity5","Pulse wave velocity6","Pulse wave velocity7",
     "Mitral E/A ratio1","Mitral E/A ratio2","Mitral E/A ratio3","Mitral E/A ratio4","Mitral E/A ratio5","Mitral E/A ratio6","Mitral E/A ratio7",
     "E/e'1","E/e'2","E/e'3","E/e'4","E/e'5","E/e'6","E/e'7",
     "Left atrial diameter1","Left atrial diameter2","Left atrial diameter3","Left atrial diameter4","Left atrial diameter5","Left atrial diameter6","Left atrial diameter7",
     "Ejection fraction1","Ejection fraction2","Ejection fraction3","Ejection fraction4","Ejection fraction5","Ejection fraction6","Ejection fraction7",
     "Fractional shortening1","Fractional shortening2","Fractional shortening3","Fractional shortening4","Fractional shortening5","Fractional shortening6","Fractional shortening7",
     "LVMI1","LVMI2","LVMI3","LVMI4","LVMI5","LVMI6","LVMI7")
rownames(phewas_con)<-phewas_con$sort
phewas_con<-phewas_con[x,]

#unadj and PCS/sex adj only
tabletext <- cbind(c("Category", "Blood pressure", NA,
                     "Blood biomarkers",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                     "Anthropometry",NA,NA,"Early atherosclerosis/","arteriosclerosis","Cardiac structure and function",NA,NA,NA,NA,NA),
                   c("Outcome (24yrs)",phewas_con$label[1],phewas_con$label[8],phewas_con$label[15],phewas_con$label[22],phewas_con$label[29],
                     phewas_con$label[36],phewas_con$label[43],phewas_con$label[50],phewas_con$label[57],phewas_con$label[64],phewas_con$label[71],
                     phewas_con$label[78],phewas_con$label[85],phewas_con$label[92],phewas_con$label[99],phewas_con$label[106],phewas_con$label[113],
                     phewas_con$label[120],phewas_con$label[127],phewas_con$label[134],phewas_con$label[141],phewas_con$label[148],phewas_con$label[155],
                     phewas_con$label[162],phewas_con$label[169]),
                   c("Sample size",phewas_con$n[1],phewas_con$n[8],phewas_con$n[15],phewas_con$n[22],phewas_con$n[29],
                     phewas_con$n[36],phewas_con$n[43],phewas_con$n[50],phewas_con$n[57],phewas_con$n[64],phewas_con$n[71],
                     phewas_con$n[78],phewas_con$n[85],phewas_con$n[92],phewas_con$n[99],phewas_con$n[106],phewas_con$n[113],
                     phewas_con$n[120],phewas_con$n[127],phewas_con$n[134],phewas_con$n[141],phewas_con$n[148],phewas_con$n[155],
                     phewas_con$n[162],phewas_con$n[169]))


forestplot(tabletext, 
           mean = cbind(c(NA,phewas_con$est[c(1,8,15,22,29,36,43,50,57,64,71,78,85,92,99,106,113,120,127,134,141,148,155,162,169)]),
                        c(NA,phewas_con$est[c(2,9,16,23,30,37,44,51,58,65,72,79,86,93,100,107,114,121,128,135,142,149,156,163,170)]),
                        c(NA,phewas_con$est[c(3,10,17,24,31,38,45,52,59,66,73,80,87,94,101,108,115,122,129,136,143,150,157,164,171)]),
                        c(NA,phewas_con$est[c(4,11,18,25,32,39,46,53,60,67,74,81,88,95,102,109,116,123,130,137,144,151,158,165,172)]),
                        c(NA,phewas_con$est[c(5,12,19,26,33,40,47,54,61,68,75,82,89,96,103,110,117,124,131,138,145,152,159,166,173)]),
                        c(NA,phewas_con$est[c(6,13,20,27,34,41,48,55,62,69,76,83,90,97,104,111,118,125,132,139,146,153,160,167,174)]),
                        c(NA,phewas_con$est[c(7,14,21,28,35,42,49,56,63,70,77,84,91,98,105,112,119,126,133,140,147,154,161,168,175)])),
           lower = cbind(c(NA,phewas_con$ci.l[c(1,8,15,22,29,36,43,50,57,64,71,78,85,92,99,106,113,120,127,134,141,148,155,162,169)]),
                         c(NA,phewas_con$ci.l[c(2,9,16,23,30,37,44,51,58,65,72,79,86,93,100,107,114,121,128,135,142,149,156,163,170)]),
                         c(NA,phewas_con$ci.l[c(3,10,17,24,31,38,45,52,59,66,73,80,87,94,101,108,115,122,129,136,143,150,157,164,171)]),
                         c(NA,phewas_con$ci.l[c(4,11,18,25,32,39,46,53,60,67,74,81,88,95,102,109,116,123,130,137,144,151,158,165,172)]),
                         c(NA,phewas_con$ci.l[c(5,12,19,26,33,40,47,54,61,68,75,82,89,96,103,110,117,124,131,138,145,152,159,166,173)]),
                         c(NA,phewas_con$ci.l[c(6,13,20,27,34,41,48,55,62,69,76,83,90,97,104,111,118,125,132,139,146,153,160,167,174)]),
                         c(NA,phewas_con$ci.l[c(7,14,21,28,35,42,49,56,63,70,77,84,91,98,105,112,119,126,133,140,147,154,161,168,175)])),
           upper = cbind(c(NA,phewas_con$ci.u[c(1,8,15,22,29,36,43,50,57,64,71,78,85,92,99,106,113,120,127,134,141,148,155,162,169)]),
                         c(NA,phewas_con$ci.u[c(2,9,16,23,30,37,44,51,58,65,72,79,86,93,100,107,114,121,128,135,142,149,156,163,170)]),
                         c(NA,phewas_con$ci.u[c(3,10,17,24,31,38,45,52,59,66,73,80,87,94,101,108,115,122,129,136,143,150,157,164,171)]),
                         c(NA,phewas_con$ci.u[c(4,11,18,25,32,39,46,53,60,67,74,81,88,95,102,109,116,123,130,137,144,151,158,165,172)]),
                         c(NA,phewas_con$ci.u[c(5,12,19,26,33,40,47,54,61,68,75,82,89,96,103,110,117,124,131,138,145,152,159,166,173)]),
                         c(NA,phewas_con$ci.u[c(6,13,20,27,34,41,48,55,62,69,76,83,90,97,104,111,118,125,132,139,146,153,160,167,174)]),
                         c(NA,phewas_con$ci.u[c(7,14,21,28,35,42,49,56,63,70,77,84,91,98,105,112,119,126,133,140,147,154,161,168,175)])),
           col = fpColors(box = c("grey90","grey75","grey60","grey45","grey30","grey15","black")),
           boxsize = 0.05,
           graphwidth=unit(6,"cm"),
           xlog=FALSE,
           txt_gp = fpTxtGp(xlab=gpar(cex=1.0),ticks=gpar(cex=0.75)),
           xlab="SD increase per SD increase in JIA PRS",
           hrzl_lines = list("2" = gpar(lty = 1),
                             "4" = gpar(lty = 1),
                             "16" = gpar(lty = 1),
                             "19" = gpar(lty = 1),
                             "21" = gpar(lty = 1)),
           xticks=c(-0.10,-0.05,0,0.05,0.10),
           clip=c(-0.10,0.10),
           lineheight = unit(1.0,"cm"),
           legend = c("1e-2", "1e-3","1e-4","1e-5","1e-6","1e-7","5e-8"),
           legend_args = fpLegend(title="P value threshold",pos=list(x=1.3,y=0.9)))

write.table(phewas_con,"",row.names=FALSE, quote=FALSE) #add file path
write.table(phewas_bin,"",row.names=FALSE, quote=FALSE) #add file path


#####Analysis age 7 years to 24 years####
#Add in dummy variables where variable not measured at specific timepoint
card$anthro_fmi_stage3_zscore_continuous<-NA
card$anthro_fmi_stage7_zscore_continuous<-NA
card$cardio_insulin_stage3_zscore_continuous<-NA
card$cardio_insulin_stage5_zscore_continuous<-NA
card$cardio_insulin_stage6_zscore_continuous<-NA
card$immuno_CRP_log_stage3_no_outliers_zscore_continuous<-NA
card$immuno_CRP_log_stage5_no_outliers_zscore_continuous<-NA
card$immuno_CRP_log_stage6_no_outliers_zscore_continuous<-NA
card$cardio_IR_stage3_zscore_continuous<-NA
card$cardio_IR_stage4_zscore_continuous<-NA
card$cardio_IR_stage5_zscore_continuous<-NA
card$cardio_IR_stage6_zscore_continuous<-NA

#####All stages
outcomes<-names(card[grepl("dbp|log|insulin|IR|waist|fmi|bmi",names(card))])
outcomes<-outcomes[grepl("dbp|insulin|IR|waist|fmi|bmi|outliers",outcomes)]
outcomes<-outcomes[grepl("zscore",outcomes)]
outcomes<-outcomes[!grepl("ish|RTA|covs|head_cir|stage1|stage2|stage0", outcomes)]
covariates1<-names(card[,c(316:325,1)]) #10 pcs, sex
key <- data.frame(outcome=outcomes,
                  var.type=ifelse(!grepl("binary",outcomes),"continuous","binary")
)

#Run PheWAS
phewas_res_unadj <- run_analysis(exposure="Pt_1e.05MHC_zscore_continuous",
                                 outcomes=outcomes,
                                 covariates = NULL,
                                 df=card)

phewas_res_adj1 <- run_analysis(exposure="Pt_1e.05MHC_zscore_continuous",
                                outcomes=outcomes,
                                covariates = covariates1,
                                df=card)

#combine regression model results
phewas_all<-NULL
phewas_all<-left_join(phewas_res_unadj,phewas_res_adj1,by=c("outcome","exposure"))

#relabel data
phewas_all$label[grepl("anthro_waist_",phewas_all$outcome)]<- "Waist circumference"
phewas_all$label[grepl("anthro_bmi_stage",phewas_all$outcome)]<- "BMI"
phewas_all$label[grepl("anthro_fmi_stage",phewas_all$outcome)]<- "Fat mass index"
phewas_all$label[grepl("cardio_dbp_stage",phewas_all$outcome)]<- "Diastolic BP"
phewas_all$label[grepl("cardio_insulin_stage",phewas_all$outcome)]<- "Blood insulin"
phewas_all$label[grepl("immuno_CRP_log_stage",phewas_all$outcome)]<- "logCRP"
phewas_all$label[grepl("cardio_IR_stage",phewas_all$outcome)]<- "Insulin resistance"

phewas_all$stage<-NA
phewas_all$stage[grep("stage9",phewas_all$outcome)]<-9
phewas_all$stage[grep("stage8",phewas_all$outcome)]<-8
phewas_all$stage[grep("stage7",phewas_all$outcome)]<-7
phewas_all$stage[grep("stage6",phewas_all$outcome)]<-6
phewas_all$stage[grep("stage5",phewas_all$outcome)]<-5
phewas_all$stage[grep("stage4",phewas_all$outcome)]<-4
phewas_all$stage[grep("stage3",phewas_all$outcome)]<-3

phewas_all$age<-NA
phewas_all$age[grep("stage9",phewas_all$outcome)]<-"24yrs"
phewas_all$age[grep("stage8",phewas_all$outcome)]<-"17yrs"
phewas_all$age[grep("stage7",phewas_all$outcome)]<-"15yrs"
phewas_all$age[grep("stage6",phewas_all$outcome)]<-"13yrs"
phewas_all$age[grep("stage5",phewas_all$outcome)]<-"11yrs"
phewas_all$age[grep("stage4",phewas_all$outcome)]<-"9yrs"
phewas_all$age[grep("stage3",phewas_all$outcome)]<-"7yrs"

phewas_bin<-phewas_all[grep("binary",phewas_all$outcome)]
phewas_con<-phewas_all[grep("continuous",phewas_all$outcome),]

#Reorder data
phewas_bin<-arrange(phewas_bin,label,stage)
phewas_con<-arrange(phewas_con,label,stage)

write.table(phewas_all,"",row.names=FALSE, quote=FALSE) #add file path

##Forest plot results##
#Diastolic
tabletext <- cbind(c("Age",phewas_con$age[grep("dbp",phewas_con$outcome)]))
forestplot(tabletext, 
           mean = cbind(c(NA,phewas_con$est.x[grep("dbp",phewas_con$outcome)]),c(NA, phewas_con$est.y[grep("dbp",phewas_con$outcome)])),
           lower = cbind(c(NA,phewas_con$ci.l.x[grep("dbp",phewas_con$outcome)]),c(NA, phewas_con$ci.l.y[grep("dbp",phewas_con$outcome)])),
           upper = cbind(c(NA,phewas_con$ci.u.x[grep("dbp",phewas_con$outcome)]),c(NA,phewas_con$ci.u.y[grep("dbp",phewas_con$outcome)])),
           col = fpColors(box = c("grey", "black")),
           boxsize = 0.1,
           graphwidth=unit(6,"cm"),
           xlog=FALSE,
           txt_gp = fpTxtGp(xlab=gpar(cex=1.0),ticks=gpar(cex=0.75)),
           xlab="SD increase per SD increase in JIA PRS",
           hrzl_lines = list("2" = gpar(lty = 1)),
           xticks=c(-0.10,-0.05,0,0.05,0.10),
           clip=c(-0.10,0.10),
           lineheight = unit(1,"cm"),
           title="Diastolic BP")


#logCRP
tabletext <- cbind(c("Age",phewas_con$age[grep("log",phewas_con$outcome)]))
forestplot(tabletext, 
           mean = cbind(c(NA,phewas_con$est.x[grep("log",phewas_con$outcome)]),c(NA, phewas_con$est.y[grep("log",phewas_con$outcome)])),
           lower = cbind(c(NA,phewas_con$ci.l.x[grep("log",phewas_con$outcome)]),c(NA, phewas_con$ci.l.y[grep("log",phewas_con$outcome)])),
           upper = cbind(c(NA,phewas_con$ci.u.x[grep("log",phewas_con$outcome)]),c(NA,phewas_con$ci.u.y[grep("log",phewas_con$outcome)])),
           col = fpColors(box = c("grey", "black")),
           boxsize = 0.1,
           graphwidth=unit(6,"cm"),
           xlog=FALSE,
           txt_gp = fpTxtGp(xlab=gpar(cex=1.0),ticks=gpar(cex=0.75)),
           xlab="SD increase per SD increase in JIA PRS",
           hrzl_lines = list("2" = gpar(lty = 1)),
           xticks=c(-0.10,-0.05,0,0.05,0.10),
           clip=c(-0.10,0.10),
           lineheight = unit(1,"cm"),
           title="log hsCRP")


#Waist circ
tabletext <- cbind(c("Age",phewas_con$age[grep("waist",phewas_con$outcome)]))
forestplot(tabletext, 
           mean = cbind(c(NA,phewas_con$est.x[grep("waist",phewas_con$outcome)]),c(NA, phewas_con$est.y[grep("waist",phewas_con$outcome)])),
           lower = cbind(c(NA,phewas_con$ci.l.x[grep("waist",phewas_con$outcome)]),c(NA, phewas_con$ci.l.y[grep("waist",phewas_con$outcome)])),
           upper = cbind(c(NA,phewas_con$ci.u.x[grep("waist",phewas_con$outcome)]),c(NA,phewas_con$ci.u.y[grep("waist",phewas_con$outcome)])),
           col = fpColors(box = c("grey", "black")),
           boxsize = 0.1,
           graphwidth=unit(6,"cm"),
           xlog=FALSE,
           txt_gp = fpTxtGp(xlab=gpar(cex=1.0),ticks=gpar(cex=0.75)),
           xlab="SD increase per SD increase in JIA PRS",
           hrzl_lines = list("2" = gpar(lty = 1)),
           xticks=c(-0.10,-0.05,0,0.05,0.10),
           clip=c(-0.10,0.10),
           lineheight = unit(1,"cm"),
           title="Waist circumference")

#FMI
tabletext <- cbind(c("Age",phewas_con$age[grep("fmi",phewas_con$outcome)]))
forestplot(tabletext, 
           mean = cbind(c(NA,phewas_con$est.x[grep("fmi",phewas_con$outcome)]),c(NA, phewas_con$est.y[grep("fmi",phewas_con$outcome)])),
           lower = cbind(c(NA,phewas_con$ci.l.x[grep("fmi",phewas_con$outcome)]),c(NA, phewas_con$ci.l.y[grep("fmi",phewas_con$outcome)])),
           upper = cbind(c(NA,phewas_con$ci.u.x[grep("fmi",phewas_con$outcome)]),c(NA,phewas_con$ci.u.y[grep("fmi",phewas_con$outcome)])),
           col = fpColors(box = c("grey", "black")),
           boxsize = 0.1,
           graphwidth=unit(6,"cm"),
           xlog=FALSE,
           txt_gp = fpTxtGp(xlab=gpar(cex=1.0),ticks=gpar(cex=0.75)),
           xlab="SD increase per SD increase in JIA PRS",
           hrzl_lines = list("2" = gpar(lty = 1)),
           xticks=c(-0.10,-0.05,0,0.05,0.10),
           clip=c(-0.10,0.10),
           lineheight = unit(1,"cm"),
           title="FMI")

#BMI
tabletext <- cbind(c("Age",phewas_con$age[grep("bmi",phewas_con$outcome)]))
forestplot(tabletext, 
           mean = cbind(c(NA,phewas_con$est.x[grep("bmi",phewas_con$outcome)]),c(NA, phewas_con$est.y[grep("bmi",phewas_con$outcome)])),
           lower = cbind(c(NA,phewas_con$ci.l.x[grep("bmi",phewas_con$outcome)]),c(NA, phewas_con$ci.l.y[grep("bmi",phewas_con$outcome)])),
           upper = cbind(c(NA,phewas_con$ci.u.x[grep("bmi",phewas_con$outcome)]),c(NA,phewas_con$ci.u.y[grep("bmi",phewas_con$outcome)])),
           col = fpColors(box = c("grey", "black")),
           boxsize = 0.1,
           graphwidth=unit(6,"cm"),
           xlog=FALSE,
           txt_gp = fpTxtGp(xlab=gpar(cex=1.0),ticks=gpar(cex=0.75)),
           xlab="SD increase per SD increase in JIA PRS",
           hrzl_lines = list("2" = gpar(lty = 1)),
           xticks=c(-0.10,-0.05,0,0.05,0.10),
           clip=c(-0.10,0.10),
           lineheight = unit(1,"cm"),
           title="BMI")

#insulin
tabletext <- cbind(c("Age",phewas_con$age[grep("insulin",phewas_con$outcome)]))
forestplot(tabletext, 
           mean = cbind(c(NA,phewas_con$est.x[grep("insulin",phewas_con$outcome)]),c(NA, phewas_con$est.y[grep("insulin",phewas_con$outcome)])),
           lower = cbind(c(NA,phewas_con$ci.l.x[grep("insulin",phewas_con$outcome)]),c(NA, phewas_con$ci.l.y[grep("insulin",phewas_con$outcome)])),
           upper = cbind(c(NA,phewas_con$ci.u.x[grep("insulin",phewas_con$outcome)]),c(NA,phewas_con$ci.u.y[grep("insulin",phewas_con$outcome)])),
           col = fpColors(box = c("grey", "black")),
           boxsize = 0.1,
           graphwidth=unit(6,"cm"),
           xlog=FALSE,
           txt_gp = fpTxtGp(xlab=gpar(cex=1.0),ticks=gpar(cex=0.75)),
           xlab="SD increase per SD increase in JIA PRS",
           hrzl_lines = list("2" = gpar(lty = 1)),
           xticks=c(-0.10,-0.05,0,0.05,0.10),
           clip=c(-0.10,0.10),
           lineheight = unit(1,"cm"),
           title="Blood insulin")

#insulin resistance
tabletext <- cbind(c("Age",phewas_con$age[grep("IR",phewas_con$outcome)]))
forestplot(tabletext, 
           mean = cbind(c(NA,phewas_con$est.x[grep("IR",phewas_con$outcome)]),c(NA, phewas_con$est.y[grep("IR",phewas_con$outcome)])),
           lower = cbind(c(NA,phewas_con$ci.l.x[grep("IR",phewas_con$outcome)]),c(NA, phewas_con$ci.l.y[grep("IR",phewas_con$outcome)])),
           upper = cbind(c(NA,phewas_con$ci.u.x[grep("IR",phewas_con$outcome)]),c(NA,phewas_con$ci.u.y[grep("IR",phewas_con$outcome)])),
           col = fpColors(box = c("grey", "black")),
           boxsize = 0.1,
           graphwidth=unit(6,"cm"),
           xlog=FALSE,
           txt_gp = fpTxtGp(xlab=gpar(cex=1.0),ticks=gpar(cex=0.75)),
           xlab="SD increase per SD increase in JIA PRS",
           hrzl_lines = list("2" = gpar(lty = 1)),
           xticks=c(-0.10,-0.05,0,0.05,0.10),
           clip=c(-0.10,0.10),
           lineheight = unit(1,"cm"),
           title="HOMA2_IR")

#####Age sensitivity analysis attended all clinics only####
cardagesens<-card[!is.na(card$covs_age_child_stage9_f24_continuous),]
cardagesens<-cardagesens[!is.na(cardagesens$covs_age_child_stage8_tf4_continuous),]
cardagesens<-cardagesens[!is.na(cardagesens$covs_age_child_stage7_tf3_continuous),]
cardagesens<-cardagesens[!is.na(cardagesens$covs_age_child_stage6_tf2_continuous),]
cardagesens<-cardagesens[!is.na(cardagesens$covs_age_child_stage5_f11_continuous),]
cardagesens<-cardagesens[!is.na(cardagesens$covs_age_child_stage4_f9_continuous),]
cardagesens<-cardagesens[!is.na(cardagesens$covs_age_child_stage3_f7_continuous),]

outcomes<-names(cardagesens[grepl("dbp|log|insulin|IR|waist|fmi|bmi",names(cardagesens))])
outcomes<-outcomes[grepl("dbp|insulin|IR|waist|fmi|bmi|outliers",outcomes)]
outcomes<-outcomes[grepl("zscore",outcomes)]
outcomes<-outcomes[!grepl("ish|RTA|covs|head_cir|stage1|stage2|stage0", outcomes)]
covariates1<-names(cardagesens[,c(316:325,1)]) #10 pcs, sex
key <- data.frame(outcome=outcomes,
                  var.type=ifelse(!grepl("binary",outcomes),"continuous","binary")
)

#Run PheWAS
phewas_res_unadj <- run_analysis(exposure="Pt_1e.05MHC_zscore_continuous",
                                 outcomes=outcomes,
                                 covariates = NULL,
                                 df=cardagesens)

phewas_res_adj1 <- run_analysis(exposure="Pt_1e.05MHC_zscore_continuous",
                                outcomes=outcomes,
                                covariates = covariates1,
                                df=cardagesens)

#combine regression model results
phewas_all<-NULL
phewas_all<-left_join(phewas_res_unadj,phewas_res_adj1,by=c("outcome","exposure"))


#relabel data
phewas_all$label<-NA
phewas_all$label[grep("anthro_waist_",phewas_all$outcome)]<- "Waist circumference"
phewas_all$label[grep("anthro_bmi_stage",phewas_all$outcome)]<- "BMI"
phewas_all$label[grep("anthro_fmi_stage",phewas_all$outcome)]<- "Fat mass index"
phewas_all$label[grep("cardio_dbp_stage",phewas_all$outcome)]<- "Diastolic BP"
phewas_all$label[grep("cardio_insulin_stage",phewas_all$outcome)]<- "Blood insulin"
phewas_all$label[grep("immuno_CRP_log_stage",phewas_all$outcome)]<- "logCRP"
phewas_all$label[grep("cardio_IR_stage",phewas_all$outcome)]<- "Insulin resistance"

phewas_all$stage<-NA
phewas_all$stage[grep("stage9",phewas_all$outcome)]<-9
phewas_all$stage[grep("stage8",phewas_all$outcome)]<-8
phewas_all$stage[grep("stage7",phewas_all$outcome)]<-7
phewas_all$stage[grep("stage6",phewas_all$outcome)]<-6
phewas_all$stage[grep("stage5",phewas_all$outcome)]<-5
phewas_all$stage[grep("stage4",phewas_all$outcome)]<-4
phewas_all$stage[grep("stage3",phewas_all$outcome)]<-3


phewas_all$age<-NA
phewas_all$age[grep("stage9",phewas_all$outcome)]<-"24yrs"
phewas_all$age[grep("stage8",phewas_all$outcome)]<-"17yrs"
phewas_all$age[grep("stage7",phewas_all$outcome)]<-"15yrs"
phewas_all$age[grep("stage6",phewas_all$outcome)]<-"13yrs"
phewas_all$age[grep("stage5",phewas_all$outcome)]<-"11yrs"
phewas_all$age[grep("stage4",phewas_all$outcome)]<-"9yrs"
phewas_all$age[grep("stage3",phewas_all$outcome)]<-"7yrs"



phewas_bin<-phewas_all[grep("binary",phewas_all$outcome)]
phewas_con<-phewas_all[grep("continuous",phewas_all$outcome),]

#Reorder data
phewas_bin<-arrange(phewas_bin,label,stage)
phewas_con<-arrange(phewas_con,label,stage)

write.table(phewas_con,"",row.names=FALSE,quote=FALSE) #add file path


