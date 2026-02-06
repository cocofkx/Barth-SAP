#-----------------------------------Barth Syndrome Survival Analysis-----------------------------------#
#Editor: Kexin Fu
#Session:
#       1.Data and Cohorts 
#       2.Descriptive analysis
#       3.Survival analysis
#       4.Death causes
#       5.Symptomatology


#-----------------------------Data and Cohorts-----------------------------#
#load packages
library(tidyverse)
library(Hmisc)
library(expss)
library(psych)
library(table1)
library(readxl)
library(MASS)
library(survival)
library(ggplot2)
library(broom)
library(patchwork) 
library(tranSurv)
library(mice)
library(binom)
library(purrr)
library(tibble)
#-------------Set up-------------#
#Upload main intake data (affetced individuals)
bf.intake <- read.csv("~/Barth_Syndrome_Data/BSFintake_known_812.csv")
#Drop one duplicates
length(unique(bf.intake$Intake_ID))
bf.intake<-bf.intake%>%filter(Patient_ID!=447)

#Upload gender data
bf_gender<-read_csv('~/Barth_Syndrome_Data/bsfdata_gender.csv')
bf_gender<-bf_gender%>%dplyr::select(Intake_ID,Gender)
#merge with intake data
bf.intake<-merge(bf.intake,bf_gender,by='Intake_ID', all.x = TRUE)

#Upload intake data (published cases)
bf.intake.na<-read_csv('~/Barth_Syndrome_Data/BSFintake_unknown_812.csv')
#indicator of known or unknown to BSF
bf.intake<-bf.intake%>%mutate(bsf_known=1)
bf.intake.na<-bf.intake.na%>%mutate(bsf_known=0)
#fill registry date, registry ID, join date, last contact year, gender information with NA
bf.intake.na<-bf.intake.na%>%mutate(Registry_Entry_Date=NA)
bf.intake.na<-bf.intake.na%>%mutate(Registry_ID=NA)
bf.intake.na<-bf.intake.na%>%mutate(DateJoined=NA)
bf.intake.na<-bf.intake.na%>%mutate(Last_Contact_Year=NA)
bf.intake.na<-bf.intake.na%>%mutate(Gender=NA)
#merge with main intake data
bf.intake.all<-rbind(bf.intake,bf.intake.na)

#Drop empty ID -->Total sample: 561
bf.intake.total<-bf.intake.all%>%filter(!is.na(bf.intake.all$Intake_ID))
table(bf.intake.total$bsf_known) 

#Recode variables: transplant, year birth, stroke, death, developed country, continent, death when wait transplant, age at first transplant, diagnosis status
bf.intake.total$YearBirth<-as.numeric(bf.intake.total$YearBirth)
bf.intake.total$YearTX1<-as.numeric(bf.intake.total$YearTX1)
bf.intake.total$tx1<-ifelse(is.na(bf.intake.total$YearTX1),'0','1')
bf.intake.total$tx2<-ifelse(is.na(bf.intake.total$YearTX2),'0','1')
bf.intake.total$tx3<-ifelse(is.na(bf.intake.total$YearTX3),'0','1')
bf.intake.total$stroke<-ifelse(is.na(bf.intake.total$Stroke)|bf.intake.total$Stroke=='','0','1')
bf.intake.total$externalfeed<-ifelse(is.na(bf.intake.total$EnteralFeeds)|bf.intake.total$EnteralFeeds=='','0','1')
#bf.intake.total$death<-ifelse(is.na(bf.intake.total$YearDeath),'0','1')
list<-c('Australia','Austria','Belgium','Canada','Czech Republic','Denmark','england','England','France','Germany',
        'Iceland','Ireland','Israel','Italy','Japan','Korea','Lithuania','Netherlands','New Zealand','Portugal',
        'Saudi Arabia','Scotland','Slovakia','Spain','Switzerland','Swizterland','United States','Russia','Poland','United Arab Emirates')
bf.intake.total$developed<-ifelse(bf.intake.total$Country%in% list,'1','0')
bf.intake.total$Continent[bf.intake.total$Continent=='']<-NA
bf.intake.total$DWWTX1[bf.intake.total$DWWTX1=='']<-NA #death when wait transplant
bf.intake.total$agetx<-bf.intake.total$YearTX1-bf.intake.total$YearBirth

#-------------Define cohorts-------------#
##Cohort details can be found in Figure 1

#Recode variables: age at death, current age, year of birth, year of death
bf.intake.total$AgeDeath[bf.intake.total$AgeDeath=='unknown'|bf.intake.total$AgeDeath=='']<-NA
bf.intake.total$CurrentAge[bf.intake.total$CurrentAge=='unknown'|bf.intake.total$CurrentAge=='']<-NA
bf.intake.total$YearBirth[bf.intake.total$YearBirth=='unknown'|bf.intake.total$YearBirth=='']<-NA
bf.intake.total$YearDeath[bf.intake.total$YearDeath=='unknown'|bf.intake.total$YearDeath=='']<-NA

#Check # of missing in time
sum(complete.cases(bf.intake.total$AgeDeath))#202 complete (death overall sample:202)
sum(complete.cases(bf.intake.total$YearDeath))#confirm 202 complete (death overall sample:202)
sum(complete.cases(bf.intake.total$YearBirth))#505 complete

#Calculate death (secondary event) time using ceiling YearDeath -Yearbirth / ceiling current age
bf.intake.total$YearDeath<-as.numeric(bf.intake.total$YearDeath)
bf.intake.total$death_ceil<-bf.intake.total$YearDeath-bf.intake.total$YearBirth + 1
bf.intake.total$CurrentAge<-as.numeric(bf.intake.total$CurrentAge)
bf.intake.total$live_ceil<-ceiling(bf.intake.total$CurrentAge)
bf.intake.total$time<-ifelse(!is.na(bf.intake.total$YearDeath),bf.intake.total$death_ceil,bf.intake.total$live_ceil)
bf.intake.total$censor<-ifelse(!is.na(bf.intake.total$YearDeath),'0','1')

#Exclude unknown current age or age at death --> 1. Overall cohort: 502
intermediate<-bf.intake.total%>%filter(!is.na(time))
table(intermediate$bsf_known)

#Recode variables: registry entry date, date of join the registry
intermediate$Registry_Entry_Date[intermediate$Registry_Entry_Date=='']<-NA
intermediate$DateJoined[intermediate$DateJoined=='']<-NA
#Exclude unregistered inidividuals --> 2. BRR cohort: 162
sample<-intermediate[!is.na(intermediate$Registry_Entry_Date)|!is.na(intermediate$DateJoined),]
table(sample$censor) #n death=19

#Calculate entry time using Registry_Year - YearBirth
sample$Registry_Year <- ifelse(is.na(sample$Registry_Entry_Date),sub("/.*", "", sample$DateJoined),sub("/.*", "", sample$Registry_Entry_Date))%>%as.numeric()
sample$entry<-sample$Registry_Year-sample$YearBirth

#Exclude post-mortem individual--> 3. BRR longitudinal survival cohort: 150
subset(sample, entry > time)%>%dplyr::select(Intake_ID,censor,time,entry,Registry_Entry_Date,DateJoined,YearBirth,YearDeath,CurrentAge)
subset(sample, entry == time)%>%dplyr::select(Intake_ID,censor,time,entry,Registry_Entry_Date,DateJoined,YearBirth,YearDeath,CurrentAge,live_ceil,death_ceil)
#since we only have birth year, some round-up issue occurs when the person entered cohort in early 2014 while still alive
#their follow-up time should be a positive number ~0.5, and we want to include them into our longitudinal BRR cohort
id<-sample$Intake_ID[sample$entry==sample$time&sample$censor=='1']
sample$time[sample$Intake_ID %in% id]<-sample$time[sample$Intake_ID %in% id]+1
subset_registry<-subset(sample, entry < time)
table(subset_registry$censor) #N death=7 (12 post-mortem cases)

#Define death cohort:202
death<-bf.intake.total%>%filter(!is.na(YearDeath))

#Upload clinical data (at least fill follow-up clinical survey once)
bf.clinical<-read.csv('/Users/mac1/Desktop/Barth_Syndrome_Data/BSFclinical_812.csv')
#merge with BRR data
bf.risk<-merge(bf.clinical,sample,by='Intake_ID', all.x = TRUE)
#Drop empty ID --> 4. Symptom Survey cohort: 125
bf.risk<-subset(bf.risk,Intake_ID!=999) 
length(unique(bf.risk$Intake_ID))
n_per_id <-table(bf.risk$Intake_ID,useNA='always')
sum(n_per_id == 1) # 88 subject have only 1 response survey
mean(n_per_id)

#Upload current symptom data
current<-read.csv('~/Barth_Syndrome_Data/current_symptom812.csv')
#Drop empty ID and missing survey age -->Current symptom Survey cohort: 125
current<-current%>%filter(!is.na(Age.At.Survey))%>%filter(Intake.ID!=999)
length(unique(current$Intake.ID))
#Catergorize age
current$AGE<-ifelse(current$Age.At.Survey<1,'<1',
                    ifelse(current$Age.At.Survey<=5&current$Age.At.Survey>=1,'1-5',
                           ifelse(current$Age.At.Survey<=10&current$Age.At.Survey>5,'5-10',
                                  ifelse(current$Age.At.Survey>10&current$Age.At.Survey<=18,'11-18',
                                         ifelse(current$Age.At.Survey>18&current$Age.At.Survey<=30,'19-30','>30')))))
#Count distinct current symptom report at each age of subjects in current symptom survey cohort: 159
current2<-current
age<-c('<1','1-5','5-10','11-18','19-30','>30')
count<-0
for (i in age) {
  distinct_age_counts <- current2 %>%filter(AGE==i)%>%
    group_by(Intake.ID) %>%
    summarise(distinct_age_count = n_distinct(Age.At.Survey))%>%
    summarise(count = sum(distinct_age_count))
  print(i)
  print(distinct_age_counts)
  count<-count+distinct_age_counts}
count


dist_fre<-current2%>%
  group_by(Intake.ID) %>%
  summarise(distinct_age_count = n_distinct(Age.At.Survey))
table(dist_fre$distinct_age_count)

#-----------------------------Descriptive analysis-----------------------------#
##Descriptive results can be found in Table 1, Table 2 and Supp. Table 2

#-------------Demographic information-------------#
#Describe intake information in Overall Cohort
table1(~ bsf_known+Gender+DiagnosisStatus+YearBirth+Continent+stroke+Continent+agetx+censor+DWWTX1,data=intermediate)
intermediate[which(intermediate$Country=='United States'),]
#Describe intake information in BRR Cohort
table1(~ bsf_known+Gender+DiagnosisStatus+YearBirth+Continent+stroke+Continent+agetx+censor+DWWTX1,data=sample)
sample[which(sample$Country=='United States'),]

#DWWTX = Deceased while waiting for Transplant For the "Diagnosis/Status" column: P = post-mortem diagnosis, D = diagnosed, S = suspected / L = living, D = deceased

#-------------Detailed clinical information-------------#
bf.risk[bf.risk == 999] <- NA

#Define a subset of every first complete record for subjects in symptom survey cohort
bf.risk$Surveytime <- strptime(bf.risk$SurveyTime, format = "%Y/%m/%d %H:%M", tz = "EST")
bf.risk<-bf.risk%>%dplyr::group_by(Intake_ID)%>%mutate(first=ifelse(Surveytime==min(Surveytime),1,0))
list<-bf.risk$Intake_ID[is.na(bf.risk$first)]
bf.risk[bf.risk$Intake_ID %in% list,]
bf.risk$first<-ifelse(bf.risk$Intake_ID %in% list,1,bf.risk$first)
bf.risk$first[bf.risk$Intake_ID=='74_1_2'&bf.risk$Registry_ID.x==3347]<-0 #all distinct intake id except for 74_1_2
bf.risk1<-bf.risk%>%
  filter(first==1)
nrow(bf.risk1) #125 records

####Age at regsitry, first symptom, first diagnosis###
bf.risk1$FirstAgegroup<-ifelse(bf.risk1$FirstAge=='Prenatal'|
                                 bf.risk1$FirstAge=='Less than 1 year'|
                                 bf.risk1$FirstAge=='At birth',1, #<1
                               ifelse(bf.risk1$FirstAge=='1 year old'|
                                        bf.risk1$FirstAge=='2 years old'|
                                        bf.risk1$FirstAge=='3 years old'|
                                        bf.risk1$FirstAge=='4 years old'|
                                        bf.risk1$FirstAge=='5 years old',2, #1-5
                                      ifelse(bf.risk1$FirstAge=='6 years old'|
                                               bf.risk1$FirstAge=='8 years old',3, #6-10
                                             ifelse(bf.risk1$FirstAge=='19 years old',4, #19-30
                                                    ifelse(bf.risk1$FirstAge=='48 years old',5,NA))))) #>30
bf.risk1$FirstAgegroup<-factor(bf.risk1$FirstAgegroup)
table1(~entry+FirstAgegroup+DiagnAge,data=bf.risk1)

###first symptom###
list<-bf.risk1$Intake_ID[is.na(bf.risk1$FirstSymptom)]
bf.risk[bf.risk$Intake_ID %in% list,c('Intake_ID','first','FirstSymptom')]
table1(~FirstSymptom,bf.risk1)
firstsymp_data<-bf.risk1%>%filter(!is.na(FirstAgegroup))
table1(~FirstSymptom|FirstAgegroup,firstsymp_data)

###family history###
#blood relatives
bf.risk <- bf.risk %>%
  mutate(across(10:33, as.numeric))
bf.risk$relative_blood<-rowSums(bf.risk[, 11:34], na.rm = TRUE)
bf.risk$relative_blood_yesno<-ifelse(bf.risk$relative_blood>0,1,0)
bf.risk%>%group_by(Intake_ID)%>% #if any response answer "yes"
  summarise(relative_blood = any(relative_blood_yesno == 1)) %>%
  summarise(count = sum(relative_blood))
bf.risk%>%group_by(Intake_ID)%>%#if all response missing -->n response=125
  summarise(relative_blood = all(is.na(relative_blood_yesno))) %>%
  summarise(count = sum(relative_blood))
#create an indicator variable for any blood relatives for survival analysis
risk_survival_familyR<-bf.risk%>%group_by(Intake_ID)%>%
  summarise(family_relative = ifelse(any(relative_blood_yesno == 1),1,ifelse(all(is.na(relative_blood_yesno) == T),NA,0)))
table(risk_survival_familyR$family_relative,useNA = 'always')

#biological mother
table(bf.risk$RelaK_Mother,useNA = 'always')
bf.risk$relative_mother<-ifelse(bf.risk$RelaK_Mother=='Biological mother is a carrier',1,ifelse(bf.risk$RelaK_Mother=='Biological mother is not a carrier',0,NA))

bf.risk$relative_mother[bf.risk$RelaK_Mother=='Biological mother is a carrier']<-1
bf.risk$relative_mother[bf.risk$RelaK_Mother=='Biological mother is not a carrier']<-2
bf.risk$relative_mother[bf.risk$RelaK_Mother=='Biological mother is an obligate carrier (not tested but has had more than one child with Barth syndrome)']<-1
bf.risk$relative_mother[bf.risk$RelaK_Mother=='Unsure']<-NA
bf.risk%>%group_by(Intake_ID)%>%
  summarise(Relative_mother = any(relative_mother==1)) %>%
  summarise(count = sum(Relative_mother,na.rm=T))
bf.risk%>%group_by(Intake_ID)%>% 
  summarise(Relative_mother = all(is.na(relative_mother))) %>%
  summarise(count = sum(Relative_mother,na.rm=T))

###Cardiac Manifestations###
#heart failure
table(bf.risk$heart_failure,useNA = 'always')
bf.risk%>%group_by(Intake_ID)%>%
  summarise(heart = any(heart_failure == 'Yes')) %>%
  summarise(count = sum(heart,na.rm = T))
bf.risk%>%group_by(Intake_ID)%>%
  summarise(heart = all(heart_failure == 'Unsure'|is.na(heart_failure))) %>%
  summarise(count = sum(heart,na.rm = T))
#create an indicator variable for heart failure for survival analysis
risk_survival_hf<-bf.risk%>%group_by(Intake_ID)%>%
  summarise(heart_failure =ifelse(any(heart_failure=='Yes'),1,ifelse(all(is.na(heart_failure)),NA,0)))
table(risk_survival_hf$heart_failure,useNA = 'always')

#cardiac arrest
table(bf.risk$cardiac_arrest,useNA = 'always')
bf.risk%>%group_by(Intake_ID)%>%
  summarise(cardiac = any(cardiac_arrest == 'Yes')) %>%
  summarise(count = sum(cardiac,na.rm = T))
bf.risk%>%group_by(Intake_ID)%>%
  summarise(cardiac = all(cardiac_arrest == 'Unsure'|is.na(cardiac_arrest))) %>%
  summarise(count = sum(cardiac,na.rm = T))
#create an indicator variable for cardiac arrest for survival analysis
risk_survival_ca<-bf.risk%>%group_by(Intake_ID)%>%
  summarise(cardiac_arrest =ifelse(any(cardiac_arrest=='Yes'),1,
                                   ifelse(all(cardiac_arrest == 'Unsure'|is.na(cardiac_arrest)),NA,0)))
table(risk_survival_ca$cardiac_arrest,useNA = 'always')

#use of devices
table(bf.risk$AICD.Age,useNA = 'always')
bf.risk$AICD.yesno<-1
bf.risk$AICD.yesno[bf.risk$AICD.Age=='']<-NA
bf.risk$AICD.yesno[bf.risk$AICD.Age=='Never had this']<-0
table(bf.risk$AICD.yesno,useNA = 'always')
table(bf.risk$AED.Age,useNA = 'always')
bf.risk$AED.yesno<-1
bf.risk$AED.yesno[bf.risk$AED.Age=='']<-NA
bf.risk$AED.yesno[bf.risk$AED.Age=='Never had this']<-0
table(bf.risk$AED.yesno,useNA = 'always')
bf.risk%>%group_by(Intake_ID)%>% #any AICD
  summarise(aicd = any(AICD.yesno == 1)) %>%
  summarise(count = sum(aicd,na.rm = T))
bf.risk%>%group_by(Intake_ID)%>%#any AED
  summarise(aed = any(AED.yesno == 1)) %>%
  summarise(count = sum(aed,na.rm = T))
bf.risk$Implantable.loop.monitor.Age[bf.risk$Implantable.loop.monitor.Age=='']<-NA
bf.risk$Berlin.Heart.Age[bf.risk$Berlin.Heart.Age=='']<-NA
bf.risk$Left.Ventricular.Assist.Device..LVAD..Age[bf.risk$Left.Ventricular.Assist.Device..LVAD..Age=='']<-NA
bf.risk$Pacemaker.Age[bf.risk$Pacemaker.Age=='']<-NA
bf.risk$Right.Ventricular.Assist.Device..RVAD..Age[bf.risk$Right.Ventricular.Assist.Device..RVAD..Age=='']<-NA
bf.risk$ECMO.Age[bf.risk$ECMO.Age=='']<-NA
bf.risk%>%group_by(Intake_ID)%>% #all missing
  summarise(na = all(is.na(AED.yesno)&is.na(AICD.yesno)&is.na(Implantable.loop.monitor.Age)&is.na(Berlin.Heart.Age)&is.na(Left.Ventricular.Assist.Device..LVAD..Age)&is.na(Pacemaker.Age)&is.na(Right.Ventricular.Assist.Device..RVAD..Age)&is.na(ECMO.Age))) %>%
  summarise(count = sum(na,na.rm = T))
risk_survival_aicdex <- bf.risk %>%
  group_by(Intake_ID) %>%
  summarise(
    aicd_aex = case_when(
      any(AICD_appro == 1, na.rm = TRUE) | any(AED_inappro == 1, na.rm = TRUE) ~ 1,
      any(AICD_neveruse == 1, na.rm = TRUE) | any(AICD_none == 1, na.rm = TRUE) ~ 0,
      TRUE ~ NA_real_
    )
  )

###Blood clot###
table(bf.risk$blood_clot,useNA = 'always')
bf.risk%>%group_by(Intake_ID)%>%
  summarise(blood = any(blood_clot == 'Yes')) %>%
  summarise(count = sum(blood,na.rm = T))
bf.risk%>%group_by(Intake_ID)%>%
  summarise(blood = all(blood_clot == 'Unsure'|is.na(blood_clot))) %>%
  summarise(count = sum(blood,na.rm = T))
#create an indicator variable for blood clot for survival analysis
risk_survival_bc<-bf.risk%>%group_by(Intake_ID)%>%
  summarise(blood_clot =ifelse(any(blood_clot=='Yes'),1,ifelse(all(is.na(blood_clot)),NA,0)))
table(risk_survival_bc$blood_clot,useNA = 'always')

###Neutropenia###
table(bf.risk$AgeNeutropenia,useNA = 'always')
bf.risk$Neutropenia<-NA
bf.risk$Neutropenia[bf.risk$AgeNeutropenia=='Never diagnosed with neutropenia']<-'never'
bf.risk$Neutropenia[bf.risk$AgeNeutropenia=='At birth']<-'<1'
bf.risk$Neutropenia[bf.risk$AgeNeutropenia=='11 - 20 days old'|
                      bf.risk$AgeNeutropenia=='2 - 3 days old'|
                      bf.risk$AgeNeutropenia=='21 - 29 days old'|
                      bf.risk$AgeNeutropenia=='4 - 7 days old'|
                      bf.risk$AgeNeutropenia=='4 - 7 months old'|
                      bf.risk$AgeNeutropenia=='8 - 10 days old'|
                      bf.risk$AgeNeutropenia=='8 - 11 months old'|
                      bf.risk$AgeNeutropenia=='1 - 3 months old'|
                      bf.risk$AgeNeutropenia=='4 - 7 months old']<-'<1'
bf.risk$Neutropenia[bf.risk$AgeNeutropenia=='1 year old'|
                      bf.risk$AgeNeutropenia=='2 years old'|
                      bf.risk$AgeNeutropenia=='3 years old'|
                      bf.risk$AgeNeutropenia=='4 years old'|
                      bf.risk$AgeNeutropenia=='5 years old']<-'1-5'
bf.risk$Neutropenia[bf.risk$AgeNeutropenia=='10 years old'|
                      bf.risk$AgeNeutropenia=='8 years old'|
                      bf.risk$AgeNeutropenia=='9 years old']<-'5-10'
bf.risk$Neutropenia[bf.risk$AgeNeutropenia=='11 years old'|bf.risk$AgeNeutropenia=='16 years old']<-'11-18'
bf.risk$Neutropenia[bf.risk$AgeNeutropenia=='19 years old']<-'19-30'
bf.risk$Neutropenia[bf.risk$AgeNeutropenia=='44 years old'|bf.risk$AgeNeutropenia=='48 years old']<-'>30'
#Define a subset of every last complete record for subjects in symptom survey cohort
bf.risk<-bf.risk%>%dplyr::group_by(Intake_ID)%>%mutate(last=ifelse(Surveytime==max(Surveytime),1,0))
list<-bf.risk$Intake_ID[is.na(bf.risk$last)]
bf.risk$last<-ifelse(bf.risk$Intake_ID %in% list,1,bf.risk$last)
bf.last<-bf.risk%>%filter(last==1)
bf.last<-bf.last%>%filter(Registry_ID.x!=3347)
table(bf.last$Neutropenia,useNA = 'always')

###GI Manifestation###
bf.risk <- bf.risk %>%
mutate(across(77:111, as.numeric))
bf.risk$gi<-rowSums(bf.risk[, 78:110], na.rm = TRUE)
bf.risk$gi_yesno<-ifelse(bf.risk$gi>0,1,ifelse(bf.risk$GIEVER_None==1,0,NA))
bf.risk%>%group_by(Intake_ID)%>%
  summarise(GI = any(gi_yesno == 1)) %>%
  summarise(count = sum(GI,na.rm = T))
bf.risk%>%group_by(Intake_ID)%>%
  summarise(GI = all(GIEVER_None == 1)) %>%
  summarise(count = sum(GI,na.rm = T))
bf.risk%>%group_by(Intake_ID)%>%
  summarise(GI = all(is.na(gi_yesno) == T)) %>%
  summarise(count = sum(GI,na.rm = T))
colSums(bf.risk[, 78:110], na.rm = TRUE)

###Feeding Manifestation###
bf.risk <- bf.risk %>%
  mutate(across(117:128, as.numeric))
bf.risk$feed<-rowSums(bf.risk[, c(118,120:126,129)], na.rm = TRUE)
# bf.risk.feed$feedsum<-rowSums(bf.risk.feed[, 123:126], na.rm = TRUE)
bf.risk$feed_yesno<-ifelse(bf.risk$feed>0,1,ifelse(bf.risk$Feed_None==1|bf.risk$Feed_past==1,0,NA))
bf.risk%>%group_by(Intake_ID)%>%
  summarise(FEED = any(feed_yesno == 1)) %>%
  summarise(count = sum(FEED,na.rm = T))
bf.risk%>%group_by(Intake_ID)%>%
  summarise(FEED = all(feed_yesno == 0)) %>%
  summarise(count = sum(FEED,na.rm = T))
bf.risk%>%group_by(Intake_ID)%>%
  summarise(FEED = all(is.na(feed_yesno) == T)) %>%
  summarise(count = sum(FEED,na.rm = T))
bf.risk%>%group_by(Intake_ID)%>%
  summarise(feed = any(Feed_odors == 1)) %>%
  summarise(count = sum(feed,na.rm = T))
bf.risk%>%group_by(Intake_ID)%>%
  summarise(feed = any(Feed_tastes == 1)) %>%
  summarise(count = sum(feed,na.rm = T))
bf.risk%>%group_by(Intake_ID)%>%
  summarise(feed = any(Feed_textures == 1)) %>%
  summarise(count = sum(feed,na.rm = T))
#create an indicator variable for feeding manifestation for survival analysis
risk_survival_feed<-bf.risk%>%group_by(Intake_ID)%>%
  summarise(feed =ifelse(any(feed_yesno==1),1,ifelse(is.na(feed_yesno),NA,0)))
table(risk_survival_feed$feed,useNA = 'always')

###Endocrine Manifestations###
bf.risk <- bf.risk %>%
  mutate(across(36:75, as.numeric))
bf.risk$endoc<-rowSums(bf.risk[, c(37:73,76)], na.rm = TRUE)
bf.risk$endoc_yesno<-ifelse(bf.risk$endoc>0,1,ifelse(bf.risk$EndocEVER_None==1,0,NA))
bf.risk%>%group_by(Intake_ID)%>%
  summarise(Endoc = any(endoc_yesno == 1)) %>%
  summarise(count = sum(Endoc,na.rm = T))
bf.risk%>%group_by(Intake_ID)%>%
  summarise(Endoc = all(EndocEVER_None == 1)) %>%
  summarise(count = sum(Endoc,na.rm = T))
bf.risk%>%group_by(Intake_ID)%>%
  summarise(Endoc = all(is.na(endoc_yesno) == T)) %>%
  summarise(count = sum(Endoc,na.rm = T))
colSums(bf.risk[, c(37:73,76)], na.rm = TRUE)

###Musculoskeletal Manifestations###
bf.risk <- bf.risk %>%
  mutate(across(130:144, as.numeric))
bf.risk$orth<-rowSums(bf.risk[, c(131:142,145)], na.rm = TRUE)
bf.risk$orth_yesno<-ifelse(bf.risk$orth>0,1,ifelse(bf.risk$EverOrthopedic_None==1,0,NA))
bf.risk%>%group_by(Intake_ID)%>%
  summarise(orth = any(orth_yesno == 1)) %>%
  summarise(count = sum(orth,na.rm = T))
bf.risk%>%group_by(Intake_ID)%>%
  summarise(orth = all(is.na(orth_yesno) == T)) %>%
  summarise(count = sum(orth,na.rm = T))
colSums(bf.risk[, 131:142], na.rm = TRUE)


#-----------------------------Survival analysis-----------------------------#
#Survival results can be found Figure 2, Table 3, Supp. Table 3 and Supp. Table 4

#-------------Overall Cox-------------#
#censoring indicator: death or heart transplant
table(intermediate$censor)
table(intermediate$tx1)
intermediate$censor2<-ifelse(!is.na(intermediate$YearDeath)|intermediate$tx1==1,0,1)
table(intermediate$censor2)
intermediate$censor<-as.numeric(intermediate$censor)
intermediate$censor2<-as.numeric(intermediate$censor2)
#event/censoring time:
intermediate$time2<-ifelse(!is.na(intermediate$YearTX1),intermediate$YearTX1-intermediate$YearBirth + 1,intermediate$time)

#cox regression with composite endpoint
coxfit<-coxph(Surv(time2, 1-censor2) ~  stroke+externalfeed+developed, data = intermediate)
summary(coxfit)

#cox regression with secondary endpoint
coxfit<-coxph(Surv(time, 1-censor) ~  tx1+stroke+externalfeed+developed, data = intermediate)
summary(coxfit)

#-------------BRR Cox-------------#
#censoring indicator: death or heart transplant
table(sample$censor)
table(intermediate$tx1)
sample$censor2<-ifelse(!is.na(sample$YearDeath)|sample$tx1==1,0,1)
table(sample$censor2)
#event/censoring time:
sample$time2<-ifelse(!is.na(sample$YearTX1),sample$YearTX1-sample$YearBirth + 1,sample$time)

sample%>%filter(entry >= time2)%>%dplyr::select(censor2,YearBirth,YearDeath,YearTX1,time2,entry)
#except 12 post-mortal cases, there are 19 post-transplant cases.

# if we get rid of all cases violate entry<=event time
#compute conditional Kendall’s tau values and the associated p-value to test independence
subset_registry2<-subset(sample, entry < time2) # 131 potential left-truncated truncated time strictly smaller than death time
(cKendall(subset_registry2$entry, subset_registry2$time2, 1-subset_registry2$censor2))
(cKendall(subset_registry2$entry, subset_registry2$time2, 1-subset_registry2$censor2, method = "IPW1"))
coxph(Surv(time2,1-as.numeric(censor2))~entry,data=sample)

table(subset_registry2$externalfeed,subset_registry2$censor2)
table(subset_registry2$tx1,subset_registry2$censor2)
table(subset_registry2$stroke,subset_registry2$censor2)

#cox regression with composite endpoint 
coxfit<-coxph(Surv(time2,1-censor2) ~ stroke+developed, data = subset_registry2)
summary(coxfit)

#same for secondary endpoint
#compute conditional Kendall’s tau values and the associated p-value to test independence
subset_registry$censor<-as.numeric(subset_registry$censor)
(cKendall(subset_registry$entry, subset_registry$time, 1-subset_registry$censor))
(cKendall(subset_registry$entry, subset_registry$time, 1-subset_registry$censor, method = "IPW1"))
coxph(Surv(time,1-as.numeric(censor))~entry,data=sample)

table(subset_registry$externalfeed,subset_registry$censor)
table(subset_registry$tx1,subset_registry$censor)
table(subset_registry$stroke,subset_registry$censor)

#cox regression with secondary endpoint 
coxfit<-coxph(Surv(time,1-as.numeric(censor)) ~ developed, data = subset_registry)
summary(coxfit)

#-------------Overall KM-------------#
# Fit KM models
fit1 <- survfit(Surv(time, 1 - censor) ~ 1, data = intermediate)
fit2 <- survfit(Surv(time2, 1 - censor2) ~ 1, data = intermediate)

# Tidy
fit1_df <- tidy(fit1) |> mutate(group = "Death")
fit2_df <- tidy(fit2) |> mutate(group = "Heart transplant or Death")

# Add starting point
fit1_df <- bind_rows(tibble(time = 0, estimate = 1, conf.low = 1, conf.high = 1, group = "Death"), fit1_df)
fit2_df <- bind_rows(tibble(time = 0, estimate = 1, conf.low = 1, conf.high = 1, group = "Heart transplant or Death"), fit2_df)

# Combine
combined_df1 <- bind_rows(fit1_df, fit2_df)

# Build plot
p1 <- ggplot(combined_df1, aes(x = time, y = estimate, color = group, fill = group)) +
  geom_step(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  labs(
    x = "Age", y = "Survival Probability", title = "Overall"
  ) +
  coord_cartesian(xlim = c(0, 50), ylim = c(0, 1)) +
  theme_minimal(base_size = 12) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10)
  ) +
  annotate("label",
           x = 20, y = 0.2,
           hjust = 0,
           vjust = 1,
           label = "N = 502\nN Deaths = 202\nN Heart transplant = 54",
           size = 2.5,
           label.size = 0.3,
           fill = "white",
           color = "black")

#Sensitivity analysis: patients w/o a confirmed diagnosis 
table(intermediate$DiagnosisStatus,useNA='always')
diagnosed<-intermediate%>%filter(DiagnosisStatus=='D/D'|DiagnosisStatus=='D/L'|DiagnosisStatus=='P/D')
nrow(diagnosed) #n=455
table(diagnosed$censor) #n death=157
table(diagnosed$censor2) #n transplant=202-157=45
fit1_dx <- survfit(Surv(time, 1 - censor) ~ 1, data = diagnosed)
fit2_dx <- survfit(Surv(time2, 1 - censor2) ~ 1, data = diagnosed)

fit1_dx_df <- tidy(fit1_dx) |> mutate(group = "Death")
fit2_dx_df <- tidy(fit2_dx) |> mutate(group = "Heart transplant or Death")

fit1_dx_df <- bind_rows(tibble(time = 0, estimate = 1, conf.low = 1, conf.high = 1, group = "Death"), fit1_dx_df)
fit2_dx_df <- bind_rows(tibble(time = 0, estimate = 1, conf.low = 1, conf.high = 1, group = "Heart transplant or Death"), fit2_dx_df)

combined_df2 <- bind_rows(fit1_dx_df, fit2_dx_df)

p2 <- ggplot(combined_df2, aes(x = time, y = estimate, color = group, fill = group)) +
  geom_step(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  labs(x = "Age", y = "Survival Probability", title = "Diagnosed cohort") +
  coord_cartesian(xlim = c(0, 50), ylim = c(0, 1)) +
  theme_minimal(base_size = 12) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10)
  )+
  annotate("label",
           x = 20, y = 0.2,
           hjust = 0,
           vjust = 1,
           label = "N = 455\nN Deaths = 157\nN Heart transplant = 45",
           size = 2.5,
           label.size = 0.3,
           fill = "white",
           color = "black")
#The CI ribbon simply stops where CI is no longer mathematically defined.
(p1 | p2)


#-------------BRR KM-------------#
subset_registry$censor<-as.numeric(subset_registry$censor)
subset_registry2$censor2<-as.numeric(subset_registry2$censor2)
# Fit truncated KM
tmax <- 25

# 1) Fit on the FULL cohorts (no subsetting by time <= 25)

fit1 <- survfit(Surv(time, 1 - censor)  ~ 1, data = subset_registry)   # Death
fit2 <- survfit(Surv(time2, 1 - censor2) ~ 1, data = subset_registry2) # HTx or Death

# 2) Tidy
fit1_df <- tidy(fit1) %>% mutate(group = "Death")
fit2_df <- tidy(fit2) %>% mutate(group = "Heart transplant or Death")

# 3) Add time 0 row (for pretty start at 1.0)
add0 <- function(df, grp) bind_rows(tibble(time=0, estimate=1, conf.low=1, conf.high=1, group=grp), df)
fit1_df <- add0(fit1_df, "Death")
fit2_df <- add0(fit2_df, "Heart transplant or Death")

combined_df2 <- bind_rows(fit1_df, fit2_df)

# 4) Clip what's PLOTTED at tmax, but keep x-axis to 50
clip_km <- function(df, tmax = 25){
  out <- df %>% filter(time <= tmax)
  if (nrow(out) == 0) return(out)
  last <- out %>% slice_tail(n = 1)
  # Extend the last horizontal step to exactly tmax (optional but looks clean)
  if (last$time < tmax) {
    last$time <- tmax
    out <- bind_rows(out, last)
  }
  out
}

plot_df <- combined_df2 %>%
  group_by(group) %>%
  group_modify(~clip_km(.x, tmax)) %>%
  ungroup()

# 5) Plot: curves stop at 25; grid/x-axis still to 50
p2 <- ggplot(plot_df, aes(x = time, y = estimate, color = group, fill = group)) +
  geom_step(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, linewidth = 0) +
  labs(x = "Age", y = "Survival Probability", title = "BRR (longitudinal survival data)") +
  coord_cartesian(xlim = c(0, 50), ylim = c(0, 1)) +
  theme_minimal(base_size = 12) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10)
  ) +
  geom_vline(xintercept = seq(0, 50, by = 5), color = "gray80", linetype = "dotted") +
  annotate("label",
           x = 25, y = 0.2,
           hjust = 0,
           vjust = 1,
           label = "N = 131\nN Deaths = 7\nN Heart transplant = 5",
           size = 3,
           label.size = 0.3,
           fill = "white",
           color = "black")


(p1 | p2) + plot_layout(guides = "collect") + theme(legend.position = "bottom")


#-------------Clinical risk factors-------------#
risk_survival<-bf.last%>%dplyr::select(Intake_ID,time,entry,censor,developed,externalfeed)
#merging and create a combined dataset for survival analysis
table(risk_survival_familyR$family_relative,useNA = 'always')
risk_survival$Neutropenia<-ifelse(is.na(bf.last$Neutropenia),NA,ifelse(bf.last$Neutropenia=='never',0,1))
risk_survival$transplant<-bf.last$tx1
table(risk_survival_hf$heart_failure,useNA = 'always')
risk_survival$stroke<-bf.last$stroke
risk_survival<-merge(risk_survival,risk_survival_familyR,by.x='Intake_ID')
risk_survival<-merge(risk_survival,risk_survival_hf,by.x='Intake_ID')
risk_survival<-merge(risk_survival,risk_survival_feed,by.x='Intake_ID')
risk_survival<-merge(risk_survival,risk_survival_ca,by.x='Intake_ID')
risk_survival<-merge(risk_survival,risk_survival_aicdex,by.x='Intake_ID')
risk_survival<-merge(risk_survival,risk_survival_bc,by.x='Intake_ID')
table(risk_survival$censor)

#exploring missing
sum(complete.cases(risk_survival))#88 complete, other 41 have missing.
1-sum(complete.cases(risk_survival$Neutropenia))/125 #25.6% missing for neutropenia
1-sum(complete.cases(risk_survival$heart_failure))/125 #12.8% missing for heart failure
1-sum(complete.cases(risk_survival$feed))/125 #14.4% missing for feeding problem
1-sum(complete.cases(risk_survival$cardiac_arrest))/125 #18.4% missing for cardiac arrest
1-sum(complete.cases(risk_survival$blood_clot))/125 #16% missing for blood clot

#impute missing data
risk_survival_impute<-cbind(risk_survival[1:6],data.frame(lapply(risk_survival[7:15], as.factor)))
imputed_data <- mice(risk_survival_impute, m = 10, method = 'pmm', seed = 500)
#cox regression with secondary endpoint 
#we exclude stroke and external feeds, perfect seperation
coxfit<-with(imputed_data,coxph(Surv(time, 1-as.numeric(censor)) ~ transplant+developed+family_relative+Neutropenia+feed+heart_failure+cardiac_arrest+aicd_aex+blood_clot))
(pooled_coxfit<-summary(pool(coxfit)))
#Small differences may exist in pooled estimates when reproducing results—even with the same data and seed

#consider the longitudinal subset
#compute conditional Kendall’s tau values and the associated p-value to test independence
coxph(Surv(time,1-as.numeric(censor))~entry,data=risk_survival)
subset_risk<-subset(risk_survival, entry < time) # 109 potential left-truncated
length(unique(subset_risk$Intake_ID))
length(unique(subset_risk$Intake_ID[subset_risk$censor==1]))
length(unique(subset_risk$Intake_ID[subset_risk$censor==0]))
subset_risk$censor<-as.numeric(subset_risk$censor)
(cKendall(subset_risk$entry, subset_risk$time, 1-subset_risk$censor))
(cKendall(subset_risk$entry, subset_risk$time, 1-subset_risk$censor, method = "IPW1"))
subset_risk_impute<-cbind(subset_risk[1:6],data.frame(lapply(subset_risk[7:15], as.factor)))
imputed_data <- mice(subset_risk_impute, m = 10, method = 'pmm', seed = 500)
#cox regression with secondary endpoint in longitudinal cohort
#we exclude stroke, external feeds, transplant, family relative, cardiac device due to perfect seperation
coxfit_adjusted<-with(imputed_data,coxph(Surv(entry,time, 1-as.numeric(censor)) ~ Neutropenia+feed+heart_failure+cardiac_arrest+blood_clot))
(pooled_coxfit<-summary(pool(coxfit_adjusted)))
#Small differences may exist in pooled estimates when reproducing results—even with the same data and seed



#-----------------------------Death causes-----------------------------#
#Results can be found in Supp. Table 6

#Upload cause of death data
cause_death<-read.csv('/Users/mac1/Desktop/Barth_Syndrome_Data/cause_death.csv')
cause_death$age<-ifelse(cause_death$AgeDeath<1,1,
                        ifelse(cause_death$AgeDeath>=1&cause_death$AgeDeath<=5,2,
                               ifelse(cause_death$AgeDeath>5&cause_death$AgeDeath<=10,3,
                                      ifelse(cause_death$AgeDeath>10&cause_death$AgeDeath<=18,4,
                                             ifelse(cause_death$AgeDeath>18&cause_death$AgeDeath<=30,5,6)))))
table1(~CauseDeath|age,data=cause_death)


#-----------------------------Symptomatology-----------------------------#
#Results can be found in Fig. 3 and Supp. Table 5

age<-c('<1','1-5','5-10','11-18','19-30','>30')
distinct_age_counts<-numeric()
for (i in 1:6) {
  distinct_age_counts[i] <- current2 %>%filter(AGE==age[i])%>%
    group_by(Intake.ID) %>%
    summarise(distinct_age_count = n_distinct(Age.At.Survey))%>%
    summarise(count = sum(distinct_age_count))}
class(distinct_age_counts)

##bino_prev function for calculating prevalance and margin in Supp.Table 5
bino_prev <- function(var){
  count <- prev <- margin <- numeric(6)
  
  for (i in 1:6) {
    yes <- current2 %>%
      filter(AGE == age[i]) %>%
      group_by(Intake.ID, Age.At.Survey) %>%
      summarise(x_yes = any({{ var }} == "Yes", na.rm = TRUE), .groups = "drop") %>%
      summarise(count = sum(x_yes), .groups = "drop")
    
    count[i] <- yes$count
    prev[i]  <- round(count[i] / distinct_age_counts[[i]] * 100, 1)
    
    ci <- binom.confint(count[i], distinct_age_counts[[i]], methods = "wilson")
    margin[i] <- round((ci$upper - ci$lower) / 2 * 100, 1)
  }
  
  overall_prev <- round(sum(count) / 159 * 100,1)
  ci <- binom.confint(sum(count), 159, methods = "wilson")
  overall_margin <- round((ci$upper - ci$lower) / 2 * 100, 1)
  
  return(list(
    count = count,
    prev = prev,
    margin = margin,
    overall_prev = overall_prev,
    overall_margin = overall_margin
  ))
}

#-------------GI Manifestations-------------#
#If "always", "very frequently",'frequently','occasionlly','rarely' -> yes 
#If ' very rarely','never' -> no 
#If 'Unsure', NA -> NA
#for each patient, among observations of each distinct age, any "yes" -> 'Yes' , all NA -> na , otherwise -> 'No'

#Abdominal Pain
current2$Abdominal_pain12<-NA
current2$Abdominal_pain12[current$Abdominal_pain12=='Never'|
                            current$Abdominal_pain12=='Very Rarely' ]<-'No'
current2$Abdominal_pain12[current$Abdominal_pain12=='Always'|
                            current$Abdominal_pain12=='Very frequently'|
                            current$Abdominal_pain12=='Frequently'|
                            current$Abdominal_pain12=='Occasionally'|
                            current$Abdominal_pain12=='Rarely']<-'Yes'
table(current2$Abdominal_pain12,useNA = 'always')

Abdominal_pain_prev<-bino_prev(Abdominal_pain12)

#Constipation
current2$Constipation12<-NA
current2$Constipation12[current$Constipation12=='Never'|
                          current$Constipation12=='Very Rarely' ]<-'No'
current2$Constipation12[current$Constipation12=='Always'|
                          current$Constipation12=='Very frequently'|
                          current$Constipation12=='Frequently'|
                          current$Constipation12=='Occasionally'|
                          current$Constipation12=='Rarely']<-'Yes'
table(current2$Constipation12,useNA = 'always')
Constipation_prev<-bino_prev(Constipation12)

#Diarrhea
current2$Diarrhea12<-NA
current2$Diarrhea12[current$Diarrhea12=='Never'|
                      current$Diarrhea12=='Very Rarely' ]<-'No'
current2$Diarrhea12[current$Diarrhea12=='Always'|
                      current$Diarrhea12=='Very frequently'|
                      current$Diarrhea12=='Frequently'|
                      current$Diarrhea12=='Occasionally'|
                      current$Diarrhea12=='Rarely']<-'Yes'
table(current2$Diarrhea12,useNA = 'always')
Diarrhea_prev<-bino_prev(Diarrhea12)

#Nausea/vomiting
current2$Nausea_vomiting12<-NA
current2$Nausea_vomiting12[current$Nausea_vomiting12=='Never'|
                             current$Nausea_vomiting12=='Very Rarely' ]<-'No'
current2$Nausea_vomiting12[current$Nausea_vomiting12=='Always'|
                             current$Nausea_vomiting12=='Very frequently'|
                             current$Nausea_vomiting12=='Frequently'|
                             current$Nausea_vomiting12=='Occasionally'|
                             current$Nausea_vomiting12=='Rarely']<-'Yes'
table(current2$Nausea_vomiting12,useNA = 'always')
Nausea_vomiting_prev<-bino_prev(Nausea_vomiting12)

#Sores on anus
current2$Sores_anus12<-NA
current2$Sores_anus12[current$Sores_anus12=='Never'|
                        current$Sores_anus12=='Very Rarely' ]<-'No'
current2$Sores_anus12[current$Sores_anus12=='Always'|
                        current$Sores_anus12=='Very frequently'|
                        current$Sores_anus12=='Frequently'|
                        current$Sores_anus12=='Occasionally'|
                        current$Sores_anus12=='Rarely']<-'Yes'
table(current2$Sores_anus12,useNA = 'always')
Sores_anus_prev<-bino_prev(Sores_anus12)

#Any GI
count <- prev <- margin <- numeric(6)
for (i in 1:6) {
  yes<-current2%>%
    filter(AGE==age[i])%>%
    group_by(Intake.ID)%>%
    group_by(Age.At.Survey,.add=T)%>%
    summarise(gi12_yes = any(Sores_anus12 == 'Yes'|
                               Nausea_vomiting12 == 'Yes'|
                               Diarrhea12 == 'Yes'|
                               Constipation12 == 'Yes'|
                               Abdominal_pain12 == 'Yes')) %>%
    summarise(count = sum(gi12_yes,na.rm = T))
    count[i] <- sum(yes$count)
    prev[i]  <- round(count[i] / distinct_age_counts[[i]] * 100, 1)
    
    ci <- binom.confint(count[i], distinct_age_counts[[i]], methods = "wilson")
    margin[i] <- round((ci$upper - ci$lower) / 2 * 100, 1)
}
overall_prev <- round(sum(count) / 159 * 100,1)
ci <- binom.confint(sum(count), 159, methods = "wilson")
overall_margin <- round((ci$upper - ci$lower) / 2 * 100, 1)
anyGI_prev<-list(count = count,
                     prev = prev,
                     margin = margin,
                     overall_prev = overall_prev,
                     overall_margin = overall_margin)

    
GI_prev<-list(anyGI_prev,Abdominal_pain_prev,Constipation_prev,Diarrhea_prev,Nausea_vomiting_prev,
              Sores_anus_prev)


#----Table----#
stay_levels <- c("Any GI", "Abdominal pain", "Constipation", "Diarrhea","Nausea vomiting","Sores anus")
age_labels <- c(
  "Age <1",
  "Age 1–5",
  "Age 6–10",
  "Age 11–18",
  "Age 19–30",
  "Age >30"
)
denominators <- c(14, 38, 27, 34, 28, 18)

GI_df <- map2_df(
  GI_prev,
  stay_levels,
  ~ tibble(
    GI = .y,
    AgeGroup     = age_labels,
    Count        = .x$count,
    Denominator  = denominators,
    Percentage   = .x$prev,
    Margin       = .x$margin,
    Category     = "GI"
  )
)

GI_table <- GI_df %>%
  mutate(Cell = sprintf("%.1f ± %.1f", Percentage, Margin)) %>%
  dplyr::select(GI, AgeGroup, Cell) %>%
  pivot_wider(names_from = AgeGroup, values_from = Cell)

GI_print_table<-knitr::kable(GI_table)

#----Heatmap----#
GI_long <- GI_df %>%
  mutate(
    AgeGroup = factor(
      AgeGroup,
      levels = c("Age <1", "Age 1–5", "Age 6–10", "Age 11–18", "Age 19–30", "Age >30"),
      labels = c("<1", "1-5", "6-10", "11-18", "19-30", ">30")
    ),
    GI = factor(GI, levels = c("Sores anus","Nausea vomiting","Diarrhea","Constipation","Abdominal pain","Any GI")
  ))

# Plot heatmap
p1<-ggplot(GI_long, aes(x = AgeGroup, y = GI, fill = Percentage)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f", Percentage)), size = 2.5) +
  scale_fill_gradient(low = "white", high = "red") +
  labs(
    title = "GI Manifestation",
    x = "Age Group",
    y = "GI",
    fill = "Prevalence (%)"
  ) +
  theme_minimal() +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 8),    
        plot.title = element_text(face = "bold", hjust = 0.5,size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0)
  )


#-------------Infection & Immune-------------#
#Abscess on anus
current2$Abscess_anus12<-NA
current2$Abscess_anus12[current$Abscess_anus12=='Never'|
                              current$Abscess_anus12=='Very Rarely' ]<-'No'
current2$Abscess_anus12[current$Abscess_anus12=='Always'|
                              current$Abscess_anus12=='Very frequently'|
                              current$Abscess_anus12=='Frequently'|
                              current$Abscess_anus12=='Occasionally'|
                              current$Abscess_anus12=='Rarely']<-'Yes'
table(current2$Abscess_anus12,useNA = 'always')
Abscess_anus_prev<-bino_prev(Abscess_anus12)

#Bronchitis
current2$Bronchitis12<-NA
current2$Bronchitis12[current$Bronchitis12=='Never'|
                        current$Bronchitis12=='Very Rarely' ]<-'No'
current2$Bronchitis12[current$Bronchitis12=='Always'|
                        current$Bronchitis12=='Very frequently'|
                        current$Bronchitis12=='Frequently'|
                        current$Bronchitis12=='Occasionally'|
                        current$Bronchitis12=='Rarely']<-'Yes'
table(current2$Bronchitis12,useNA = 'always')
Bronchitis_prev<-bino_prev(Bronchitis12)

#fever
current2$Fever12<-NA
current2$Fever12[current$Fever12=='Never'|
                   current$Fever12=='Very Rarely' ]<-'No'
current2$Fever12[current$Fever12=='Always'|
                   current$Fever12=='Very frequently'|
                   current$Fever12=='Frequently'|
                   current$Fever12=='Occasionally'|
                   current$Fever12=='Rarely']<-'Yes'
table(current2$Fever12,useNA = 'always')
fever_prev<-bino_prev(Fever12)

#ear infection
current2$Ear_infections12<-NA
current2$Ear_infections12[current$Ear_infections12=='Never'|
                            current$Ear_infections12=='Very Rarely' ]<-'No'
current2$Ear_infections12[current$Ear_infections12=='Always'|
                            current$Ear_infections12=='Very frequently'|
                            current$Ear_infections12=='Frequently'|
                            current$Ear_infections12=='Occasionally'|
                            current$Ear_infections12=='Rarely']<-'Yes'
table(current2$Ear_infections12,useNA = 'always')
ear_infections_prev<-bino_prev(Ear_infections12)

#Gingivitis
current2$Gingivitis_12<-NA
current2$Gingivitis_12[current$Gingivitis_12=='Never'|
                         current$Gingivitis_12=='Very Rarely' ]<-'No'
current2$Gingivitis_12[current$Gingivitis_12=='Always'|
                         current$Gingivitis_12=='Very frequently'|
                         current$Gingivitis_12=='Frequently'|
                         current$Gingivitis_12=='Occasionally'|
                         current$Gingivitis_12=='Rarely']<-'Yes'
table(current2$Gingivitis_12,useNA = 'always')
Gingivitis_prev<-bino_prev(Gingivitis_12)

#Mouth throat ulcers
current2$Mouth_throat_ulcers12<-NA
current2$Mouth_throat_ulcers12[current$Mouth_throat_ulcers12=='Never'|
                                 current$Mouth_throat_ulcers12=='Very Rarely' ]<-'No'
current2$Mouth_throat_ulcers12[current$Mouth_throat_ulcers12=='Always'|
                                 current$Mouth_throat_ulcers12=='Very frequently'|
                                 current$Mouth_throat_ulcers12=='Frequently'|
                                 current$Mouth_throat_ulcers12=='Occasionally'|
                                 current$Mouth_throat_ulcers12=='Rarely']<-'Yes'
table(current2$Mouth_throat_ulcers12,useNA = 'always')
Mouth_throat_ulcers_prev<-bino_prev(Mouth_throat_ulcers12)

#Pneumonia
current2$Pneumonia12<-NA
current2$Pneumonia12[current$Pneumonia12=='Never'|
                       current$Pneumonia12=='Very Rarely' ]<-'No'
current2$Pneumonia12[current$Pneumonia12=='Always'|
                       current$Pneumonia12=='Very frequently'|
                       current$Pneumonia12=='Frequently'|
                       current$Pneumonia12=='Occasionally'|
                       current$Pneumonia12=='Rarely']<-'Yes'
table(current2$Pneumonia12,useNA = 'always')
Pneumonia_prev<-bino_prev(Pneumonia12)

#Sinus infection
current2$Sinus12<-NA
current2$Sinus12[current$Sinus12=='Never'|
                   current$Sinus12=='Very Rarely' ]<-'No'
current2$Sinus12[current$Sinus12=='Always'|
                   current$Sinus12=='Very frequently'|
                   current$Sinus12=='Frequently'|
                   current$Sinus12=='Occasionally'|
                   current$Sinus12=='Rarely']<-'Yes'
table(current2$Sinus12,useNA = 'always')
Sinus_prev<-bino_prev(Sinus12)

#Skin infection
current2$Skin_infections12<-NA
current2$Skin_infections12[current$Skin_infections12=='Never'|
                             current$Skin_infections12=='Very Rarely' ]<-'No'
current2$Skin_infections12[current$Skin_infections12=='Always'|
                             current$Skin_infections12=='Very frequently'|
                             current$Skin_infections12=='Frequently'|
                             current$Skin_infections12=='Occasionally'|
                             current$Skin_infections12=='Rarely']<-'Yes'
table(current2$Skin_infections12,useNA = 'always')
Skin_infections_prev<-bino_prev(Skin_infections12)

#Sore glands
current2$Sore_glands12<-NA
current2$Sore_glands12[current$Sore_glands12=='Never'|
                         current$Sore_glands12=='Very Rarely' ]<-'No'
current2$Sore_glands12[current$Sore_glands12=='Always'|
                         current$Sore_glands12=='Very frequently'|
                         current$Sore_glands12=='Frequently'|
                         current$Sore_glands12=='Occasionally'|
                         current$Sore_glands12=='Rarely']<-'Yes'
table(current2$Sore_glands12,useNA = 'always')
Sore_glands_prev<-bino_prev(Sore_glands12)


#Throat or tonsil infection
current2$Throat_tonsil_infections12<-NA
current2$Throat_tonsil_infections12[current$Throat_tonsil_infections12=='Never'|
                                      current$Throat_tonsil_infections12=='Very Rarely' ]<-'No'
current2$Throat_tonsil_infections12[current$Throat_tonsil_infections12=='Always'|
                                      current$Throat_tonsil_infections12=='Very frequently'|
                                      current$Throat_tonsil_infections12=='Frequently'|
                                      current$Throat_tonsil_infections12=='Occasionally'|
                                      current$Throat_tonsil_infections12=='Rarely']<-'Yes'
table(current2$Throat_tonsil_infections12,useNA = 'always')
Throat_tonsil_infections_prev<-bino_prev(Throat_tonsil_infections12)

#Tooth abscess
current2$Tooth_abscess12<-NA
current2$Tooth_abscess12[current$Tooth_abscess12=='Never'|
                           current$Tooth_abscess12=='Very Rarely' ]<-'No'
current2$Tooth_abscess12[current$Tooth_abscess12=='Always'|
                           current$Tooth_abscess12=='Very frequently'|
                           current$Tooth_abscess12=='Frequently'|
                           current$Tooth_abscess12=='Occasionally'|
                           current$Tooth_abscess12=='Rarely']<-'Yes'
table(current2$Tooth_abscess12,useNA = 'always')
Tooth_abscess_prev<-bino_prev(Tooth_abscess12)

#Urinary tract infection
current2$Urinary_tract_infection12<-NA
current2$Urinary_tract_infection12[current$Urinary_tract_infection12=='Never'|
                                     current$Urinary_tract_infection12=='Very Rarely' ]<-'No'
current2$Urinary_tract_infection12[current$Urinary_tract_infection12=='Always'|
                                     current$Urinary_tract_infection12=='Very frequently'|
                                     current$Urinary_tract_infection12=='Frequently'|
                                     current$Urinary_tract_infection12=='Occasionally'|
                                     current$Urinary_tract_infection12=='Rarely']<-'Yes'
table(current2$Urinary_tract_infection12,useNA = 'always')
Urinary_tract_infection_prev<-bino_prev(Urinary_tract_infection12)

#Any infection
count <- prev <- margin <- numeric(6)
for (i in 1:6) {
  yes<-current2%>%
    filter(AGE==age[i])%>%
    group_by(Intake.ID)%>%
    group_by(Age.At.Survey,.add=T)%>%
    summarise(infection12_yes = any(Urinary_tract_infection12 == 'Yes'|
                                      Tooth_abscess12 == 'Yes'|
                                      Throat_tonsil_infections12 == 'Yes'|
                                      Sore_glands12 == 'Yes'|
                                      Skin_infections12 == 'Yes'|
                                      Sinus12 == 'Yes'|
                                      Pneumonia12 == 'Yes'|
                                      Mouth_throat_ulcers12 == 'Yes'|
                                      Gingivitis_12 == 'Yes'|
                                      Ear_infections12 == 'Yes'|
                                      Fever12 == 'Yes'|
                                      Bronchitis12 == 'Yes'|
                                      Abscess_anus12 == 'Yes')) %>%
    summarise(count = sum(infection12_yes,na.rm = T))
  count[i] <- sum(yes$count)
  prev[i]  <- round(count[i] / distinct_age_counts[[i]] * 100, 1)
  
  ci <- binom.confint(count[i], distinct_age_counts[[i]], methods = "wilson")
  margin[i] <- round((ci$upper - ci$lower) / 2 * 100, 1)
}
overall_prev <- round(sum(count) / 159 * 100,1)
ci <- binom.confint(sum(count), 159, methods = "wilson")
overall_margin <- round((ci$upper - ci$lower) / 2 * 100, 1)
anyinfection_prev<-list(count = count,
                 prev = prev,
                 margin = margin,
                 overall_prev = overall_prev,
                 overall_margin = overall_margin)

infection_prev<-list(anyinfection_prev,Abscess_anus_prev,Bronchitis_prev,
                     ear_infections_prev,fever_prev,Gingivitis_prev,Mouth_throat_ulcers_prev,
                     Pneumonia_prev,Sinus_prev,Skin_infections_prev,Sore_glands_prev,
                     Throat_tonsil_infections_prev,Tooth_abscess_prev,Urinary_tract_infection_prev)


#-----Table----#
stay_levels <- c("Any", "Abscess on anus", "Bronchitis", "Ear infection","Fever"," Gingivitis",
                 "Mouth/throat ulcers","Pneumonia","Sinus infection","Skin infection",
                 "Sore glands","Throat or tonsil infection","Tooth abscess","Urinary tract infection")


infection_df <- map2_df(
  infection_prev,
  stay_levels,
  ~ tibble(
    infection = .y,
    AgeGroup     = age_labels,
    Count        = .x$count,
    Denominator  = denominators,
    Percentage   = .x$prev,
    Margin       = .x$margin,
    Category     = "infection"
  )
)

infection_table <- infection_df %>%
  mutate(Cell = sprintf("%.1f ± %.1f", Percentage, Margin)) %>%
  dplyr::select(infection, AgeGroup, Cell) %>%
  pivot_wider(names_from = AgeGroup, values_from = Cell)

infection_print_table<-knitr::kable(infection_table)

#-----Heatmap-----#
infection_long <- infection_df %>%
  mutate(
    AgeGroup = factor(
      AgeGroup,
      levels = c("Age <1", "Age 1–5", "Age 6–10", "Age 11–18", "Age 19–30", "Age >30"),
      labels = c("<1", "1-5", "6-10", "11-18", "19-30", ">30")
    ),
    infection = factor(infection, levels = rev(stay_levels)
    ))

# Plot heatmap
p2<-ggplot(infection_long, aes(x = AgeGroup, y = infection, fill = Percentage)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f", Percentage)), size = 2.5) +
  scale_fill_gradient(low = "white", high = "red") +
  labs(
    title = "Infection & Immune Manifestation",
    x = "Age Group",
    y = "infection",
    fill = "Prevalence (%)"
  ) +
  theme_minimal() +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 8),    
        plot.title = element_text(face = "bold", hjust = 0.5,size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0)
  )


(p1 | p2)
ggsave("combined_plot12.png", width = 8, height = 3, dpi = 400)

#-------------Musculoskeletal-------------#
#Bone / joint pain    
current2$Bone_jointpain12<-NA
current2$Bone_jointpain12[current$Bone_jointpain12=='Never'|
                            current$Bone_jointpain12=='Very Rarely' ]<-'No'
current2$Bone_jointpain12[current$Bone_jointpain12=='Always'|
                            current$Bone_jointpain12=='Very frequently'|
                            current$Bone_jointpain12=='Frequently'|
                            current$Bone_jointpain12=='Occasionally'|
                            current$Bone_jointpain12=='Rarely']<-'Yes'
table(current2$Bone_jointpain12,useNA = 'always')
Bone_jointpain_prev<-bino_prev(Bone_jointpain12)
                     
#Excessive fatigue / tired
current2$Excessive.fatigue_tired12<-NA
current2$Excessive.fatigue_tired12[current$Excessive.fatigue_tired12=='Never'|
                                     current$Excessive.fatigue_tired12=='Very Rarely' ]<-'No'
current2$Excessive.fatigue_tired12[current$Excessive.fatigue_tired12=='Always'|
                                     current$Excessive.fatigue_tired12=='Very frequently'|
                                     current$Excessive.fatigue_tired12=='Frequently'|
                                     current$Excessive.fatigue_tired12=='Occasionally'|
                                     current$Excessive.fatigue_tired12=='Rarely']<-'Yes'
table(current2$Excessive.fatigue_tired12,useNA = 'always')
Excessive.fatigue_tired_prev<-bino_prev(Excessive.fatigue_tired12)

#Any Musculoskeletal
for (i in 1:6) {
  yes<-current2%>%
    filter(AGE==age[i])%>%
    group_by(Intake.ID)%>%
    group_by(Age.At.Survey,.add=T)%>%
    summarise(Musculoskeletal12_yes = any(Excessive.fatigue_tired12 == 'Yes'|
                                            Bone_jointpain12 == 'Yes')) %>%
    summarise(count = sum(Musculoskeletal12_yes,na.rm = T))
  count[i] <- sum(yes$count)
  prev[i]  <- round(count[i] / distinct_age_counts[[i]] * 100, 1)
  
  ci <- binom.confint(count[i], distinct_age_counts[[i]], methods = "wilson")
  margin[i] <- round((ci$upper - ci$lower) / 2 * 100, 1)
}
overall_prev <- round(sum(count) / 159 * 100,1)
ci <- binom.confint(sum(count), 159, methods = "wilson")
overall_margin <- round((ci$upper - ci$lower) / 2 * 100, 1)
anymuscul_prev<-list(count = count,
                        prev = prev,
                        margin = margin,
                        overall_prev = overall_prev,
                        overall_margin = overall_margin)

muscul_prev<-list(anymuscul_prev,Bone_jointpain_prev,Excessive.fatigue_tired_prev)                 


#-----Table-----#
stay_levels <- c("Any", "Bone / joint pain", "Fatigue / tired")

muscul_df <- map2_df(
  muscul_prev,
  stay_levels,
  ~ tibble(
    muscul = .y,
    AgeGroup     = age_labels,
    Count        = .x$count,
    Denominator  = denominators,
    Percentage   = .x$prev,
    Margin       = .x$margin,
    Category     = "muscul"
  )
)

muscul_table <- muscul_df %>%
  mutate(Cell = sprintf("%.1f ± %.1f", Percentage, Margin)) %>%
  dplyr::select(muscul, AgeGroup, Cell) %>%
  pivot_wider(names_from = AgeGroup, values_from = Cell)

muscul_print_table<-knitr::kable(muscul_table)

#-----Heatmap-----#
muscul_long <- muscul_df %>%
  mutate(
    AgeGroup = factor(
      AgeGroup,
      levels = c("Age <1", "Age 1–5", "Age 6–10", "Age 11–18", "Age 19–30", "Age >30"),
      labels = c("<1", "1-5", "6-10", "11-18", "19-30", ">30")
    ),
    muscul = factor(muscul, levels = rev(stay_levels)
    ))

# Plot heatmap
p3<-ggplot(muscul_long, aes(x = AgeGroup, y = muscul, fill = Percentage)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f", Percentage)), size = 2.5) +
  scale_fill_gradient(low = "white", high = "red") +
  labs(
    title = "Musculoskeletal Manifestation",
    x = "Age Group",
    y = "infection",
    fill = "Prevalence (%)"
  ) +
  theme_minimal() +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 8),    
        plot.title = element_text(face = "bold", hjust = 0.5,size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0)
  )                


#-------------Neurological-------------#
#Depression
current2$Depression12<-NA
current2$Depression12[current$Depression12=='Never'|
                        current$Depression12=='Very Rarely' ]<-'No'
current2$Depression12[current$Depression12=='Always'|
                        current$Depression12=='Very frequently'|
                        current$Depression12=='Frequently'|
                        current$Depression12=='Occasionally'|
                        current$Depression12=='Rarely']<-'Yes'
table(current2$Depression12,useNA = 'always')
Depression_prev<-bino_prev(Depression12)

#Headaches
current2$Headaches12<-NA
current2$Headaches12[current$Headaches12=='Never'|
                       current$Headaches12=='Very Rarely' ]<-'No'
current2$Headaches12[current$Headaches12=='Always'|
                       current$Headaches12=='Very frequently'|
                       current$Headaches12=='Frequently'|
                       current$Headaches12=='Occasionally'|
                       current$Headaches12=='Rarely']<-'Yes'
table(current2$Headaches12,useNA = 'always')
Headaches_prev<-bino_prev(Headaches12)

#Any neurological
for (i in 1:6) {
  yes<-current2%>%
    filter(AGE==age[i])%>%
    group_by(Intake.ID)%>%
    group_by(Age.At.Survey,.add=T)%>%
    summarise(Neurological12_yes = any(Headaches12 == 'Yes'|
                                         Depression12 == 'Yes')) %>%
    summarise(count = sum(Neurological12_yes,na.rm = T))
  count[i] <- sum(yes$count)
  prev[i]  <- round(count[i] / distinct_age_counts[[i]] * 100, 1)
  
  ci <- binom.confint(count[i], distinct_age_counts[[i]], methods = "wilson")
  margin[i] <- round((ci$upper - ci$lower) / 2 * 100, 1)
}
overall_prev <- round(sum(count) / 159 * 100,1)
ci <- binom.confint(sum(count), 159, methods = "wilson")
overall_margin <- round((ci$upper - ci$lower) / 2 * 100, 1)
anyneuro_prev<-list(count = count,
                     prev = prev,
                     margin = margin,
                     overall_prev = overall_prev,
                     overall_margin = overall_margin)

neuro_prev<-list(anyneuro_prev,Depression_prev,Headaches_prev)         


#-----Table-----#
stay_levels <- c("Any", "Depression", "Headaches")

neuro_df <- map2_df(
  neuro_prev,
  stay_levels,
  ~ tibble(
    neuro = .y,
    AgeGroup     = age_labels,
    Count        = .x$count,
    Denominator  = denominators,
    Percentage   = .x$prev,
    Margin       = .x$margin,
    Category     = "neuro"
  )
)

neuro_table <- neuro_df %>%
  mutate(Cell = sprintf("%.1f ± %.1f", Percentage, Margin)) %>%
  dplyr::select(neuro, AgeGroup, Cell) %>%
  pivot_wider(names_from = AgeGroup, values_from = Cell)

neuro_print_table<-knitr::kable(neuro_table)

#-----Heatmap-----#
neuro_long <- neuro_df %>%
  mutate(
    AgeGroup = factor(
      AgeGroup,
      levels = c("Age <1", "Age 1–5", "Age 6–10", "Age 11–18", "Age 19–30", "Age >30"),
      labels = c("<1", "1-5", "6-10", "11-18", "19-30", ">30")
    ),
    neuro = factor(neuro, levels = rev(stay_levels)
    ))

# Plot heatmap
p4<-ggplot(neuro_long, aes(x = AgeGroup, y = neuro, fill = Percentage)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f", Percentage)), size = 2.5) +
  scale_fill_gradient(low = "white", high = "red") +
  labs(
    title = "Neurological Manifestation",
    x = "Age Group",
    y = "Neurological",
    fill = "Prevalence (%)"
  ) +
  theme_minimal() +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 8),    
        plot.title = element_text(face = "bold", hjust = 0.5,size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0)
  )  

(p3 | p4)
ggsave("combined_plot34.png", width = 8, height = 3, dpi = 400)

#-------------Neutropenia-------------#
table(current2$neutropenia_current)
#use the last observation of each distinct age of each subject
current2$Surveytime <- strptime(current2$Survey.Time, format = "%Y/%m/%d %H:%M", tz = "EST")
current2<-current2%>%dplyr::group_by(Intake.ID)%>%group_by(Age.At.Survey,.add=T)%>%mutate(last=ifelse(Surveytime==max(Surveytime),1,0))
current2%>%dplyr::select(Intake.ID,Survey.Time,last)
current2$last[is.na(current2$last)]<-1
count <- prev <- margin <- numeric(6)
bino_netro_prev<-function(x){
  for (i in 1:6) {
    yes<-current2%>%
      filter(AGE==age[i])%>%
      filter(last==1)%>%#use the last observation of each distinct age of each subject
      filter(neutropenia_current==x)%>%
      group_by(Intake.ID)%>%
      summarise(distinct_age_count = n_distinct(Age.At.Survey))%>%
      summarise(count = sum(distinct_age_count))
      count[i] <- yes$count
      prev[i]  <- round(count[i] / distinct_age_counts[[i]] * 100, 1)
    
      ci <- binom.confint(count[i], distinct_age_counts[[i]], methods = "wilson")
      margin[i] <- round((ci$upper - ci$lower) / 2 * 100, 1)
  }
  
    overall_prev <- round(sum(count) / 159 * 100,1)
    ci <- binom.confint(sum(count), 159, methods = "wilson")
    overall_margin <- round((ci$upper - ci$lower) / 2 * 100, 1)
  
  return(list(
    count = count,
    prev = prev,
    margin = margin,
    overall_prev = overall_prev,
    overall_margin = overall_margin
  ))
  }

#Chronic
chronic_netro_prev<-bino_netro_prev('Chronic (always neutropenic/ constantly less than 1500/mm3)')

#Cyclic
cyclic_netro_prev<-bino_netro_prev('Cyclic (regular and predictable cycle with neutrophil count less than 1500/mm3)')

#Intermittent
interm_netro_prev<-bino_netro_prev('Intermittent (unpredictable cycle/ with a neutrophil count sometimes less than 1500/mm3 or 1.5 x10*9/L)')

#Neutropenia a problem in the past but no longer a problem
past_netro_prev<-bino_netro_prev('Neutropenia a problem in the past but no longer a problem')

neutro_prev<-list(chronic_netro_prev,cyclic_netro_prev,interm_netro_prev,past_netro_prev)


#-----Table-----#
stay_levels <- c("Chronic", " Cyclic","Intermittent","Past but no now")

neutro_df <- map2_df(
  neutro_prev,
  stay_levels,
  ~ tibble(
    neutro = .y,
    AgeGroup     = age_labels,
    Count        = .x$count,
    Denominator  = denominators,
    Percentage   = .x$prev,
    Margin       = .x$margin,
    Category     = "neutro"
  )
)

neutro_table <- neutro_df %>%
  mutate(Cell = sprintf("%.1f ± %.1f", Percentage, Margin)) %>%
  dplyr::select(neutro, AgeGroup, Cell) %>%
  pivot_wider(names_from = AgeGroup, values_from = Cell)

neutro_print_table<-knitr::kable(neutro_table)

#-----Heatmap-----#
neutro_long <- neutro_df %>%
  mutate(
    AgeGroup = factor(
      AgeGroup,
      levels = c("Age <1", "Age 1–5", "Age 6–10", "Age 11–18", "Age 19–30", "Age >30"),
      labels = c("<1", "1-5", "6-10", "11-18", "19-30", ">30")
    ),
    neutro = factor(neutro, levels = stay_levels
    ))

# Plot heatmap
p5<-ggplot(neutro_long, aes(x = AgeGroup, y = neutro, fill = Percentage)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f", Percentage)), size = 2.5) +
  scale_fill_gradient(low = "white", high = "red") +
  labs(
    title = "Neutropenia",
    x = "Age Group",
    y = "Neutropenia",
    fill = "Prevalence (%)"
  ) +
  theme_minimal() +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 8),    
        plot.title = element_text(face = "bold", hjust = 0.5,size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0)
  ) 


#-------------Cardiovascular symptoms-------------#

count <- prev <- margin <- numeric(6)
bino_cardio_prev<-function(var,x){
  for (i in 1:6) {
    yes<-current2%>%
      filter(AGE==age[i])%>%
      filter(last==1)%>%#use the last observation of each distinct age of each subject
      filter({{ var }}==x)%>%
      group_by(Intake.ID)%>%
      summarise(distinct_age_count = n_distinct(Age.At.Survey))%>%
      summarise(count = sum(distinct_age_count))
    count[i] <- yes$count
    prev[i]  <- round(count[i] / distinct_age_counts[[i]] * 100, 1)
    
    ci <- binom.confint(count[i], distinct_age_counts[[i]], methods = "wilson")
    margin[i] <- round((ci$upper - ci$lower) / 2 * 100, 1)
  }
  
  overall_prev <- round(sum(count) / 159 * 100,1)
  ci <- binom.confint(sum(count), 159, methods = "wilson")
  overall_margin <- round((ci$upper - ci$lower) / 2 * 100, 1)
  
  return(list(
    count = count,
    prev = prev,
    margin = margin,
    overall_prev = overall_prev,
    overall_margin = overall_margin
  ))
}
#Ejection fraction
table(current2$Ejection_Fraction_current)
current2$Ejection_Fraction_current<-NA
current2$Ejection_Fraction_current[current$Ejection_Fraction_current=='71% - 75%'|
                                     current$Ejection_Fraction_current=='55% - 70%' ]<-'Normal'
current2$Ejection_Fraction_current[current$Ejection_Fraction_current=='40% - 54%'|
                                     current$Ejection_Fraction_current=='35% - 39%'|
                                     current$Ejection_Fraction_current=='30% - 34%'|
                                     current$Ejection_Fraction_current=='24% - 29%'|
                                     current$Ejection_Fraction_current=='15% - 23%'|
                                     current$Ejection_Fraction_current=='10% - 14%'  ]<-'Abnormal'
table(current2$Ejection_Fraction_current,useNA = 'always')
normal_cardio_prev<-bino_cardio_prev(Ejection_Fraction_current,'Normal')
ab_cardio_prev<-bino_cardio_prev(Ejection_Fraction_current,'Abnormal')

#Shortening fraction
table(current2$Shortening_Fraction_current)
current2$Shortening_Fraction_current<-NA
current2$Shortening_Fraction_current[current$Shortening_Fraction_current=='71% - 75%'|
                                       current$Shortening_Fraction_current=='55% - 70%' ]<-'Normal'
current2$Shortening_Fraction_current[current$Shortening_Fraction_current=='40% - 54%'|
                                       current$Shortening_Fraction_current=='35% - 39%'|
                                       current$Shortening_Fraction_current=='30% - 34%'|
                                       current$Shortening_Fraction_current=='24% - 29%'|
                                       current$Shortening_Fraction_current=='15% - 23%'|
                                       current$Shortening_Fraction_current=='10% - 14%'  |
                                       current$Shortening_Fraction_current=='Less than 10%']<-'Abnormal'
table(current2$Shortening_Fraction_current,useNA = 'always')
Shorten_cardio_prev<-bino_cardio_prev(Shortening_Fraction_current,'Abnormal')

#Racing heart/palpitation
table(current2$racing_heart12)
current2$racing_heart12[current$racing_heart12=='Very frequently'|
                          current$racing_heart12=='Frequently' |
                          current$racing_heart12=='Always'|
                          current$racing_heart12=='Occasionally'|
                          current$racing_heart12=='Rarely']<-'Yes'

current2$racing_heart12[current$racing_heart12=='Very Rarely'|
                          current$racing_heart12=='Never']<-'No'

current2$racing_heart12[current$racing_heart12=='Unsure'|
                          current$racing_heart12=='999' ]<-NA

table(current2$racing_heart12,useNA='always')
palp_cardio_prev<-bino_cardio_prev(racing_heart12,'Yes')

cardio_prev<-list(ab_cardio_prev,palp_cardio_prev)


#-----Table-----#
stay_levels <- c("Abnormal Ejection (<55%)","Palpitation")

cardio_df <- map2_df(
  cardio_prev,
  stay_levels,
  ~ tibble(
    cardio = .y,
    AgeGroup     = age_labels,
    Count        = .x$count,
    Denominator  = denominators,
    Percentage   = .x$prev,
    Margin       = .x$margin,
    Category     = "cardio"
  )
)

cardio_table <- cardio_df %>%
  mutate(Cell = sprintf("%.1f ± %.1f", Percentage, Margin)) %>%
  dplyr::select(cardio, AgeGroup, Cell) %>%
  pivot_wider(names_from = AgeGroup, values_from = Cell)

cardio_print_table<-knitr::kable(cardio_table)

#-----Heatmap-----#
cardio_long <- cardio_df %>%
  mutate(
    AgeGroup = factor(
      AgeGroup,
      levels = c("Age <1", "Age 1–5", "Age 6–10", "Age 11–18", "Age 19–30", "Age >30"),
      labels = c("<1", "1-5", "6-10", "11-18", "19-30", ">30")
    ),
    neutro = factor(cardio, levels = rev(stay_levels)
    ))

# Plot heatmap
p6<-ggplot(cardio_long, aes(x = AgeGroup, y = neutro, fill = Percentage)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f", Percentage)), size = 2.5) +
  scale_fill_gradient(low = "white", high = "red") +
  labs(
    title = "Cardiovascular Manifestations",
    x = "Age Group",
    y = "Cardiovascular",
    fill = "Prevalence (%)"
  ) +
  theme_minimal() +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 8),    
        plot.title = element_text(face = "bold", hjust = 0.5,size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0)
  ) 


(p5 | p6)
ggsave("combined_plot56.png", width = 8, height = 3, dpi = 400)

###Nights in the hospital last year###
table(current2$hospital12,useNA = 'always')
current2$hospital12[current2$hospital12=='999'|current2$hospital12=='Unsure']<-NA
current2$hospital12[current2$hospital12=='100 or more']<-101
current2$hospital12<-as.numeric(current2$hospital12)
current2$hospital12[current2$hospital12>0&current2$hospital12<=30]<-1
current2$hospital12[current2$hospital12>30&current2$hospital12<=100]<-2
current2$hospital12[current2$hospital12==101]<-3
table(current2$hospital12,useNA = 'always')

bino_hostay_prev<-function(x){
    for (i in 1:6) {
      yes<-current2%>%
        filter(AGE==age[i])%>%
        filter(last==1)%>%#use the last observation of each distinct age of each subject
        filter(hospital12==x)%>%
        group_by(Intake.ID)%>%
        summarise(distinct_age_count = n_distinct(Age.At.Survey))%>%
        summarise(count = sum(distinct_age_count))
      count[i] <- yes$count
      prev[i]  <- round(count[i] / distinct_age_counts[[i]] * 100, 1)
      
      ci <- binom.confint(count[i], distinct_age_counts[[i]], methods = "wilson")
      margin[i] <- round((ci$upper - ci$lower) / 2 * 100, 1)
    }
    
    overall_prev <- round(sum(count) / 159 * 100,1)
    ci <- binom.confint(sum(count), 159, methods = "wilson")
    overall_margin <- round((ci$upper - ci$lower) / 2 * 100, 1)
    
    return(list(
      count = count,
      prev = prev,
      margin = margin,
      overall_prev = overall_prev,
      overall_margin = overall_margin
    ))
  }
hostay0_prev<-bino_hostay_prev(0)
hostay1_prev<-bino_hostay_prev(1)
hostay2_prev<-bino_hostay_prev(2)
hostay3_prev<-bino_hostay_prev(3)

hostay_prev<-list(hostay0_prev,hostay1_prev,hostay2_prev,hostay3_prev)
hostay_prev
#-----Table-----#
stay_levels <- c("0", "1–30", "30–100", ">100")
age_labels <- c(
  "Age <1",
  "Age 1–5",
  "Age 6–10",
  "Age 11–18",
  "Age 19–30",
  "Age >30"
)
denominators <- c(14, 38, 27, 34, 28, 18)

hospital_df <- map2_df(
  hostay_prev,
  stay_levels,
  ~ tibble(
    HospitalStay = .y,
    AgeGroup     = age_labels,
    Count        = .x$count,
    Denominator  = denominators,
    Percentage   = .x$prev,
    Margin       = .x$margin,
    Category     = "Hospital Stay"
  )
)

hospital_table <- hospital_df %>%
  mutate(Cell = sprintf("%.1f ± %.1f", Percentage, Margin)) %>%
  dplyr::select(HospitalStay, AgeGroup, Cell) %>%
  pivot_wider(names_from = AgeGroup, values_from = Cell)

hospital_print_table<-knitr::kable(hospital_table)


#-----Heatmap-----#
hospital_long <- hospital_long %>%
  mutate(
    AgeGroup = factor(AgeGroup, levels = c("Age <1", "Age 1-5", "Age 6-10", "Age 11-18", "Age 19-30", "Age >30"),
                      labels = c("<1", "1-5", "6-10", "11-18", "19-30", ">30")),
    HospitalStay = factor(HospitalStay, levels = c("0", "1–30", "30–100", ">100"))
  )


# Plot heatmap
ggplot(hospital_long, aes(x = AgeGroup, y = HospitalStay, fill = Percentage)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f", Percentage)), size = 2.5) +
  scale_fill_gradient(low = "white", high = "red") +
  labs(
    title = "Nights in the hospital last year",
    x = "Age Group",
    y = "Hospital Nights",
    fill = "Prevalence (%)"
  ) +
  theme_minimal() +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 8),    
        plot.title = element_text(face = "bold", hjust = 0.5,size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0)
  )


###Healthcare Usage table: Specialists visited in the last year (mean rate)###
#clean "More than" (impute with 13)
symptom <- current2 %>% dplyr::select(Intake.ID,Age.At.Survey, immunologist_allergist12,endocrinologist12,metabolic_specialist12,vision_specialist12,gastroenterologist12,hepatologist12,cardiologist12,electrophysiologist12,hematologist12,neurologist12,nutritionist_dietition12,orthopedic_specialist12)%>%
  mutate(across(where(is.character), ~ gsub("More than", "12", .))) %>%
  mutate(across(where(is.character), ~ gsub("999", NA, .))) %>%
  mutate(across(where(is.character), ~ gsub("Unsure", NA, .))) %>%
  mutate(across(where(is.character), as.numeric))
symptom$AGE<-current2$AGE
df <- symptom

specialist_cols <- colnames(df)[grepl("12$", colnames(df)) & colnames(df) != "Age.At.Survey"]

#Collapse each (Intake.ID, Age.At.Survey) into a single "individual" row using max
df_individual <- df %>%
  group_by(Intake.ID, Age.At.Survey) %>%
  summarise(
    AGE = first(AGE),
    across(all_of(specialist_cols), ~ if (all(is.na(.))) NA_real_ else max(.x, na.rm = TRUE)), .groups = "drop")

#Add row for total visits across all specialists
df_individual$total_visits <- rowSums(df_individual[specialist_cols], na.rm = TRUE)

#Reshape to long format
df_long <- df_individual %>%
  pivot_longer(cols = c("total_visits", all_of(specialist_cols)), names_to = "specialist", values_to = "count")

#Function to compute Poisson rate + 95% CI
poisson_ci <- function(counts) {
  total <- sum(counts, na.rm = TRUE)
  n <- sum(!is.na(counts))
  pt <- poisson.test(total, T = n)
  sprintf("%.1f ± %.1f", pt$estimate, (pt$conf.int[2]-pt$conf.int[1])/2)
}

# Group by AGE and specialist to get table
summary_by_age <- df_long %>%
  group_by(AGE, specialist) %>%
  summarise(result = poisson_ci(count), .groups = "drop")
table_by_age <- summary_by_age %>%
  pivot_wider(names_from = AGE, values_from = result)
overall_col <- df_long %>%
  group_by(specialist) %>%
  summarise(Overall = poisson_ci(count), .groups = "drop")

final_table <- left_join(table_by_age, overall_col, by = "specialist")
Healthcare_table <- final_table %>%
  mutate(specialist = ifelse(specialist == "total_visits", "Overall", specialist),
         ) %>%
  arrange(factor(specialist, levels = c("Overall", setdiff(final_table$specialist, "Overall"))))



                     
