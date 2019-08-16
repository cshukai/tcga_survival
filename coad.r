library(RTCGA)
library(RTCGA.clinical)
library(survival)
library(survminer)
#################################STROMA#################################

stromaIdx=colnames(COAD.clinical)[grep(x=colnames(COAD.clinical),pattern="stroma")] # find where stroma -related colnames are

stroma_avg=NULL
for(i in 1:nrow(COAD.clinical)){
  reps=0
  percents=0
  this_stroma=-1
  for(j in 1:length(stromaIdx)){
    this_percent=COAD.clinical[i,stromaIdx[j]]
    if(!is.na(this_percent)){
      percents=percents+as.numeric(this_percent)
      reps=reps+1
    }
  }
  if(reps>0){
    this_stroma=percents/reps
  }
  this_row=cbind(COAD.clinical[i,"patient.bcr_patient_barcode"],this_stroma)
  stroma_avg=rbind(stroma_avg,this_row)
  
}
stroma_avg=stroma_avg[-which(stroma_avg[,2]=="-1"),]
colnames(stroma_avg)=c("barcode","str_avg")
#########################TUMOR ######################
tumorCols=colnames(COAD.clinical)[intersect(grep(x=colnames(COAD.clinical),pattern="percent_tumor_cells"),grep(x=colnames(COAD.clinical),pattern="slide"))]

tumor_avg=NULL
for(i in 1:nrow(COAD.clinical)){
  reps=0
  percents=0
  this_tumor=-1
  for(j in 1:length(tumorCols)){
    this_percent=COAD.clinical[i,tumorCols[j]]
    if(!is.na(this_percent)){
      percents=percents+as.numeric(this_percent)
      reps=reps+1
    }
  }
  if(reps>0){
    this_tumor=percents/reps
  }
  this_row=cbind(COAD.clinical[i,"patient.bcr_patient_barcode"],this_tumor)
  tumor_avg=rbind(tumor_avg,this_row)
  
}
tumor_avg=tumor_avg[-which(tumor_avg[,2]=="-1"),]
colnames(tumor_avg)=c("barcode","tumor_avg")
############################ tumor/stroma for survival####################################
tum_str=merge(tumor_avg,stroma_avg)
st_ratio=as.numeric(as.character(tum_str[,3]))/as.numeric(as.character(tum_str[,2]))
tum_str=cbind(tum_str,st_ratio)

clin=survivalTCGA(COAD.clinical)
clin[,"bcr_patient_barcode"]=tolower(clin[,"bcr_patient_barcode"])
colnames(clin)[2]="barcode"

sur_tum_str=merge(clin,tum_str)
group =ifelse(sur_tum_str$st_ratio>median(sur_tum_str$st_ratio),'high','low')

sfit =survfit(Surv(times, patient.vital_status)~group, data=sur_tum_str)
summary(sfit)
surv_pvalue(sfit)$pval
ggsurvplot(sfit, conf.int=F, pval=TRUE)

########################other image features  #################################
features=c("percent_granulocyte_infiltration", "percent_inflam_infiltration", "percent_lymphocyte_infiltration", "percent_monocyte_infiltration", "percent_necrosis", "percent_neutrophil_infiltration")

for(i in 1:length(features)){
  interestingCols=colnames(COAD.clinical)[intersect(grep(x=colnames(COAD.clinical),pattern=features[i]),grep(x=colnames(COAD.clinical),pattern="slide"))]
  avg=NULL
  for(j in 1:nrow(COAD.clinical)){
    reps=0
    percents=0
    this_avg=-1
    for(k in 1:length(interestingCols)){
      this_percent=COAD.clinical[j,interestingCols[k]]
      if(!is.na(this_percent)){
        percents=percents+as.numeric(this_percent)
        reps=reps+1
      }
    }
    if(reps>0){
      this_avg=percents/reps
    }
    this_row=cbind(COAD.clinical[j,"patient.bcr_patient_barcode"],this_avg)
    avg=rbind(avg,this_row)
    
  }
  avg=avg[-which(avg[,2]=="-1"),]
  
  
}
#"patient.samples.sample-2.portions.portion.slides.slide.percent_necrosis"
#"patient.samples.sample.portions.portion.slides.slide.percent_eosinophil_infiltration"[1]
#[1]https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4848192/
[220] "patient.samples.sample.portions.portion.slides.slide.percent_granulocyte_infiltration"
[221] "patient.samples.sample.portions.portion.slides.slide.percent_inflam_infiltration"
[222] "patient.samples.sample.portions.portion.slides.slide.percent_lymphocyte_infiltration"
[223] "patient.samples.sample.portions.portion.slides.slide.percent_monocyte_infiltration"
[224] "patient.samples.sample.portions.portion.slides.slide.percent_necrosis"
[225] "patient.samples.sample.portions.portion.slides.slide.percent_neutrophil_infiltration"

save.image("tumorstroma.RData")
