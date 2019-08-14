library(RTCGA)
library(RTCGA.clinical)

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
############################ x=colnames(COAD.clinical),pattern="stroma"tumor/stroma for survival####################################
clin=survivalTCGA(COAD.clinical)

save.image("tumorstroma.RData")
