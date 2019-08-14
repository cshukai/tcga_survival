library(RTCGA)
library(RTCGA.clinical)

#################################EDA#################################
stromaIdx=colnames(COAD.clinical)[grep(x=colnames(COAD.clinical),pattern="stroma")] # find where stroma -related colnames are

#################################STROMA#################################
# 
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
# stroma percentage survival
clin=survivalTCGA(COAD.clinical)
