library(RTCGA)
library(RTCGA.clinical)
library(RTCGA.mRNA)

clin = survivalTCGA(COAD.clinical)
# parse barcode
COAD.mRNA
d=merge(clin,COAD.mRNA,by="bcr_patient_barcode")
