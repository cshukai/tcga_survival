library(RTCGA)
library(RTCGA.clinical)

clin=survivalTCGA(COAD.clinical)

colnames(COAD.clinical)[grep(x=colnames(COAD.clinical),pattern="stroma")] # find where stroma -related colnames are

