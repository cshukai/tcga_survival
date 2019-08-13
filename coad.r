library(RTCGA)
library(RTCGA.clinical)

#################################EDA#################################
stromaIdx=colnames(COAD.clinical)[grep(x=colnames(COAD.clinical),pattern="stroma")] # find where stroma -related colnames are

#################################STROMA#################################
# normalization
for(i in 1:nrow(COAD.clinical)){
}
# stroma percentage survival
clin=survivalTCGA(COAD.clinical)




