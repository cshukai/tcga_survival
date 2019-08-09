setwd("../Downloads/")
d=read.csv("COADREAD_70Vars_ColumnWise.csv",header=T)
d=d[which(d["admin.disease_code"]=="read"),]
table(d[,"patient.history_of_neoadjuvant_treatment"])


 
