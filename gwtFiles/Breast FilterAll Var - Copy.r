##PREPARES SEER BREAST CACER DATA FOR CLUSTER ANALYSIS
##SELECTS FOR BASIC PROGNOSTIC FACTORS FOR BREAST CANCER
##REVIEWED DEC 12, 2010
##CREATES MULTIPLE DATA SETS WITH DIFFERENT VARIABLES
##NOTE-Variables must be same as output under "write table"
#READS OM "BreastDataSet" LOCATED IN FOLDER "DATA SETS"
#MODIFIED 3/16/2012
#MODIFIED 7/1/2013
#########################################################
#########################################################
dataa<-read.table("BreastDataSet.txt",sep = "\t",header=TRUE);
dim(dataa);

data1=dataa;
##########################################################

#Reads in 9 variabiles
data<-dataa[,c("size","nodes","Grade","agedx","race","ER","vs","time","delta")];
dim(data)

#########################################################
#Eliminates non Grade cases
data1=data[data$Grade>0 & data$Grade<5,];
remove(data);

##########################################################
#Selects all cases less than 10 cm in size
data2=data1[data1$size>0 & data1$size<100,];

dim(data2)
##########################################################
#Selects all patients with less than or equal to 20 nodes
data4=data2[data2$nodes<=20,];
dim(data4)
remove(data1);

#Selects for hormone receptor results
#Projesterone not used in this script
data5=data4[data4$ER>0 & data4$ER<3,];
dim(data5)
data5[2:30,]

#####################################################
#Converts alive to 0
data5[data5$vs==1,]$vs<-0;

#Converts deaths from #4 to #1 for calcuation
data5[data5$vs==4,]$vs<-1;
dim(data5)

remove(data4);
#######################################################3
#Histology, set for infiltrating ductal and lobular
#data5[data5$histo==500,]$histo<-1; 
#data5[data5$histo==140,]$histo<-1;
#data5[data5$histo==520,]$histo<-2;
#dim(data5)
########################################################
#Recodes sizes according to 3 T categories for breast cancer.
#Needs to be changed for other sites.

data4=data5;
data4[data4$size>0 & data4$size<=20,]$size<-1;
data4[data4$size>20 & data4$size<=50,] $size<-2;
data4[data4$size>50 & data4$size<=100,]$size<-3;

remove(data5);
 #######################################################
#Recodes nodes
data4[data4$nodes==0,]$nodes<-0;
#Need 0 here first, then turn to 1 later
data4[data4$nodes>0 & data4$nodes<4,]$nodes<-2;
data4[data4$nodes>3 & data4$nodes<=10,]$nodes<-3;
data4[data4$nodes>10 & data4$nodes<=20,]$nodes<-4;
data4[data4$nodes==0,]$nodes<-1;
#Puts nodes for #0 to 1, avoid conflict with 1-4 group

###########################################################
#Converts grade 4 to grade 3
data4[data4$Grade==4,]$Grade<-3;

############################################################
#Converts age to less than or more than 50.
data4[data4$agedx>50,]$agedx<-2;
data4[data4$agedx!=2,]$agedx<-1;

############################################################
#Converts to specific cause of death	
data4[data4$delta!=46,]$delta<-0;
data4[data4$delta==46,]$delta<-1;

##############################################################
data4[10:25,];

##########################################################
#CREATS MULTIPLE DATA SETS WITH ANY NUMBER OF VARIABLES
dim(data4)
##############################################################
#Saves major data set and all variables
data5=data4[,c("size","nodes","Grade","race","agedx","ER","delta")];
	dim(data5)
write.table(data5, file="test5var.txt", sep="\t", col.names=TRUE, row.names=FALSE)

################################################################
data4[,c("time","Grade","size","nodes","ER","agedx","delta")];
data5=data4;
	dim(data5)
write.table(data5, file="test5var.txt", sep="\t", col.names=TRUE, row.names=FALSE)

########################################################

data6=data4[,c("time","size","nodes","ER","delta")];
write.table(data6, file="test3var.txt", sep="\t", col.names=TRUE, row.names=FALSE)
  dim(data6)
remove(data6)

##############################################################
data7=data4[,c("Grade","time","size","nodes","ER","agedx","delta")];
write.table(data7, file="testvar.txt", sep="\t", col.names=TRUE, row.names=FALSE)
  dim(data7)
  
 ############################################################ 
#Removes all data
rm(list=ls())
objects()

