#REVISED CLUSTER ANALYSIS FOR CANCER
#CREATED BY jOYCE, MAY 28, 2010
#Modified June 15, 2010. Part 4
#LAST SECTION DELETED, AUG. 25, 2010
#TESTED 8/3/2013

library(splines)          
library(survival)
library(cluster)
library(stats)

######################################################################
# 		PART 0: CHOOSE DATA AND ENTER IMPORTANT PARAMETERS                 
#
######################################################################

data<-read.table("c1.txt",sep = "\t",header=T)
attach(data)

vars<-cbind(Grade, size, nodes, ER)
### MUST INDICATE NUMBER OF VARIABLES
numVar<-4
elimSize<-100

######################################################################
# 		PART 1: CREATE A WORKING FILE WITH USABLE COMBOS                 
#
######################################################################

## create a unique identifying code for each group combination

index<-seq(1,numVar,1)
mat<-matrix(0,nrow(data),numVar)
for (i in 1:numVar){
	mat[,i]=vars[,index[i]]*(10^(max(index)-i))
}
code<-rowSums(mat)
mydata<-cbind(data,code)


## create the working dataset 

new<-as.matrix(rowsum(rep(1,nrow(mydata)),mydata$code))
new2<-new[new[,1]>elimSize,] 
#new2 contains a list of combo codes and num of patients in each combo

codenames<-as.numeric(rownames(as.matrix(new2)))

ind<-match(mydata$code,codenames)
indtemp<-ifelse(!is.na(match(mydata$code,codenames)),1,0)

mydata<-cbind(mydata,ind,indtemp)

testData<-mydata[mydata$indtemp==1,] # testData has only the usable combos

numGroup<-nrow(as.matrix(new2)) 


######################################################################
# 		PART 2: CLUSTER EMSEMBLING                 
#
######################################################################

## change the next 5 lines, if necessary, before running this section
n=numGroup
N<-1000
k1<-2
k2<-numGroup-1
Alg<-"complete" # doublecheck the method to be used


# myChisq is the initial dissimilarity matrix 
# (larger the value, the less similar the two survival curves)

myChisq<-matrix(0, nrow = numGroup, ncol=numGroup)
for (i in 1:(numGroup-1)){
	for (j in (i+1):numGroup){	
		testsample<-0
		temp<-0
		testsample<-testData[testData$ind==i | testData$ind==j,]
		temp<-survdiff(Surv(time, delta) ~ ind, data=testsample)
		myChisq[i,j]<-temp$chisq
		myChisq[j,i]<-temp$chisq
	}
}

testChisq<-0
for (i in 1:(numGroup-1))
	for (j in (i+1):numGroup)
		{testChisq[(i-1)*numGroup-i*(i-1)/2+j-i]<-myChisq[i,j]}

kpam.1<-function(x,clusnum)
	{
	classes<-pam(x, clusnum, diss=T)$cluster # vector of group assignment
	dissimilarity<-c()
	for (i in 1:n)
		{ 
		dissimilarity<-cbind(dissimilarity, classes!=classes[i])
		}
	dissimilarity
	}

## learning dissimiarity and hierachical clustering
ensemble<-function(x,N,k1,k2,Alg,n){
	dissimilarity<-matrix(rep(FALSE, n^2), ncol=n)
	
	for (i in 1:N) # number of times two groups assigned to two different clusters
		{
		dissimilarity<-dissimilarity + kpam.1(x,sample(k1:k2,1,replace=T))
		}

	dissimilarity<-dissimilarity / N 
	
	dis<-dist(dissimilarity) # create a distance matrix in order to use the structures
	
	k<-1
	for (i in 1:(n-1))
		{
		for (j in (i+1):n)
			{
			dis[k]<-dissimilarity[j,i]
			k      <- k + 1
			}
		}

	tempClust<-hclust(dis,method=Alg)
	tempClust
}

myClust<-ensemble(testChisq,N,k1,k2,Alg,n)
#Stop here, then start PART 3 
######################################################################
# 		PART 3: DENDROGRAMS                 
#
######################################################################
## Select number one or two depending if a 5 or 10 year survival rate is desired

##----------------------------
## 1
## approximate 5-year survival rate of each combo

surv5yr<-matrix(0,numGroup,1)
for (i in 1:numGroup){
	obj<-0
	obj<-summary(survfit(Surv(time,delta)~code,data=testData[testData$ind==i,]))
	surv5yr[i,]<-round( max(obj$surv[obj$time==60-min(abs(obj$time-60)) | obj$time==60+min(abs(obj$time-60))]),2)}

##---------------------------
## 2
## approximate 10-year survival rate of each combo

surv10yr<-matrix(0,numGroup,1)
for (i in 1:numGroup){
	obj<-0
	obj<-summary(survfit(Surv(time,delta)~code,data=testData[testData$ind==i,]))
	surv10yr[i,]<-round( max(obj$surv[obj$time==120-min(abs(obj$time-120)) | obj$time==120+min(abs(obj$time-120))]),2)
}

###################################################################
#STOP TEMPORARILY HERE. SELECT ONE OF THE FOLLWOING 2 POSSIBILITIES
#AFTER SELECTION, CONTINUE IN PART 4


## creating code name labels for the dendrogram
##newlab1<-paste(100*surv5yr,": " ,codenames,sep="")
newlab2<-paste(100*surv10yr,": " ,codenames,sep="")


##SELECT ONLY ONE OF THE FOLLOWING 4 POSSIBILITIES AS A LABEL
## THE FOLLOWING LINE PRODUCES DENDROGRAM WITH COMBO CODES AS LABELS
plclust(myClust,labels=as.character(Clusters))	

## PLOT THE ORIGINAL DENDROGRAM WITH "5-YR" SURVIVAL RATES AS LABELS (named "Clusters")
##plclust(myClust,labels=as.character(100*surv5yr))	

## PLOT THE ORIGINAL DENDROGRAM WITH "10-YR" SURVIVAL RATES AS LABELS (named "Clusters")
plclust(myClust,labels=as.character(100*surv10yr))	

## PLOT THE ORIGINAL DENDROGRAM WITH CODENAMES AND "5-YR" SURVIVAL RATES AS LABELS (named "Clusters")
##plclust(myClust,labels=as.character(newlab1))	

## PLOT THE ORIGINAL DENDROGRAM WITH CODENAMES AND "10-YR" SURVIVAL RATES AS LABELS (named "Clusters")
plclust(myClust,labels=as.character(newlab2))
	

###########################################################################
#CREATES VECTOR OF COMBOS (NOT ALWAYS NECESSARY TO INCLUDE), OR ELSE
#GO TO PARTs 7 and 8 TO CREATE OUTPUT FOR ALL COMGINATIONS
###########################################################################

## CREATES A MATRIX OF CODENAMES AND SURVIVAL
ordercombo=cbind(codenames, 100*surv10yr);	
ordersurvival=cbind(100*surv10yr, codenames);
os=ordersurvival; 
oc=ordercombo;

######################################################################
# 		PART 4: PLOT SURVIVAL CURVES FOR THE COMBOS               
#
######################################################################

## enter combo codes of interest in the following line, in no particular order
plotGroups<-cbind(1112,2111, 3112)

## no need to change the code below
## the following code generates survival curves for the groups specified above

numPlot<-sum(!is.na(plotGroups))
plotData<-testData[as.numeric(!is.na(match(testData$code,sort(plotGroups)    )))==1,]
test<-survfit(Surv(time,delta)~code, data=plotData)
logtest<-(survdiff(Surv(time,delta)~code,data=plotData))
logtest # log-rank test results (comparing the groups specified above)
pval=round(1-pchisq(logtest$chisq,(numPlot-1)),7) # significance level

plot(test, xmax=120, ylim=c(.5,1.0),#main=cbind("logrank=",round(logtest$chisq,4),";","pval = ",pval),
	xlab="Survival Time in  Months",
	ylab="Proportion Surviving",cex=3, col=1:numPlot, lty=1:numPlot, lwd=4, mark.time=FALSE) 
#legend("center", "bottom", as.character(sort(plotGroups)), col=1:numPlot, lty = 1:numPlot, xjust=1, yjust=1)

codenames


######################################################################
# 		PART 5: PLOT HAZARD RATES FOR THE COMBOS                
#
######################################################################

library(muhaz)
data1<-testData[testData$code==11121,]
data2<-testData[testData$code==21122,]
data3<-testData[testData$code==31221,]

Hazard1<-muhaz(data1$time, data1$delta, max.time=60)
Hazard2<-muhaz(data2$time, data2$delta, max.time=60)
Hazard3<-muhaz(data3$time, data3$delta, max.time=60)

#in the line below, adjust the ylim based on the appropriate maximum y-value for the graph
plot(Hazard1, ylim=c(0, 0.035), xlab="Time in Months", col=1, lty=1:3,lwd=4)
lines.muhaz(Hazard2, col=2, lty=2,lwd=4)
lines.muhaz(Hazard3, col=4, lty=4,lwd=5)
legend(locator(1), c("T2, N0, G1, ER+, 2662 cases ", "T2, N0, G2, ER+, 10898 cases", "T3, N3, G3, ER-,  323 cases"), col=c(1,2,4), lty=c(1,2,4), lwd=2, xjust=1, yjust=1)
## note that the legend needs to be changed
##once the graph appears left click where you would like the upper righthand corner of the legend to appear


######################################################################
# 		PART 6: CREATE CLUSTER GROUPS                
#
######################################################################

## CUT DOWN DENDROGRAM TO MAKE FEWWER LEAVES
## specify number of clusters (boxes) to be draw on the dendrogram
## change the following codes
plclust(myClust,labels=as.character(codenames))	
assign<-rect.hclust(myClust, k=numclust, border="red") # draw boxes

testData<-cbind(testData,cluster=0)
for (i in 1:numclust) {
testData$cluster[!is.na(match(testData$code, codenames[assign[[i]]]))] <- i
}
numclust=8 #Produces number of clusters

## no need t
## results checking
## print out combos in each cluster (as indicated by the dendrogram)
## rows are combo codes and columns are cluster labels
table(testData$code,testData$cluster)


################################################################
## PLOT SURVIVAL CURVES FOR CLUSTERS (one or more on the same plot)

#################################################################

## enter cluster numbers of interest in the following line, in no particular order
plotGroups<-cbind(1,2,3,4,5,6,7,8)

## no need to change the code below
## the following code generates survival curves for the groups specified above

numPlot<-sum(!is.na(plotGroups))

plotData<-testData[as.numeric(!is.na(match(testData$cluster,sort(plotGroups)    )))==1,]

test<-survfit(Surv(time,delta)~cluster, data=plotData)
logtest<-(survdiff(Surv(time,delta)~cluster,data=plotData))
logtest # log-rank test results (comparing the groups specified above)
pval=round(1-pchisq(logtest$chisq,(numPlot-1)),7) # significance level
#Does not print logrank results
plot(test, xmax=120, ylim=c(.5, 1.0), #main=cbind("logrank=",round(logtest$chisq,4),";","pval = ",pval),
	xlab="Survival Time in  Months",
	ylab="Proportion Surviving",cex=2, col=1:numPlot, lwd=4, lty=1:numPlot, mark.time=FALSE) 
#legend("center", "bottom", as.character(sort(plotGroups)), col=1:numPlot, lty = 1:numPlot, xjust=1, yjust=1)


##----------------------------------
## PLOT HAZARD RATES FOR CLUSTERS (one of more on the same plot)
## DATA<-lines of testData should equal number of plot groups

library(muhaz)
data1<-testData[testData$cluster==1,]
data2<-testData[testData$cluster==2,]
data3<-testData[testData$cluster==3,]
data4<-testData[testData$cluster--4,]

Hazard1<-muhaz(data1$time, data1$delta, max.time=120)
Hazard2<-muhaz(data2$time, data2$delta, max.time=120)
Hazard3<-muhaz(data3$time, data3$delta, max.time=120)
Hazard4<-muhaz(data4$time, data4$delta, max.time=120)

plot(Hazard1, ylim=c(0, 0.045), xlab="Time in Months", col=1, lty=1:3,lwd=4)
lines.muhaz(Hazard2, col=2, lty=2,lwd=4)
lines.muhaz(Hazard3, col=4, lty=4,lwd=5)
lines.muhaz(Hazard4, col=5, lty=5,lwd=6)
#legend(120, .03, c("T2, N0, G1, ER+, 2662 cases ", "T2, N0, G2, ER+, 10898 cases", "T3, N3, G3, ER-,  323 cases"), col=c(1,2,4), lty=c(1,2,4), lwd=2, xjust=1, yjust=1)


######################################################################
# 		PART 7: 
#		GENERATE A SPREADSHEET OF SURVIVAL RATES FOR EACH COMBO
#
######################################################################

survmat<-matrix(0,max(time)+1,numGroup)
timelist<-seq(0,max(time),1)
rownames(survmat)<-timelist
colnames(survmat)<-paste("combo#",codenames)


for (i in 1:numGroup) {
obj<-summary(survfit(Surv(time,delta)~code,data=testData[testData$ind==i,]))
survmat[as.vector(match(obj$time,timelist)),i] <- round(obj$surv,3)
}

survmatgraph<-rbind(rep(1,numGroup),survmat)
for (j in 1:ncol(survmatgraph)){
for (i in 2:nrow(survmatgraph)) {
	if (survmatgraph[i,j]==0) {survmatgraph[i,j]=survmatgraph[i-1,j]}
}}

# the following matrix (survmat) contains survival rates for each combination
survmat

# the following matrix (survmatgraph) is used for graphing purposes only
# zero survival rate is replaced by teh previous survival rate
survmatgraph



######################################################################
# 		PART 8: 
#		GENERATE A SPREADSHEET OF HAZARD RATES FOR EACH COMBO
#
######################################################################

library(muhaz)

temphaz <-muhaz(testData[testData$ind==1,]$time, testData[testData$ind==1,]$delta, max.time=60)
temphaz$pin$n.est.grid
hazmat<-matrix(0,temphaz$pin$n.est.grid,numGroup)
timelist<-seq(1,temphaz$pin$n.est.grid,1)
rownames(hazmat)<-temphaz$est.grid
colnames(hazmat)<-paste("combo#",codenames)
for (i in 1:numGroup) {
hazmat[,i]<-muhaz(testData[testData$ind==i,]$time, testData[testData$ind==i,]$delta, max.time=60)$haz.est
}

# the following matrix (hazmat) contains hazard rates for each combination
# row represents estimated event times in months (different from survival rates)
hazmat

#SAVE OUTPUT AS TEXT FILE, NEED TO CONVERT TO EXCEL
#Save file to make chart. Must select name of file
write.table(survmatgraph, file="SurFourVariables.txt", sep="\t", col.names=T,row.names=F)

#Save file to make spreadsheet of hazard data.
write.table(hazmat, file="HazFourVariables.txt", sep="\t", col.names=T,row.names=F)

codenames
ordercombo
ordersurvival=os[order(100*surv10yr),]


#Stop here, next 3 functions eliminates all variables in program
###################################################################

length(codenames)
 rm(list=ls())
 objects()

 


#Stop Here
#---------------------------------------------------------------------------
