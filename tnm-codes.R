## The purpose of this script is to investiage distributions of TNM-codes in our data
## Use the prospective bc data from Run 1, with outliers removed. 


################################
## Set paths
################################
setwd("/home/hbo010/breastCancer/scripts")
wd <- getwd()
setwd("..")
dataPath <- file.path(getwd(), "data")
#resultPath <- file.path(getwd(),"results","04-03-2013") ## folder to store results
resultPath <- file.path(getwd(),"results","16-04-2013") ## folder to store results
setwd(wd)

dataIn <-  "18-04-2013" ## location of the main data file

################################
## Import libraries
################################
require(glmnet)
require(lme4)
require(ggplot2)
#install.packages("pROC")
require(pROC)
require(limma)
require(lumi)
## source("http://bioconductor.org/biocLite.R")
## biocLite("illuminaHumanv4.db")
require(illuminaHumanv4.db)

#############################
#### IMPORT DATA ############
#############################
files <- file.path(dataPath, dataIn, "prospective-bc-358cc-pairs-Outliers-removed-Run175-background-corrected-filtered.Rdata")
if(file.exists(files))
{
  print(paste("File exists. Loading data", files, sep=" "))
  load(files)
} else {
  print("File does not exist.")
}



############################
## Set some basic variables:
############################
n <- nrow(exprs) ## number of individuals
p <- ncol(exprs) ## number of probes
case_id <- seq(1,n,by=2) ## Cases are always on odd rows, controls on even rows.
ctrl_id <- seq(2,n,by=2)
## A small check:
table(background[case_id,"Case_ctrl"])
table(background[ctrl_id,"Case_ctrl"])

## Use the log2 of probe intensities:
log2.exprs <- log2(exprs)
d.log2.exprs <- log2.exprs[case_id,] - log2.exprs[ctrl_id,] ## case-control log2-differenece


##################################
## TNM codes
###################################

## Want to make a variable with the TNM code for each case 
## at the first bc diagnosis after blood sample. 

### Pick out names of spesific columns:
invasName <- grep("invas", colnames(background), value=T) #find all column names with invasiveX
insitName <- grep("insit", colnames(background), value=T) # find all column names with insituX
allT <- grep("P_T", colnames(background), value=T)
allN <- grep("P_N", colnames(background), value=T)
allM <- grep("P_M", colnames(background), value=T)

## use names to pick out the submatrix with informaiton:
invas <- background[,invasName]
insit <- background[,insitName]
size <- as.matrix(background[,allT])
node <- as.matrix(background[,allN])
meta <- as.matrix(background[,allM])

## There is probably a punching error in meta, leading to 
## some individuals having meta = " 0" instead of meta = "0". 
meta[grep(" ", meta)] 
## set these to "0": 
meta[grep(" ", meta)] <- "0"


## Make matrix with invasive1:6 and insitu1:6 logical. This will be used to pick out the right 
## size/node/meta informaiton: 
invasLogic <- sapply(invas,is.na) ## if invasive, then FALSE
insituLogic <- sapply(insit,is.na) ## if insitu, then FALSE

## coherce the logical information from insitu and invasive, and make cancer diagnosis == TRUE.
cancerDiag <- !(invasLogic & insituLogic) 

## use the above logical matrix:
nrOfBcDiagnosis <- apply(cancerDiag,1,sum) #count the nr of breast cancer diagnosis for each ind.
diagIndex <- which(nrOfBcDiagnosis!=0) ## pick out the bc cases

## Investigate the variables: 
apply(meta, 2, function(x) sum(x != "", na.rm=T))
apply(size, 2, function(x) sum(x != "", na.rm=T))
apply(node, 2, function(x) sum(x != "", na.rm=T))
apply(cancerDiag, 2, function(x) sum(x == T, na.rm=T))

## Investigate the invasive variable, to get an idea
## of how many individuals we have with multiple diagnoses
## of breast cancer, and with e.g. bc as the second cancer
## diagnosis. 
hege <- apply(invas,1,sum,na.rm=T) ## number of invasive bc diangoses
invas[which(hege >1),] ## multiple bc diagnoses
invas[which(is.na(invas[,1]) & invas[,2]==1),] ## bc as the second diag.
## From the code above, we learn that the meta/size/node information
## is only present for the breast cancer diagnoses. 
## Further, we see that the invas variable gives an 1 on the breast cancer cases, 
## but can e.g. be NA on invas1, and 1 on invas2. 
## In order to pick out the right meta/size/node information, we thus need to be careful. 
## ---

## SIZE: 
Psize <- data.frame(matrix(rep(NA,length(cancerDiag)),ncol=ncol(cancerDiag))) ## make empty data.frame
rownames(Psize) <- rownames(background) ## set rownames equal to labnr
colnames(Psize) <- paste("P_T",1:ncol(cancerDiag),sep="") ## set column names
## pick out all breast cancer diagnosis dates. Remember that diagdat contain all cancer diagnosis dates, 
## and thus we need to pick out the ones that correspond to bc diagnosis dates. 
## "diagIndex" corresponds to the row number of the bc cases.
for (ind in diagIndex){
  Psize[ind,seq_len(nrOfBcDiagnosis[ind])] <- size[ind,cancerDiag[ind,]]
}

## NODE:
Pnode <- data.frame(matrix(rep(NA,length(cancerDiag)),ncol=ncol(cancerDiag))) ## make empty data.frame
rownames(Pnode) <- rownames(background) ## set rownames equal to labnr
colnames(Pnode) <- paste("P_N",1:ncol(cancerDiag),sep="") ## set column names
## pick out all breast cancer diagnosis dates. Remember that diagdat contain all cancer diagnosis dates, 
## and thus we need to pick out the ones that correspond to bc diagnosis dates. 
for (ind in diagIndex){
  Pnode[ind,seq_len(nrOfBcDiagnosis[ind])] <- node[ind,cancerDiag[ind,]]
}

## META: 
Pmeta <- data.frame(matrix(rep(NA,length(cancerDiag)),ncol=ncol(cancerDiag))) ## make empty data.frame
rownames(Pmeta) <- rownames(background) ## set rownames equal to labnr
colnames(Pmeta) <- paste("P_M",1:ncol(cancerDiag),sep="") ## set column names
## pick out all breast cancer diagnosis dates. Remember that diagdat contain all cancer diagnosis dates, 
## and thus we need to pick out the ones that correspond to bc diagnosis dates. 
for (ind in diagIndex){
  Pmeta[ind,seq_len(nrOfBcDiagnosis[ind])] <- as.character(meta[ind,cancerDiag[ind,]])
}

## DIAGNOSIS DATE:
allDiag <- grep("DIAGNOSE", colnames(background), value=T)
diag <- as.matrix(background[,allDiag])

bcDiagDate <- data.frame(matrix(rep(NA,length(cancerDiag)),ncol=ncol(cancerDiag))) ## make empty data.frame
rownames(bcDiagDate) <- rownames(background) ## set rownames equal to labnr
colnames(bcDiagDate) <- paste("BcDiag",1:ncol(cancerDiag),sep="") ## set column names
## pick out breast cancer diagnosis dates. Remember that diagdat contain all cancer diagnosis dates, 
## and thus we need to pick out the ones that correspond to bc diagnosis dates. 
for (ind in diagIndex){
  bcDiagDate[ind,seq_len(nrOfBcDiagnosis[ind])] <- diag[ind,cancerDiag[ind,]]    
}

## We will only work with the first breast cancer diagnosis after blood sample. 
## We thus need to pick out the right date and the right TNM and stage information:
firstDiagDate <- as.Date(bcDiagDate[,1], format="%Y-%m-%d")
followUp <- firstDiagDate - background[,"BPROVEDATO"] ## Calculate follow-up time in days.
all(followUp[case_id] == background[case_id,"follow_up_time"]) ## Test if we have the same dates as calculated already.
all(followUp > 0, na.rm=T) ## Test to see that we only have positive follow-up times.
## From the above, we see that all follow up times are positive, i.e. the first diagnosis
## recored here is the one we should consider as the "first diagnosis after blood sample". 

## Check if we have individuals with two diagnosis at the same date. 
## (If we do, we should use the most severe diagnosis)
any(bcDiagDate[,1] == bcDiagDate[,2]) # logical test. We obtain "TRUE" and continue to search:
multiple.diag <- rownames(bcDiagDate)[which(bcDiagDate[,1] == bcDiagDate[,2])] ## labnr of cases with multiple diagnosis. 
any(bcDiagDate[,1] == bcDiagDate[,3]) ## just to check if there are individuals with three diag on the same date. (none)

## We observe that there are 11 cases with two diagnosis at the same date. 
## We should always use the most severe diagnosis, so check if it differs:
apply(Psize[multiple.diag,1:2], 1, unique)
apply(Pnode[multiple.diag,1:2], 1, unique) ## all are equal
apply(Pmeta[multiple.diag,1:2], 1, unique)
background[multiple.diag,c("STADIUM_B1","STADIUM_B2")] ## Check stage.

## Four individuals have a different value for Psize (labnr 146590, 109382, 118129, and 138221).
## Two individuals should change stage from 0 to 100 and 999 (labnr 118129 and 138221). 
## For labnr 146590 we also have two different choices of Pmeta: Y or X. 
## We choose to keep the size and meta information of the first bc diagnosis for labnr 146590 and 109382, 
## since the first listed is the most severe. 
## For labnr 118129 and 138221 we change pTNM and stage information to the second diagnosis, which 
## is the most severe. 

## For our purpose, we will thus only use the first column of Psize, Pnode, and Pmeta 
## exept for labnr 118129 and 138221.
pTNM_size_post_blood <- Psize[,1]
names(pTNM_size_post_blood) <- rownames(background)
pTNM_node_post_blood <- Pnode[,1]
names(pTNM_node_post_blood) <- rownames(background)
pTNM_meta_post_blood <- Pmeta[,1]
names(pTNM_meta_post_blood) <- rownames(background)

pTNM_size_post_blood[c("118129","138221")] <- Psize[c("118129","138221"),2]
pTNM_node_post_blood[c("118129","138221")] <- Pnode[c("118129","138221"),2]
pTNM_meta_post_blood[c("118129","138221")] <- Pmeta[c("118129","138221"),2]

## Change stage for labnr "118129","138221":
background[c("118129","138221"),"stage_post_blood"] <- background[c("118129","138221"),"STADIUM_B2"]

## merge information in size to make it more easy to use: 
pTNM <- cbind(background[,"stage_post_blood"],pTNM_size_post_blood, pTNM_node_post_blood, pTNM_meta_post_blood)[case_id,]
colnames(pTNM)[1] <- "stage"
rownames(pTNM) <- background[case_id,"labnr"]

## examine the values to see if the coding looks reasonable:
s0 <- pTNM[pTNM[,"stage"] == 0,]
s1 <- pTNM[pTNM[,"stage"] == 100,]
s2 <- pTNM[pTNM[,"stage"] > 100 & pTNM[,"stage"] < 300,]
s3 <- pTNM[pTNM[,"stage"] > 300,]

apply(s2,2, table)


## make the codes a bit easier to access by re-coding them: 
tnm <- pTNM
tnm <- cbind(tnm,tnm[,"pTNM_size_post_blood"],tnm[,"pTNM_node_post_blood"])
colnames(tnm) <- c("stage", "T", "N", "M","ET","EN")
tnm[grep("1",tnm[,"ET"]),"ET"] <- "1"
tnm[grep("1",tnm[,"EN"]),"EN"] <- "1"


#### RECODING OF TNM IN CORPORATION WITH NICOLLE #####
### RECODING OF T: ###
## If T is X or 4, then the "new" code ET (Eiliv's T) is "missing". 
mis <- (tnm[,"T"]=="X" | tnm[,"T"]=="4")
tnm[which(mis),"ET"] <- "missing"

## If T is (X, Y, 4, or missing) and stage is 1 (stadium_B = 100), 
## then ET = 1. (14 individuals)
oo <- (tnm[,"stage"] == 100) & (tnm[,"T"] == "X" | tnm[,"T"] == "Y" | tnm[,"T"] == "4" | tnm[,"T"] == "missing")
tnm[which(oo),"ET"] <- 1


### RECODING OF N: ###
## If stage is 1 (stadium_B=100) then EN = 0. 
## (Four individuals with N="Y" will be changed to "N"=0):
tnm[(tnm[,"stage"] == 100),"EN"] <- 0

## If stage is 2 (stadium_B=200-299), and ET is 0 or 1, then EN = 1.
## (This is true for all values of "EN" to begin with. Nothing will be changed.)
ii <- which((tnm[,"stage"] == "200" | tnm[,"stage"] == "210" | tnm[,"stage"] == "220") & (tnm[,"ET"] == 0 | tnm[,"ET"] == 1))
tnm[ii,"EN"] <- 1

## If "EN" is X or Y, then set "EN" to missing. 
mis <- tnm[,"EN"] == "X" | tnm[,"EN"] == "Y"
tnm[mis,"EN"] <- "missing"


## Make a cross-table with the recoded and original values of 
## "ET" and "T" and "EN" and "N". Compare with the cross-table
## that Nicolle sent on June 21 2013. 
table(tnm[,"T"],tnm[,"ET"])
table(tnm[,"N"],tnm[,"EN"])



# ## Pick out different subclasses of stage and summarize: 
# ## T2N0M0:
# t2n0m0 <- (tnm[(tnm[,"T"] == 2 & tnm[,"N"] == 0 & tnm[,"M"] == 0),])
# nrow(t2n0m0)
# which((tnm[,"T"] == 2 & tnm[,"N"] == 0 & tnm[,"M"] == 0))
# table(years_to_diag[which((tnm[,"T"] == 2 & tnm[,"N"] == 0 & tnm[,"M"] == 0))]) ## distribution of years to diag for t2n0m0
# 
# 
# length(which((tnm[,"T"] == 3 & tnm[,"N"] == 0 & tnm[,"M"] == 0)))
# length(which((tnm[,"T"] == 4 & tnm[,"N"] == 0 & tnm[,"M"] == 0)))
# 
# table(years_to_diag[which((tnm[,"T"] == 3 & tnm[,"N"] == 0 & tnm[,"M"] == 0))])
# table(years_to_diag[which((tnm[,"T"] == 4 & tnm[,"N"] == 0 & tnm[,"M"] == 0))])

## Save information to file: 
## Add labnr: 
tnm <- cbind(background[case_id,"labnr"],tnm)
colnames(tnm)[1] <- "labnr"

files = file.path(dataPath,"24-06-2013/tnm-information",paste("tnm-information",n/2,"cc-pairs.Rdata", sep="-"))
save(tnm,file=files)


