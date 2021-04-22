# load required libraries
library("glmnet")
library("ggplot2")
library("gtools")
library("gridExtra")
library("igraph")
library("Rmpfr")
library("survival")
library("survminer")

# set the seed
set.seed(12345)

# load results
load(file="final_results/clinical_information.RData")
load(file="final_results/model_based_features.RData")
load(file="final_results/results.RData")

# preprocessing
clinical_information_cancers = NULL
for(cancer in names(clinical_information)) {
    clinical_information_cancers = rbind(clinical_information_cancers,clinical_information[[cancer]][,c("AGE_DIAGNOSIS","OS_STATUS","OS_DAYS","PFI_STATUS","PFI_DAYS")])
}
clinical_information = clinical_information_cancers
PATIENTS = rownames(clinical_information)
OS_STATUS = as.character(clinical_information[,"OS_STATUS"])
OS_YEARS = as.numeric(clinical_information[,"OS_DAYS"])
OS_STATUS[which(OS_YEARS<30)] = NA
OS_YEARS[which(OS_YEARS<30)] = NA
OS_YEARS = OS_YEARS/365
OS_STATUS[which(OS_YEARS>5)] = "ALIVE"
OS_YEARS[which(OS_YEARS>5)] = 5
OS_STATUS[which(as.numeric(clinical_information[,"AGE_DIAGNOSIS"])>80)] = NA
OS_YEARS[which(as.numeric(clinical_information[,"AGE_DIAGNOSIS"])>80)] = NA
OS_STATUS[which(OS_STATUS=="ALIVE")] = 0
OS_STATUS[which(OS_STATUS=="DEAD")] = 1
OS_STATUS = as.numeric(OS_STATUS)
PFI_STATUS = as.character(clinical_information[,"PFI_STATUS"])
PFI_YEARS = as.numeric(clinical_information[,"PFI_DAYS"])
PFI_STATUS[which(PFI_YEARS<30)] = NA
PFI_YEARS[which(PFI_YEARS<30)] = NA
PFI_YEARS = PFI_YEARS/365
PFI_STATUS[which(PFI_YEARS>5)] = "TUMOR_FREE"
PFI_YEARS[which(PFI_YEARS>5)] = 5
PFI_STATUS[which(as.numeric(clinical_information[,"AGE_DIAGNOSIS"])>80)] = NA
PFI_YEARS[which(as.numeric(clinical_information[,"AGE_DIAGNOSIS"])>80)] = NA
PFI_STATUS[which(PFI_STATUS=="TUMOR_FREE")] = 0
PFI_STATUS[which(PFI_STATUS=="WITH_TUMOR")] = 1
PFI_STATUS = as.numeric(PFI_STATUS)

# perform analysis
dev.new(width=15,height=10)
survival_analysis = list()
survival_analysis[["OS"]] = list()
survival_analysis[["PFI"]] = list()
for(i in names(model_based_features)) {

    # create data structures
    curr_OS_YEARS = NULL
    curr_OS_STATUS = NULL
    curr_PFI_YEARS = NULL
    curr_PFI_STATUS = NULL
    for(j in rownames(model_based_features[[i]])) {
        curr_OS_YEARS = c(curr_OS_YEARS,OS_YEARS[which(PATIENTS==j)])
        curr_OS_STATUS = c(curr_OS_STATUS,OS_STATUS[which(PATIENTS==j)])
        curr_PFI_YEARS = c(curr_PFI_YEARS,PFI_YEARS[which(PATIENTS==j)])
        curr_PFI_STATUS = c(curr_PFI_STATUS,PFI_STATUS[which(PATIENTS==j)])
    }
    curr_lambda = as.numeric(results[[i]][["patients_progression"]][rownames(model_based_features[[i]]),])
    curr_data = cbind(curr_OS_STATUS,curr_OS_YEARS,curr_lambda,model_based_features[[i]])
    rownames(curr_data) = 1:nrow(curr_data)
    colnames(curr_data)[1:3] = c("Status","Times","Lambda")
    os_data = data.frame(curr_data)
    curr_data = cbind(curr_PFI_STATUS,curr_PFI_YEARS,curr_lambda,model_based_features[[i]])
    rownames(curr_data) = 1:nrow(curr_data)
    colnames(curr_data)[1:3] = c("Status","Times","Lambda")
    pfi_data = data.frame(curr_data)

    # perform analysis for OS survival
    if(length(which(is.na(os_data$Status)|is.na(os_data$Times)))>0) {
        os_data = os_data[-which(is.na(os_data$Status)|is.na(os_data$Times)),]
    }
    y = Surv(as.numeric(os_data$Times),as.numeric(os_data$Status))
    lasso_cov = os_data[,colnames(os_data)[3:ncol(os_data)]]
    string_test = paste0("~",paste0("lasso_cov$",colnames(lasso_cov),collapse="+"))
    x = model.matrix(as.formula(string_test))
    cv.fit = cv.glmnet(x,y,family="cox",maxit=10000,alpha=1)
    survival_analysis[["OS"]][[i]][["data"]] = os_data
    survival_analysis[["OS"]][[i]][["fit"]] = cv.fit
    print(plot(cv.fit,main=paste0(i," - OS Survival")))

    # perform analysis for PF interval
    if(length(which(is.na(pfi_data$Status)|is.na(pfi_data$Times)))>0) {
        pfi_data = pfi_data[-which(is.na(pfi_data$Status)|is.na(pfi_data$Times)),]
    }
    y = Surv(as.numeric(pfi_data$Times),as.numeric(pfi_data$Status))
    lasso_cov = pfi_data[,colnames(pfi_data)[3:ncol(pfi_data)]]
    string_test = paste0("~",paste0("lasso_cov$",colnames(lasso_cov),collapse="+"))
    x = model.matrix(as.formula(string_test))
    cv.fit = cv.glmnet(x,y,family="cox",maxit=10000,alpha=1)
    survival_analysis[["PFI"]][[i]][["data"]] = pfi_data
    survival_analysis[["PFI"]][[i]][["fit"]] = cv.fit
    print(plot(cv.fit,main=paste0(i," - PF Interval")))

}

# save results
save(survival_analysis,file="final_results/survival_analysis.RData")
