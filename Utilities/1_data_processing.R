# load the required libraries and scripts
library("igraph")
source("R/utils.R")

# perform data processing
results = list()
clinical_information_data = list()
corrected_genotypes = list()
observed_genotypes = list()
for(f in list.files("results_files/")) {
    results[[gsub("_aic.txt","",f)]] = list()
    load(paste0("RData/",gsub("_aic.txt","",f),"/data.RData"))
    load(paste0("full_RData/",gsub("_aic.txt","",f),"/clinical_information.RData"))
    load(paste0("full_RData/",gsub("_aic.txt","",f),"/point_mutations.RData"))
    hesbcn = import.hesbcn(paste0("results_files/",f),genes=as.character(data$annotations[,"event"]))
    results[[gsub("_aic.txt","",f)]][["hesbcn"]] = hesbcn
    predictability = compute.predictability(hesbcn)
    results[[gsub("_aic.txt","",f)]][["predictability"]] = predictability
    valid_genotypes = valid.genotypes(hesbcn,predictability)
    results[[gsub("_aic.txt","",f)]][["valid_genotypes"]] = valid_genotypes
    patients_progression = patients.progression(valid_genotypes,data$genotypes,hesbcn$epsilon)
    results[[gsub("_aic.txt","",f)]][["patients_progression"]] = patients_progression
    clinical_information = clinical_information[which(clinical_information[,"PATIENT_ID"]%in%rownames(patients_progression)),]
    rownames(clinical_information) = clinical_information[,"PATIENT_ID"]
    clinical_information[,"PATIENT_ID"] = NULL
    clinical_information_tmp = array(NA,c(nrow(patients_progression),ncol(clinical_information)))
    rownames(clinical_information_tmp) = rownames(patients_progression)
    colnames(clinical_information_tmp) = colnames(clinical_information)
    clinical_information_tmp[rownames(clinical_information),colnames(clinical_information)] = as.matrix(clinical_information[rownames(clinical_information),colnames(clinical_information)])
    clinical_information = clinical_information_tmp
    clinical_information_data[[gsub("_aic.txt","",f)]] = clinical_information
    corrected_genotypes[[gsub("_aic.txt","",f)]] = corrected.genotypes(valid_genotypes,data$genotypes,hesbcn$epsilon)
    curr_observed_genotypes = data$genotypes
    curr_observed_genotypes = curr_observed_genotypes[rownames(corrected_genotypes[[gsub("_aic.txt","",f)]]),]
    colnames(curr_observed_genotypes) = colnames(corrected_genotypes[[gsub("_aic.txt","",f)]])
    observed_genotypes[[gsub("_aic.txt","",f)]] = curr_observed_genotypes
}
clinical_information = clinical_information_data

# save results
save(clinical_information,file="final_results/clinical_information.RData")
save(corrected_genotypes,file="final_results/corrected_genotypes.RData")
save(observed_genotypes,file="final_results/observed_genotypes.RData")
save(results,file="final_results/results.RData")
