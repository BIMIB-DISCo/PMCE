# load data
load(file="final_results/corrected_genotypes.RData")
load(file="final_results/results.RData")

# consider each cancer type
model_based_features = list()
for(cancer in names(corrected_genotypes)) {
    
    # consider current cancer type
    curr_genotype = corrected_genotypes[[cancer]]
    curr_adj_matrix = results[[cancer]][["hesbcn"]][["adjacency_matrix"]]
    
    # consider each patient
    arcs = NULL
    for(i in 1:nrow(curr_genotype)) {
        if(length(which(curr_genotype[i,]!="*"))>0) {
            patient_adj_matrix = curr_adj_matrix[c("Root",names(which(curr_genotype[i,]!="*"))),c("Root",names(which(curr_genotype[i,]!="*")))]
            curr_arcs = NULL
            for(j in colnames(patient_adj_matrix)[2:ncol(patient_adj_matrix)]) {
                if(results[[cancer]][["hesbcn"]][["parent_set"]][[j]]=="Single") {
                    curr_arc = paste0(names(which(curr_adj_matrix[,j]==1))," to ",j)
                }
                else {
                    curr_arc = paste0("OR(",paste0(names(which(curr_adj_matrix[,j]==1)),collapse=","),")"," to ",j)
                }
                curr_arcs = cbind(rownames(curr_genotype)[i],curr_arc)
            }
            arcs = rbind(arcs,curr_arcs)
        }
    }
    
    # build arcs matrix
    arcs_matrix = array(0,c(nrow(curr_genotype),length(unique(arcs[,2]))))
    rownames(arcs_matrix) = rownames(curr_genotype)
    colnames(arcs_matrix) = sort(unique(unique(arcs[,2])))
    for(i in 1:nrow(arcs)) {
        arcs_matrix[arcs[i,1],arcs[i,2]] = 1
    }
    model_based_features[[cancer]] = data.frame(arcs_matrix)
    
}

# save results
save(model_based_features,file="final_results/model_based_features.RData")
