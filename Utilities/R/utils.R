# import the inferred HESBCN to be used in the analyses
"import.hesbcn" <- function( file, genes = NULL ) {

    # read results from file
    results = read.table(file=file,header=FALSE,sep="\n",stringsAsFactors=FALSE,check.names=FALSE)

    # build adjacency matrix
    start = 2
    end = grep("accepted",results[,1])-1
    adjacency_matrix = NULL
    for(entry in results[start:end,1]) {
        adjacency_matrix = rbind(adjacency_matrix,as.numeric(strsplit(entry,split=",")[[1]]))
    }
    rownames(adjacency_matrix) = paste0("Gene_",1:nrow(adjacency_matrix))
    colnames(adjacency_matrix) = paste0("Gene_",1:ncol(adjacency_matrix))

    # get lambda values and parent set types
    lambdas = results[(grep("Best Lambdas",results[,1])+1),1]
    lambdas = strsplit(lambdas,split=" ")[[1]][1:(length(strsplit(lambdas,split=" ")[[1]])-1)]
    lambdas = as.numeric(lambdas)
    parent_set = strsplit(results[(grep("best theta types",results[,1])+1),1],split=",")[[1]]
    for(i in 1:length(parent_set)) {
        if(parent_set[i]=="-1") {
            parent_set[i] = "Single"
        }
        else if(parent_set[i]=="0") {
            parent_set[i] = "AND"
        }
        else if(parent_set[i]=="1") {
            parent_set[i] = "OR"
        }
        else if(parent_set[i]=="2") {
            parent_set[i] = "XOR"
        }
    }
    names(parent_set) = paste0("Gene_",1:nrow(adjacency_matrix))

    # make final results and normalize lambda values
    adjacency_matrix = cbind(rep(0,nrow(adjacency_matrix)),adjacency_matrix)
    adjacency_matrix = rbind(rep(0,ncol(adjacency_matrix)),adjacency_matrix)
    rownames(adjacency_matrix)[1] = "Root"
    colnames(adjacency_matrix)[1] = "Root"
    adjacency_matrix["Root",as.numeric(which(colSums(adjacency_matrix)==0))[-1]] = 1
    lambdas_matrix = adjacency_matrix
    for(i in 2:nrow(lambdas_matrix)) {
        curr_in_lambda = lambdas[(i-1)]
        # we assume equally distributed lambdas among components of a single formula
        lambdas_matrix[as.numeric(which(adjacency_matrix[,i]==1)),i] = curr_in_lambda / length(which(adjacency_matrix[,i]==1))
    }

    # set genes names
    if(!is.null(genes)) {
        rownames(adjacency_matrix) = c("Root",genes)
        colnames(adjacency_matrix) = c("Root",genes)
        rownames(lambdas_matrix) = c("Root",genes)
        colnames(lambdas_matrix) = c("Root",genes)
        names(parent_set) = genes
    }

    # estimated epsilon
    epsilon = as.numeric(gsub("Best epsilon = ","",results[nrow(results),1]))
    
    # return results
    hesbcn = list(adjacency_matrix=adjacency_matrix,lambdas_matrix=lambdas_matrix,parent_set=parent_set,epsilon=epsilon)
    return(hesbcn)

}

# compute predictability of the HESBCN
"compute.predictability" <- function( hesbcn ) {

    # compute needed statistics to calculate predictability from graph
    results = list()
    graph = graph_from_adjacency_matrix(hesbcn$adjacency_matrix)
    paths_predictability = NULL
    number_path = 0
    for(i in 1:nrow(hesbcn$adjacency_matrix)) {
        # consider all leaves
        if(all(hesbcn$adjacency_matrix[i,]==0)) {
            leaf = rownames(hesbcn$adjacency_matrix)[i] # consider each leaf
            paths = all_simple_paths(graph=graph,from="Root",to=leaf,mode="out")
            paths_nodes = list()
            transitions = list()
            lambdas = list()
            for(j in 1:length(paths)) {
                curr_path = names(paths[[j]]) # consider each simple path
                number_path = number_path + 1
                curr_nodes = NULL
                curr_transitions = NULL
                curr_lambdas = NULL
                curr_path_predictability = 1
                for(k in 2:length(curr_path)) { # consider each node visited in this path
                    curr_node = curr_path[k]
                    curr_nodes = c(curr_nodes,curr_node)
                    curr_in_lambda = sum(hesbcn$lambdas_matrix[,curr_node])
                    curr_out_lambda = sum(hesbcn$lambdas_matrix[curr_node,])
                    if(k==length(curr_path)) {
                        curr_out_lambda = 1
                    }
                    curr_transitions = c(curr_transitions,(curr_in_lambda/curr_out_lambda))
                    curr_lambdas = c(curr_lambdas,curr_in_lambda)
                    curr_path_predictability = curr_path_predictability * (curr_in_lambda/curr_out_lambda)
                }
                # save current results
                paths_nodes[[j]] = curr_nodes
                transitions[[j]] = curr_transitions
                lambdas[[j]] = curr_lambdas
                paths_predictability = c(paths_predictability,curr_path_predictability)
            }
            results[[paste0("Evolution_",(length(results)+1))]] = list(leaf_node=leaf,simple_paths=paths_nodes,transitions_vector=transitions,lambdas=lambdas)
        }
    }
    paths_predictability = paths_predictability / sum(paths_predictability) # normalize path predictabilities to sum to 1
    normalized_paths_predictability = paths_predictability

    # compute graph metrics
    graph_predictability = 0.0
    for(path in normalized_paths_predictability) {
        graph_predictability = graph_predictability + (path*log(path))
    }
    graph_entropy = - graph_predictability
    max_entropy = - log((1/number_path))
    graph_predictability = 1 - (graph_entropy/max_entropy)

    # return results
    predictability = list(statistics=results,paths_probabilities=normalized_paths_predictability,graph_entropy=graph_entropy,number_path=number_path,graph_predictability=graph_predictability)
    return(predictability)

}

# enumerate valid genotypes given the HESBCN
"valid.genotypes" <- function( hesbcn, predictability ) {

    # consider all valid genotypes
    valid_genotypes = NULL
    template_genotype = rep("*",(nrow(hesbcn$adjacency_matrix)-1))
    names(template_genotype) = rownames(hesbcn$adjacency_matrix)[2:nrow(hesbcn$adjacency_matrix)]
    for(i in 1:length(predictability$statistics)) { # consider each independent evolution
        for(j in 1:length(predictability$statistics[[i]][["simple_paths"]])) { # consider each path within the current evolution
            curr_genotype = template_genotype # start from empty/template genotype
            path_time = 0.00
            for(k in 1:length(predictability$statistics[[i]][["simple_paths"]][[j]])) {
                curr_gene = predictability$statistics[[i]][["simple_paths"]][[j]][[k]]
                curr_genotype[curr_gene] = 1
                if(hesbcn$parent_set[curr_gene]=="AND") {
                    curr_parents_extra_nodes = names(which(hesbcn$adjacency_matrix[,curr_gene]==1)) # consider all parents of current node
                    curr_genotype[curr_parents_extra_nodes] = 1 # AND --> all parents are 1
                }
                if(hesbcn$parent_set[curr_gene]=="XOR") {
                    curr_parents_extra_nodes = names(which(hesbcn$adjacency_matrix[,curr_gene]==1)) # consider all parents of current node
                    curr_parents_extra_nodes = curr_parents_extra_nodes[which(!curr_parents_extra_nodes==predictability$statistics[[i]][["simple_paths"]][[j]][(k-1)])]
                    curr_genotype[curr_parents_extra_nodes] = 0 # XOR --> only one parent is 1
                }
                path_time = path_time + (1/predictability$statistics[[i]][["lambdas"]][[j]][[k]])
                valid_genotypes = rbind(valid_genotypes,c(curr_genotype,path_time))
            }
        }
    }
    rownames(valid_genotypes) = paste0("Genotype_",1:nrow(valid_genotypes))
    colnames(valid_genotypes)[ncol(valid_genotypes)] = "Progression_Time"
    
    # return results
    return(valid_genotypes)

}

# estimate patients' corrected genotypes from the HESBCN
"corrected.genotypes" <- function( genotypes, patients, epsilon ) {

    # estimate patients' corrected genotypes
    corrected_genotypes = NULL
    for(i in 1:nrow(patients)) {
        # compute likelihood of patients' attachment to corrected genotypes
        patients_attachments_likelihoods = NULL
        curr_patient = as.character(patients[i,])
        for(j in 1:nrow(genotypes)) {
            curr_genotype = as.character(genotypes[j,1:(ncol(genotypes)-1)])
            consistent = length(which(curr_patient[which(curr_genotype!="*")]==curr_genotype[which(curr_genotype!="*")]))
            inconsistent = length(which(curr_patient[which(curr_genotype!="*")]!=curr_genotype[which(curr_genotype!="*")]))
            wild_card = which(curr_genotype=="*")
            wild_card_p = 1
            if(length(wild_card)>0) {
                for(k in wild_card) {
                    if(as.numeric(curr_patient[k])==0) {
                        curr_prob = length(which(patients[,k]==0))/nrow(patients)
                    }
                    if(as.numeric(curr_patient[k])==1) {
                        curr_prob = length(which(patients[,k]==1))/nrow(patients)
                    }
                    wild_card_p = wild_card_p * ((curr_prob*(1-epsilon))+((1-curr_prob)*epsilon))
                }
            }
            curr_attachment_likelihood = ((1-epsilon)^consistent) * (epsilon^inconsistent) * wild_card_p
            patients_attachments_likelihoods = c(patients_attachments_likelihoods,curr_attachment_likelihood)
        }
        max_likelihood_attachment = which(patients_attachments_likelihoods==max(patients_attachments_likelihoods))[1]
        curr_corrected_genotype = as.character(genotypes[max_likelihood_attachment,1:(ncol(genotypes)-1)])
        corrected_genotypes = rbind(corrected_genotypes,curr_corrected_genotype)
    }
    rownames(corrected_genotypes) = rownames(patients)
    colnames(corrected_genotypes) = colnames(genotypes)[1:(ncol(genotypes)-1)]
    
    # return results
    return(corrected_genotypes)

}

# estimate patients' progression relative to cancer timeline (evolution)
"patients.progression" <- function( genotypes, patients, epsilon ) {

    # estimate patients' progression time
    patients_progression = NULL
    progression_time = as.numeric(genotypes[,"Progression_Time"])
    for(i in 1:nrow(patients)) {
        # compute likelihood of patients' attachment to corrected genotypes
        patients_attachments_likelihoods = NULL
        curr_patient = as.character(patients[i,])
        for(j in 1:nrow(genotypes)) {
            curr_genotype = as.character(genotypes[j,1:(ncol(genotypes)-1)])
            consistent = length(which(curr_patient[which(curr_genotype!="*")]==curr_genotype[which(curr_genotype!="*")]))
            inconsistent = length(which(curr_patient[which(curr_genotype!="*")]!=curr_genotype[which(curr_genotype!="*")]))
            wild_card = which(curr_genotype=="*")
            wild_card_p = 1
            if(length(wild_card)>0) {
                for(k in wild_card) {
                    if(as.numeric(curr_patient[k])==0) {
                        curr_prob = length(which(patients[,k]==0))/nrow(patients)
                    }
                    if(as.numeric(curr_patient[k])==1) {
                        curr_prob = length(which(patients[,k]==1))/nrow(patients)
                    }
                    wild_card_p = wild_card_p * ((curr_prob*(1-epsilon))+((1-curr_prob)*epsilon))
                }
            }
            curr_attachment_likelihood = ((1-epsilon)^consistent) * (epsilon^inconsistent) * wild_card_p
            patients_attachments_likelihoods = c(patients_attachments_likelihoods,curr_attachment_likelihood)
        }
        patients_attachments_likelihoods = patients_attachments_likelihoods/sum(patients_attachments_likelihoods) # normalize likelihoods of attachments
        # estimate progression time given attachments' likelihood
        curr_progression_time = 0.00
        for(k in 1:length(patients_attachments_likelihoods)) {
            curr_progression_time = curr_progression_time + progression_time[k]*patients_attachments_likelihoods[k]
        }
        patients_progression = rbind(patients_progression,curr_progression_time)
    }
    rownames(patients_progression) = rownames(patients)
    colnames(patients_progression) = "Progression_Time"
    
    # return results
    return(patients_progression)

}
