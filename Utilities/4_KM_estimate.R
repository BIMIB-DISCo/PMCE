# load required libraries
library("glmnet")
library("pracma")
library("survival")
library("survminer")

# load data
load(file="final_results/survival_analysis.RData")

# perform analysis
for(i in names(survival_analysis[["OS"]])) {
    
    # OS
    os_data = survival_analysis[["OS"]][[i]][["data"]]
    os_fit = survival_analysis[["OS"]][[i]][["fit"]]
    beta_temp = matrix(NA,nrow=dim(os_data)[1],ncol=1)

    Coefficients <- as.numeric(coef(os_fit,s=os_fit$lambda.min))[2:length(as.numeric(coef(os_fit,s=os_fit$lambda.min)))]
    for(k in 1:nrow(os_data)) {
        x = Coefficients
        y = os_data[k,3:length(os_data[k,])]
        beta_temp[k] = dot(x,as.numeric(y))
    }
    beta = as.numeric(beta_temp)
    os_data$beta = beta
    os_data$cluster = "High risk"
    os_data$cluster[os_data$beta<=0] = "Low risk"
    os_data$cluster[os_data$beta==0] = "Medium risk"
    pdf(paste("plots_km/KM_",i,"_OS.pdf",sep=""),width=6,height=6)
    print(ggsurvplot(survfit(Surv(as.numeric(os_data$Times),as.numeric(os_data$Status))~cluster,data=os_data), 
        xlab = "years", 
        ylab = "Overall survival probability OS",mark.time=TRUE,pval=TRUE,ggtheme=theme_bw(), 
        title = i, 
        font.main = 18, 
        font.x = 18, 
        font.y = 18,
        font.caption = 18, 
        font.legend = 18, 
        font.tickslab = 18))
    dev.off()
}
