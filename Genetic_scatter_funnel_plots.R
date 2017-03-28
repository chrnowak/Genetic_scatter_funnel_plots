################## ################## ################## ################## ################## ################## 
################## March 2017: Chris, genetic scatter and funnel plots      ################## ################## 
################## ################## ################## ################## ################## ################## 

## Proviso: instrument strength defined "as genetic association with exposure divided by standard error of 
## genetic association with outcome"; see: Bowden et al. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4849733/ 
## Funnel plots as instrument strength against "genetic association with O divided by genetic association with E)"

#### data set with columns
## rsid       markername  
## beta_out   beta snp-outcome
## se_out     se snp-outcome
## p_out      (not required, p snp-outcome)
## for option to plot only SNP above p-threshold for outcome association
## beta_exp   beta snp-exp
## se_exp    se snp-exp

library(ggplot2)
library(gridExtra)
library(grid)

## simulate dataset of 20 SNPs, 10 with GWAS-sign association with outcome
set.seed(123)
test_set <- data.frame(paste(rep("rs", 20), 1:20, sep=""), runif(20, 0.1, 1), runif(20, 0.01, 0.1),
                       c(runif(10, 1e-20, 5e-08), runif(10, 1e-04, 0.05)), runif(20, -1, 1), runif(20, 0.01, 0.1))
names(test_set) <- c("rsid", "beta_out", "se_out",  "p_out", "beta_exp" , "se_exp" )

############################
## scatter plot function  ##
############################

plot_scatter <- function(data, title , x_label, y_label, p_cutoff){
  if(p_cutoff != F){
    data <- data[which(data$p_out < p_cutoff), ]
  }
  df <- data.frame(x = data$beta_exp, xmin = data$beta_exp - (1.96*data$se_exp), xmax = data$beta_exp + (1.96*data$se_exp),
                   y = data$beta_out, ymin = data$beta_out - (1.96*data$se_out), ymax = data$beta_out + (1.96*data$se_out),
                   snp_name = data$rsid)
  p <- ggplot(data = df,aes(x = x,y = y)) + geom_point() + 
    geom_errorbar(aes(ymin = ymin, ymax = ymax)) + 
    geom_errorbarh(aes(xmin = xmin, xmax = xmax)) +
    geom_smooth(method = 'lm', formula = y~x, se = F,col = "red", fullrange = T) +
    geom_hline(aes(yintercept = 0)) +
    geom_vline(aes(xintercept = 0)) +
    ggtitle(title) + 
    theme(axis.text = element_text(size=10), axis.title = element_text(size = 12)) +
    labs(y = y_label, x = x_label) 
  return(p)
}

############################
## funnel plot function   ##
############################

plot_funnel <- function(data, title, vertical = F , egger = F, p_cutoff = F){
  if(p_cutoff != F){
    data <- data[which(data$p_out < p_cutoff),]
  }
  if(vertical == F){
    vertical <- 0
    } 
  if(egger == F){
    egger <- 0
  }
  dfunnel <- data.frame(y = data$beta_exp / data$se_out, x = data$beta_out / data$beta_exp,
                        xmin = (data$beta_out / data$beta_exp) - (1.96*(data$se_out / data$beta_exp)), 
                        xmax = (data$beta_out / data$beta_exp) + (1.96*(data$se_out / data$beta_exp)))
  p2 <- ggplot(data = dfunnel,aes(x = x, y = abs(y))) + 
    geom_point(size = 1) + 
    geom_errorbarh(aes(xmin = xmin, xmax = xmax)) +
    geom_vline(xintercept = 0) + 
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = vertical, col = "blue") +
    geom_vline(xintercept = egger, col = "red") +
    ggtitle(title) + 
    labs(y = "IV strength",x = "IV estimate")  +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12)) 
  return(p2)
}

### Plot! ###
plot_scatter(test_set, title="exp effect on outcome","snp-exposure","snp-outcome", p_cutoff=F)
plot_scatter(test_set, title="exp effect on outcome","snp-exposure","snp-outcome", p_cutoff= 5e-08)# plot only GWAS-sign snps with outcome

plot_funnel(test_set , "funnel")
plot_funnel(test_set , "funnel only sign", vertical = -1, egger = -2, p_cutoff = 5e-08) 
  ## e.g., vertical = causal estimate in 2SLS, egger = egger estimate

p1 <- plot_scatter(test_set, title="exp effect on outcome","snp-exposure","snp-outcome", p_cutoff=F)
p2 <- plot_funnel(test_set , "funnel")
grid.arrange(p1, p2, ncol = 2)

### PDF export ###
# pdf(w = 12, h = 5.5, paste("exposure", "outcome", "pdf", sep = "."))
print(grid.arrange(p1, p2, ncol = 2))
# dev.off()


#########################################################################################################################