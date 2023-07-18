source('C:\\MDICC-main\\NetworkFusion.R')
library(Rcpp)
library(parallel)
library(Matrix)
library(survival)
library("survminer")
setwd("C:\\MDICC-main\\")
# system("R CMD SHLIB projsplx_R.c")
dyn.load("projsplx_R.dll")

# connect anaconda environment
library(reticulate)
use_virtualenv("C:\\MDICC-main\\mdicc_env")
# the path of python 
Sys.setenv(RETICULATE_PYTHON="C:/MDICC-main/mdicc_env/Scripts/python.exe") 
use_python("C:/MDICC-main/mdicc_env/Scripts/python.exe", required = TRUE)
py_config()
py_available()
source_python("LocalAffinityMatrix.py")
source_python("score.py")
source_python("label.py")

# read data
setwd("C:\\MDICC-main\\data_SKCM\\data\\skcm") # data path
list <- list.files()
data <- data.frame()
data1 <- list()
X <- list()
for(i in list){
  path <- i
  data <- read.csv(file = path, header = TRUE)
  #rownames(data) <- data$X
  data <- data[-1]
  data11 <- as.matrix(data)
  data1[[i]] = scale(data11, center=TRUE, scale=TRUE) 
  data2 = t(data1[[i]])
  d1 = dist(data2)
  d1 = as.matrix(d1)
  X[[i]] <- d1
}

all_i = list()
# parameter setting
k1 = 18 # 18 the neighbor of affinity matrix
k2 = 42 # 
k3 = 2  # number of cluster
c  = 3  # c = k3(c>2) or c = 3(c=2)

# calculate affinity matrix
aff_MDICC = list()

for(i in 1:3){
  a = as.matrix(X[[i]])
  xxx = testaff(a,k1)
  aff_MDICC[[i]] = xxx
}

# Get proposed affinity matrices based on CCA
survival_path <- "C:\\MDICC-main\\data_SKCM\\data\\survival.csv"
aff_files <- c("A12_SKCM.csv", "A21_SKCM.csv", "A20_SKCM.csv", "A02_SKCM.csv", "A01_SKCM.csv", "A10_SKCM.csv")
aff_proposed <- vector("list", length(aff_files))
# Read affinity files
for (i in 1:length(aff_files)) {
  file_path <- paste("C:\\MDICC-main\\data_SKCM\\data\\", aff_files[i], sep = "")
  aff_proposed[[i]] <- as.matrix(read.csv(file = file_path))[, -1]
}

# Read survival data
surv_data = read.csv(file = survival_path)

# Convert variables to appropriate types
surv_data$Survival = as.numeric(surv_data$Survival)
surv_data$Death = as.logical(surv_data$Death)

# Vector to store the -10log values
logrank_values = numeric()

# List to store all p-values
best_p = list(0,0,0)

# Initialize values for iterating
score_list = list()
#cluster_range = 2:5
cluster_range = 3:5
fused_list = list()
for(cluster in cluster_range){
  k3 = cluster
  c  = cluster
  temp_best_p = 0
  if(cluster == 2){
    c = 3
    
    # Set path for class/clustering data
    label_path <- "C:\\MDICC-main\\data_BRCA\\data\\label.csv"
    class_name <- 'label1'
    
    # Get the fused matrix from MDICC
    MDICC_fusion = MDICC(aff_MDICC,c = c,k = 42, fusion_number = 1)
    MDICC_S = as.matrix(MDICC_fusion)
    
    # Compute ARI and NMI 
    MDICC_score = MDICCscore(MDICC_S,k3,label_path,class_name)
    
    # Fuse the proposed affinity matrices
    proposed_fusion = MDICC(aff_proposed,c = c,k = 2,fusion_number = 2)
    proposed_S = as.matrix(proposed_fusion)
    proposed_score = MDICCscore(proposed_S,k3,label_path,class_name)
    
    # Scaling the fused matrices
    final_aff = list()
    final_aff[[1]] = MDICC_S
    final_aff[[2]] = proposed_S
    scaled_final_aff <- list()  # Create an empty list to store the scaled matrices
    for (i in 1:2) {
      scaled_final_aff[[i]] <- testaff(final_aff[[i]], k1)
    }
    
    # Final fusion step, determine k2 in an exhaustive way
    for(i in 2:100){
      final_fuse = MDICC(scaled_final_aff,c = c,k = i, fusion_number = 3)
      fuse_S = as.matrix(final_fuse)
      temp_score = MDICCscore(fuse_S,k3,label_path,class_name)
      score_list[[i]] = temp_score
    }
  }
  else{
    
    # Get the fused matrix from MDICC
    MDICC_fusion = MDICC(aff_MDICC,c = c,k = 42, fusion_number = 1)
    MDICC_S = as.matrix(MDICC_fusion)

    # Fuse the proposed affinity matrices
    proposed_fusion = MDICC(aff_proposed,c = c,k = 42,fusion_number = 2)
    proposed_S = as.matrix(proposed_fusion)

    # Scaling the fused matrices
    final_aff = list()
    final_aff[[1]] = MDICC_S
    final_aff[[2]] = proposed_S
    scaled_final_aff <- list()  # Create an empty list to store the scaled matrices
    for (i in 1:2) {
      scaled_final_aff[[i]] <- testaff(final_aff[[i]], k1)
    }
    
    # Final fusion step, determine k2 in an exhaustive way
    for(i in 2:100){
      final_fuse = MDICC(scaled_final_aff,c = c,k = i, fusion_number = 3)
      fuse_S = as.matrix(final_fuse)
        
      # Get the labels from clustering the fused affinity matrix S
      clustering = MDICClabel(fuse_S, cluster)
        
      # Create and initialize the needs for computing the cox model
      surv_obj = Surv(surv_data$Survival, surv_data$Death)
      clustering_labels = as.factor(clustering)
      data_cluster = data.frame(Survival = surv_obj, Clustering = clustering_labels)
      
      # Perform the log-rank test
      logrank_test = survdiff(surv_obj ~ clustering, data = data_cluster)
      
      # Calculate the -10log value
      p_value = -log10(logrank_test$pvalue)
      
      # If the p value is the biggest so far, save it
      if(p_value > temp_best_p){
        temp_best_p = p_value
        best_p[[(cluster - 2)]] = p_value
      }
    }
  }
}
print(best_p)
