# Functions ------------------------

# This function was created to selected to minimal sampling area observed in each mesophotic zone
# Function created to selected the minimal sampling area 
choosing_transec <- function(DATA, min.samp.area=msa){
  resu <- unlist (lapply (split (DATA, DATA$ID_Locality_Mesophotic_Zone), function(x){
    x   <- droplevels(x)
    y   <- split (x,x$Mesophotic2)
    tab <- rep(1:length(y),1000)
    transec <- numeric();  Area=0; i=1
    while (Area<msa){ ## Define the limite for each Mesophotic zone
      pool <- y[[ tab[i] ]][,"transect_id"]; pool
      pool <- pool[!pool%in%transec];pool
      if (length(pool)==0)  {i <- i+1} else{   
        transec[i] <- sample(pool,1,replace=F)
        rm(pool)
        i <- i+1
        Area <- sum (x$area [x$transect_id%in%transec]);
      }
    }
    transec
  }))
  resu=resu[!is.na(resu)]
  names(resu)=NULL
  return(resu)
}

# This function was created to estimated the beta diversity in each mesophothic zone 

beta_diversity <- function (x) {
  # Transform abundance data in presence/absence data
  x [x>0] <- 1
  
  # Estimate beta diversity per each location 
  # BERMUDA -----------------
  Bermuda <- beta.pair (x [1:2, 1:ncol(x)], index.family = "jaccard")
  Beta_Bermuda <- matrix (NA, ncol = 3, nrow=3)
  
  for (i in 1:3) {
    Beta_Bermuda [i,] <- array(Bermuda[[i]])
  } 
  
  rm (Bermuda)
  Beta_Bermuda <- data.frame (Beta_Bermuda)
  Beta_Bermuda$Location  <- "Bermuda"
  Beta_Bermuda$Component <- c("beta.jtu", "beta.jne", "beta.jac")
  
  # CURACAO -----------------
  Curacao <- beta.pair (x [4:6, 1:ncol(x)], index.family = "jaccard")
  Beta_Curacao <- matrix (NA, ncol = 3, nrow=3)
  
  for (i in 1:3) {
    Beta_Curacao [i,] <- array(Curacao[[i]])
  } 
  
  rm (Curacao)
  Beta_Curacao <- data.frame (Beta_Curacao)
  Beta_Curacao$Location  <- "Curacao"
  Beta_Curacao$Component <- c("beta.jtu", "beta.jne", "beta.jac")
  
  # NORONHA -----------------
  Noronha <- beta.pair (x [7:9, 1:ncol(x)], index.family = "jaccard")
  Beta_Noronha <- matrix (NA, ncol = 3, nrow=3)
  
  for (i in 1:3) {
    Beta_Noronha [i,] <- array(Noronha[[i]])
  } 
  
  rm (Noronha)
  Beta_Noronha <- data.frame (Beta_Noronha)
  Beta_Noronha$Location  <- "Noronha"
  Beta_Noronha$Component <- c("beta.jtu", "beta.jne", "beta.jac")
  
  # GUAM -----------------
  Guam <- beta.pair (x [10:12, 1:ncol(x)], index.family = "jaccard")
  Beta_Guam <- matrix (NA, ncol = 3, nrow=3)
  
  for (i in 1:3) {
    Beta_Guam [i,] <- array(Guam[[i]])
  } 
  
  rm (Guam)
  Beta_Guam <- data.frame (Beta_Guam)
  Beta_Guam$Location  <- "Guam"
  Beta_Guam$Component <- c("beta.jtu", "beta.jne", "beta.jac")
  
  # HAWAII -----------------
  Hawaii <- beta.pair (x [13:15, 1:ncol(x)], index.family = "jaccard")
  Beta_Hawaii <- matrix (NA, ncol = 3, nrow=3)
  
  for (i in 1:3) {
    Beta_Hawaii [i,] <- array(Hawaii[[i]])
  } 
  
  rm (Hawaii)
  Beta_Hawaii <- data.frame (Beta_Hawaii)
  Beta_Hawaii$Location  <- "Hawaii"
  Beta_Hawaii$Component <- c("beta.jtu", "beta.jne", "beta.jac")
  
  # MARSHALL -----------------
  Marshall <- beta.pair (x [16:18, 1:ncol(x)], index.family = "jaccard")
  Beta_Marshall <- matrix (NA, ncol = 3, nrow=3)
  
  for (i in 1:3) {
    Beta_Marshall [i,] <- array(Marshall[[i]])
  } 
  
  rm (Marshall)
  Beta_Marshall <- data.frame (Beta_Marshall)
  Beta_Marshall$Location  <- "Marshall"
  Beta_Marshall$Component <- c("beta.jtu", "beta.jne", "beta.jac")
  
  # MOOREA -----------------
  Moorea <- beta.pair (x [19:21, 1:ncol(x)], index.family = "jaccard")
  Beta_Moorea <- matrix (NA, ncol = 3, nrow=3)
  
  for (i in 1:3) {
    Beta_Moorea [i,] <- array(Moorea[[i]])
  } 
  
  rm (Moorea)
  Beta_Moorea <- data.frame (Beta_Moorea)
  Beta_Moorea$Location  <- "Moorea"
  Beta_Moorea$Component <- c("beta.jtu", "beta.jne", "beta.jac")
  
  # PALAU -----------------
  Palau <- beta.pair (x [22:24, 1:ncol(x)], index.family = "jaccard")
  Beta_Palau <- matrix (NA, ncol = 3, nrow=3)
  
  for (i in 1:3) {
    Beta_Palau [i,] <- array(Palau[[i]])
  } 
  
  rm (Palau)
  Beta_Palau <- data.frame (Beta_Palau)
  Beta_Palau$Location  <- "Palau"
  Beta_Palau$Component <- c("beta.jtu", "beta.jne", "beta.jac")
  
  # PHILIPPINES -----------------
  Philippines <- beta.pair (x [25:27, 1:ncol(x)], index.family = "jaccard")
  Beta_Philippines <- matrix (NA, ncol = 3, nrow=3)
  
  for (i in 1:3) {
    Beta_Philippines [i,] <- array(Philippines[[i]])
  } 
  
  rm (Philippines)
  Beta_Philippines <- data.frame (Beta_Philippines)
  Beta_Philippines$Location  <- "Philippines"
  Beta_Philippines$Component <- c("beta.jtu", "beta.jne", "beta.jac")
  
  # POHNPEI -----------------
  Pohnpei <- beta.pair (x [28:30, 1:ncol(x)], index.family = "jaccard")
  Beta_Pohnpei <- matrix (NA, ncol = 3, nrow=3)
  
  for (i in 1:3) {
    Beta_Pohnpei [i,] <- array(Pohnpei[[i]])
  } 
  
  rm (Pohnpei)
  Beta_Pohnpei <- data.frame (Beta_Pohnpei)
  Beta_Pohnpei$Location  <- "Pohnpei"
  Beta_Pohnpei$Component <- c("beta.jtu", "beta.jne", "beta.jac")
  
  # ST PAUL ROCKS -----------------
  SPaulRocks <- beta.pair (x [31:33, 1:ncol(x)], index.family = "jaccard")
  Beta_SPaulRocks <- matrix (NA, ncol = 3, nrow=3)
  
  for (i in 1:3) {
    Beta_SPaulRocks [i,] <- array(SPaulRocks[[i]])
  } 
  
  rm (SPaulRocks)
  Beta_SPaulRocks <- data.frame (Beta_Pohnpei)
  Beta_SPaulRocks$Location  <- "SPaulRocks"
  Beta_SPaulRocks$Component <- c("beta.jtu", "beta.jne", "beta.jac")
  
  # TAHITI -----------------
  Tahiti <- beta.pair (x [34:36, 1:ncol(x)], index.family = "jaccard")
  Beta_Tahiti <- matrix (NA, ncol = 3, nrow=3)
  
  for (i in 1:3) {
    Beta_Tahiti [i,] <- array(Tahiti[[i]])
  } 
  
  rm (Tahiti)
  Beta_Tahiti <- data.frame (Beta_Pohnpei)
  Beta_Tahiti$Location  <- "Tahiti"
  Beta_Tahiti$Component <- c("beta.jtu", "beta.jne", "beta.jac")
  
  Global.Beta <- rbind (Beta_Bermuda, Beta_Curacao, Beta_Noronha, Beta_Guam, Beta_Hawaii, 
                        Beta_Marshall, Beta_Moorea, Beta_Palau, Beta_Philippines, Beta_Pohnpei,
                        Beta_SPaulRocks, Beta_Tahiti)
  colnames(Global.Beta)[1:3] <- c("Shallow_vs_Deeper","Shallow_vs_Upper","Upper_vs_Deeper")
  Global.Beta <- melt (Global.Beta, id=c("Location","Component"))
  Global.Beta <- cast (Location + variable ~ Component, value = "value", data=Global.Beta)
  colnames (Global.Beta) [2] <- "Zone.Comparison"
  Global.Beta
}
