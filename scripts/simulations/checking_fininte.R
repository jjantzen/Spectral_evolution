


str(bm_data_sim)

i <- 1

bm_data_sim[[i]]

list_of_tests <- as.data.frame(matrix(ncol = 3, nrow = 100))
colnames(list_of_tests) <- c("nrep", "finite", "missing")


for (i in 1:length(bm_data_sim)){
  list_of_tests$nrep[i] <- i
  list_values_inf <- c()
  list_values_mis <- c()
  for (j in 1:nrow(bm_data_sim[[i]]$trait)){
    is_it_finite <- is.finite(bm_data_sim[[i]]$trait[j,])
    if (FALSE %in% is_it_finite){
      list_values_inf <- c(list_values_inf, FALSE)
    } else if (TRUE %in% is_it_finite){
      list_values_inf <- c(list_values_inf, TRUE)
    }
    is_it_missing <- is.na(bm_data_sim[[i]]$trait[j,])
    if (TRUE %in% is_it_missing) {
      list_values_mis <- c(list_values_mis, TRUE)
    } else if (FALSE %in% is_it_missing){
      list_values_mis <- c(list_values_mis, FALSE)
    }
  list_of_tests$finite[i] <- unique(list_values_inf)
  #list_of_tests$finite[i] <- "FALSE"
  list_of_tests$missing[i] <- unique(list_values_mis)
  #list_of_tests$missing[i] <- "FALSE"
  }
  
  #unique(is_it_finite[i]) #want true
  #unique(is_it_missing[i]) #want false
  #list_of_tests$nrep[i] <- i
  #list_of_tests$finite[i] <- list(unique(is_it_finite[2,]))
  #list_of_tests$missing[i] <- list(unique(is_it_missing))
}

j
i

list_of_tests

#all finite and no missing so not sure how it returns that error 

#do I need all positive numbers?

list_of_negs <- as.data.frame(matrix(ncol = 2, nrow = 100))
colnames(list_of_negs) <- c("nrep", "neg")

for (i in 1:length(bm_data_sim)){
  list_of_negs$nrep[i] <- i
  list_values_neg <- c()
  for (j in 1:nrow(bm_data_sim[[i]]$trait)){
    is_it_negative <- bm_data_sim[[i]]$trait[j,] < 0 
    if (TRUE %in% is_it_negative){
      list_values_neg <- c(list_values_neg, TRUE)
    } else {
      list_values_neg <- c(list_values_neg, FALSE)
    }
  }
  #unique(is_it_negative) #want true
  list_of_negs$neg[i] <- list(unique(list_values_neg))
}

list_of_negs


str(is_it_negative)
bm_data_sim[[2]]$trait[1:10]
i
is_it_negative[1:10]
list_of_negs
