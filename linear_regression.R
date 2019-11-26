#Use the genes that are significant in linear regression
linear_regression <- function(all_gene){
  p_v <- c()
  for (i in 1:(ncol(all_gene)-2)) {
    fit <- lm(ending~all_gene[,i],data =all_gene )
    a <- summary(fit)
    p_value <- a$coefficients[2,4]
    p_v <- c(p_v,p_value)
  }
  ad <- ifelse(p_v<0.05,1,0)
  sig <-cbind(colnames(gene),p_v,ad)
  sig_gene <- subset(sig,ad==1)
  colnames(sig_gene) <- c('gene','p','exist')
  return(sig_gene)
}
