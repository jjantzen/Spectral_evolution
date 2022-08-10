model

covR <- model$sigma$Pinv
if (mode == "corr") 
  covR <- cov2cor(covR)
eig <- eigen(covR)
values <- eig$values
U <- eig$vectors
resids <- model$residuals
S <- resids %*% U

str(S)
axes <- c(1, 2)

tot <- sum(values)
valX <- round(values[axes[1]] * 100/tot, digits = 2)
valY <- round(values[axes[2]] * 100/tot, digits = 2)
xlabel <- paste("PC", axes[1], " (", valX, " %)", sep = "")
ylabel <- paste("PC", axes[2], " (", valY, " %)", sep = "")


str(S[,axes[2]])
df_S <- as.data.frame(S[,axes[1:2]])
rownames(df_S)
colnames(df_S)
df_S$V1
ggplot(data = df_S, aes(x = V1, y = V2))+
  geom_point()
class(S)

rownames(df_S)
