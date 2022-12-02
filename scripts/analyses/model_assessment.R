#Model assessment 

library(mvMORPH)
library(cowplot)


#load models

#example one
example_model <- readRDS("./analysis/myc_models/transformed/myc_model_trans_BM_iteration8.rds")
#example_model <- lp_model_BM_iteration12
#make format for plot
attach(mtcars)

jpeg("./output/model_assessment/assessment_myc_trans_BM_8_best.jpg", width = 8, height = 6, res = 500, units = "in")
par(mfrow=c(2,2))

#plot residuals vs independent variable
plot(example_model$variables$Y, example_model$residuals)

#Distribution of residuals (normal)
qplot <- qqnorm(example_model$residuals)
qqline(example_model$residuals, col = "red") 

#density plot
plot(density(example_model$residuals))

#Homoskedasticity (residuals vs fitted values)
plot(fitted(example_model), example_model$residuals)
abline(0,0, col = "red")
dev.off()

#get standardized residuals - how?
#par(mfrow = c(2,2))


#Influential cases


#Model stability (look at rare combinations)



#How do parameters relate to min/max of preset bounds
example_model$start_values
example_model$param


#Full-null models

