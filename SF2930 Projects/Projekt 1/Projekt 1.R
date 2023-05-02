#project1
getwd() #used to se where my document are stored
library(MPV)
table <- read.csv('bodyfatmen.csv')

#View(table)
#skapar y(density) = full model

table.model <- lm(density ~ age + weight + height + neck + chest + 
                    abdomen + hip + thigh + knee + ankle + biceps + 
                    + forearm + wrist,data = table)
summary(table.model)
#stars on neck, abdomen, forearm and wrist
anova(table.model)

# Normal probability plot of residuals
qqnorm(rstudent(table.model))
qqline(rstudent(table.model))
#Good normal probability plot, 

# Residuals vs. fitted values
plot(table.model$fitted.values, table.model$residuals)
### Horizontal band, so satisfactory distribution

# Residuals vs. regressor variables
par(mfrow = c(3, 3))
plot(table$age, table.model$residuals)
plot(table$weight, table.model$residuals)
plot(table$height, table.model$residuals)
plot(table$neck, table.model$residuals)

plot(table$chest, table.model$residuals)
plot(table$abdomen, table.model$residuals)
plot(table$hip, table.model$residuals)
par(mfrow = c(3, 3))
plot(table$thigh, table.model$residuals)
plot(table$knee, table.model$residuals)
plot(table$ankle, table.model$residuals)
plot(table$biceps, table.model$residuals)
plot(table$forearm, table.model$residuals)

plot(table$wrist, table.model$residuals)
### We can't say much, except for??

# Partial residual plots
library(car)  # avPlots
avPlots(table.model)
### can't see any linear relations in knee, chest, height, use this with caution

# Studentized residuals
library(MASS)  # studres
par(mfrow = c(1, 2))
plot(studres(table.model))
plot(rstudent(table.model))
#eventuellt avvikande punkter i hörnen



#from exc 4
#Q1
par(mfrow = c(2, 3))
plot(table.model, which = 1:5)
#39  avviker 


#Q2
n <- nrow(table)
p <- ncol(table)
leverage.cutoff <- 2 * p / n  # MPV p. 213
cooks.cutoff <- qf(0.5, p, n - p, lower.tail = FALSE)  # MPV p. 215
dfbetas.cutoff <- 2 / sqrt(n)  # MPV p. 218
dffits.cutoff <- 2 * sqrt(p / n)  # MPV p. 219
studres.cutoff <- qt(0.05 / 2, n - p, lower.tail = FALSE)  # MPV p. 135

### leverage points
model.hat <- hatvalues(table.model)
model.hat[model.hat > leverage.cutoff]
#influetial points candiates: new ones

### Cook's distance
table.cooks <- cooks.distance(table.model)
table.cooks[table.cooks > cooks.cutoff]
#non or wrong


### DFBETAS NOTHING with newfullmodel
model.dfbetas <- dfbetas(table.model)
model.dfbetas[abs(model.dfbetas[, 1]) > dfbetas.cutoff, 1]  # beta_0
#28 31 32 39 52 79 105 108 136 167 178 200 217 227
model.dfbetas[abs(model.dfbetas[, 2]) > dfbetas.cutoff, 2]  # beta_1
# 3 22 28 32 39 60 65 73 83 167 168 171 200 220
model.dfbetas[abs(model.dfbetas[, 3]) > dfbetas.cutoff, 3]  # etc
# 28 31 32 33 39 52 77 79 105 108 136 178 200 203 217 227
model.dfbetas[abs(model.dfbetas[, 4]) > dfbetas.cutoff, 4]
# 25 31 39 52 105 108 124 136 152 167 168 178 200 203 212 217 219 227 236
model.dfbetas[abs(model.dfbetas[, 5]) > dfbetas.cutoff, 5]
model.dfbetas[abs(model.dfbetas[, 6]) > dfbetas.cutoff, 6]
model.dfbetas[abs(model.dfbetas[, 7]) > dfbetas.cutoff, 7]
model.dfbetas[abs(model.dfbetas[, 8]) > dfbetas.cutoff, 8]
model.dfbetas[abs(model.dfbetas[, 9]) > dfbetas.cutoff, 9]
model.dfbetas[abs(model.dfbetas[, 10]) > dfbetas.cutoff, 10]
model.dfbetas[abs(model.dfbetas[, 11]) > dfbetas.cutoff, 11]
model.dfbetas[abs(model.dfbetas[, 12]) > dfbetas.cutoff, 12]
model.dfbetas[abs(model.dfbetas[, 13]) > dfbetas.cutoff, 13]
model.dfbetas[abs(model.dfbetas[, 14]) > dfbetas.cutoff, 14]

### DFFITS
model.dffits <- dffits(table.model)
model.dffits[abs(model.dffits) > dffits.cutoff]
# 3 31 39 52 78 79 83 124 136 171 178 293 212 217 227 246

# Built-in functions
par(mfrow = c(2, 3))
plot(table.model, 1:5)
# 39 83 217 cook's
#shall we remove these or not

### Removing outliers
newtable <- table[-c(39), ]
#View(mynewtable)
newfullmodel <- lm(density ~ . , data = newtable)
summary(newfullmodel) # R^2 = 0.73
summary(table.model)

#########transaformation for response, not motivated to use box cox since constant variance.
##box cox on response
library(MASS)

#fit linear regression model
modelbox <- newfullmodel

#find optimal lambda for Box-Cox transformation 
bc <- boxcox(density ~ . , data = newtable)
(lambda <- bc$x[which.max(bc$y)])

[1] -2

#fit new linear regression model using the Box-Cox transformation
new_model <- lm(((density^-2)/lambda) ~ ., data = newtable)
#summary(new_model)
#summary(newfullmodel)
##########################

#Q3
X <- model.matrix(newfullmodel)
WtW <- cor(X[,-1])
WtW
#many strong correlations - chest vs weight, weight vs abdomen, wrist vs neck etc


#Variance inflation factor and condition number
library(car)  # vif
vif(table.model)  # cutoff value is 10 (MPV p. 118)
#weight, chest, abdomen, hip are all above 10
solve(WtW)  # = (W^T W)^{-1}. Compare diagonal values with VIFs

### Condition number, cutoff value is 100 (MPV p.298)
model.eigen <- eigen(WtW)
max(model.eigen$values) / min(model.eigen$values)
### There is multicollinearity, not servesince condition number < 1000



####Backward elimination 
backwardmodel <- lm(density ~ .,data = newtable)
summary(backwardmodel)
step(backwardmodel, direction ="backward")
# y = B0 + B(age, height, neck,chest, abdomen, hip, thigh, forearm, wrist)
#############

###forward selection 
forwardmodel <- lm(density ~ 1, data = newtable)
head(table)
summary(forwardmodel)
step(forwardmodel, direction = "forward", scope = formula(backwardmodel))
#here backwardmodel = full model
# y = B0 + B(weight, abdomen, neck, biceps,  forearm, wrist)
##########

#### Stepwise elimination 
stepwise <- lm(density ~ 1, data = newtable)
head(table)
summary(stepwise)
step(stepwise, direction = "both", scope = formula(backwardmodel))
##########

#### Mallow's Cp #ruunning the one with lowest ACI
library(olsrr)
fullmodel1 <-  lm(density ~ .,data = newtable)
fullmodel2 <- lm(density ~ .,data = table)

#backwards choice, without potential outlier
model1 <- lm(density ~ age+ height+ neck + abdomen+ hip+
               thigh+ forearm + wrist, data = newtable)

#forwards and stepwise choice without potential outlier
model2 <-lm(density ~abdomen + weight + wrist+ biceps, data = newtable)


model3 <-   lm(density ~ age+ height+ neck + abdomen+ hip+
                 thigh+ forearm + wrist, data = table)

model4 <- lm(density ~abdomen + weight + wrist+ biceps, data = table)


ols_mallows_cp(model1, fullmodel1) # mallows CPs = 21.51341
ols_mallows_cp(model2, fullmodel1) # mallows CPs = 22.97961
ols_mallows_cp(model3, fullmodel2) # mallows CPs = 8.33443
ols_mallows_cp(model4, fullmodel2) # mallows CPs = 11.79142
summary(model1) #adjusted R_squared= 0.7368
summary(model2) #adjusted R_squared= 0.7314
summary(model3) #adjusted R_squared= 0.7317
summary(model4) #adjusted R_squared= 0.7325

#### PRESS ,small value is desired 
PRESS(model1)  #0.02358523
PRESS(model2)  #0.02383791
PRESS(model3)  #0.02495297
PRESS(model4)  #0.02524917

#### All possible regressions
newdata <- c(table)
newdata_2 <- c(newtable)
model_allpossiblewithoutlier <- lm(density ~. , data = newdata)
model_allpossiblewithoutoutlier <- lm(density ~. , data = newdata_2)

allpossible_withoutlier <- ols_step_all_possible(model_allpossiblewithoutlier)
allpossible_without_outlier <- ols_step_all_possible(model_allpossiblewithoutoutlier)



#with outlier lowest cp and highest adjr
allpossible_withoutlier[which.max(allpossible_withoutlier$adjr),]
#age + weight + neck + abdomen + hip + thigh + biceps + forearm + wrist

allpossible_withoutlier[which.min(allpossible_withoutlier$cp),]
#age + weight + neck + abdomen + hip + thigh + forearm + wrist



############# without otuliers!!!!
#without outlier lowest cp and highest adjr
allpossible_without_outlier[which.max(allpossible_withoutlier$adjr),]
#age + height + neck + chest + abdomen + hip + thigh  + forearm + wrist

allpossible_without_outlier[which.min(allpossible_without_outlier$cp),]
#age + weight + neck + abdomen  + thigh + forearm + wrist




########################

#### leaps for all pos regression
library(leaps)
y <- newtable$density
x_table <- newtable[,-1]
x <- x_table
leaps(x, newtable$density)
#####

###### PRESS 
library(stringr)
predictors = allpossible_without_outlier$predictors
number_iterations = length(predictors)
press_vector = c(1:length(predictors), nrows = length(predictors))
press_vector
for (i in 1:number_iterations)
{  
  input = predictors[i]
  input  
  input = gsub(" ", ",", predictors[i])
  
  #  if (commas == 0)
  # {
  #  predictors[i] = paste(predictors[i],last = ",")
  # predictors[i] = gsub(" ", "", predictors[i])
  #}  
  input = paste(input,last = ",")
  predictors[i] = input
  predictors[i]
  predictors[i] = gsub(" ", "", predictors[i])  
  commas = str_count(input, ",")
  commas
  input = predictors[i]
  input
  predictors[i]
  
  
  
  commas_location = str_locate_all(input, ",")
  
  commas_location
  size_matrix = number_iterations*13
  size_matrix
  index_vector = matrix(1:size_matrix, nrow = number_iterations, byrow=TRUE)
  #print(index_vector)
  letter_index = 1  
  for (j in 1:commas)
  {
    #j = 1
    input
    next_commas = str_locate_all(input, ",")[[1]]
    next_commas[j]
    
    test_input = substr(input, letter_index, next_commas[j])
    test_input
    length_test_input = nchar(test_input)
    length_test_input
    search_str = substr(test_input, 1,length_test_input-1)
    search_str
    
    
    length = nchar(test_input)
    index = which(search_str == names(table))
    index
    
    index_vector[i, j] = index
    letter_index = letter_index + length
    
    
  }
  #  index_vector[i,] = c(index_vector[i][1:length(input)])
  
  to_insert = colnames(newtable[index_vector[i,][1:commas]])
  outcome = "density"
  
  formula_to_insert = as.formula(paste(outcome, paste(to_insert, collapse = " + "), sep = " ~ "))
  formula_to_insert
  i
  test_model = lm(formula_to_insert, data = newtable)
  press_vector[i] = PRESS(test_model)
  
}

press_vector

which.min(press_vector)





 

### PRESS for our models with outliers
P_model1 <- lm(density ~ age + weight + neck + abdomen + hip + thigh + biceps + forearm + wrist, data = table)
P_model2 <- lm(density ~ age + weight + neck + abdomen + hip + thigh + forearm + wrist, data = table)
P_model3 <- lm(density ~ age + weight + neck + abdomen + hip + thigh + forearm, data = table)
P_model4 <- lm(density ~ age + weight + neck + abdomen + hip + thigh + biceps + forearm + wrist, data = table)

### PRESS for our models without outliers
P_model5 <- lm(density ~ age + height + neck + chest + abdomen + hip + thigh  + forearm + wrist, data = newtable)
P_model6 <- lm(density ~ age + weight + neck + abdomen  + thigh + forearm + wrist, data = newtable)
P_model7 <- lm(density ~ age + weight + neck + abdomen + hip + thigh + ankle + biceps + forearm + wrist, data = newtable)
P_model8 <- lm(density ~ weight + neck + abdomen + biceps + forearm + wrist, data = newtable)

p1 <- PRESS(P_model1) #0.02457851
p2 <- PRESS(P_model2) #0.02449232
p3 <- PRESS(P_model3) #0.02495787
p4 <- PRESS(P_model4) #0.02457851

p5 <- PRESS(P_model5) #0.02364599
p6 <- PRESS(P_model6) #0.02356691
p7 <- PRESS(P_model7) #0.02391456
p8 <- PRESS(P_model8) #0.02380427



#Ranking: 1. 
# Perfrom a thorugh analysis of the "best" models

# Built-in functions with outlier
par(mfrow = c(2, 3))
plot(P_model1, 1:5)
summary(P_model1)
anova(P_model1)

par(mfrow = c(2, 3))
plot(P_model2, 1:5)
summary(P_model2)
anova(P_model2)

par(mfrow = c(2, 3))
plot(P_model3, 1:5)
summary(P_model3)
anova(P_model3)

par(mfrow = c(2, 3))
plot(P_model4, 1:5)
summary(P_model4)
anova(P_model4)

#without outliers
par(mfrow = c(2, 3))
plot(P_model5, 1:5)
summary(P_model5)
anova(P_model5)

par(mfrow = c(2, 3))
plot(P_model6, 1:5)
summary(P_model6)
anova(P_model6)

par(mfrow = c(2, 3))
plot(P_model7, 1:5)
summary(P_model7)
anova(P_model7)

par(mfrow = c(2, 3))
plot(P_model8, 1:5)
summary(P_model8)
anova(P_model8)



##### traning and test data:
library(ISLR)


# (a) Split data into training and test sets
#set.seed(689)
n <- nrow(newtable)
train <- sample(1:n, n / 2) # split in half
density.train <- newtable[train, ] # 124 obs of 14 variables
density.test <- newtable[-train, ] #123 obs of 14 variables

# (b) Fit least squares linear model on training set
density.model <- lm(density ~ ., data = density.train)
### Report test error
model.pred <- predict(density.model, density.test)
model.MSE <- mean((density.test$density - model.pred) ^ 2)
model.MSE  # MSE = 0.0001080797

# (c) Ridge regression
library(glmnet)  # elastic net
View(newtable)
x.train <- data.matrix(density.train[,-1]) #Everything except the first column(density)
y.train <- data.matrix(density.train[,"density"])
density.ridge <- cv.glmnet(x.train, y.train, alpha=0)  # alpha = 0 means ridge regression
density.ridge$lambda.1se  # one standard deviation away from minimizing lambda = 0.004529917 what does this number mean?
coef(density.ridge)  # coefficients for lambda.1se, What does this number means?

### Report test error
x.test <- data.matrix(density.test[,-1])
ridge.pred <- predict(density.ridge, s = density.ridge$lambda.1se, newx = x.test)
ridge.MSE <- mean((density.test$density - ridge.pred)^2)
ridge.MSE # = 0.0001300618 what does this number mean, this one in larger than MSE

# (d) Lasso
density.lasso <- cv.glmnet(x.train, y.train, alpha=1)  # alpha = 1 means lasso
density.lasso$lambda.1se # =  0.001320495 what does this mean?

### Report test error
lasso.pred <- predict(density.lasso, s =density.lasso$lambda.1se, newx = x.test)
lasso.MSE <- mean((density.test$density - lasso.pred)^2)
lasso.MSE # = 0.0001148993 what does this mean?

### Report coefficient estimates
coefs.lasso <- coef(density.lasso)  # coefficients for lambda.1se
coefs.lasso[rowSums(coefs.lasso) != 0,]  # nonzero coefficients, height, abdomen, ankle & wrist, what does this mean

# (e) Principal component regression
library(pls)  # pcr
density.pcr <- pcr(density ~ ., data = newtable, scale = TRUE, validation = "CV")
coef(density.pcr)
summary(density.pcr)
validationplot(density.pcr, val.type = "MSEP")   # MSE vs. number of components
### Stabilization at m = 5, where is our stabilization?
### Report test error
pcr.pred <- predict(density.pcr, x.test, ncomp = 12)
pcr.MSE <- mean((density.test$density - pcr.pred)^2)
pcr.MSE


# (g) Compare the methods
model.MSE
ridge.MSE
lasso.MSE
pcr.MSE

### Isaac comments, The least-squares model minimizes the MSE by definition. Ridge and lasso
### have comparable predicting power, with PCR falling behind.

#cvmodel1 = age + height + neck + chest + abdomen + hip + thigh + forearm + wrist
#cvmodel2 = age + weight + neck + abdomen + thigh + forearm + wrist

### trying CV 
library(caret)

# specify the cross-validation method
ctrl <- trainControl(method = "cv", number = 5, savePredictions = TRUE)

# fit a regression model and use 5-fold CV to evaluate performance
model1 <- train(density~ age + height + neck + chest + abdomen + hip + thigh + forearm + wrist, data = newtable, method = "lm", trControl = ctrl, metric = "Rsquared")
model2 <-  train(density~ age + weight + neck + abdomen + thigh + forearm + wrist, data = newtable, method = "lm", trControl = ctrl, metric = "Rsquared")
model3 <- train(density~ age + height + neck + chest + abdomen + hip + thigh + biceps + forearm + wrist, data = newtable, method = "lm", trControl = ctrl, metric = "Rsquared")
model4 <- train(density~ age + weight + abdomen + thigh + forearm + wrist, data = newtable, method = "lm", trControl = ctrl, metric = "Rsquared")
model5 <- train(density~ age + weight + neck + abdomen + hip + thigh + forearm + wrist, data = newtable, method = "lm", trControl = ctrl, metric = "Rsquared")

n <- nrow(newtable)
train <- sample(1:n, n / 2)
traningdata <- newtable[train, ]
testdata <- newtable[-train, ]

View(newtable)
library(glmnet)  # elastic net
x.train <- data.matrix(traningdata[,-1])
y.train <- data.matrix(traningdata[,"density"])


# (d) Lasso
lasso <- cv.glmnet(x.train, y.train, alpha=1)  # alpha = 1 means lasso
lasso$lambda.1se

### Report test error
lasso.pred <- predict(lasso, s = lasso$lambda.1se, newx = x.test)
lasso.MSE <- mean((testdata$density - lasso.pred)^2)
lasso.MSE

### Report coefficient estimates
coefs.lasso <- coef(lasso)  # coefficients for lambda.1se
coefs.lasso[rowSums(coefs.lasso) != 0,]  # nonzero coefficients



library(boot)
library(ggplot2)
library(boot)

#define function to calculate fitted regression coefficients
coef_function <- function(formula, data, indices) {
  d <- data[indices,] #allows boot to select sample
  fit <- lm(formula, data=d) #fit regression model
  return(coef(fit)) #return coefficient estimates of model
}

#perform bootstrapping with 2000 replications
reps <- boot(data=newtable, statistic=coef_function, R=2000, formula=density~age + abdomen + height)
reps <- boot(data=newtable, statistic=coef_function, R=2000, formula=density~age + height + neck + chest + abdomen + hip + thigh + forearm + wrist)

#view results of boostrapping
reps





boot(data = newtable, statistic = coef_function, R = 2000, formula = density ~ 
       age + abdomen + height)

plot(reps)

boot.ci(boot.out = reps, type = "bca", index = 1)

library(boot)
library(ggplot2)

r_squared <- function(formula, data, indices) {
  val <- data[indices,] # selecting sample with boot 
  fit <- lm(formula, data=val)
  return(summary(fit)$r.square)
} 
# Performing 1500 replications with boot 
output <- boot(data=newtable, statistic=r_squared, 
               R=2000, formula=density~ age + height + neck + chest + abdomen + hip + thigh + forearm + wrist)
# Plotting the output
output 
plot(output)
# Obtaining a confidence interval of 95%
boot.ci(output, type="bca")






model1 <- lm(density~ age + height + neck + chest + abdomen + hip + thigh + forearm + wrist, data = newtable)
model2 <-  lm(density~ age + weight + neck + abdomen + thigh + forearm + wrist, data = newtable)
model3 <- lm(density~ age + height + neck + chest + abdomen + hip + thigh + biceps + forearm + wrist, data = newtable)
model4 <- lm(density~ age + weight + abdomen + thigh + forearm + wrist, data = newtable)
model5 <- lm(density~ age + weight + neck + abdomen + hip + thigh + forearm + wrist, data = newtable)
model6 <- lm(density ~ age + abdomen + height, data = newtable)
summary(model1)
summary(model2)
summary(model3)
summary(model4)
summary(model5)
summary(model6)

# show a summary about the cross-validation
print(model1)
print(model2)
print(model3)
print(model4)
print(model5)


Bestmodel = lm(density ~ age +  height + neck +  chest + abdomen + hip +  thigh + forearm + wrist,data=newtable)
#Variance inflation factor and condition number for the best model to check for multicoll
library(car)  # vif
vif(Bestmodel)  # cutoff value is 10 (MPV p. 118)
#weight, chest, abdomen, hip are all above 10
solve(WtW)  # = (W^T W)^{-1}. Compare diagonal values with VIFs












# shows some properties of the different folds
model1$resample

# lists the data points in each fold
model1$pred

# lists the data points in the first fold
subset(model1$pred, Resample == "Fold1") 




