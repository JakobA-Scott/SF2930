#install.packages("ggplot2")
#install.packages("foreach")
#install.packages("openxlsx")

source("multiplot.R")

library(ggplot2)
library(foreach)
library(openxlsx)

######################### Section 1: Read data #########################

# This is where your GLM data is read form Cancellation.csv into a table in R
# Note that the folder in which you have Cancellation.csv must be set as the working directory
###### You do not need to change anything in this section. The data will be sorted in a table named glmdata

glmdata <- read.table("Cancellation.csv", header=TRUE, sep=";", dec=",", quote="")

######################### Section 2: Create groups & aggregate data #########################

# Now you need to modify your data so that you can perform a GLM analysis

# First, any continuous variable needs to be grouped into discrete groups
# The code below groups the variable Number of Persons, from you table glmdata, into six groups, and stores this in a new column, NoPGroup
# It also groups the variable Activity into three groups.


###### This is only an example. You need to create your own groups, with breaks that suit your data
###### For example visualize the data for the different variables and refer to the information in the Appendix group the Activity codes on similar types of business.
###### Remember to think about risk homogeneous and stable groups, for example having groups with zero # of claims or claims cost and almost no Duration makes it hard to determine its risk.
###### You might also want to group other variables from glmdata, in a similar manner


glmdata$NoPGroup <- cut(glmdata$Number.of.Persons,
                       breaks = c(-Inf, 2, 5, 15, 50, 100, Inf),
                       labels = c("01_<2", "02_2_5", "03_5-15", "04_15-50", "05_50-100", "06_>=100"),
                       right = FALSE)



Industry <- c('B','D', 'E', 'U')
Service <- c('G','I','S','T')
glmdata$ActivityGroup <- ifelse(glmdata$Activity %in% Industry, "Industry",
                                ifelse(glmdata$Activity %in% Service,"Service",
                                       "Other" ))

# Secondly, we want to aggregate the data.
# That is, instead of having one row per company&year, we want one row for each existing combination of variables
# This code aggregates columns for Duration, # of Claims & ClaimCost of glmdata, by two variables: NoPGroup and ActivityGroup
# Tha aggregated data is stored in a new table, glmdata2
##### You need to consider if there are any other variables you want to aggregate by, and modify the code accordingly

glmdata2 <- aggregate(glmdata[c("Duration", "NumberOfClaims", "ClaimCost")],by=list(NoP_group = glmdata$NoPGroup,
                                                                                    Activity_group = as.factor(glmdata$ActivityGroup)
), FUN=sum, na.rm=TRUE)

# We then do some preparation for the output the GLM function will give.
# This piece of code creates a new table, glmdata3, with a row per variable and group, and with data on the total duration corresponding to this group.
##### You need ot modify the code to take into account any changes in variables you're using

glmdata3 <-
  data.frame(rating.factor =
               c(rep("NoPGroup", nlevels(glmdata2$NoP_group)),
                 rep("ActivityGroup", nlevels(glmdata2$Activity_group))),
             class =
               c(levels(glmdata2$NoP_group),
                 levels(glmdata2$Activity_group)),
             stringsAsFactors = FALSE)

new.cols <-
  foreach (rating.factor = c("NoP_group", "Activity_group"),
           .combine = rbind) %do%
  {
    nclaims <- tapply(glmdata2$NumberOfClaims, glmdata2[[rating.factor]], sum)
    sums <- tapply(glmdata2$Duration, glmdata2[[rating.factor]], sum)
    n.levels <- nlevels(glmdata2[[rating.factor]])
    contrasts(glmdata2[[rating.factor]]) <-
      contr.treatment(n.levels)[rank(-sums, ties.method = "first"), ]
    data.frame(duration = sums, n.claims = nclaims)
  }
glmdata3 <- cbind(glmdata3, new.cols)
rm(new.cols)

######################### Section 3: GLM analysis #########################

# Now we get to the fun part - the GLM analysis. It is performed using R's built in GLM function

# First, we model the claims frequency.
# The first part of this performs a GLM analysis, with glmdata2 as the data source modelling NoOfClaims, by the Duration. It looks at three variables: weight_group, Climate, and ActivityCode.
##### This is where you can modify the model by adding or removing variables

model.frequency <-
  glm(NumberOfClaims ~ NoP_group + Activity_group + offset(log(Duration)),
      data = glmdata2, family = poisson)

# Then we save the coefficients resulting from the GLM analysis in an array
##### You should not need to modify this part of the code

rels <- coef(model.frequency)
rels <- exp(rels[1] + rels[-1])/exp(rels[1])

# Finally, we attach the coefficients to the already prepared table glmdata3, in a column named rels.frequency
# There is no good way of doing this automatically, so we need to do some manual tricks
# This code creates a vector with 6 positions consisting of the integer 1, and then positions number 1-5 in the rels array.
# Then it attaches this to rows 1-6 of glmdata3, sorted from highest to lowest duration, since the GLM data is on this form.
# In other words, the code takes the GLM coefficients for the six NoP groups and saves those in glmdata3, in the rows corresponding to those groups.
# After that, it does the same thing for the rest of the GLM coefficients, belonging to climate and activity code variables.

##### You need to modify this code to suit your set of variables and groups, to make sure each GLM coefficient is saved in the correct place.

##### You need to modify this code to fit your variables
variableLevels <- c(nlevels(glmdata2[["NoP_group"]]),
                    nlevels(glmdata2[["Activity_group"]]))

#You do not need to modify this part
cs <- cumsum(variableLevels)
cs_rels <- cs
for(i in 1:length(variableLevels)){
  cs_rels[i] <- cs[i]-i
}

# The following code needs to be used at two different places so we put it in a function.
##### This part needs to be modified if you change which variables are included in the model, but not if you change the groups inside a variable
attachRels <- function(rels_vec, vector, cs, cs_rels) {
  c(c(1, rels_vec[ 1 : cs_rels[1] ])[rank(-vector[ 1 : cs[1] ], ties.method = "first")],
    c(1, rels_vec[ (cs_rels[1]+1) : cs_rels[2] ])[rank(-vector[ (cs[1]+1) : cs[2] ], ties.method = "first")])
}




# Use the function created above, you do not need to modify this part
glmdata3$rels.frequency <- attachRels(rels, glmdata3$duration, cs, cs_rels)

# We then do the same thing again, now modelling severity instead of claim frequency.
# That means that, in this part, we want to look at the average claim. So first, we calculate the average claim for each row in glmdata2
##### You should not need to change anything in this piece of code.

glmdata2$avgclaim=glmdata2$ClaimCost/glmdata2$NumberOfClaims

# Then we do the same thing as we did when modelling claims frequency, but we look at average claim;
# A GLM analysis is run, the coefficients stored, and saved in a new column, named rels.severity, glmdata3
##### You need to modify this part of the code in the same way as you did for the GLM model for frequency.  Add or remove variables
##### Remember that, according to the project instructions, you need to use the same variables for the severity as for the frequency.

model.severity <-
  glm(avgclaim ~ NoP_group + Activity_group ,
      data = glmdata2[glmdata2$avgclaim>0,], family = Gamma("log"), weight=NumberOfClaims)

# You do not need to change this part
rels <- coef(model.severity)
rels <- exp( rels[1] + rels[-1] ) / exp( rels[1] )
glmdata3$rels.severity <- attachRels(rels, glmdata3$duration, cs, cs_rels)

# Finally, the final risk factor is calculated, as the product of the frequency and severity factors.
##### You should not have to modify this coed.
##### Congratulations! You now have a model for the risk!
glmdata3$rels.risk <- with(glmdata3, rels.frequency*rels.severity)
######################### Section 4: Plotting #########################

# In this section, the results from the GLM are plotted.

# First, long variable names need to be cut, to fit into the plots.
# This row of code cuts away everything except for the first letter for variable names belonging to activity codes.
##### If you have long variable names, modify here to cut them.
glmdata3[glmdata3$rating.factor == "Activity",2] <- substr(glmdata3$class,1,1)[7:29]


# Then the results are plotted. This code plots the GLM factors for frequency, severity, and total risk, for the three variables Weight, Climate, and Activity code.
##### If you have changed what variables are included in your model, add, remove, or modify sections of this code to plot them.
##### This is also where you can make changes to change the look of your plots, if you would like to.

p1 <- ggplot(subset(glmdata3, rating.factor=="NoPGroup"), aes(x=class, y=rels.frequency)) +
  geom_point(colour="blue") + geom_line(aes(group=1), colour="blue") + ggtitle("NoP: frequency factors") +
  geom_text(aes(label=paste(round(rels.frequency,2))), nudge_y=1) +theme(axis.text.x = element_text(angle = 30, hjust = 1))

p2 <- ggplot(subset(glmdata3, rating.factor=="NoPGroup"), aes(x=class, y=rels.severity)) +
  geom_point(colour="blue") + geom_line(aes(group=1), colour="blue") + ggtitle("NoP: severity factors") +
  geom_text(aes(label=paste(round(rels.severity,2))), nudge_y=0.5)+theme(axis.text.x = element_text(angle = 30, hjust = 1))

p3 <- ggplot(subset(glmdata3, rating.factor=="NoPGroup"), aes(x=class, y=rels.risk)) +
  geom_point(colour="blue") + geom_line(aes(group=1), colour="blue") + ggtitle("NoP: risk factors") +
  geom_text(aes(label=paste(round(rels.risk,2))), nudge_y=1.6)+theme(axis.text.x = element_text(angle = 30, hjust = 1))

p4 <- ggplot(subset(glmdata3, rating.factor=="ActivityGroup"), aes(x=class, y=rels.frequency)) +
  geom_point(colour="blue") + geom_line(aes(group=1), colour="blue") + ggtitle("Activity: frequency factors") +
  geom_text(aes(label=paste(round(rels.frequency,2))), nudge_y=0.05)

p5 <- ggplot(subset(glmdata3, rating.factor=="ActivityGroup"), aes(x=class, y=rels.severity)) +
  geom_point(colour="blue") + geom_line(aes(group=1), colour="blue") + ggtitle("Activity: severity factors") +
  geom_text(aes(label=paste(round(rels.severity,2))), nudge_y=0.1)

p6 <- ggplot(subset(glmdata3, rating.factor=="ActivityGroup"), aes(x=class, y=rels.risk)) +
  geom_point(colour="blue") + geom_line(aes(group=1), colour="blue") + ggtitle("Activity: risk factors") +
  geom_text(aes(label=paste(round(rels.risk,2))), nudge_y=0.1)



multiplot(p1,p2,p3,p4,p5,p6, cols=2)



######################### Section 5: Leveling #########################

#Here is a simple example on how you can estimate the base factor for step 2.
#The base factor relies heavily on how you estimate the expected claim cost for next year
#In this example we choose the simplest possible method and predicts it should be the same as last year

###### However it is always important to look into the data and analyze if there are any outliers or other strange outcomes that you need to consider.
############ Do you think 2021 will be similar to 2020, and how similar?
############ How much should older years affect?
############ Are there any trends over the years?

###Estimate Claim cost next year


agg_vars = c('ClaimCost','Duration','NumberOfClaims')
Claims_Per_Year = aggregate(glmdata[agg_vars], by=list(Category=glmdata$RiskYear), FUN=sum)
plot(Claims_Per_Year$Category,Claims_Per_Year$ClaimCost)

### Here are some additional metrics to determine expected claim cost

# Claims_Per_Year$ClaimCostperDuration = Claims_Per_Year$ClaimCost / Claims_Per_Year$Duration
# Claims_Per_Year$NumOfClaimsperDuration = Claims_Per_Year$NumberOfClaims / Claims_Per_Year$Duration
# Claims_Per_Year

Expected_ClaimCost = Claims_Per_Year$ClaimCost[Claims_Per_Year$Category==2020]

Ratio = 0.9

Expected_Premium = Expected_ClaimCost/Ratio


### Estimate total risk for all the active policies

#Now we have an example of calculating the total risk for policies in 2020.
#This is also the most simple method to just choose risk factor based on the group
#If you want to further differentiate within groups you can try to make the risk factor continuous for numerical variables
#For example so that if you a group for 15-50 NoP, then it might be an idea to have different risks if the value is 15 or 50.

#Choose only those policies that are active
Active = subset(glmdata, glmdata$RiskYear == 2020)

#Set up for looping and merging
rating_factors = unique(glmdata3$rating.factor)
keep_vars = c("class","rels.risk")

for (r in rating_factors){

  #choose the factors to compare
  factors = subset(glmdata3[keep_vars], glmdata3$rating.factor == r)

  #Add correct values to each active policy
  Active = merge(Active, factors, by.x=r, by.y='class')

  #Adjust names by variable
  names(Active)[names(Active) == 'rels.risk'] <- paste("Risk",r,sep = "_")

}

#Calculate total risk per policy
Active[is.na(Active)] <- 1
Active$Risk <- apply(Active[grep("Risk_",names(Active))], 1, prod)

#Calculate total risk for all active policies
Total_Risk_Last_Year = sum(Active$Risk)

### Calculate Base factor needed on each current risk to reach desired premium
Base_factor = Expected_Premium/Total_Risk_Last_Year

Base_factor



######################### Section 6: Export factors to Excel #########################

#As a last step, the risk factors are exported to excel.
# The dopcument will be saved in the folder set as your working directory.

write.xlsx(glmdata3, "glmfactors.xlsx")


