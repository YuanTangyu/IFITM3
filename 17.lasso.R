library(glmnet)
library(survival)
library(tidyverse)

dta=readxl::read_excel("") %>%   
  as.data.frame()
rownames(dta) <- dta$ID 

dta <- dta[-1] %>% t() %>% as.data.frame()

ch <- c(rep(1,218),rep(0,35)) %>% as.data.frame() %>%  
  set_names('ch')

set.seed(1)
x <- as.matrix(dta)
y <- model.matrix(~ch-1,ch) %>% data.matrix() %>% as.factor()
cvla <- glmnet(x,y,family = 'binomial')

cv.fit <- cv.glmnet(x,y,family='binomial')

plot(cvla,xvar='lambda',label=T)
plot(cv.fit)

cv.fit$lambda.min
cv.fit$lambda.1se

total=cv.fit$lambda.min+cv.fit$lambda.1se

coef=coef(cvla,s=total)
index=which(coef!=0)
actcoef=coef[index]
LassoGene=row.names(coef)[index]
genecoef=cbind(Gene=LassoGene,coef=actcoef)
genecoef

coef.min=coef(cvla,s=cv.fit$lambda.min)
index.min=which(coef.min!=0)
actcoef.min=coef.min[index.min]
LassoGene.min=row.names(coef.min)[index.min]
genecoef.min=cbind(Gene=LassoGene.min,coef=actcoef.min)
genecoef.min
