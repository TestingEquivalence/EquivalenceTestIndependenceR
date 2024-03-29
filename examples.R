source("asymptotic_test.R")
source("bootstrap_test.R")

# This project contains asymptotic and bootstrap tests,
# which are based on the article:
# Vladimir Ostrovski "New Equivalence Tests for Approximate Independence in Contingency Tables".
# Stats 2019, 2, 239–246; doi:10.3390/stats2020018
# http://www.mdpi.com/journal/stats

# The project contains equivalence tests for approximate row column independence
# in two way contingency tables. 
# Two different test statistics are used:
# 
# 1) First test statistic is scaled Euclidian distance between the contingency table
#    and the product measure of the marginal distributions.
#    The asymptotic test for the test statistic is the function "asymptotic_test_absolute".
#    The bootstrap test is the function "bootstrap_test_absolute".
# 2) Second test statistic is scaled Euclidian norm of the relative deviation 
#    between the contingency table and the product measure of the marginal distributions.
#    The asymptotic test for the test statistic is the function "asymptotic_test_relative".
#    The bootstrap test is the function "bootstrap_test_relative".
# All tests return the minimum tolerance parameter epsilon,
# for which the approximate independence can be shown.

# The real data examples, which are considered in the paper, are presented below.
# -------------------------------------------------------------------------------

# significance level
alpha=0.05

#table of tests results
results=matrix(data=NA, nrow=4,ncol=3)
colnames(results)=c("Nitren","Color","Children")
rownames(results)=c("asymptotic absolute","bootstrap absolute",
                    "asymptotic relative","bootstrap relative")

# Nitren data set
# -------------------------------------------------------------------------------
# Contingency table relating gender and treatment outcome on nitrendipine mono-therapy
# in patients suffering from mild arterial hypertension

nitren_data  = matrix(data=c(9, 13, 13, 48,24, 18, 20, 72),
                      nrow=2, ncol=4, byrow=TRUE)

results[1,1]=asymptotic_test_absolute(nitren_data,alpha)
results[2,1]=bootstrap_test_absolute(nitren_data,alpha)
results[3,1]=asymptotic_test_relative(nitren_data,alpha)
results[4,1]=bootstrap_test_relative(nitren_data,alpha)

# Cross-classification of eye color and hair color
# -------------------------------------------------------------------------------

eye_hair_color_data=matrix(
  data=c(68, 119, 26, 7,20, 84, 17, 94,15, 54, 14, 10,5, 29, 14, 16),
  nrow=4,ncol=4, byrow=TRUE)

results[1,2]=asymptotic_test_absolute(eye_hair_color_data,alpha)
results[2,2]=bootstrap_test_absolute(eye_hair_color_data,alpha)
results[3,2]=asymptotic_test_relative(eye_hair_color_data,alpha)
results[4,2]=bootstrap_test_relative(eye_hair_color_data,alpha)

# Cross-classification of number of children by annual income
# -------------------------------------------------------------------------------

children_income_data=matrix(
  data=c(2161, 3577, 2184, 1636,
         2755, 5081, 2222, 1052,
         936, 1753, 640, 306,
         225, 419, 96, 38,
         39, 98, 31, 14),
  nrow=5, ncol=4, byrow = TRUE)


results[1,3]=asymptotic_test_absolute(children_income_data,alpha)
results[2,3]=bootstrap_test_absolute(children_income_data,alpha)
results[3,3]=asymptotic_test_relative(children_income_data,alpha)
results[4,3]=bootstrap_test_relative(children_income_data,alpha)
