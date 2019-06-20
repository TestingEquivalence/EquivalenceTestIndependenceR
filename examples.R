source("asymptotic_test.R")
alpha=0.05

nitren_data  = matrix(data=c(9, 13, 13, 48,24, 18, 20, 72),
                      nrow=2, ncol=4, byrow=TRUE)

eye_hair_color_data=matrix(
  data=c(68, 119, 26, 7,20, 84, 17, 94,15, 54, 14, 10,5, 29, 14, 16),
  nrow=4,ncol=4, byrow=TRUE)

children_income_data=matrix(
  data=c(2161, 3577, 2184, 1636,
         2755, 5081, 2222, 1052,
         936, 1753, 640, 306,
         225, 419, 96, 38,
         39, 98, 31, 14),
  nrow=5, ncol=4, byrow = TRUE)

results=matrix(data=NA, nrow=4,ncol=3)
colnames(results)=c("Nitren","Color","Children")
rownames(results)=c("asymptotic absolute","bootstrap absolute",
                    "asymptotic relative","bootstrap relative")


results[1,1]=asymptotic_test_absolute(nitren_data,alpha)
results[1,2]=asymptotic_test_absolute(eye_hair_color_data,alpha)
results[1,3]=asymptotic_test_absolute(children_income_data,alpha)
 
results[3,1]=asymptotic_test_relative(nitren_data,alpha)
results[3,2]=asymptotic_test_relative(eye_hair_color_data,alpha)
results[3,3]=asymptotic_test_relative(children_income_data,alpha)

results[2,1]=bootstrap_test_absolute(nitren_data,alpha)
results[2,2]=bootstrap_test_absolute(eye_hair_color_data,alpha)
results[2,3]=bootstrap_test_absolute(children_income_data,alpha)

results[4,1]=bootstrap_test_relative(nitren_data,alpha)
results[4,2]=bootstrap_test_relative(eye_hair_color_data,alpha)
results[4,3]=bootstrap_test_relative(children_income_data,alpha)
