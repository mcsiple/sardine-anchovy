## playing with tests between surrogate(s) and observed WMRs

# Load all the necessary functions and pre-treated data
source(here::here("Corrs/NullModel.R"))

# Perform tests with full null accross each combination of datasource, region, variable
compares.RBF <- RBF %>% distinct(region, datasource, variable) # get unique comparisons

# make table results from Wilcoxon-Mann-Whitney test comparing observed vs null WMRs
test.table <- vector()
for(s in 1:nrow(compares.RBF)){ 
  giantnull <- get_large_null(dat = RBF,dsource = compares.RBF$datasource[s],reg = compares.RBF$region[s],var = compares.RBF$variable[s],nsims=50)
  test <- test_wmr(obs = get_obs(dat = RBF,dsource = compares.RBF$datasource[s],reg = compares.RBF$region[s],var = compares.RBF$variable[s]), null.combined = giantnull)
  test$region <- compares.RBF$region[s]
  test$datasource <- compares.RBF$datasource[s]
  test$variable <- compares.RBF$variable[s]
  test.table <- rbind(test.table,test)
}

test.table.out <- test.table %>% select(region, variable, datasource, period, U, N, p = p.value, d = diff, CI.L, CI.U, r = r_p)

test.table.out$period <- recode(test.table.out$period, less.than.5 = "< 5 yr", five.ten = "5-10 yr", ten.plus="10+ yr")
test.table.out$p <- round(test.table.out$p, 3)
test.table.out$p[test.table.out$p < 0.001] <- "< 0.001"
test.table.out$d <- round(test.table.out$d, 3)
test.table.out$CI.L <- round(test.table.out$CI.L, 3)
test.table.out$CI.U <- round(test.table.out$CI.U, 3)
test.table.out$r <- round(test.table.out$d, 2)

test.table.out <- arrange(test.table.out, region, datasource, variable)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# checking for consistency between test_wmr and test_wmr_sub (full null vs subsampled null equal in length to obs wmr)
null.combined <- get_large_null(dat = RBF,dsource = "Barange",reg = "California",var = "ssb",nsims=50)
full <- test_wmr(obs = get_obs(dat = RBF,dsource = "Barange",reg = "California",var = "ssb"),
                 null.combined = null.combined)
test1 = test5 = test10 = diff1 = diff5 = diff10 = CL.L1 = CL.L5 = CL.L10 = vector()
for(i in 1:100){
  blah = test_wmr_sub(obs = get_obs(dat = RBF,dsource = "Barange",reg = "California",var = "ssb"),
                      null.combined = null.combined, n.factor = 1)[[1]]
  test1 <- c(test1, blah[1,2])
  test5 <- c(test5, blah[2,2])
  test10 <- c(test10, blah[3,2])
  diff1 <- c(diff1, blah[1,3])
  diff5 <- c(diff5, blah[2,3])
  diff10 <- c(diff10, blah[3,3])
  CL.L1 <- c(CL.L1, blah[1,4])
  CL.L5 <- c(CL.L5, blah[2,4])
  CL.L10 <- c(CL.L10, blah[3,4])
}
hist(test1)
hist(test5)
hist(test10)
hist(diff1)
abline(v=full[1,3], col = "red")
hist(diff5)
abline(v=full[2,3], col = "red")
hist(diff10)
abline(v=full[3,3], col = "red")
hist(CL.L1)
abline(v=full[1,4], col = "red")
hist(CL.L5)
abline(v=full[2,4], col = "red")
hist(CL.L10)
abline(v=full[3,4], col = "red")

# results appear fairly consistent among individual surrogate tests and full null test
# diff equivalent or within 0.005 that from full null test
# lower CL typically =< 0.01 that from full null test

# could also try to get threshold values by percentile method of bootstrapped medians but seems unnecessary 
# unless we really want test of diff in medians vs medians of diffs, but estimated differences are likely similar
# (also, this test is already a form of permutation test)