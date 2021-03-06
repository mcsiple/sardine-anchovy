basedir <- "C:/Users/siplem"

source(here::here("R/WMR/NullModel.R"))



# playing with tests between surrogate(s) and observed WMRs
null.combined <- get_large_null(dat = RBF,dsource = "Barange",reg = "California",var = "ssb",nsims=50)
full <- test_wmr(obs = get_obs(dat = RBF,dsource = "Barange",reg = "California",var = "ssb"),
                 null.combined = null.combined)
full

# Try with one of the negative ones
nc <- get_large_null(dat = RBF,dsource = "Barange",reg = "Humboldt",var = "landings",nsims=50)
negobs <- get_obs(dat = RBF,dsource = "Barange",reg = "Humboldt",var = "landings")
plot(1:length(negobs$std_anchovy),negobs$std_anchovy,type='l')
lines(1:length(negobs$std_sardine),negobs$std_sardine,col='red')
test2 <- test_wmr(obs = get_obs(dat = RBF,dsource = "Barange",reg = "Humboldt",var = "landings"),
                 null.combined = nc)
test2

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