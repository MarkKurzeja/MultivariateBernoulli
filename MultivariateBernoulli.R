rm(list = ls())

library(plyr)
library(dplyr)
library(testit)
library(matrixStats)
library(mcmc)

gen_power_series <- function(elements, up_to_order) {
  # Function that generates all of the elements of the power series of elements
  # up to and including order up_to_order.
  # Example: elements = 1,2,3, up_to_order = 2:
  # {1,2,3,12,13,23}
  stopifnot(up_to_order <= length(elements))
  unlist(lapply(1:up_to_order, function(i) {
    combn(elements, i, simplify = F)
  }),
  recursive = F)
}

test_gen_power_series <- function() {
  assert("Testing gen_power_series", {
    (gen_power_series(c(1, 2, 3), 2)[[1]] == 1)
    (gen_power_series(c(1, 2, 3), 2)[[2]] == 2)
    (gen_power_series(c(1, 2, 3), 2)[[3]] == 3)
    (gen_power_series(c(1, 2, 3), 2)[[4]] == c(1, 2))
    (gen_power_series(c(1, 2, 3), 2)[[5]] == c(1, 3))
    (gen_power_series(c(1, 2, 3), 2)[[6]] == c(2, 3))
    (length(gen_power_series(c(1, 2, 3), 3)) == 2 ^ 3 - 1)
    (length(gen_power_series(c(1, 2, 3, 4), 4)) == 2 ^ 4 - 1)
    (gen_power_series(c(1, 2, 3, 4, 5), 5)[[1]] == 1)
    (gen_power_series(c(1, 2, 3, 4, 5), 5)[[5]] == 5)
    (gen_power_series(c(1, 2, 3, 4, 5), 5)[[31]] == c(1, 2, 3, 4, 5))
  })
}; test_gen_power_series()

f_positions <- function(vector_pos) {
  # Takes in a vector of positive integers >= 1 of length <= N where n is the
  # number of Bernoulli variables in the joint distribution of y and returns the
  # position in the range of 1-2^(n-1) that corresponds to the position in the
  # f-vector.
  # Example:
  # (1,3,4) -> binary(001101) ->
  # 0 * 2^5 +  0 * 2^4 + 1 * 2^3 + 1 * 2^2 + 0 * 2^1 + 1 * 2^0 ->
  # 8 + 4 + 1 -> 13 in the f-vector
  # NOTE, the algorithm stores the f-positions like so:
  # (1, 2, 1-2, 3, 1-3, 2-3, 1-2-3, 4, 1-4, 2-4, ....) so that even if you add
  # more elements to the array at position n + 1, you will have stable positions
  # for the elements 1:n which is desirable since many of the f-parameters will
  # be used in many joint distributions.
  stopifnot(length(vector_pos) == length(unique(vector_pos)))
  stopifnot(sum(vector_pos) != 0)
  sum(2 ^ (vector_pos - 1))
}

test_f_positions <- function() {
  assert("f_positions should be returning the corresponding input in the high dim vector", {
    (f_positions(c(1)) == 1)
    (f_positions(c(1, 2)) == 3)
    (f_positions(c(2)) == 2)
    (f_positions(c(3)) == 4)
    (f_positions(c(1, 2, 3)) == 7)
    (f_positions(c(1, 3)) == 5)
    (f_positions(c(1, 2, 3, 4)) == 15)
  })
}; test_f_positions()


b_function <- function(y_vector, indexing_vector) {
  # Return 1 if all of the elements that are referenced in the indexing vector
  # are 1 in the y_vector, else return 0
  return(prod(y_vector[indexing_vector]))
}

test_b_function <- function() {
  assert("b-function should return the indexing vector of the y-vectors", {
    (b_function(c(1, 1, 1), c(1, 2, 3)) == 1)
    (b_function(c(1, 0, 1), c(1, 2)) == 0)
    (b_function(c(1), c(1)) == 1)
    (b_function(c(1, 0, 0, 0, 1), c(1)) == 1)
    (b_function(c(1, 0, 0, 0, 1), c(1, 3)) == 0)
  })
}; test_b_function()


propto_log_likli_multi_bern <- function(y, f) {
  # Compute a value that is proportional to -ln(p(y|f)) for a multivariate
  # Bernoulli distribution - this function is basically summing over the values
  # of f in such a way that only the non-zero interactions are considered.
  # TODO(markkurzeja): Implement this function more efficiently through the
  # following transformation:
  # y = (1, 1, 0, 0, 1) -> let k = (1,2,5) where k is just giving the non-zero
  # entries. Then you take all of the combinations you find that only a few of
  # them can be non-zero at any given time. Especially for computations where
  # the vector y is non-dense, this is an efficient computation since there are
  # few elements in the power-set that will be non-zero. You are taking advantage
  # of the fact that this sparsity is computable.
  stopifnot(length(f) == 2^length(y) - 1)
  all_combins = gen_power_series(seq_along(y), length(y))
  sum(sapply(all_combins, function(combination) {
    if (b_function(y, combination) == 0) {
      return(0)
    }
    return(f[f_positions(combination)])
  }))
}


test_propto_log_likli_multi_bern <- function() {
  assert("Testing_log_multi_bern", {
    (propto_log_likli_multi_bern(c(0, 0, 0), 1:7) == 0)
    (propto_log_likli_multi_bern(c(1, 0, 0), 1:7) == 1)
    (propto_log_likli_multi_bern(c(0, 1, 0), 1:7) == 2)
    (propto_log_likli_multi_bern(c(0, 0, 1), 1:7) == 4)
    (propto_log_likli_multi_bern(c(1, 1, 0), 1:7) == 6)
    (propto_log_likli_multi_bern(c(1, 0, 1), 1:7) == 10)
    (propto_log_likli_multi_bern(c(0, 1, 1), 1:7) == 12)
    (propto_log_likli_multi_bern(c(1, 1, 1), 1:7) == 28)
    (propto_log_likli_multi_bern(c(0, 0, 0), c(1, 1, 1, 1, 1, 1, 1)) == 0)
    (propto_log_likli_multi_bern(c(1, 0, 0), c(1, 1, 1, 1, 1, 1, 1)) == 1)
    (propto_log_likli_multi_bern(c(0, 1, 0), c(1, 1, 1, 1, 1, 1, 1)) == 1)
    (propto_log_likli_multi_bern(c(0, 0, 1), c(1, 1, 1, 1, 1, 1, 1)) == 1)
    (propto_log_likli_multi_bern(c(1, 1, 0), c(1, 1, 1, 1, 1, 1, 1)) == 3)
    (propto_log_likli_multi_bern(c(1, 1, 1), c(1, 1, 1, 1, 1, 1, 1)) == 7)
  })
}; test_propto_log_likli_multi_bern()


number_to_bin_vector <- function(decimal_number, length_of_bin_vector) {
  # Turn a number into its binary vector equivalent:
  # 4 -> c(1,0,0), 15 -> c(1,1,1,1)
  # decimal_number is the number to convert.
  # length_of_bin_vector is the maximum length of the binary vector.
  stopifnot(decimal_number < 2^length_of_bin_vector)
  result = rep(0, length_of_bin_vector)
  number_left = decimal_number
  for(i in (length_of_bin_vector-1):0) {
    if(number_left >= 2^i) {
      number_left = number_left - 2^i
      result[length_of_bin_vector - i] = 1
    }
  }
  return(result)
}

test_number_to_bin_vector <- function() {
  assert("Testing test_number_to_bin_vector function", {
    (number_to_bin_vector(0, 1) == c(0))
    (number_to_bin_vector(0, 4) == c(0,0,0,0))
    (number_to_bin_vector(1, 1) == c(1))
    (number_to_bin_vector(1, 3) == c(0,0,1))
    (number_to_bin_vector(2, 2) == c(1,0))
    (number_to_bin_vector(5, 6) == c(0,0,0,1,0,1))
    (number_to_bin_vector(15, 6) == c(0,0,1,1,1,1))
    (number_to_bin_vector(32, 6) == c(1,0,0,0,0,0))
  })

}; test_number_to_bin_vector()

almost_equal <- function(x,y, tol = 1e-5) {
  # A simple function that takes in two numbers, x & y and returns true if they
  # are within tol of each other - this function is necessary for almost
  # equality due to the issues with floats.
  return(abs(x - y) < tol)
}

log_normalizing_constant_multi_bern <- function(size_of_bin_vec, f) {
  # Computes the denominator for the log probability function of the
  # multivariate Bernoulli. size_of_bin_vec is the size of the y-vector that
  # contains the joint probabilities.
  result = rep(0, 2^(size_of_bin_vec))
  for(i in 0:(2^(size_of_bin_vec) - 1)) {
    y_temp = number_to_bin_vector(i, size_of_bin_vec)
    result[i + 1] = propto_log_likli_multi_bern(y_temp, f)
  }
  return(matrixStats::logSumExp(result))
}

test_log_normalizing_constant_multi_bern <- function() {
  assert("Testing the by-hand cases to tie out to the spreadsheet.", {
    (almost_equal(log_normalizing_constant_multi_bern(3, rep(-1,7)), 0.8126671116))
    (almost_equal(log_normalizing_constant_multi_bern(3, rep(0,7)), 2.079441542))
    (almost_equal(log_normalizing_constant_multi_bern(3, 1:7), 28))
    (almost_equal(log_normalizing_constant_multi_bern(2, c(0,0,0)), 1.386294361))
    (almost_equal(log_normalizing_constant_multi_bern(2, c(1,1,1)), 3.27797837))
    (almost_equal(log_normalizing_constant_multi_bern(2, c(1,2,3)), 6.027160139))
  })

}; test_log_normalizing_constant_multi_bern()

log_pmultibern <- function(y, f) {
  # A function that implements the joint distribution of the multivariate
  # Bernoulli distribution in the log form.
  propto_log_likli_multi_bern(y, f) - 
    log_normalizing_constant_multi_bern(length(y), f)
}

test_log_pmultibern <- function() {
  assert("Testing to see that the uniform cases are correct.", {
    (almost_equal(log_pmultibern(c(0,0), rep(0,3)), log(1/4)))
    (almost_equal(log_pmultibern(c(1,0), rep(0,3)), log(1/4)))
    (almost_equal(log_pmultibern(c(0,1), rep(0,3)), log(1/4)))
    (almost_equal(log_pmultibern(c(1,1), rep(0,3)), log(1/4)))
    (almost_equal(log_pmultibern(c(0,0,0), rep(0,7)), log(1/8)))
    (almost_equal(log_pmultibern(c(1,0,0), rep(0,7)), log(1/8)))
    (almost_equal(log_pmultibern(c(0,1,1), rep(0,7)), log(1/8)))
    (almost_equal(log_pmultibern(c(1,1,1), rep(0,7)), log(1/8)))
  })
  assert("Testing to see that the by-hand examples are correct.", {
    (almost_equal(log_pmultibern(c(0,0), 1:3), -6.027160139))
    (almost_equal(log_pmultibern(c(1,0), 1:3), -5.027160139))
    (almost_equal(log_pmultibern(c(0,1), 1:3), -4.027160139))
    (almost_equal(log_pmultibern(c(1,1), 1:3), -0.027160139))
    # I know that it seems very improbable that logsumexp actually produced
    # exactly 28 as its answer, but that is the case. When you are working with
    # really large exponents, it is the case that logsumexp almost works like
    # softmax and since the y = (1,1,1) vector has an non-normalized likelihood
    # of 28 and the rest are much smaller (the next largest is 12), we find that
    # the 28 term dominates and we see a normalizing constant of 28.
    (almost_equal(log_pmultibern(c(0,0,0), 1:7), -28.000000))
    (almost_equal(log_pmultibern(c(1,0,0), 1:7), -27.000000))
    (almost_equal(log_pmultibern(c(0,1,1), 1:7), -16.000000))
    (almost_equal(log_pmultibern(c(1,1,1), 1:7), 0.000000))
  })
}; test_log_pmultibern()


pmultibern <- function(y,f) {
  # Return the joint probability of a vector y \in {0,1}^N under the
  # Multivariate Bernoulli model with natural parameters f(Powerset(y)).
  return(exp(log_pmultibern(y,f)))
}

test_pmultibern <- function() {
  assert("Testing that all zeros always yields a valid symmetric density.", {
    (almost_equal(pmultibern(c(0,0), c(0,0,0)), 1/4))
    (almost_equal(pmultibern(c(0,1), c(0,0,0)), 1/4))
    (almost_equal(pmultibern(c(1,0), c(0,0,0)), 1/4))
    (almost_equal(pmultibern(c(1,1), c(0,0,0)), 1/4))
    (almost_equal(pmultibern(c(0,0,0), c(0,0,0,0,0,0,0)), 1/8))
    (almost_equal(pmultibern(c(1,0,1), c(0,0,0,0,0,0,0)), 1/8))
    (almost_equal(pmultibern(c(1,1,1), c(0,0,0,0,0,0,0)), 1/8))
  })
  assert("Testing some of the normal use cases.", {
    (almost_equal(pmultibern(c(0,0), c(1,1,1)), 0.037704))
    (almost_equal(pmultibern(c(1,0), c(1,1,1)), 0.1024911))
    (almost_equal(pmultibern(c(0,1), c(1,1,1)), 0.1024911))
    (almost_equal(pmultibern(c(1,1), c(1,1,1)), 0.7573132))
    RANDOM_F_SEQ = seq(from = -0.5, to = 0.7, by = 0.2) 
    (almost_equal(pmultibern(c(0,0,0), RANDOM_F_SEQ), 0.1230384837))
    (almost_equal(pmultibern(c(0,1,0), RANDOM_F_SEQ), 0.09114915058))
    (almost_equal(pmultibern(c(1,1,0), RANDOM_F_SEQ), 0.05002371446))
    (almost_equal(pmultibern(c(1,1,1), RANDOM_F_SEQ), 0.24776907))
  })
}; test_pmultibern()


posterior_update_f_multibern_indep_normal <- function(Y, prior_f) {
  # Returns the posterior probability of seeing a particular prior_f given the
  # data in Y. Y is an nxk dimensional matrix where n is the number of
  # observations and k is the dimension of the vector. The length of prior_f is
  # equal to 2^k - 1. prior_f is assumed to come from a multi-variate normal
  # with independent marginals (i.e. the covariance matrix is diagional), and
  # further, the variance of each of the priors is one.
  # We implement this function using the log of the posterior for numerical
  # stability but return the result in the normal units.
  prior_prob = sum(dnorm(prior_f, mean = 0, sd = 1, log = TRUE))
  # cat("prior_prob: ", prior_prob, "\n")
  data_prob = sapply(1:nrow(Y), function(i) {
    # cat("data_prob: ", log_pmultibern(Y[i,], prior_f), "\n")
    log_pmultibern(Y[i,], prior_f)
  })
  return(exp(prior_prob + sum(data_prob)))
}



posterior_update_f_multibern_indep_normal(matrix(c(0), ncol = 1), 0)



test_posterior_update_f_multibern_indep_normal <- function() {
  assert("Testing the by-hand examples", {
    (almost_equal(posterior_update_f_multibern_indep_normal(matrix(c(0), ncol = 1), c(0)), 0.1994711402))
    (almost_equal(posterior_update_f_multibern_indep_normal(matrix(c(1), ncol = 1), c(0)), 0.1994711402))
    (almost_equal(posterior_update_f_multibern_indep_normal(matrix(c(0), ncol = 1), c(1)), 0.06507595058))
    (almost_equal(posterior_update_f_multibern_indep_normal(matrix(c(1), ncol = 1), c(1)), 0.1768947739))
    (almost_equal(posterior_update_f_multibern_indep_normal(matrix(c(0,1), ncol = 1), c(1)), 0.04757433194))
    (almost_equal(posterior_update_f_multibern_indep_normal(matrix(c(1,0), ncol = 1), c(1)), 0.04757433194))
    RANDOM_F_SEQ = seq(from = -0.5, to = 0.7, by = 0.2) 
    (almost_equal(
      posterior_update_f_multibern_indep_normal(
        matrix(
          c(1,1,1, 0,1,0), 
          ncol = 3, 
          byrow = T), 
        RANDOM_F_SEQ), 2.0033875e-5, tol = 1e-9))
  })
}; test_posterior_update_f_multibern_indep_normal()




extract_counts_from_matrix <- function(data) {
  # This function takes in a matrix of data and counts the occurrence of each of
  # the elements in the matrix data - for instance, if it sees 3 counts of the
  # (1,1,1) vector in a 3 dimensional dataset, then it will return 3 in position
  # 7 of the counts vector (since position 3 corresponds to f^123 in the
  # f-vector nomenclature. 

  
  # index_multipliers = 2^(rev(seq_len(ncol(data))) - 1)
  # browser()
  index_multipliers = 2^(seq_len(ncol(data)) - 1)
  index_multipliers_matrix = matrix(index_multipliers, ncol = 1)
  indicies = data %*% index_multipliers_matrix %>% as.numeric() 
  
  counts = rep(0, 2^ncol(data) - 1)
  for(i in indicies) {
    if(i == 0) {
      next
    }
    counts[i] <- counts[i] + 1
  }
  
  return(list(counts = counts, n = nrow(data)))
}


test_extract_counts_from_matrix <- function() {
  fake_data = matrix(c(1,1,1,
                       0,1,1,
                       0,1,1,
                       1,1,1,
                       0,0,0,
                       1,1,1, 
                       1,0,0,
                       0,0,1), byrow = T, ncol = 3)
  assert("Testing the construction case", {
    (extract_counts_from_matrix(fake_data)$n == 8) 
    (extract_counts_from_matrix(fake_data)$counts == c(1,0,0,1,0,2,3))
  })
  fake_data = matrix(c(0,0,
                       0,1,
                       1,0,
                       1,1), byrow = T, ncol = 2)
  assert("Testing the bivariate case fully", {
    (extract_counts_from_matrix(fake_data)$n == 4)
    (extract_counts_from_matrix(fake_data)$counts == c(1,1,1))
  })
}; test_extract_counts_from_matrix()







S_function <- function(indicies, f_vector) {
  indicies <- c(1,3,5) 
  all_combins <- gen_power_series(indicies,3)
  f_pos <- sapply(all_combins, f_positions)
  sum(f_vector[f_pos])  
}

assert("Need to test S function", (F == T))







# Now we just messing around :)


k = c(1,5, 1)
n = 100





log_liklihood_multibern_normal_prior <- function(f) {
  assert("F-vector only implemented for bivariate case", {
    (length(f) == 3)
  })
  
  k[1] * f[1] + k[2] * f[2] + # Univariate terms 
    k[3] * f[3] - # bivariate terms
    n * log(1 + exp(f[1]) + exp(f[2]) + exp(f[1] + f[2] + f[3])) - # normalizing terms 
    sum(f^2 / 2) # normal prior
}

log_liklihood_multibern_normal_prior(c(0,0,0))



out <- mcmc::metrop(log_liklihood_multibern_normal_prior, c(0,0,0), 1e5)
cat("We had ", out$accept , " acceptance rate\n")

results <- out$batch

f1 <- sample(results[,1], size = 1000)
f2 <- sample(results[,2], size = 1000)
f12 <- sample(results[,3], size = 1000)

hist(f1, freq = F)
curve(dnorm, add = T)

hist(f2, freq = F)
curve(dnorm, add = T)

hist(f12, freq = F)
curve(dnorm, add = T)


p00 = {1 / (1 + exp(f1) + exp(f2) + exp(f1 + f2 + f12))}; hist(p00)
p10 = {exp(f1) / (1 + exp(f1) + exp(f2) + exp(f1 + f2 + f12))} %T>% hist()
p01 = {exp(f2) / (1 + exp(f1) + exp(f2) + exp(f1 + f2 + f12))} %T>% hist()
p11 = {exp(f1 + f2 + f12) / (1 + exp(f1) + exp(f2) + exp(f1 + f2 + f12))} %T>% hist()










