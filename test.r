#+ include=FALSE
remove(list=ls())
library(MASS)
library(tidyverse)
library(survival)
library(rootSolve)

#+ echo=FALSE
IVHR.multi <- function(Time, Status, X, IV, Z = NULL) {
    # X is the exposure
    # IV is the instrumental variable
    # Z is an optional matrix of covariates
    n <- length(Time)
    ord <- order(-Time)
    Time <- Time[ord]
    Status <- Status[ord]
    X <- as.matrix(as.matrix(X)[ord, ])
    IV <- as.matrix(as.matrix(IV)[ord, ])
    if (dim(X)[2] != dim(IV)[2]) stop("X and IV must have equal number of columns.")
    if (dim(X)[1] != dim(IV)[1]) stop("X and IV must have equal number of rows.")
    if (is.null(Z)) Z <- matrix(nrow = n, ncol = 0)
    Z <- as.matrix(Z)[ord, ]
    XX <- cbind(X, Z)
    W <- cbind(IV, Z)
    S.W1 <- matrix(nrow = n, ncol = dim(W)[2])
    Est.Equat <- function(beta) {
        HR <- exp(XX %*% beta)
        S.0 <- cumsum(HR)
        for (i.col in 1:dim(W)[2]) {
            S.W1[, i.col] <- cumsum(W[, i.col] * HR)
        }
        colSums(Status * (W - S.W1 / S.0))
    }
    out.solution <- multiroot(Est.Equat, start = rep(0, dim(W)[2]))
    no.root.found <- is.na(out.solution$estim.precis) | out.solution$estim.precis >= 0.00001
    if (no.root.found) {
        beta.hat <- rep(NA, dim(W)[2])
        Sandwich <- matrix(NA, nrow = dim(W)[2], ncol = dim(W)[2])
    }
    if (!no.root.found) {
        beta.hat <- out.solution$root
        # Variance
        HR <- exp(XX %*% beta.hat)
        S.0 <- cumsum(HR)
        W1 <- S.X1 <- matrix(nrow = n, ncol = dim(W)[2])
        for (i.col in 1:dim(W)[2]) {
            S.W1[, i.col] <- cumsum(W[, i.col] * HR)
            S.X1[, i.col] <- cumsum(XX[, i.col] * HR)
        }
        S.W1X1 <- array(dim = c(n, dim(W)[2], dim(W)[2]))
        Var.E.E <- Deriv <- matrix(nrow = dim(W)[2], ncol = dim(W)[2])
        for (i.col in 1:dim(W)[2]) {
            for (j.col in 1:dim(W)[2]) {
                S.W1X1[, i.col, j.col] <- cumsum(W[, i.col] * XX[, j.col] * HR)
                Var.E.E[i.col, j.col] <- sum(Status * (W[, i.col] - S.W1[, i.col] / S.0) * (W[, j.col] - S.W1[, j.col] / S.0))
                Deriv[i.col, j.col] <- -sum(Status * (S.W1X1[, i.col, j.col] / S.0 - (S.W1[, i.col] * S.X1[, j.col]) / S.0^2))
            }
        }
        Inv.Deriv <- ginv(Deriv)
        Sandwich <- t(Inv.Deriv) %*% Var.E.E %*% Inv.Deriv
    }
    list(Est.log.HR = beta.hat, SE = diag(Sandwich)^0.5, Sandwich = Sandwich)
}

#' read and clean data
dat <- read.csv("./NewTestTable9Dec2022.csv", sep = ";") %>% as_tibble()
dat %>% ggplot(aes(DSSTDY)) + geom_histogram(binwidth = 5)
dat %>%
    count(Death, DSSTDY) %>%
    drop_na()  %>%
    ggplot() +
    geom_col(aes(DSSTDY, n, fill = factor(Death)), position = position_dodge2(,'single'))
dat_cleaned <- dat %>%
    select(DSSTDY, Death, TestResult, AGE, SEX, p_total) %>%
    mutate(
        TestResult = if_else(TestResult == "POS", 1, 0),
        SEX = if_else(SEX == "M", 1, 0),
        AGE = as.numeric(AGE)
    ) %>%
    drop_na()

head(dat_cleaned)
dat_cleaned %>% count(Death)
dat_cleaned %>% count(TestResult)

#' Test run
IVHR.multi(
    Time = dat_cleaned$DSSTDY,
    Status = dat_cleaned$Death,
    X = dat_cleaned$TestResult,
    IV = dat_cleaned$p_total,
    Z = dat_cleaned[, c("AGE", "SEX")]
)
#' Not working, due to `rootSolve` failed to find solutions
#' 
#' A simpler data test with simulated IV
set.seed(123)
IVHR.multi(
    Time = dat_cleaned$DSSTDY,
    Status = dat_cleaned$Death,
    X = dat_cleaned$TestResult,
    IV = rnorm(nrow(dat_cleaned)),
    # Z = dat_cleaned[, c("AGE", "SEX")]
)
set.seed(1234)
IVHR.multi(
    Time = dat_cleaned$DSSTDY,
    Status = dat_cleaned$Death,
    X = dat_cleaned$TestResult,
    IV = rnorm(nrow(dat_cleaned)),
    # Z = dat_cleaned[, c("AGE", "SEX")]
)
#' Sometime work sometime not
#' 
#' Need to:
#' 
#' - Review of variables
#' - A more robust estimation algoritm than `rootSolve` (probably numerical or Baysian methods)