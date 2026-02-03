# Model median inference
median_inf_g <- function(Vtheta, lambda, beta, alpha, t, xmean = NULL) {
  xbar <- c(1, xmean)
  if (any(!is.finite(lambda)))  {
    stop("lambda contains non-finite values.")
  }
  if (any(!is.finite(beta))) {
    stop("beta contains non-finite values.")
  }
  if (any(!is.finite(alpha))) {
    stop("alpha contains non-finite values.")
  }
  if (!is.matrix(beta)) {
    stop("beta must be a tt x (K+1) matrix.")
  }
  if (ncol(beta) != length(xbar)) {
    stop("ncol(beta) must equal length(c(1, xmean)).")
  }
  tt <- length(lambda)
  lambdat <- lambda[t]
  xit <- c(((lambdat * t(beta[t, ]) %*% xbar) + 1) ^ (1 / lambdat))

  lambda00 <- numeric(tt)
  Dltl <- lambda00
  Dltl[t] <- (lambdat ^ (-2) * xit * (1 - lambdat * log(xit) -
                                        xit ^ (-lambdat)))
  Dltb1 <- lambda00
  Dltb1[t] <- 1
  Dltb <- xbar %x% Dltb1 * xit ^ (1 - lambdat)
  Dlt <- c(Dltl, Dltb, numeric(length(alpha)))
  SE <- as.numeric(sqrt(t(Dlt) %*% Vtheta %*% Dlt))
  return(list(median = xit, SE = SE))
}

# Inference on model median difference
median_diff_inf <- function(med_g1, med_g2, level = 0.95){
  d  <- med_g1$median - med_g2$median
  se <- sqrt(med_g1$SE^2 + med_g2$SE^2)

  z <- stats::qnorm(1 - (1 - level)/2)
  lower <- d - z * se
  upper <- d + z * se

  zval <- d / se
  pval <- 2 * (1 - stats::pnorm(abs(zval)))

  list(meddiff = d, SE = se, lower.cl = lower, upper.cl = upper,
       z.value = zval, p.value = pval, level = level)
}

# Inference on probability based treatment effect
prob_eff_inf <- function(theta1, theta2, vtheta1, vtheta2, t,
                         xmean = NULL, smaller, robust = TRUE) {
  xbar <- c(1, xmean)
  tt <- length(theta1$lambda)
  Kp1 <- ncol(theta1$beta)
  na <- length(theta1$alpha)
  if (length(theta2$lambda) != tt) {
    stop("theta1/theta2 must have same dimension.", call. = FALSE)
  }
  if (na != tt * (tt + 1) / 2) {
    stop("alpha length mismatch.", call. = FALSE)
  }
  if (smaller) {
    l1 <- theta1$lambda[t]
    l2 <- theta2$lambda[t]
    m1 <- sum(xbar * theta1$beta[t, ])
    m2 <- sum(xbar * theta2$beta[t, ])
    s1 <- sqrt(theta1$Sigma[t, t])
    s2 <- sqrt(theta2$Sigma[t, t])
    if (robust) {
      Vt <- as.matrix(Matrix::bdiag(vtheta1$V.rob, vtheta2$V.rob))
    } else {
      Vt <- as.matrix(Matrix::bdiag(vtheta1$V.mod, vtheta2$V.mod))
    }
  } else {
    l1 <- theta2$lambda[t]
    l2 <- theta1$lambda[t]
    m1 <- sum(xbar * theta2$beta[t, ])
    m2 <- sum(xbar * theta1$beta[t, ])
    s1 <- sqrt(theta2$Sigma[t, t])
    s2 <- sqrt(theta1$Sigma[t, t])
    if (robust) {
      Vt <- as.matrix(Matrix::bdiag(vtheta2$V.rob, vtheta1$V.rob))
    } else {
      Vt <- as.matrix(Matrix::bdiag(vtheta2$V.mod, vtheta1$V.mod))
    }
  }
  med1 <- (l1 * m1 + 1) ^ (1 / l1)
  med2 <- (l2 * m2 + 1) ^ (1 / l2)
  up1 <- (l1 * (m1 + s1 * qnorm(0.999)) + 1) ^ (1 / l1)
  if (is.nan(up1)) {
    up1 <- Inf
  }
  up2 <- (l2 * (m2 + s2 * qnorm(0.999)) + 1) ^ (1 / l2)
  if (is.nan(up2)) {
    up2 <- Inf
  }
  up <- max(min(med1 * 100, up1), min(med2 * 100, up2))
  prob_integrand <- function(y, l1, m1, s1, l2, m2, s2) {
    z1 <- (y ^ l1 - 1) / l1
    z2 <- (y ^ l2 - 1) / l2
    res <- pnorm(z1, m1, s1) * y ^ (l2 - 1) * dnorm(z2, m2, s2)
    return(res)
  }
  prob <- integrate(prob_integrand, 1e-8, up, l1 = l1, m1 = m1, s1 = s1,
                    l2 = l2, m2 = m2, s2 = s2)$value

  dl10 <- function(y, l1, m1, s1, l2, m2, s2) {
    z1 <- (y ^ l1 - 1) / l1
    z2 <- (y ^ l2 - 1) / l2
    dz1 <- (y ^ l1 * (l1 * log(y) - 1) + 1) / l1 ^ 2
    res <- dz1 * dnorm(z1, m1, s1) *
      y ^ (l2 - 1) * dnorm(z2, m2, s2)
    return(res)
  }

  dl1 <- integrate(dl10, 1e-8, up, l1 = l1, m1 = m1, s1 = s1,
                   l2 = l2, m2 = m2, s2 = s2)$value

  dl20 <- function(y, l1, m1, s1, l2, m2, s2) {
    ly <- log(y)
    yl2 <- y ^ l2
    z2 <- (yl2 - 1) / l2
    dz2 <- (yl2 * (l2 * ly - 1) + 1) / l2 ^ 2
    res <- (ly - dz2 * (z2 - m2) / s2 ^ 2) *
      prob_integrand(y, l1, m1, s1, l2, m2, s2)
    return(res)
  }

  dl2 <- integrate(dl20, 1e-8, up, l1 = l1, m1 = m1, s1 = s1,
                   l2 = l2, m2 = m2, s2 = s2)$value

  db10 <- function(y, l1, m1, s1, l2, m2, s2) {
    z1 <- (y ^ l1 - 1) / l1
    z2 <- (y ^ l2 - 1) / l2
    res <- dnorm(z1, m1, s1) * y ^ (l2 - 1) * dnorm(z2, m2, s2)
    return(res)
  }
  db11 <- integrate(db10, 1e-8, up, l1 = l1, m1 = m1, s1 = s1,
                    l2 = l2, m2 = m2, s2 = s2)$value

  db1 <- -xbar * db11

  db20 <-  function(y, l1, m1, s1, l2, m2, s2) {
    z2 <- (y ^ l2 - 1) / l2
    res <- (z2 - m2) / s2 *
      prob_integrand(y, l1, m1, s1, l2, m2, s2)
    return(res)
  }
  db20 <- integrate(db20, 1e-8, up, l1 = l1, m1 = m1, s1 = s1,
                    l2 = l2, m2 = m2, s2 = s2)$value

  db2 <- xbar / s2 * db20

  ds10 <- function(y, l1, m1, s1, l2, m2, s2) {
    z1 <- (y ^ l1 - 1) / l1
    z2 <- (y ^ l2 - 1) / l2
    res <- (z1 - m1) / s1  * dnorm(z1, m1, s1) *
      y ^ (l2 - 1) * dnorm(z2, m2, s2)
    return(res)
  }
  da1 <- -1 / 2 / s1 * integrate(ds10, 1e-8, up, l1 = l1, m1 = m1, s1 = s1,
                                 l2 = l2, m2 = m2, s2 = s2)$value

  ds20 <-  function(y, l1, m1, s1, l2, m2, s2) {
    z2 <- (y ^ l2 - 1) / l2
    res <- ((z2 - m2) / s2) ^ 2 *  prob_integrand(y, l1, m1, s1, l2, m2, s2)
    return(res)
  }
  da2 <- 1 / s2 ^ 2 / 2 * (integrate(ds20, 1e-8, up, l1 = l1, m1 = m1, s1 = s1,
                                     l2 = l2, m2 = m2, s2 = s2)$value - prob)
  npar <- tt + tt * Kp1 + na
  dt1 <- numeric(npar)
  dt2 <- numeric(npar)
  dt1[t] <- dl1
  dt2[t] <- dl2
  for (i in seq_len(Kp1)) {
    idx_b <- tt + (i - 1) * tt + t
    dt1[idx_b] <- db1[i]
    dt2[idx_b] <- db2[i]
  }
  idx_a <- tt + tt * Kp1 + t * (t + 1) / 2
  dt1[idx_a] <- da1
  dt2[idx_a] <- da2
  dt <- (1 / prob + 1 / (1 - prob)) * c(dt1, dt2)
  se <- sqrt(t(dt) %*% Vt %*% dt)
  lgtprob <- log(prob / (1 - prob))
  z.value <- lgtprob / se
  p.value <- (1 - pnorm(abs(z.value))) * 2
  upp0 <- log(prob / (1 - prob)) + 1.96 * se
  low0 <- log(prob / (1 - prob)) - 1.96 * se
  upper.cl <-1 / (1 + exp(-upp0))
  lower.cl <- 1 / (1 + exp(-low0))
  return(data.frame(prob = prob, se.logit = se, lower.cl = lower.cl,
                    upper.cl = upper.cl, z.value = z.value,
                    p.value = p.value))
}

