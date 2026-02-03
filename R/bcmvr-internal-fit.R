# Estimating parameters for the BCMVR model
bcmvr_fit_g <- function(y, x = NULL){
  if (any(y[!is.na(y)] <= 0)) {
    stop("y must be positive for Box-Cox transform.")
  }
  if (!is.null(x) && anyNA(x)) {
    stop("Missing value is included in x.")
  }
  tr <- function(X) {
    sum(diag(X))
  }
  dec2bin <- function(num, digit = 0){
    if(num <= 0 && digit <= 0){
      return(NULL)
    }else{
      return(append(Recall(num %/% 2,digit - 1), num %% 2))
    }
  }
  SgmC <- function(alpha, tt){
    Sgm0 <- matrix(0, tt, tt)
    Sgm0[upper.tri(Sgm0, diag = TRUE)] <- alpha
    Sgm <- Sgm0 + t(Sgm0)
    diag(Sgm) <- diag(Sgm) / 2
    return(Sgm)
  }
  alpha_index <- function(tt) {
    which(upper.tri(matrix(0, tt, tt), diag = TRUE), arr.ind = TRUE)
  }

  misscalc <- function(y) {
    n <- nrow(y)
    tt <- ncol(y)
    dp <- numeric(n)
    tt2 <- 2 ^ tt - 1
    colf <- list()
    for (j in seq_len(tt2)){
      mv <- dec2bin(j - 1, tt)
      dp[colSums(apply(is.na(y), 1, '==', mv)) == tt] <- j
      colf[[j]] <- seq_len(tt)[mv == 0]
    }

    ndp <- numeric(tt2)
    for (j in seq_len(tt2)) {
      ndp[j] <- sum(dp == j)
    }
    return(list(colf = colf, dp = dp, ndp = ndp))
  }
  chol_spd_ridge <- function(S) {
    p <- ncol(S)
    d <- diag(S)
    scale0 <- if (all(is.finite(d)) && any(d > 0)) {
      median(d[d > 0])
    } else {
      1
    }

    cholS <- try(chol(S), silent = TRUE)
    if (!inherits(cholS, "try-error")) {
      return(list(chol = cholS, eps = 0))
    }
    for (mult in c(1e-8, 1e-6, 1e-4)) {
      eps <- mult * scale0
      cholS <- try(chol(S + diag(eps, p)), silent = TRUE)
      if (!inherits(cholS, "try-error")) return(list(chol = cholS, eps = eps))
    }
    NULL
  }

  mvn_logdens_sum_from_chol <- function(r, cholS) {
    n <- nrow(r)
    p <- ncol(r)
    v <- forwardsolve(t(cholS), t(r))
    quad <- sum(v * v)
    logdet <- 2 * sum(log(diag(cholS)))
    -0.5 * (n * (p * log(2*pi) + logdet) + quad)
  }

  mvreg_xa <- function(z, x = NULL, ini = 0, missres){
    mvreglik <- function (z, X, missres, beta, alpha){
      colf <- missres$colf
      dp <- missres$dp
      ndp <- missres$ndp
      tt <- ncol(z)
      tt2 <- 2 ^ tt - 1
      Sgm <- SgmC(alpha, tt)
      ln <- 0
      for  (ti in seq_len(tt2)){
        if (ndp[ti] != 0){
          obs_cols <- colf[[ti]]
          rti <- z[dp == ti, obs_cols, drop = FALSE] -
            X[dp == ti, , drop = FALSE] %*% t(beta[obs_cols, , drop = FALSE])

          S <- Sgm[obs_cols, obs_cols, drop = FALSE]
          chol_sub <- chol_spd_ridge(S)
          if (is.null(chol_sub)) return(Inf)

          ll <- mvn_logdens_sum_from_chol(rti, chol_sub$chol)
          ln <- ln + ll
        }
      }
      return(-ln)  # negative log-likelihood
    }

    dmvreglik <- function (z, X, missres, beta, alpha){
      colf <- missres$colf
      dp <- missres$dp
      ndp <- missres$ndp
      nb <- length(c(beta))
      na <- length(alpha)
      Kp1 <- ncol(beta)
      tt <- ncol(z)
      tt2 <- 2 ^ tt - 1
      dlnb <- numeric(nb)
      dlns <- numeric(na)
      ddlnb <- matrix(0, nb, nb)
      ddlns <- matrix(0, na, na)
      ddlnbs <- matrix(0, nb, na)
      Sgm <- SgmC(alpha, tt)
      idx <- alpha_index(tt)  # na x 2
      dS <- vector("list", na)

      for (l1 in seq_len(na)) {
        i <- idx[l1, 1]
        j <- idx[l1, 2]
        M <- matrix(0, tt, tt)

        if (i == j) {
          M[i, i] <- 1
        } else {
          M[i, j] <- 1
          M[j, i] <- 1
        }
        dS[[l1]] <- M
      }

      for  (ti in seq_len(tt2)){
        if (ndp[ti]!=0){
          rti <- z[dp == ti, colf[[ti]], drop = FALSE] -
            X[dp == ti, , drop = FALSE] %*% t(beta[colf[[ti]], , drop = FALSE])
          nt <- nrow(rti)
          S <- Sgm[colf[[ti]], colf[[ti]], drop = FALSE]

          chol_sub <- chol_spd_ridge(S)
          if (is.null(chol_sub)) {
            # fail the whole derivative evaluation so outer step is rejected
            return(list(rep(Inf, nb + na), matrix(Inf, nb + na, nb + na)))
          }
          cholS <- chol_sub$chol
          iS <- chol2inv(cholS)
          sSr <- backsolve(cholS, t(rti), transpose = TRUE)
          sxb <- list()
          xjb <- list()
          bloc1 <- rep(seq_len(tt), times = Kp1)
          bloc2 <- rep(seq_len(Kp1), each  = tt)
          for (jb in seq_len(nb)){
            if (sum(colf[[ti]] == bloc1[jb]) == 1) {
              xti <- matrix(0, nt, tt)
              xti[, bloc1[jb]] <- X[dp == ti, bloc2[jb]]
              xti <- xti[, colf[[ti]], drop = FALSE]
              xjb[[jb]] <- xti
              sxb0 <- backsolve(cholS, t(xti), transpose = TRUE)
              dlnb[jb] <- dlnb[jb] + sum(sxb0 * sSr)
              sxb[[jb]] <- sxb0
            } else {
              sxb[[jb]] <- NA
            }
          }
          sxs <- list()
          rAA <- list()
          AA <- list()
          dSti <- list()
          seqi <- list()
          for (l in seq_len(na)){
            dS0 <- dS[[l]][colf[[ti]], colf[[ti]], drop = FALSE]
            pobs <- length(colf[[ti]])

            if (all(dS0 == 0)) {
              dSti[[l]] <- dS0
              AA[[l]]   <- matrix(0, pobs, pobs)
              seqi[[l]] <- seq_len(pobs)
              rAA[[l]]  <- matrix(0, pobs, pobs)   # これ重要（後段の zA を 0 にする）
              sxs[[l]]  <- matrix(0, pobs, nt)     # zA と同じ次元にしておく
              next
            }
            dSti[[l]] <- dS0
            AA0 <- -iS %*% dS0 %*% iS
            AA[[l]] <- AA0
            num <- seq_len(nrow(AA0))
            indep_num <- c()
            for(rk in seq_len(qr(AA0)$rank)){
              for(s in num){
                num0 <- c(indep_num, s)
                rk0 <- qr(AA0[, num0])$rank
                if(rk0 == rk){
                  indep_num <- num0
                  break
                }
              }
            }
            dep_num <- num[-indep_num]

            rti0 <- rti
            seql <- c(indep_num,dep_num)
            if (length(dep_num) != 0) {
              AA0 <- AA0[seql, seql, drop = FALSE]
              rti0 <- rti[, seql, drop = FALSE]
            }

            qr0 <- qr(AA0)
            seqi[[l]] <- seql
            rAA[[l]] <- qr.R(qr0)
            qx0 <- t(qr.Q(qr0)) %*% t(rti0)
            rx0 <- qr.R(qr0) %*% t(rti0)
            dlns[l] <- dlns[l]  -0.5 * nt * tr(iS %*% dS0) -
              0.5 * sum(qx0 * rx0)
            sxs[[l]] <- qx0
          }
          for (j1 in seq_len(nb)){
            for (j2 in j1:nb){
              if (!is.na(sxb[[j1]][1]) & !is.na(sxb[[j2]][1])){
                ddlnb[j1,j2] <- ddlnb[j1,j2] - sum(sxb[[j1]] * sxb[[j2]])
              }
            }
          }
          for (l1 in seq_len(na)){
            for (l2 in l1:na){
              if (!is.na(sxs[[l1]][1]) & !is.na(sxs[[l2]][1])){
                AA2 <- 2 * iS %*% (dSti[[l1]] %*% iS %*% dSti[[l2]]) %*% iS
                rkAA2 <- qr(AA2)$rank
                if (rkAA2 == 0L) {
                  next
                }
                num <- seq_len(nrow(AA2))
                indep_num <- c()
                for(rk in seq_len(rkAA2)){
                  for(s in num){
                    num0 <- c(indep_num,s)
                    rk0 <- qr(AA2[,num0])$rank
                    if(rk0 == rk){
                      indep_num <- num0
                      break
                    }
                  }
                }
                dep_num <- num[-indep_num]

                rti0 <- rti
                seql <- c(indep_num,dep_num)
                if (length(dep_num) != 0) {
                  AA2 <- AA2[seql, seql, drop = FALSE]
                  rti0 <- rti[, seql, drop = FALSE]
                }
                qr0 <- qr(AA2)
                qx0 <- t(qr.Q(qr0)) %*% t(rti0)
                rx0 <- qr.R(qr0) %*% t(rti0)
                ddlns[l1, l2] <- ddlns[l1, l2] -
                  0.5 * nt * tr(AA[[l1]] %*% dSti[[l2]]) -
                  0.5 * sum(qx0*rx0)
              }
            }
          }
          for (j in seq_len(nb)){
            for (l in seq_len(na)){
              if (!is.na(sxb[[j]][1]) & !is.na(sxs[[l]][1])){
                zA <- rAA[[l]] %*% t(xjb[[j]][, seqi[[l]], drop = FALSE])
                ddlnbs[j, l] <- ddlnbs[j, l] + sum(zA * sxs[[l]])
              }
            }
          }
        }
      }
      for (j1 in seq_len(nb)){
        for (j2 in seq_len(j1)){
          ddlnb[j1, j2] <- ddlnb[j2, j1]
        }
      }
      for (l1 in seq_len(na)){
        for (l2 in seq_len(l1)){
          ddlns[l1, l2] <- ddlns[l2, l1]
        }
      }
      dln <- -c(dlnb, dlns)
      ddln <- rbind(cbind(ddlnb, ddlnbs), cbind(t(ddlnbs), ddlns))
      ddln <- -ddln
      return(list(dln,ddln))
    }
    if (ini == 1){
      alpc <- 1e-9
      alpb <- 0.1
      maxt <- 1000
      maxc <- 1000
    } else {
      alpc <- 1e-7
      alpb <- 0.2
      maxt <- 100
      maxc <- 100
    }
    n <- nrow(z)
    tt <- ncol(z)
    if (is.null(x)) {
      X <- matrix(1, nrow(z), 1)
      res0 <- lm(z ~ 1, na.action = na.exclude)
      beta0 <- as.matrix(t(res0$coefficients))
      Sgm <- cov(res0$residuals)
    } else {
      X <- cbind(1, x)
      res0 <- lm(z ~ x, na.action = na.exclude)
      beta0 <- as.matrix(t(res0$coefficients))
      if (anyNA(beta0)) {
        for (j in seq_len(tt)) {
          if (anyNA(beta0[j, ])) {
            fitj <- lm(z[, j] ~ x, na.action = na.exclude)
            bj <- fitj$coefficients
            bj[is.na(bj)] <- 0

            na_pos <- is.na(beta0[j, ])
            beta0[j, na_pos] <- bj[na_pos]
          }
        }
        r0 <- z - X %*% t(beta0)
        Sgm <- cov(r0, use = "pairwise.complete.obs")
      } else {
        Sgm <- cov(res0$residuals)
      }
    }
    alpha0 <- Sgm[upper.tri(Sgm, diag = TRUE)]
    nb <- length(beta0)
    na <- length(alpha0)
    nbs <- nb + na
    crt <- 1e-7
    ln0 <- mvreglik(z, X, missres, beta0, alpha0)
    dln0 <- dmvreglik(z, X, missres, beta0, alpha0)
    dk <- dln0[[1]]
    Hk <- dln0[[2]]
    if (min(eigen(Hk)$values)<0) {
      cntdg <- 0
      Hkdg <- abs(diag(diag(Hk))) * 0.01
      while (min(eigen(Hk)$values) < 0 & cntdg < 100) {
        Hk <- Hk + Hkdg
        cntdg <- cntdg + 1
      }
    }
    dd <- try(qr.solve(Hk, dk), silent = TRUE)
    if (inherits(dd, "try-error")){
      dd <- dk
    }
    dhd <- t(dk) %*% dd
    flg <- 1
    theta <- c(beta0, alpha0)
    aa <- c()
    count <- 0
    dhd <- 1000
    time1 <- proc.time()
    time2 <- time1
    while (flg == 1){
      alp <- 1
      flgc <- 1
      count <- count + 1
      while (flgc == 1){
        thetap <- theta - alp * dd
        beta1 <- matrix(thetap[seq_len(nb)], tt, ncol(X))
        alpha1 <- thetap[(nb + 1):nbs]
        Sgm <- SgmC(alpha1, tt)
        egn <- eigen(Sgm)$values
        if (sum(egn < 0) > 0) {
          alp <- alp * alpb
        } else {
          ln1 <- mvreglik(z, X, missres, beta1, alpha1)
          if (!is.finite(ln1) || ln1 > ln0 || ln1 > ln0 + 0.45 * alp * dhd) {
            alp <- alp * alpb
          } else {
            flgc <- 0
          }
        }
        if (alp < alpc){
          flgc <- 0
        }
      }
      if (alp > alpc){
        dln0 <- dmvreglik(z, X, missres, beta1, alpha1)
        dk <- dln0[[1]]
        Hk <- dln0[[2]]
        cntdg <- 0
        if (min(eigen(Hk)$values) < 0) {
          Hkdg <- matrix(1, nbs, nbs) + diag(sign(diag(Hk))) * 0.01
          while(min(eigen(Hk)$values) < 0 & cntdg < 100){
            Hk <- Hk * Hkdg
            cntdg <- cntdg + 1
          }
        }
        dd <- try(qr.solve(Hk, dk), silent = TRUE)
        if (inherits(dd, "try-error")){
          dd <- dk
        }
        dhd <- t(dk) %*% dd
        time2 <- proc.time() - time1
        if (dhd / abs(ln1) < crt | count > maxc | alp < alpc |
            time2[3] > maxt) {
          flg <- 0
        }
        theta <- thetap
        beta0 <- beta1
        alpha0 <- alpha1
        ln0 <- ln1
      } else {
        flg <- 0
      }
    }
    if (count > maxc | (alp < alpc & dhd / abs(ln1) > crt * 10) |
        time2[3] > maxt) {
      ln0 <- Inf
    }
    return(list(beta = beta0, alpha = alpha0, lik = -ln0, dl = -dk, ddl = -Hk))
  }

  lik_ldif_x <- function(lambda, y, x = NULL, ini = 0, missres){
    z <- c()
    tt <- ncol(y)
    yl <- 0
    n <- nrow(y)
    for (i in seq_len(tt)){
      z <- cbind(z, (y[, i] ^ lambda[i] - 1) / lambda[i])
      yl <- yl + sum((lambda[i] - 1) * log(y[, i]), na.rm = TRUE)
    }
    lik0 <- try(mvreg_xa(z, x, ini, missres)$lik, silent = TRUE)
    if (inherits(lik0, "try-error")){
      lik0 <- -Inf
    }
    lik <- lik0 + yl
    if (lik < -1e+9) lik <- -Inf
    return(-lik)
  }

  dl <- function(lambda, y, x = NULL, ln = NULL, ini = 0, missres) {
    tt <- ncol(y)
    g <- numeric(tt)

    for (i in seq_len(tt)) {
      dd <- 1e-3
      ddi <- numeric(tt)
      ddi[i] <- dd

      f1 <- lik_ldif_x(lambda + ddi, y, x, ini = ini, missres = missres)
      f2 <- lik_ldif_x(lambda - ddi, y, x, ini = ini, missres = missres)

      if (!is.finite(f1) || !is.finite(f2)) {
        f0 <- if (is.null(ln)) {
          lik_ldif_x(lambda, y, x, ini = ini, missres = missres)
        } else {
          ln
        }
        f1f <- lik_ldif_x(lambda + ddi, y, x, ini = ini, missres = missres)
        g[i] <- (f1f - f0) / dd
      } else {
        g[i] <- (f1 - f2) / (2 * dd)
      }
    }
    return(g)
  }
  qnr <- function(lambda, y, x = NULL, missres){
    tt <- ncol(y)
    Bk <- diag(1, tt)
    ln0 <- lik_ldif_x(lambda, y, x, ini = 1, missres)
    dk <- dl(lambda, y, x, ini = 1, missres = missres)
    flg <- 1
    count <- 0
    cnt2 <- 0
    repeat{
      alp <- 1
      flgc <- 1
      repeat{
        qrschk <- try(qr.solve(Bk, dk), silent = TRUE)
        if (inherits(qrschk, "try-error")) {
          lambdap <- lambda - dk * alp
        } else {
          lambdap <- lambda - as.numeric(qrschk) * alp
        }

        ln1 <- try(lik_ldif_x(lambdap, y, x, missres = missres), silent = TRUE)

        if (inherits(ln1, "try-error") || !is.finite(ln1)) {
          ln1 <- Inf
        }
        if (ln0 < ln1) {
          alp <- alp * 0.2
        } else {
          flgc <- 0
        }
        if (flgc == 0 | alp < 1e-7) break
      }
      sk <- lambdap - lambda
      Bksk <- Bk %*% sk

      dkp <- dl(lambdap, y = y, x = x, missres = missres)
      yk <- dkp - dk
      Bk <- Bk - (Bksk %*% t(Bksk)) / (sum(Bksk * sk)) + (yk %*% t(yk)) /
        (sum(sk * yk))
      if (!is.finite(ln1) || !is.finite(sum(dkp ^ 2)) ||
          (sum(dkp ^ 2) / abs(ln1)) < 1e-7) {
        flg <- 0
      }
      if (alp < 1e-5 & cnt2 > 0) {
        flg <- 0
      }
      if (alp < 1e-5 & cnt2 == 0) {
        cnt2 <- cnt2 + 1
        Bk <- diag(1, tt)
      }
      if (alp > 1e-5) {
        cnt2 <- 0
      }
      lambda <- lambdap
      dk <- dkp
      ln0 <- ln1
      count <- count + 1
      if (count > 50) {
        flg <- 0
      }
      if (flg==0) {
        break
      }
    }
    if ((alp < 1e-7 & sum(dkp ^ 2) / abs(ln1) > 1e-4) | count > 50) {
      err <- TRUE
    } else {
      err <- FALSE
    }
    return(list(lambda = lambda, lik = ln0, err = err))
  }
  if (!is.null(x)){
    x <- as.matrix(x)
  }
  tt <- ncol(y)
  lambda <- numeric(tt)
  for (i in seq_len(tt))  {
    if (!is.null(x)){
      dft <- data.frame(y = y[, i])
      formini <- formula(paste0("y ~ ", paste(paste0("x_", seq_len(ncol(x))),
                                              collapse = " + ")))
      for (k in seq_len(ncol(x))) {
        dft[[paste0("x_", k)]] <- x[, k]
      }
      lambda[i] <- bcmixed::bcreg(formini, dft)$lambda
    } else {
      lambda[i] <- bcmixed::bct.v(y[!is.na(y[, i]), i])$lambda
    }
  }
  missres <- misscalc(y)
  res <- try(qnr(lambda, y, x, missres = missres), silent = FALSE)
  if (inherits(res, "try-error")){
    lambda <- NA
    beta <- NA
    alpha <- NA
    lik <- NA
    err <- TRUE
    rsd <- NA
    Sgm <- NA
    z <- NA
  } else {
    lambda <- res$lambda
    lik <- res$lik
    err <- res$err
    z <- c()
    for (i in seq_len(tt)){
      z <- cbind(z, (y[, i] ^ lambda[i] - 1) / lambda[i])
    }
    n <- nrow(z)
    resbs <- mvreg_xa(z, x, missres = missres)
    beta <- resbs$beta
    alpha <- resbs$alpha
    Sgm <- SgmC(alpha, tt)
  }
  return(list(lambda = lambda, beta = beta, alpha = alpha, Sigma = Sgm,
              lik = -lik, err = err))
}


# Asymptotic variance-covariance matrix estimation for MLEs of parameters for
# the BCMVR model
bcmvr_inf_g <- function (y, x = NULL, lambda, beta, alpha){
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
  if (any(y[!is.na(y)] <= 0)) {
    stop("y must be positive for Box-Cox transform.")
  }
  if (!is.null(x) && anyNA(x)) {
    stop("Missing value is included in x.")
  }
  tr <- function(X) {
    sum(diag(X))
  }
  dec2bin <- function(num, digit = 0){
    if(num <= 0 && digit <= 0){
      return(NULL)
    }else{
      return(append(Recall(num %/% 2,digit - 1), num %% 2))
    }
  }
  SgmC <- function(alpha, tt){
    Sgm0 <- matrix(0, tt, tt)
    Sgm0[upper.tri(Sgm0, diag = TRUE)] <- alpha
    Sgm <- Sgm0 + t(Sgm0)
    diag(Sgm) <- diag(Sgm) / 2
    return(Sgm)
  }
  alpha_index <- function(tt) {
    which(upper.tri(matrix(0, tt, tt), diag = TRUE), arr.ind = TRUE)
  }
  n <- nrow(y)
  tt <- ncol(y)
  if (length(lambda) != tt) {
    stop("length(lambda) must equal ncol(y).")
  }
  if (nrow(beta) != tt) {
    stop("nrow(beta) must equal ncol(y).")
  }
  if (is.null(x)) {
    X <- matrix(1, n, 1)
  } else {
    X <- cbind(1, as.matrix(x))
  }
  if (ncol(beta) != ncol(X)) {
    stop("ncol(beta) must equal ncol(cbind(1,x)).")
  }
  Kp1 <- ncol(X)
  nb <- length(beta)
  na <- length(alpha)
  na0 <- tt * (tt + 1) / 2
  if (na != na0) {
    stop("length(alpha) must be tt*(tt+1)/2.")
  }
  dp <- numeric(n)
  tt2 <- 2 ^ tt - 1
  nba <- nb + na
  nsf <- c()
  maxc <- 1000
  for (j1 in seq_len(tt)){
    for (j2 in seq_len(j1)){
      nsf <- rbind(nsf, t(c(j1, j2)))
    }
  }
  colf <- list()
  for (j in seq_len(tt2)){
    mv <- dec2bin(j - 1, tt)
    dp[colSums(apply(is.na(y), 1, '==', mv)) == tt] <- j
    colf[[j]] <- (seq_len(tt))[mv == 0]
  }
  Ncc <- sum(colSums(apply(is.na(y), 1, '==', dec2bin(0, tt))) == tt)
  ndp <- numeric(tt2)
  for (j in seq_len(tt2)) {
    ndp[j] <- sum(dp == j)
  }
  dlnl1 <- numeric(tt)
  dlnl2 <- numeric(tt)
  dlnl <- numeric(tt)
  dlnb <- numeric(nb)
  dlna <- numeric(na)
  Hl <- matrix(0, tt, tt)
  Hb <- matrix(0, nb, nb)
  Ha <- matrix(0, na, na)
  Hlb <- matrix(0, tt, nb)
  Hla <- matrix(0, tt, na)
  Hba <- matrix(0, nb, na)
  Jl <- matrix(0, tt, tt)
  Jb <- matrix(0, nb, nb)
  Ja <- matrix(0, na, na)
  Jlb <- matrix(0, tt, nb)
  Jla <- matrix(0, tt, na)
  Jba <- matrix(0, nb, na)
  z <- y
  dz <- y
  ddz <- y
  ly <- log(y)
  lly <- y
  Sgm <- SgmC(alpha, tt)
  # symmetry check (should be exact, but allow tiny noise)
  asym <- max(abs(Sgm - t(Sgm)))
  if (asym > 1e-10 * max(1, max(abs(Sgm)))) {
    stop("Sgm is not symmetric (beyond numerical tolerance).")
  }

  d <- diag(Sgm)
  if (any(!is.finite(d)) || any(d <= 0)) {
    stop("diag(Sgm) must be positive.")
  }

  # SPD check
  ch <- try(chol(Sgm), silent = TRUE)
  if (inherits(ch, "try-error")) {
    stop("Sgm is not positive definite.")
  }

  idx <- alpha_index(tt)  # na x 2
  dS <- vector("list", na)

  for (l1 in seq_len(na)) {
    i <- idx[l1, 1]
    j <- idx[l1, 2]
    M <- matrix(0, tt, tt)

    if (i == j) {
      M[i, i] <- 1
    } else {
      M[i, j] <- 1
      M[j, i] <- 1
    }
    dS[[l1]] <- M
  }

  for (jl in seq_len(tt)){
    yjl <- y[, jl, drop = FALSE]
    ljl <- lambda[jl]
    lyjl <- ly[, jl, drop = FALSE]

    yljl <- yjl ^ ljl
    z[, jl] <- (yljl - 1) / ljl
    dz[, jl] <- (yljl * (ljl * lyjl - 1) + 1) / ljl ^ 2
    ddz[, jl] <- 1 / ljl * yljl * lyjl ^ 2 - 2 / ljl ^ 2 * yljl * lyjl +
      2 / ljl ^ 3 * (yljl - 1)
    lly[, jl] <- (ljl - 1) * lyjl
  }
  # ln <- 0
  for  (j in seq_len(tt2)){
    if (ndp[j] != 0){
      zj <- z[dp == j, colf[[j]], drop = FALSE]
      rj <- zj - X[dp == j, , drop = FALSE] %*%
        t(beta[colf[[j]], , drop = FALSE])
      nj <- ndp[j]
      lyj <- ly[dp == j, , drop = FALSE]

      S <- Sgm[colf[[j]], colf[[j]], drop = FALSE]
      chS <- try(chol(S), silent = TRUE)
      if (inherits(chS, "try-error")) {
        stop("Submatrix of Sgm is not PD for an observed pattern.")
      }
      cholS <- chS
      iS <- try(chol2inv(cholS), silent = TRUE)
      if (inherits(iS, "try-error")) {
        iS <- try(qr.solve(S, diag(nrow(S))), silent = TRUE)
        if (inherits(iS, "try-error")) {
          iS <- MASS::ginv(S)
        }
      }
      sSr <- backsolve(cholS, t(rj), transpose = TRUE)

      sxb <- list()
      sdz <- list()
      sddz <- list()
      xjb <- list()
      dzjl <- list()
      cs_zSr <- list()
      cs_xSr <- list()
      for (jl in seq_len(tt)){
        if (sum(colf[[j]] == jl) == 1){
          dzj <- matrix(0, nj, tt)
          dzj[, jl] <- dz[dp == j, jl, drop = FALSE]
          dzj <- dzj[, colf[[j]], drop = FALSE]
          dzjl[[jl]] <- dzj
          ddzj <- matrix(0, nj, tt)
          ddzj[, jl] <- ddz[dp == j, jl, drop = FALSE]
          ddzj <- ddzj[, colf[[j]], drop = FALSE]
          sdz[[jl]] <- backsolve(cholS, t(dzj), transpose = TRUE)
          sddz[[jl]] <- backsolve(cholS, t(ddzj), transpose = TRUE)
          cs_zSr[[jl]] <- colSums(sdz[[jl]] * sSr)
          dlnl[jl] <- dlnl[jl] + sum(lyj[, jl]) - sum(cs_zSr[[jl]])
        } else {
          dzjl[[jl]] <- NA
        }
      }
      bloc1 <- rep(seq_len(tt), times = Kp1)
      bloc2 <- rep(seq_len(Kp1), each  = tt)
      for (jb in seq_len(nb)){
        if (sum(colf[[j]] == bloc1[jb]) == 1) {
          xj <- matrix(0, nj, tt)
          xj[, bloc1[jb]] <- X[dp == j, bloc2[jb]]
          xj <- xj[, colf[[j]], drop = FALSE]
          xjb[[jb]] <- xj
          sxb0 <- backsolve(cholS, t(xj), transpose = TRUE)
          dlnb[jb] <- dlnb[jb] + sum(sxb0 * sSr)
          sxb[[jb]] <- sxb0
          cs_xSr[[jb]] <- colSums(sxb0 * sSr)
        } else {
          sxb[[jb]] <- NA
        }
      }
      qAr <- list()
      cs_rAr <- list()
      rAA <- list()
      AA <- list()
      tris <- numeric(na)
      sxs <- list()
      rAA <- list()
      AA <- list()
      dSj <- list()
      for (ja in seq_len(na)){
        j1 <- nsf[ja, 1]
        j2 <- nsf[ja, 2]
        if (sum(colf[[j]] == j1) == 1 & sum(colf[[j]] == j2) == 1){
          dS0 <- dS[[ja]][colf[[j]], colf[[j]], drop = FALSE]
          dSj[[ja]] <- dS0
          AA0 <- -iS %*% dS0 %*% iS
          AA[[ja]] <- AA0
          num <- seq_len(nrow(AA0))
          indep_num <- c()
          for(rk in seq_len(qr(AA0)$rank)){
            for(s in num){
              num0 <- c(indep_num, s)
              rk0 <- qr(AA0[, num0])$rank
              if(rk0 == rk){
                indep_num <- num0
                break
              }
            }
          }
          dep_num <- num[-indep_num]
          seql <- c(indep_num,dep_num)
          qr0 <- qr(AA0)
          rAA[[ja]] <- qr.R(qr0)
          tris[ja] <- tr(iS %*% dS0)
          qx0 <- t(qr.Q(qr0)) %*% t(rj)
          rx0 <- qr.R(qr0) %*% t(rj)
          qAr[[ja]] <- qx0
          cs_rAr[[ja]] <- colSums(qx0 * rx0)
          dlna[ja] <- dlna[ja]  -0.5 * nj * tris[ja]  - 0.5 * sum(qx0 * rx0)
        } else {
          qAr[[ja]] <- NA
        }
      }
      for (jl1 in seq_len(tt)){
        for (jl2 in jl1:tt){
          if (!is.na(dzjl[[jl1]][1]) & !is.na(dzjl[[jl2]][1])){
            if (jl1 == jl2) {
              Hl[jl1, jl2] <- Hl[jl1, jl2] - sum(sddz[[jl1]] * sSr) -
                sum(sdz[[jl1]] ^ 2)
            } else {
              Hl[jl1, jl2] <- Hl[jl1, jl2] - sum(sdz[[jl1]] * sdz[[jl2]])
            }
            Jl[jl1, jl2] <- Jl[jl1, jl2] + sum(lyj[, jl1]*lyj[, jl2]) -
              sum(lyj[, jl1] * cs_zSr[[jl2]]) -
              sum(lyj[, jl2] * cs_zSr[[jl1]]) +
              sum(cs_zSr[[jl1]] * cs_zSr[[jl2]])
          }
        }
      }

      for (jb1 in seq_len(nb)){
        for (jb2 in jb1:nb){
          if (!is.na(sxb[[jb1]][1]) & !is.na(sxb[[jb2]][1])){
            Hb[jb1, jb2] <- Hb[jb1, jb2] - sum(sxb[[jb1]] * sxb[[jb2]])
            Jb[jb1, jb2] <- Jb[jb1, jb2] + sum(cs_xSr[[jb1]] * cs_xSr[[jb2]])
          }
        }
      }

      for (ja1 in seq_len(na)){
        for (ja2 in ja1:na){
          if (!is.na(qAr[[ja1]][1]) & !is.na(qAr[[ja2]][1])){
            AA2 <- iS %*% (dSj[[ja1]] %*% iS %*% dSj[[ja2]] +
                             dSj[[ja2]] %*% iS %*% dSj[[ja1]]) %*% iS
            qr0 <- qr(AA2)
            qx0 <- t(qr.Q(qr0)) %*% t(rj)
            rx0 <- qr.R(qr0) %*% t(rj)
            Ha[ja1, ja2] <- Ha[ja1, ja2] - 1 / 2 * nj *
              tr(AA[[ja1]] %*% dSj[[ja2]]) - 1 / 2 * sum(qx0 * rx0)
            Ja[ja1, ja2] <- Ja[ja1, ja2] + 0.25 * nj * tris[ja1] * tris[ja2] +
              0.25 * sum(tris[ja1] * cs_rAr[[ja2]]) +
              0.25 * sum(tris[ja2] * cs_rAr[[ja1]]) +
              0.25 * sum(cs_rAr[[ja1]] * cs_rAr[[ja2]])
          }
        }
      }

      for (jl in seq_len(tt)){
        for (jb in seq_len(nb)){
          if (!is.na(dzjl[[jl]][1]) & !is.na(sxb[[jb]][1])){
            Hlb[jl, jb] <- Hlb[jl, jb] + sum(sdz[[jl]] * sxb[[jb]])
            Jlb[jl, jb] <- Jlb[jl, jb] + sum(lyj[, jl] * cs_xSr[[jb]]) -
              sum(cs_zSr[[jl]] * cs_xSr[[jb]])
          }
        }
      }

      for (jl in seq_len(tt)){
        for (ja in seq_len(na)){
          if (!is.na(dzjl[[jl]][1]) & !is.na(qAr[[ja]][1])){
            dzrA <- rAA[[ja]] %*% t(dzjl[[jl]])
            Hla[jl, ja] <- Hla[jl, ja] - sum(dzrA * qAr[[ja]])
            Jla[jl, ja] <- Jla[jl, ja] - 0.5 * sum(lyj[, jl] * tris[ja]) -
              0.5 * sum(lyj[, jl] * cs_rAr[[ja]]) +
              0.5 * sum(cs_zSr[[jl]] * tris[ja]) +
              0.5 * sum(cs_zSr[[jl]] * cs_rAr[[ja]])
          }
        }
      }

      for (jb in seq_len(nb)){
        for (ja in seq_len(na)){
          if (!is.na(sxb[[jb]][1]) & !is.na(qAr[[ja]][1])){
            xrA <- rAA[[ja]] %*% t(xjb[[jb]])
            Hba[jb, ja] <- Hba[jb, ja] + sum(xrA * qAr[[ja]])
            Jba[jb, ja] <- Jba[jb, ja] - 0.5 * sum(tris[ja] * cs_xSr[[jb]]) -
              0.5 * sum(cs_xSr[[jb]] * cs_rAr[[ja]])
          }
        }
      }
    }
  }
  for (jl1 in seq_len(tt)){
    for (jl2 in seq_len(jl1)){
      Hl[jl1, jl2] <- Hl[jl2, jl1]
      Jl[jl1, jl2] <- Jl[jl2, jl1]
    }
  }
  for (jb1 in seq_len(nb)){
    for (jb2 in seq_len(jb1)){
      Hb[jb1, jb2] <- Hb[jb2, jb1]
      Jb[jb1, jb2] <- Jb[jb2, jb1]
    }
  }
  for (ja1 in seq_len(na)){
    for (ja2 in seq_len(ja1)){
      Ha[ja1, ja2] <- Ha[ja2, ja1]
      Ja[ja1, ja2] <- Ja[ja2, ja1]
    }
  }
  H <- rbind(cbind(Hl, Hlb, Hla), cbind(t(Hlb), Hb, Hba),
             cbind(t(Hla), t(Hba), Ha))
  J <- rbind(cbind(Jl, Jlb, Jla), cbind(t(Jlb), Jb, Jba),
             cbind(t(Jla), t(Jba), Ja))
  Vm <- try(qr.solve(-H, diag(nrow(H))), silent = TRUE)
  if (inherits(Vm, "try-error")) {
    Vm <- MASS::ginv(-H)
  }
  Vr <- Vm %*% J %*% Vm
  return(list(V.mod = Vm, V.rob = Vr))
}
