#' Treatment effect inference using a Box--Cox multivariate regression (BCMVR)
#' model for longitudinal data in randomized controlled trials
#'
#' Fit Box--Cox multivariate regression models separately for each group and
#' perform prespecified group comparisons based on (i) median differences and
#' (ii) probability-based treatment effect measures.
#'
#' This function is the main user-facing interface of the \pkg{bcmvr} package.
#' Internally, it fits group-specific BCMVR models, estimates variance--
#' covariance matrices, and computes inferential summaries for:
#' \itemize{
#'   \item group-specific marginal medians (with Wald-type confidence
#'     intervals),
#'   \item median differences between prespecified group pairs (with CI and
#'     p-values),
#'   \item probability-based treatment effect measures comparing full
#'     marginal distributions.
#' }
#'
#' @param data A data frame in long format.
#' @param outcome Name of the outcome variable (must be positive).
#' @param id Name of the subject identifier.
#' @param group Name of the group variable.
#' @param time Name of the time variable.
#' @param cov Optional character vector of baseline covariate names
#'   (time-invariant within subject). Covariates that are constant within any
#'   group are dropped from the analysis with a warning, to ensure comparability
#'   of covariate-adjusted effects across groups.
#'
#' @param pairs A list of length-2 vectors specifying ordered group pairs
#'   \code{c(g1, g2)} for comparison, where \code{g1} is the test group and
#'   \code{g2} is the control group. The order matters. Each element must be
#'   either group labels (as they appear in \code{group}) or their positions in
#'   the sorted unique group labels.
#'
#' @param smaller Logical. Must be specified. \code{TRUE} if a smaller outcome
#'   value is clinically favorable, \code{FALSE} otherwise. This argument
#'   determines the direction of the probability-based effect measure.
#'
#' @param t_eval Optional integer vector specifying time points at which
#'   effects are evaluated. If \code{NULL}, all time points are used.
#'
#' @param robust Logical. If \code{TRUE} (default), robust sandwich variance
#'   estimators are used for inference. If \code{FALSE}, model-based variance
#'   estimators are used.
#'
#' @param level Confidence level for Wald-type confidence intervals.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{call}{The matched function call.}
#'   \item{datalist}{A list returned by the internal data reshaping routine,
#'     including group-wise outcome matrices and covariate matrices.}
#'   \item{glst}{Sorted unique group labels used in the analysis.}
#'   \item{tt}{Number of time points (columns of the wide outcome matrix).}
#'   \item{t_eval}{Time points at which effects are evaluated.}
#'   \item{xmean}{Pooled covariate means across all subjects (or \code{NULL} if
#'     \code{cov = NULL}).}
#'   \item{pairs}{A data frame of ordered group pairs actually used for
#'     comparisons.}
#'   \item{fit}{A list of fitted group-specific BCMVR model objects.}
#'   \item{vcov}{A list of variance--covariance estimates (model-based and
#'     robust) for each group.}
#'   \item{median}{A data frame of group-specific marginal medians with
#'     standard errors and Wald-type confidence intervals.}
#'   \item{meddiff}{A data frame of median differences for specified ordered
#'     pairs, computed as \code{median(g1) - median(g2)}, with standard errors,
#'     confidence intervals, and Wald-type p-values under independent groups.}
#'   \item{prob}{A data frame of probability-based effect measures. It reports
#'     \eqn{P(Y_{g1,t} < Y_{g2,t})} on the original outcome scale, with
#'     inference performed on the logit scale. Under \code{smaller = TRUE},
#'     this can be interpreted as the probability that the test group
#'     (\code{g1}) has a more favorable outcome than the control group
#'     (\code{g2}).}
#' }
#'
#' @examples
#' data(aidscd4, package = "bcmixed")
#' set.seed(1)
#'
#' # sample ~140 subjects (IDs), keep all rows for those IDs
#' id_keep <- sample(unique(aidscd4$id), size = 140, replace = FALSE)
#' dat <- subset(aidscd4, id %in% id_keep)
#'
#' # (optional) make treatment a factor so pairs are readable/stable
#' dat$treatment <- factor(dat$treatment)
#'
#' # define ordered pairs: g1=test, g2=control (here: compare 2/3/4 vs 1)
#' pairs <- list(
#'   c("2", "1"),
#'   c("3", "1"),
#'   c("4", "1")
#'   )
#'
#' # CD4 is typically "larger is better" -> smaller = FALSE
#' res <- bcmvr(
#'   dat,
#'   outcome = cd4,
#'   id      = id,
#'   group   = treatment,
#'   time    = weekc,
#'   cov     = c("age", "cd4.bl"),
#'   pairs   = pairs,
#'   smaller = FALSE
#'   )
#'
#' res$median
#' res$meddiff
#' res$prob
#'
#' @references
#' Maruo K et al. Inference on treatment effects for longitudinal outcomes with
#' skewed distributions. Submitting (available from author).
#'
#' @importFrom stats cov dnorm formula integrate lm median na.exclude pnorm qnorm setNames
#'
#' @export
bcmvr <- function(data, outcome, id, group, time, cov = NULL,
                  pairs,
                  smaller,
                  t_eval = NULL,
                  robust = TRUE,
                  level = 0.95) {

  if (is.null(pairs)) {
    stop("pairs must be provided (direction matters).", call. = FALSE)
  }
  if (missing(smaller)) {
    stop("Argument 'smaller' must be specified: TRUE if smaller outcome is favorable, FALSE otherwise.",
         call. = FALSE)
  }
  # --- drop covariates that are constant in ANY group (global drop) ---
  drop_any_group_const <- function(xg_list, glst) {
    mats <- xg_list[!vapply(xg_list, is.null, logical(1))]
    if (length(mats) == 0) {
      return(list(keep = NULL, drop = NULL, by_group = NULL))
    }

    # column names must be consistent across groups
    cn <- colnames(mats[[1]])
    for (k in seq_along(mats)) {
      if (!identical(colnames(mats[[k]]), cn)) {
        stop("Covariate columns are not aligned across groups.", call. = FALSE)
      }
    }

    is_const_group <- matrix(FALSE, nrow = length(cn), ncol = length(mats),
                             dimnames = list(cn, glst[!vapply(xg_list, is.null,
                                                              logical(1))]))

    for (g in seq_along(mats)) {
      Xg <- mats[[g]]
      is_const_group[, g] <- vapply(seq_len(ncol(Xg)), function(j) {
        v <- Xg[, j]
        v <- v[is.finite(v)]
        (length(v) == 0) || (length(unique(v)) <= 1)
      }, logical(1))
    }

    drop <- rownames(is_const_group)[apply(is_const_group, 1, any)]
    keep <- setdiff(cn, drop)

    by_group <- lapply(seq_len(ncol(is_const_group)), function(j) {
      rownames(is_const_group)[is_const_group[, j]]
    })
    names(by_group) <- colnames(is_const_group)

    list(keep = keep, drop = drop, by_group = by_group)
  }

  datalist <- data_reshape(data, {{ outcome }}, {{ id }}, {{ group }},
                           {{ time }}, cov)

  glst <- datalist$glst

  chk <- drop_any_group_const(datalist$xg, glst)

  if (!is.null(chk$drop) && length(chk$drop) > 0) {
    # warn with details
    msg <- paste0(
      "Dropping covariate(s) that are constant within at least one group: ",
      paste(chk$drop, collapse = ", "), ".\n",
      "Constant-by-group: ",
      paste(
        vapply(names(chk$by_group), function(g) {
          if (length(chk$by_group[[g]]) == 0) return(paste0(g, ": none"))
          paste0(g, ": ", paste(chk$by_group[[g]], collapse = ", "))
        }, character(1)),
        collapse = " | "
      )
    )
    warning(msg, call. = FALSE)

    # drop consistently from every group's xg
    for (k in seq_along(datalist$xg)) {
      if (!is.null(datalist$xg[[k]]) && ncol(datalist$xg[[k]]) > 0) {
        datalist$xg[[k]] <- datalist$xg[[k]][, chk$keep, drop = FALSE]
      }
    }
  }

  # pooled xmean using remaining covariates
  xmean <- (function(xg_list) {
    mats <- xg_list[!vapply(xg_list, is.null, logical(1))]
    if (length(mats) == 0) return(NULL)
    Xall <- do.call(rbind, mats)
    if (ncol(Xall) == 0) return(NULL)
    colMeans(Xall, na.rm = TRUE)
  })(datalist$xg)

  G <- length(glst)
  if (G < 2) {
    stop("Need at least 2 groups.", call. = FALSE)
  }

  tt <- ncol(datalist$yg[[1]])
  if (is.null(t_eval)) t_eval <- seq_len(tt)
  if (any(!t_eval %in% seq_len(tt))) {
    stop("t_eval contains out-of-range values. Must be within 1:tt.",
         call. = FALSE)
  }

  zcrit <- stats::qnorm(1 - (1 - level)/2)

  # fit + inf per group
  thetag  <- vector("list", G)
  vthetag <- vector("list", G)
  names(thetag) <- names(vthetag) <- glst

  for (g0 in seq_len(G)) {
    thetag[[g0]] <- bcmvr_fit_g(datalist$yg[[g0]], datalist$xg[[g0]])
    vthetag[[g0]] <- bcmvr_inf_g(datalist$yg[[g0]], datalist$xg[[g0]],
                                 thetag[[g0]]$lambda, thetag[[g0]]$beta,
                                 thetag[[g0]]$alpha)
  }

  # map group label -> index
  g_index <- setNames(seq_len(G), glst)

  # normalize pairs input to data.frame(g1,g2)
  pair_df <- do.call(rbind, lapply(pairs, function(p) {
    if (length(p) != 2) {
      stop("Each element of pairs must have length 2.", call. = FALSE)
    }
    if (is.numeric(p)) {
      p <- glst[p]
    }
    data.frame(g1 = as.character(p[1]), g2 = as.character(p[2]),
               stringsAsFactors = FALSE)
  }))

  if (any(!pair_df$g1 %in% glst) || any(!pair_df$g2 %in% glst)) {
    stop("pairs contains unknown group labels. Available: ",
         paste(glst, collapse = ", "), call. = FALSE)
  }
  if (any(pair_df$g1 == pair_df$g2)) {
    stop("pairs contains g1 == g2, which is not allowed.", call. = FALSE)
  }

  # --- 1) median for each group x time (computed once) + CI ---
  median_rows <- list()
  for (g0 in seq_len(G)) {
    for (t0 in t_eval) {
      m <- median_inf_g(vthetag[[g0]]$V.rob,
                        thetag[[g0]]$lambda, thetag[[g0]]$beta,
                        thetag[[g0]]$alpha, t0, xmean = xmean)

      lower <- m$median - zcrit * m$SE
      upper <- m$median + zcrit * m$SE

      median_rows[[length(median_rows) + 1]] <- data.frame(
        g = glst[g0], t = t0,
        median = m$median, SE = m$SE,
        lower.cl = lower, upper.cl = upper,
        stringsAsFactors = FALSE
      )
    }
  }
  median_df <- do.call(rbind, median_rows)

  # quick lookup: (g,t) -> (median,SE)
  key_gt <- paste(median_df$g, median_df$t, sep = "||")
  get_median_row <- function(g, t) {
    ii <- match(paste(g, t, sep = "||"), key_gt)
    if (is.na(ii)) stop("Internal error: median lookup failed.", call. = FALSE)
    list(median = median_df$median[ii], SE = median_df$SE[ii])
  }

  # --- 2) median difference + prob for each pair x time ---
  meddiff_rows <- list()
  prob_rows <- list()

  for (r in seq_len(nrow(pair_df))) {
    g1 <- pair_df$g1[r]
    g2 <- pair_df$g2[r]
    i1 <- g_index[[g1]]
    i2 <- g_index[[g2]]

    for (t0 in t_eval) {
      m1 <- get_median_row(g1, t0)
      m2 <- get_median_row(g2, t0)

      md <- median_diff_inf(m1, m2, level = level)

      meddiff_rows[[length(meddiff_rows) + 1]] <- data.frame(
        g1 = g1, g2 = g2, t = t0,
        meddiff = md$meddiff, SE = md$SE,
        lower.cl = md$lower.cl, upper.cl = md$upper.cl,
        z.value = md$z.value, p.value = md$p.value,
        stringsAsFactors = FALSE
      )

      prb <- prob_eff_inf(thetag[[i1]], thetag[[i2]],
                          vthetag[[i1]], vthetag[[i2]],
                          t = t0, xmean = xmean, smaller = smaller,
                          robust = robust)

      prob_rows[[length(prob_rows) + 1]] <- cbind(
        data.frame(g1 = g1, g2 = g2, t = t0, stringsAsFactors = FALSE),
        prb
      )
    }
  }

  list(
    call = match.call(),
    datalist = datalist,
    glst = glst,
    tt = tt,
    t_eval = t_eval,
    xmean = xmean,
    pairs = pair_df,
    fit = thetag,
    vcov = vthetag,
    median = median_df,
    meddiff = do.call(rbind, meddiff_rows),
    prob = do.call(rbind, prob_rows)
  )
}
