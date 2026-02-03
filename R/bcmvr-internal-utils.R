#' @importFrom rlang .data
# Reshape long format data to data matrix
data_reshape <- function(data, outcome, id, group, time, cov = NULL) {

  id_col    <- rlang::as_name(rlang::ensym(id))
  group_col <- rlang::as_name(rlang::ensym(group))
  time_col  <- rlang::as_name(rlang::ensym(time))

  # cov can be NULL / character(0)
  if (is.null(cov)) cov <- character(0)
  if (!is.character(cov)) {
    stop("cov must be NULL or a character vector of column names.",
         call. = FALSE)
  }

  # group list
  glst <- data |>
    dplyr::pull({{ group }}) |>
    as.character() |>
    unique() |>
    sort()

  G <- length(glst)

  yg  <- vector("list", G)
  xg  <- vector("list", G)
  ids <- vector("list", G)
  names(yg) <- names(xg) <- names(ids) <- glst

  for (gg in seq_len(G)) {

    dat_g <- data |> dplyr::filter({{ group }} == glst[gg])

    # --- 0) duplicate check for (id, time) ---
    dup <- dat_g |>
      dplyr::count({{ id }}, {{ time }}) |>
      dplyr::filter(.data$n > 1)
    if (nrow(dup) > 0) {
      stop(
        paste0("Duplicate rows detected for some (id, time) pairs in group=",
               glst[gg], "."),
        call. = FALSE
      )
    }

    # --- 1) time order (to fix column order of wide y) ---
    time_vals <- dat_g |> dplyr::pull({{ time }})
    time_order <- if (is.factor(time_vals)) {
      levels(time_vals)
    } else if (is.numeric(time_vals)) {
      sort(unique(time_vals))
    } else {
      sort(unique(as.character(time_vals)))
    }
    time_order_chr <- as.character(time_order)

    # --- 2) check cov time-invariant within id ---
    if (length(cov) > 0) {
      miss_cov <- setdiff(cov, names(dat_g))
      if (length(miss_cov) > 0) {
        stop(
          paste0("cov contains unknown column(s): ",
                 paste(miss_cov, collapse = ", ")),
          call. = FALSE
        )
      }

      chk <- dat_g |>
        dplyr::group_by({{ id }}) |>
        dplyr::summarise(
          dplyr::across(dplyr::all_of(cov), ~ dplyr::n_distinct(.x,
                                                                na.rm = TRUE)),
          .groups = "drop"
        )

      bad <- chk[rowSums(chk[, cov, drop = FALSE] > 1, na.rm = TRUE) > 0, ,
                 drop = FALSE]
      if (nrow(bad) > 0) {
        bad_ids <- utils::head(bad[[id_col]], 10)
        stop(
          paste0(
            "cov is not time-invariant within some id(s) in group=", glst[gg],
            ".\n",
            "First few problematic ids: ", paste(bad_ids, collapse = ", ")
          ),
          call. = FALSE
        )
      }
    }

    # ---- 3) wide y ----
    yg_df <- dat_g |>
      tidyr::pivot_wider(
        id_cols     = {{ id }},
        names_from  = {{ time }},
        values_from = {{ outcome }}
      )

    # enforce column order of time points
    y_cols <- intersect(time_order_chr, names(yg_df))
    yg_df  <- yg_df[, c(id_col, y_cols), drop = FALSE]

    # drop all-missing outcome rows
    y_only <- yg_df[, setdiff(names(yg_df), id_col), drop = FALSE]
    keep <- rowSums(!is.na(as.matrix(y_only))) > 0
    yg_df <- yg_df[keep, , drop = FALSE]

    keep_id <- yg_df[[id_col]]

    # ---- 4) baseline covariates matrix (or NULL) ----
    if (length(cov) == 0) {
      xg_df <- yg_df[, id_col, drop = FALSE]
      xg_df$..dummy.. <- 0  # placeholder (dropped below)
      xg_df <- xg_df[, id_col, drop = FALSE]
      xg[[gg]] <- NULL
    } else {
      xg_df <- dat_g |>
        dplyr::filter(.data[[id_col]] %in% keep_id) |>
        dplyr::select({{ id }}, dplyr::all_of(cov)) |>
        dplyr::group_by({{ id }}) |>
        dplyr::summarise(
          dplyr::across(dplyr::all_of(cov), ~ dplyr::first(.x)),
          .groups = "drop"
        )

      # align order with yg_df
      xg_df <- xg_df[match(keep_id, xg_df[[id_col]]), , drop = FALSE]

      stopifnot(identical(yg_df[[id_col]], xg_df[[id_col]]))

      xg[[gg]] <- as.matrix(xg_df[, setdiff(names(xg_df), id_col),
                                  drop = FALSE])
    }

    ids[[gg]] <- keep_id
    yg[[gg]]  <- as.matrix(yg_df[, setdiff(names(yg_df), id_col), drop = FALSE])
  }

  list(yg = yg, xg = xg, ids = ids, glst = glst)
}
