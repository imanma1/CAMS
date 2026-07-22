
# ================================================================
# Adaptive local-mapping CAMS using IPCW
# ================================================================
#
# This file ADDS a local calibration method to the existing CAMS
# implementation. Source the original CAMS/model file first so that
# sc_prob() is available. This file defines its own lower-bound helper
# and therefore does not overwrite cams() or est_alpha_ipcw().
#
# Main entry point:
#
#   cams_local_ipcw(...)
#
# Validity groups R are fixed as:
#
#   R0 = {X1 = 0},  R1 = {X1 = 1}.
#
# The learned B cells are decision cells only. The optimization imposes
# ONE IPCW risk constraint per R group, not one constraint per B cell.
#
# Partition learning:
#   1. Stratify-split the supplied calibration sample into audit and
#      final-calibration samples.
#   2. Start with one B cell inside each R group.
#   3. Propose shallow axis-aligned binary splits.
#   4. Retain a split only when cross-validated audit utility improves
#      under the same R-level IPCW risk constraint.
#   5. Freeze the partition before final calibration.
#
# No AIPCW estimator is used in this file.
# ================================================================


# ------------------------------------------------
# Small utilities
# ------------------------------------------------

local_empty_df <- function(columns) {
  out <- as.data.frame(
    setNames(
      replicate(length(columns), logical(0), simplify = FALSE),
      columns
    ),
    stringsAsFactors = FALSE
  )
  out
}


local_bind_rows <- function(x) {
  x <- Filter(
    function(z) {
      is.data.frame(z) && nrow(z) > 0L
    },
    x
  )

  if (length(x) == 0L) {
    return(data.frame())
  }

  all_names <- unique(unlist(lapply(x, names)))

  x <- lapply(
    x,
    function(z) {
      missing_names <- setdiff(all_names, names(z))
      for (nm in missing_names) {
        z[[nm]] <- NA
      }
      z[, all_names, drop = FALSE]
    }
  )

  do.call(rbind, x)
}


local_logit <- function(p, floor = 1e-6) {
  p <- pmin(
    pmax(as.numeric(p), floor),
    1 - floor
  )
  log(p / (1 - p))
}


local_event_weight_ess <- function(event, weight) {
  w <- as.numeric(event) * as.numeric(weight)
  w <- w[is.finite(w) & w > 0]

  if (length(w) == 0L) {
    return(0)
  }

  denominator <- sum(w^2)

  if (!is.finite(denominator) || denominator <= 0) {
    return(0)
  }

  (sum(w)^2) / denominator
}


local_predict_bound_matrix <- function(mdl, newdata, v_list) {
  newdata <- as.data.frame(newdata)
  v_list <- as.numeric(v_list)

  n <- nrow(newdata)
  m <- length(v_list)

  if (n == 0L) {
    return(
      matrix(
        numeric(0),
        nrow = 0L,
        ncol = m
      )
    )
  }

  if (m == 0L) {
    return(
      matrix(
        numeric(0),
        nrow = n,
        ncol = 0L
      )
    )
  }

  predicted <- predict(
    mdl,
    newdata = newdata,
    type = "quantile",
    p = v_list
  )

  if (m == 1L) {
    return(
      matrix(
        as.numeric(predicted),
        nrow = n,
        ncol = 1L
      )
    )
  }

  predicted <- as.matrix(predicted)

  if (
    nrow(predicted) == n &&
    ncol(predicted) == m
  ) {
    return(predicted)
  }

  if (
    nrow(predicted) == m &&
    ncol(predicted) == n
  ) {
    return(t(predicted))
  }

  if (length(predicted) == n * m) {
    return(
      matrix(
        as.numeric(predicted),
        nrow = n,
        ncol = m
      )
    )
  }

  stop(
    sprintf(
      paste0(
        "Unexpected quantile-prediction dimensions: ",
        "n = %d, number of v values = %d, returned length = %d."
      ),
      n,
      m,
      length(predicted)
    )
  )
}


local_bound <- function(mdl, newdata, v) {
  newdata <- as.data.frame(newdata)

  if (nrow(newdata) == 0L) {
    return(numeric(0))
  }

  if (
    length(v) == 0L ||
    !is.finite(v) ||
    v <= 0
  ) {
    return(rep(0, nrow(newdata)))
  }

  as.numeric(
    predict(
      mdl,
      newdata = newdata,
      type = "quantile",
      p = v
    )
  )
}


# ------------------------------------------------
# Stratified calibration split
# ------------------------------------------------

local_stratified_calibration_split <- function(
    data_calib,
    audit_fraction = 0.40,
    split_seed = 1L,
    group_name = "X1") {

  if (!(group_name %in% names(data_calib))) {
    stop(sprintf("The calibration data do not contain %s.", group_name))
  }

  if (
    !is.finite(audit_fraction) ||
    audit_fraction <= 0 ||
    audit_fraction >= 1
  ) {
    stop("audit_fraction must be strictly between 0 and 1.")
  }

  set.seed(split_seed)

  audit_idx <- integer(0)

  group_values <- sort(unique(data_calib[[group_name]]))

  for (group_value in group_values) {
    idx <- which(data_calib[[group_name]] == group_value)
    n_group <- length(idx)

    if (n_group <= 1L) {
      next
    }

    n_audit <- floor(audit_fraction * n_group)
    n_audit <- max(1L, min(n_group - 1L, n_audit))

    audit_idx <- c(
      audit_idx,
      sample(idx, size = n_audit, replace = FALSE)
    )
  }

  audit_idx <- sort(unique(audit_idx))
  final_idx <- setdiff(seq_len(nrow(data_calib)), audit_idx)

  if (
    length(audit_idx) == 0L ||
    length(final_idx) == 0L
  ) {
    stop(
      paste0(
        "The stratified split produced an empty audit or ",
        "final-calibration sample."
      )
    )
  }

  list(
    audit = data_calib[audit_idx, , drop = FALSE],
    final_calibration = data_calib[final_idx, , drop = FALSE],
    audit_original_row = audit_idx,
    final_original_row = final_idx
  )
}


# ------------------------------------------------
# Feature construction for the shallow tree
# ------------------------------------------------

local_build_split_features <- function(
    mdl,
    mdl0,
    data,
    xnames,
    alpha,
    split_features = c(
      "event_lp",
      "logit_G_qalpha",
      "X2",
      "X3",
      "X4"
    ),
    probability_floor = 1e-6) {

  data <- as.data.frame(data)

  missing_x <- setdiff(xnames, names(data))

  if (length(missing_x) > 0L) {
    stop(
      sprintf(
        "Missing covariates for local mapping: %s",
        paste(missing_x, collapse = ", ")
      )
    )
  }

  x <- data[, xnames, drop = FALSE]

  output <- data.frame(
    row.names = seq_len(nrow(data))
  )

  for (feature_name in split_features) {

    if (identical(feature_name, "event_lp")) {

      output[[feature_name]] <- as.numeric(
        predict(
          mdl,
          newdata = x,
          type = "lp"
        )
      )

    } else if (identical(feature_name, "logit_G_qalpha")) {

      q_alpha <- local_bound(
        mdl = mdl,
        newdata = x,
        v = alpha
      )

      G_qalpha <- sc_prob(
        mdl0 = mdl0,
        data = data,
        xnames = xnames,
        t = q_alpha
      )

      output[[feature_name]] <- local_logit(
        G_qalpha,
        floor = probability_floor
      )

    } else {

      if (!(feature_name %in% names(data))) {
        stop(
          sprintf(
            "Requested split feature %s is not available.",
            feature_name
          )
        )
      }

      output[[feature_name]] <- as.numeric(
        data[[feature_name]]
      )
    }
  }

  for (feature_name in names(output)) {
    bad <- !is.finite(output[[feature_name]])

    if (any(bad)) {
      finite_values <- output[[feature_name]][!bad]

      replacement <- if (length(finite_values) == 0L) {
        0
      } else {
        stats::median(finite_values)
      }

      output[[feature_name]][bad] <- replacement
    }
  }

  output
}


# ------------------------------------------------
# IPCW cache for audit or final calibration
# ------------------------------------------------

local_make_ipcw_cache <- function(
    mdl,
    mdl0,
    data,
    xnames,
    v_list,
    probability_floor) {

  required_names <- c(
    xnames,
    "censored_T",
    "event",
    "X1"
  )

  missing_names <- setdiff(required_names, names(data))

  if (length(missing_names) > 0L) {
    stop(
      sprintf(
        "IPCW cache is missing columns: %s",
        paste(missing_names, collapse = ", ")
      )
    )
  }

  x <- data[, xnames, drop = FALSE]

  Y <- as.numeric(data$censored_T)
  event <- as.numeric(data$event)

  G_y_raw <- as.numeric(
    sc_prob(
      mdl0 = mdl0,
      data = data,
      xnames = xnames,
      t = Y
    )
  )

  G_y <- pmax(
    G_y_raw,
    probability_floor
  )

  G_y[!is.finite(G_y)] <- probability_floor

  list(
    data = data,
    x = x,
    Y = Y,
    event = event,
    G_y_raw = G_y_raw,
    G_y = G_y,
    weight = 1 / G_y,
    probability_floor = probability_floor,
    bounds = local_predict_bound_matrix(
      mdl = mdl,
      newdata = x,
      v_list = v_list
    ),
    v_list = as.numeric(v_list)
  )
}


# ------------------------------------------------
# Cell-support diagnostics and cell limits
# ------------------------------------------------

local_support_summary <- function(cache, row_idx) {
  row_idx <- as.integer(row_idx)
  row_idx <- row_idx[
    row_idx >= 1L &
      row_idx <= length(cache$Y)
  ]

  n <- length(row_idx)

  if (n == 0L) {
    return(
      list(
        n = 0L,
        events = 0L,
        censoring_rate = NA_real_,
        event_weight_ess = 0,
        fraction_G_truncated = NA_real_
      )
    )
  }

  list(
    n = n,
    events = sum(cache$event[row_idx] == 1),
    censoring_rate = mean(cache$event[row_idx] == 0),
    event_weight_ess = local_event_weight_ess(
      event = cache$event[row_idx],
      weight = cache$weight[row_idx]
    ),
    fraction_G_truncated = mean(
      !is.finite(cache$G_y_raw[row_idx]) |
        cache$G_y_raw[row_idx] < cache$probability_floor
    )
  )
}


local_compute_cell_limit <- function(
    cache,
    group_rows,
    min_absolute_n = 40L,
    min_group_fraction = 0.15,
    min_events = 30L,
    min_event_ess = 30,
    global_cap = 3L) {

  group_support <- local_support_summary(
    cache = cache,
    row_idx = group_rows
  )

  required_n <- max(
    as.integer(min_absolute_n),
    ceiling(
      min_group_fraction *
        group_support$n
    )
  )

  max_by_n <- if (required_n <= 0L) {
    1L
  } else {
    floor(group_support$n / required_n)
  }

  max_by_events <- if (min_events <= 0L) {
    global_cap
  } else {
    floor(group_support$events / min_events)
  }

  max_by_ess <- if (min_event_ess <= 0) {
    global_cap
  } else {
    floor(group_support$event_weight_ess / min_event_ess)
  }

  K_max <- max(
    1L,
    min(
      as.integer(global_cap),
      as.integer(max_by_n),
      as.integer(max_by_events),
      as.integer(max_by_ess)
    )
  )

  list(
    n_group = group_support$n,
    observed_events = group_support$events,
    event_weight_ess = group_support$event_weight_ess,
    censoring_rate = group_support$censoring_rate,
    fraction_G_truncated = group_support$fraction_G_truncated,
    required_min_raw_n = required_n,
    min_group_fraction = min_group_fraction,
    min_events = min_events,
    min_event_ess = min_event_ess,
    max_by_raw_n = max_by_n,
    max_by_events = max_by_events,
    max_by_event_ess = max_by_ess,
    global_cap = global_cap,
    K_max = K_max
  )
}


local_child_is_supported <- function(
    support,
    required_n,
    min_events,
    min_event_ess) {

  isTRUE(
    support$n >= required_n &&
      support$events >= min_events &&
      support$event_weight_ess >= min_event_ess
  )
}


# ------------------------------------------------
# Risk and utility curves for a fixed partition
# ------------------------------------------------

local_compute_cell_curves <- function(
    cache,
    row_idx,
    leaf_assignment,
    expected_leaves = NULL) {

  row_idx <- as.integer(row_idx)
  leaf_assignment <- as.integer(leaf_assignment)

  if (length(row_idx) != length(leaf_assignment)) {
    stop("row_idx and leaf_assignment must have the same length.")
  }

  if (is.null(expected_leaves)) {
    expected_leaves <- sort(unique(leaf_assignment))
  }

  expected_leaves <- as.integer(expected_leaves)
  n_group <- length(row_idx)

  output <- vector(
    "list",
    length(expected_leaves)
  )

  names(output) <- as.character(expected_leaves)

  for (j in seq_along(expected_leaves)) {
    leaf_id <- expected_leaves[j]
    local_pos <- which(leaf_assignment == leaf_id)
    idx <- row_idx[local_pos]
    n_cell <- length(idx)

    if (n_cell == 0L) {
      output[[j]] <- list(
        leaf_id = leaf_id,
        n = 0L,
        p = 0,
        risk = rep(0, length(cache$v_list)),
        utility = rep(0, length(cache$v_list))
      )
      next
    }

    bound_matrix <- cache$bounds[
      idx,
      ,
      drop = FALSE
    ]

    indicator_matrix <- sweep(
      bound_matrix,
      MARGIN = 1L,
      STATS = cache$Y[idx],
      FUN = ">"
    )

    weighted_event <- (
      cache$event[idx] *
        cache$weight[idx]
    )

    weighted_loss <- sweep(
      indicator_matrix,
      MARGIN = 1L,
      STATS = weighted_event,
      FUN = "*"
    )

    risk <- colSums(weighted_loss) / n_cell
    risk[!is.finite(risk)] <- 1
    risk <- cummax(risk)

    utility <- colMeans(bound_matrix)
    utility[!is.finite(utility)] <- 0
    utility <- cummax(utility)

    output[[j]] <- list(
      leaf_id = leaf_id,
      n = n_cell,
      p = if (n_group == 0L) 0 else n_cell / n_group,
      risk = risk,
      utility = utility
    )
  }

  output
}


# ------------------------------------------------
# Exact small-K optimization under one R constraint
# ------------------------------------------------

local_optimize_r_map <- function(
    cell_curves,
    v_list,
    risk_budget) {

  K <- length(cell_curves)

  if (K == 0L) {
    return(
      list(
        selected_option = integer(0),
        selected_v = numeric(0),
        estimated_risk = 0,
        estimated_utility = 0,
        feasible = TRUE
      )
    )
  }

  if (K > 3L) {
    stop(
      paste0(
        "local_optimize_r_map() supports at most three B cells. ",
        "Keep global_cap <= 3."
      )
    )
  }

  risk_budget <- max(
    0,
    as.numeric(risk_budget)
  )

  v_augmented <- c(0, as.numeric(v_list))

  contribution_risk <- vector("list", K)
  contribution_utility <- vector("list", K)

  for (k in seq_len(K)) {
    cell <- cell_curves[[k]]

    risk <- c(0, cummax(as.numeric(cell$risk)))
    utility <- c(0, cummax(as.numeric(cell$utility)))

    contribution_risk[[k]] <- cell$p * risk
    contribution_utility[[k]] <- cell$p * utility
  }

  best_option <- rep(1L, K)
  best_utility <- 0
  best_risk <- 0

  if (K == 1L) {

    feasible <- which(
      contribution_risk[[1]] <= risk_budget + 1e-12
    )

    if (length(feasible) > 0L) {
      utilities <- contribution_utility[[1]][feasible]
      j <- feasible[which.max(utilities)]

      best_option[1] <- j
      best_utility <- contribution_utility[[1]][j]
      best_risk <- contribution_risk[[1]][j]
    }

  } else if (K == 2L) {

    risk_1 <- contribution_risk[[1]]
    risk_2 <- contribution_risk[[2]]
    utility_1 <- contribution_utility[[1]]
    utility_2 <- contribution_utility[[2]]

    for (j1 in seq_along(risk_1)) {
      remaining <- risk_budget - risk_1[j1]

      if (remaining < -1e-12) {
        next
      }

      j2 <- findInterval(
        remaining + 1e-12,
        risk_2
      )

      if (j2 < 1L) {
        next
      }

      candidate_utility <- (
        utility_1[j1] +
          utility_2[j2]
      )

      candidate_risk <- (
        risk_1[j1] +
          risk_2[j2]
      )

      if (
        candidate_utility > best_utility + 1e-12 ||
        (
          abs(candidate_utility - best_utility) <= 1e-12 &&
          candidate_risk < best_risk
        )
      ) {
        best_option <- c(j1, j2)
        best_utility <- candidate_utility
        best_risk <- candidate_risk
      }
    }

  } else if (K == 3L) {

    risk_1 <- contribution_risk[[1]]
    risk_2 <- contribution_risk[[2]]
    risk_3 <- contribution_risk[[3]]

    utility_1 <- contribution_utility[[1]]
    utility_2 <- contribution_utility[[2]]
    utility_3 <- contribution_utility[[3]]

    for (j1 in seq_along(risk_1)) {

      remaining_after_1 <- (
        risk_budget -
          risk_1[j1] -
          risk_2
      )

      valid_j2 <- which(
        remaining_after_1 >= -1e-12
      )

      if (length(valid_j2) == 0L) {
        next
      }

      j3 <- findInterval(
        remaining_after_1[valid_j2] + 1e-12,
        risk_3
      )

      valid_pair <- j3 >= 1L

      if (!any(valid_pair)) {
        next
      }

      valid_j2 <- valid_j2[valid_pair]
      j3 <- j3[valid_pair]

      candidate_utility <- (
        utility_1[j1] +
          utility_2[valid_j2] +
          utility_3[j3]
      )

      best_local_position <- which.max(
        candidate_utility
      )

      j2_best <- valid_j2[best_local_position]
      j3_best <- j3[best_local_position]

      candidate_best_utility <- candidate_utility[
        best_local_position
      ]

      candidate_best_risk <- (
        risk_1[j1] +
          risk_2[j2_best] +
          risk_3[j3_best]
      )

      if (
        candidate_best_utility > best_utility + 1e-12 ||
        (
          abs(candidate_best_utility - best_utility) <= 1e-12 &&
          candidate_best_risk < best_risk
        )
      ) {
        best_option <- c(
          j1,
          j2_best,
          j3_best
        )
        best_utility <- candidate_best_utility
        best_risk <- candidate_best_risk
      }
    }
  }

  selected_v <- v_augmented[best_option]
  names(selected_v) <- names(cell_curves)

  names(best_option) <- names(cell_curves)

  list(
    selected_option = best_option,
    selected_v = selected_v,
    estimated_risk = best_risk,
    estimated_utility = best_utility,
    feasible = best_risk <= risk_budget + 1e-10
  )
}


local_evaluate_map <- function(
    cache,
    row_idx,
    leaf_assignment,
    selected_option) {

  row_idx <- as.integer(row_idx)
  leaf_assignment <- as.integer(leaf_assignment)

  if (length(row_idx) == 0L) {
    return(
      list(
        risk = NA_real_,
        utility = NA_real_,
        bounds = numeric(0)
      )
    )
  }

  option_by_row <- selected_option[
    as.character(leaf_assignment)
  ]

  if (any(is.na(option_by_row))) {
    stop("A validation row was assigned to an unknown B leaf.")
  }

  lower_bound <- numeric(length(row_idx))

  positive_option <- option_by_row > 1L

  if (any(positive_option)) {
    matrix_column <- option_by_row[positive_option] - 1L

    lower_bound[positive_option] <- cache$bounds[
      cbind(
        row_idx[positive_option],
        matrix_column
      )
    ]
  }

  weighted_loss <- (
    cache$event[row_idx] *
      cache$weight[row_idx] *
      as.numeric(cache$Y[row_idx] < lower_bound)
  )

  list(
    risk = mean(weighted_loss),
    utility = mean(lower_bound),
    bounds = lower_bound
  )
}


# ------------------------------------------------
# Cross-validation folds
# ------------------------------------------------

local_make_folds <- function(n, K, seed) {
  if (n <= 0L) {
    return(integer(0))
  }

  K <- min(
    as.integer(K),
    n
  )

  K <- max(2L, K)

  set.seed(seed)

  sample(
    rep(
      seq_len(K),
      length.out = n
    )
  )
}


# ------------------------------------------------
# Candidate-split evaluation
# ------------------------------------------------

local_evaluate_candidate_split_cv <- function(
    cache,
    group_rows,
    current_assignment,
    proposed_assignment,
    current_leaves,
    proposed_leaves,
    alpha,
    folds,
    required_n,
    min_events,
    min_event_ess,
    audit_risk_tolerance = 0.03) {

  fold_values <- sort(unique(folds))
  fold_log <- vector("list", length(fold_values))

  minimum_train_n <- max(
    10L,
    floor(0.50 * required_n)
  )

  minimum_train_events <- max(
    5L,
    floor(0.50 * min_events)
  )

  minimum_train_ess <- max(
    5,
    0.50 * min_event_ess
  )

  for (j in seq_along(fold_values)) {
    fold_value <- fold_values[j]

    train_pos <- which(folds != fold_value)
    valid_pos <- which(folds == fold_value)

    train_rows <- group_rows[train_pos]
    valid_rows <- group_rows[valid_pos]

    proposed_train_assignment <- proposed_assignment[train_pos]
    proposed_valid_assignment <- proposed_assignment[valid_pos]

    train_supported <- TRUE

    for (leaf_id in proposed_leaves) {
      leaf_train_rows <- train_rows[
        proposed_train_assignment == leaf_id
      ]

      support <- local_support_summary(
        cache = cache,
        row_idx = leaf_train_rows
      )

      if (
        support$n < minimum_train_n ||
        support$events < minimum_train_events ||
        support$event_weight_ess < minimum_train_ess
      ) {
        train_supported <- FALSE
        break
      }
    }

    validation_supported <- all(
      vapply(
        proposed_leaves,
        function(leaf_id) {
          sum(proposed_valid_assignment == leaf_id) >= 2L
        },
        logical(1)
      )
    )

    if (
      !train_supported ||
      !validation_supported
    ) {
      fold_log[[j]] <- data.frame(
        fold = fold_value,
        valid_fold = FALSE,
        current_risk = NA_real_,
        proposed_risk = NA_real_,
        current_utility = NA_real_,
        proposed_utility = NA_real_,
        absolute_gain = NA_real_,
        relative_gain = NA_real_,
        proposed_risk_ok = FALSE,
        stringsAsFactors = FALSE
      )
      next
    }

    current_curves <- local_compute_cell_curves(
      cache = cache,
      row_idx = train_rows,
      leaf_assignment = current_assignment[train_pos],
      expected_leaves = current_leaves
    )

    proposed_curves <- local_compute_cell_curves(
      cache = cache,
      row_idx = train_rows,
      leaf_assignment = proposed_train_assignment,
      expected_leaves = proposed_leaves
    )

    current_map <- local_optimize_r_map(
      cell_curves = current_curves,
      v_list = cache$v_list,
      risk_budget = alpha
    )

    proposed_map <- local_optimize_r_map(
      cell_curves = proposed_curves,
      v_list = cache$v_list,
      risk_budget = alpha
    )

    current_eval <- local_evaluate_map(
      cache = cache,
      row_idx = valid_rows,
      leaf_assignment = current_assignment[valid_pos],
      selected_option = current_map$selected_option
    )

    proposed_eval <- local_evaluate_map(
      cache = cache,
      row_idx = valid_rows,
      leaf_assignment = proposed_valid_assignment,
      selected_option = proposed_map$selected_option
    )

    absolute_gain <- (
      proposed_eval$utility -
        current_eval$utility
    )

    relative_gain <- absolute_gain / max(
      abs(current_eval$utility),
      1e-8
    )

    proposed_risk_ok <- isTRUE(
      proposed_eval$risk <=
        alpha + audit_risk_tolerance
    )

    fold_log[[j]] <- data.frame(
      fold = fold_value,
      valid_fold = TRUE,
      current_risk = current_eval$risk,
      proposed_risk = proposed_eval$risk,
      current_utility = current_eval$utility,
      proposed_utility = proposed_eval$utility,
      absolute_gain = absolute_gain,
      relative_gain = relative_gain,
      proposed_risk_ok = proposed_risk_ok,
      stringsAsFactors = FALSE
    )
  }

  fold_log <- local_bind_rows(fold_log)
  valid_log <- fold_log[
    fold_log$valid_fold,
    ,
    drop = FALSE
  ]

  if (nrow(valid_log) == 0L) {
    return(
      list(
        summary = data.frame(
          valid_folds = 0L,
          positive_gain_folds = 0L,
          risk_ok_folds = 0L,
          mean_absolute_gain = NA_real_,
          mean_relative_gain = NA_real_,
          mean_current_risk = NA_real_,
          mean_proposed_risk = NA_real_,
          mean_current_utility = NA_real_,
          mean_proposed_utility = NA_real_,
          stringsAsFactors = FALSE
        ),
        folds = fold_log
      )
    )
  }

  list(
    summary = data.frame(
      valid_folds = nrow(valid_log),
      positive_gain_folds = sum(
        valid_log$absolute_gain > 0
      ),
      risk_ok_folds = sum(
        valid_log$proposed_risk_ok
      ),
      mean_absolute_gain = mean(
        valid_log$absolute_gain
      ),
      mean_relative_gain = mean(
        valid_log$relative_gain
      ),
      mean_current_risk = mean(
        valid_log$current_risk
      ),
      mean_proposed_risk = mean(
        valid_log$proposed_risk
      ),
      mean_current_utility = mean(
        valid_log$current_utility
      ),
      mean_proposed_utility = mean(
        valid_log$proposed_utility
      ),
      stringsAsFactors = FALSE
    ),
    folds = fold_log
  )
}


# ------------------------------------------------
# Learn one shallow adaptive partition inside one R
# ------------------------------------------------

local_learn_group_partition <- function(
    group_value,
    audit_cache,
    audit_features,
    alpha,
    cell_limit,
    split_quantiles = seq(0.20, 0.80, by = 0.10),
    cv_folds = 4L,
    cv_seed = 1L,
    min_relative_gain = 0.01,
    min_positive_fold_fraction = 0.75,
    min_risk_ok_fold_fraction = 0.75,
    audit_risk_tolerance = 0.03) {

  group_rows <- which(
    audit_cache$data$X1 == group_value
  )

  n_group <- length(group_rows)

  current_assignment <- rep(
    1L,
    n_group
  )

  current_leaves <- 1L
  next_leaf_id <- 2L

  rules <- data.frame(
    step = integer(0),
    parent_leaf = integer(0),
    left_leaf = integer(0),
    right_leaf = integer(0),
    feature = character(0),
    threshold = numeric(0),
    stringsAsFactors = FALSE
  )

  candidate_logs <- list()
  accepted_fold_logs <- list()

  folds <- local_make_folds(
    n = n_group,
    K = cv_folds,
    seed = cv_seed + as.integer(group_value)
  )

  candidate_counter <- 0L
  step <- 0L

  while (
    length(current_leaves) <
      cell_limit$K_max
  ) {

    step <- step + 1L

    step_candidates <- list()
    step_candidate_assignments <- list()
    step_candidate_rules <- list()
    step_fold_logs <- list()

    for (parent_leaf in current_leaves) {

      parent_pos <- which(
        current_assignment == parent_leaf
      )

      if (length(parent_pos) == 0L) {
        next
      }

      parent_rows <- group_rows[parent_pos]

      for (feature_name in names(audit_features)) {

        feature_values <- audit_features[
          parent_rows,
          feature_name
        ]

        thresholds <- unique(
          as.numeric(
            stats::quantile(
              feature_values,
              probs = split_quantiles,
              na.rm = TRUE,
              names = FALSE,
              type = 8
            )
          )
        )

        thresholds <- thresholds[
          is.finite(thresholds)
        ]

        for (threshold in thresholds) {

          candidate_counter <- candidate_counter + 1L
          candidate_id <- candidate_counter

          left_pos <- parent_pos[
            audit_features[
              parent_rows,
              feature_name
            ] <= threshold
          ]

          right_pos <- setdiff(
            parent_pos,
            left_pos
          )

          left_rows <- group_rows[left_pos]
          right_rows <- group_rows[right_pos]

          left_support <- local_support_summary(
            cache = audit_cache,
            row_idx = left_rows
          )

          right_support <- local_support_summary(
            cache = audit_cache,
            row_idx = right_rows
          )

          support_pass <- (
            local_child_is_supported(
              support = left_support,
              required_n = cell_limit$required_min_raw_n,
              min_events = cell_limit$min_events,
              min_event_ess = cell_limit$min_event_ess
            ) &&
              local_child_is_supported(
                support = right_support,
                required_n = cell_limit$required_min_raw_n,
                min_events = cell_limit$min_events,
                min_event_ess = cell_limit$min_event_ess
              )
          )

          candidate_log <- data.frame(
            candidate_id = candidate_id,
            step = step,
            group = group_value,
            parent_leaf = parent_leaf,
            feature = feature_name,
            threshold = threshold,
            left_n = left_support$n,
            right_n = right_support$n,
            left_events = left_support$events,
            right_events = right_support$events,
            left_event_weight_ess = left_support$event_weight_ess,
            right_event_weight_ess = right_support$event_weight_ess,
            support_pass = support_pass,
            valid_folds = NA_integer_,
            positive_gain_folds = NA_integer_,
            risk_ok_folds = NA_integer_,
            mean_absolute_gain = NA_real_,
            mean_relative_gain = NA_real_,
            mean_current_risk = NA_real_,
            mean_proposed_risk = NA_real_,
            mean_current_utility = NA_real_,
            mean_proposed_utility = NA_real_,
            passes_cv_rule = FALSE,
            accepted = FALSE,
            stringsAsFactors = FALSE
          )

          if (!support_pass) {
            step_candidates[[length(step_candidates) + 1L]] <- candidate_log
            next
          }

          proposed_assignment <- current_assignment

          proposed_assignment[
            left_pos
          ] <- next_leaf_id

          proposed_assignment[
            right_pos
          ] <- next_leaf_id + 1L

          proposed_leaves <- sort(
            c(
              setdiff(
                current_leaves,
                parent_leaf
              ),
              next_leaf_id,
              next_leaf_id + 1L
            )
          )

          cv_result <- local_evaluate_candidate_split_cv(
            cache = audit_cache,
            group_rows = group_rows,
            current_assignment = current_assignment,
            proposed_assignment = proposed_assignment,
            current_leaves = current_leaves,
            proposed_leaves = proposed_leaves,
            alpha = alpha,
            folds = folds,
            required_n = cell_limit$required_min_raw_n,
            min_events = cell_limit$min_events,
            min_event_ess = cell_limit$min_event_ess,
            audit_risk_tolerance = audit_risk_tolerance
          )

          cv_summary <- cv_result$summary

          minimum_valid_folds <- ceiling(
            0.75 *
              length(unique(folds))
          )

          minimum_positive_folds <- ceiling(
            min_positive_fold_fraction *
              cv_summary$valid_folds
          )

          minimum_risk_ok_folds <- ceiling(
            min_risk_ok_fold_fraction *
              cv_summary$valid_folds
          )

          passes_cv_rule <- isTRUE(
            cv_summary$valid_folds >= minimum_valid_folds &&
              cv_summary$positive_gain_folds >= minimum_positive_folds &&
              cv_summary$risk_ok_folds >= minimum_risk_ok_folds &&
              is.finite(cv_summary$mean_relative_gain) &&
              cv_summary$mean_relative_gain >= min_relative_gain
          )

          candidate_log$valid_folds <- cv_summary$valid_folds
          candidate_log$positive_gain_folds <- cv_summary$positive_gain_folds
          candidate_log$risk_ok_folds <- cv_summary$risk_ok_folds
          candidate_log$mean_absolute_gain <- cv_summary$mean_absolute_gain
          candidate_log$mean_relative_gain <- cv_summary$mean_relative_gain
          candidate_log$mean_current_risk <- cv_summary$mean_current_risk
          candidate_log$mean_proposed_risk <- cv_summary$mean_proposed_risk
          candidate_log$mean_current_utility <- cv_summary$mean_current_utility
          candidate_log$mean_proposed_utility <- cv_summary$mean_proposed_utility
          candidate_log$passes_cv_rule <- passes_cv_rule

          cv_result$folds$candidate_id <- candidate_id
          cv_result$folds$step <- step
          cv_result$folds$group <- group_value
          cv_result$folds$parent_leaf <- parent_leaf
          cv_result$folds$feature <- feature_name
          cv_result$folds$threshold <- threshold

          candidate_index <- length(step_candidates) + 1L

          step_candidates[[candidate_index]] <- candidate_log
          step_candidate_assignments[[as.character(candidate_id)]] <- proposed_assignment
          step_candidate_rules[[as.character(candidate_id)]] <- data.frame(
            step = step,
            parent_leaf = parent_leaf,
            left_leaf = next_leaf_id,
            right_leaf = next_leaf_id + 1L,
            feature = feature_name,
            threshold = threshold,
            stringsAsFactors = FALSE
          )
          step_fold_logs[[as.character(candidate_id)]] <- cv_result$folds
        }
      }
    }

    step_candidate_log <- local_bind_rows(
      step_candidates
    )

    if (nrow(step_candidate_log) == 0L) {
      break
    }

    passing <- step_candidate_log[
      step_candidate_log$passes_cv_rule,
      ,
      drop = FALSE
    ]

    if (nrow(passing) == 0L) {
      candidate_logs[[length(candidate_logs) + 1L]] <- step_candidate_log
      break
    }

    ordering <- order(
      -passing$mean_relative_gain,
      -passing$mean_absolute_gain,
      passing$candidate_id
    )

    accepted_candidate <- passing[
      ordering[1],
      ,
      drop = FALSE
    ]

    accepted_id <- as.character(
      accepted_candidate$candidate_id
    )

    step_candidate_log$accepted[
      step_candidate_log$candidate_id ==
        accepted_candidate$candidate_id
    ] <- TRUE

    candidate_logs[[length(candidate_logs) + 1L]] <- step_candidate_log

    current_assignment <- step_candidate_assignments[[accepted_id]]

    accepted_rule <- step_candidate_rules[[accepted_id]]

    rules <- rbind(
      rules,
      accepted_rule
    )

    accepted_fold_logs[[length(accepted_fold_logs) + 1L]] <-
      step_fold_logs[[accepted_id]]

    current_leaves <- sort(
      unique(current_assignment)
    )

    next_leaf_id <- next_leaf_id + 2L
  }

  candidate_log <- local_bind_rows(
    candidate_logs
  )

  accepted_fold_log <- local_bind_rows(
    accepted_fold_logs
  )

  leaf_label_map <- data.frame(
    group = group_value,
    leaf_id = sort(unique(current_assignment)),
    cell = paste0(
      "B",
      group_value,
      "_",
      seq_along(sort(unique(current_assignment)))
    ),
    stringsAsFactors = FALSE
  )

  list(
    group = group_value,
    rules = rules,
    final_assignment = current_assignment,
    final_leaves = sort(unique(current_assignment)),
    leaf_label_map = leaf_label_map,
    candidate_log = candidate_log,
    accepted_fold_log = accepted_fold_log,
    cell_limit = cell_limit,
    audit_group_rows = group_rows
  )
}


# ------------------------------------------------
# Learn the full partition for R0 and R1
# ------------------------------------------------

local_learn_partition <- function(
    mdl,
    mdl0,
    audit_data,
    xnames,
    alpha,
    audit_v_list = seq(0.01, 0.99, by = 0.01),
    split_features = c(
      "event_lp",
      "logit_G_qalpha",
      "X2",
      "X3",
      "X4"
    ),
    split_quantiles = seq(0.20, 0.80, by = 0.10),
    cv_folds = 4L,
    cv_seed = 1L,
    min_absolute_n = 40L,
    min_group_fraction = 0.15,
    min_events = 30L,
    min_event_ess = 30,
    global_cap = 3L,
    min_relative_gain = 0.01,
    min_positive_fold_fraction = 0.75,
    min_risk_ok_fold_fraction = 0.75,
    audit_risk_tolerance = 0.03,
    audit_probability_floor = NULL) {

  if (global_cap > 3L) {
    stop(
      paste0(
        "This implementation uses an exact optimizer for at most ",
        "three cells. Set global_cap <= 3."
      )
    )
  }

  if (is.null(audit_probability_floor)) {
    audit_probability_floor <- 1 / log(
      max(nrow(audit_data), 3)
    )
  }

  audit_cache <- local_make_ipcw_cache(
    mdl = mdl,
    mdl0 = mdl0,
    data = audit_data,
    xnames = xnames,
    v_list = audit_v_list,
    probability_floor = audit_probability_floor
  )

  audit_features <- local_build_split_features(
    mdl = mdl,
    mdl0 = mdl0,
    data = audit_data,
    xnames = xnames,
    alpha = alpha,
    split_features = split_features
  )

  group_models <- list()
  limit_rows <- list()

  for (group_value in c(0L, 1L)) {

    group_rows <- which(
      audit_data$X1 == group_value
    )

    cell_limit <- local_compute_cell_limit(
      cache = audit_cache,
      group_rows = group_rows,
      min_absolute_n = min_absolute_n,
      min_group_fraction = min_group_fraction,
      min_events = min_events,
      min_event_ess = min_event_ess,
      global_cap = global_cap
    )

    limit_rows[[length(limit_rows) + 1L]] <- data.frame(
      group = group_value,
      n_group = cell_limit$n_group,
      observed_events = cell_limit$observed_events,
      event_weight_ess = cell_limit$event_weight_ess,
      censoring_rate = cell_limit$censoring_rate,
      fraction_G_truncated = cell_limit$fraction_G_truncated,
      required_min_raw_n = cell_limit$required_min_raw_n,
      min_group_fraction = cell_limit$min_group_fraction,
      min_events = cell_limit$min_events,
      min_event_ess = cell_limit$min_event_ess,
      max_by_raw_n = cell_limit$max_by_raw_n,
      max_by_events = cell_limit$max_by_events,
      max_by_event_ess = cell_limit$max_by_event_ess,
      global_cap = cell_limit$global_cap,
      K_max = cell_limit$K_max,
      stringsAsFactors = FALSE
    )

    group_models[[as.character(group_value)]] <- local_learn_group_partition(
      group_value = group_value,
      audit_cache = audit_cache,
      audit_features = audit_features,
      alpha = alpha,
      cell_limit = cell_limit,
      split_quantiles = split_quantiles,
      cv_folds = cv_folds,
      cv_seed = cv_seed,
      min_relative_gain = min_relative_gain,
      min_positive_fold_fraction = min_positive_fold_fraction,
      min_risk_ok_fold_fraction = min_risk_ok_fold_fraction,
      audit_risk_tolerance = audit_risk_tolerance
    )
  }

  structure(
    list(
      group_models = group_models,
      split_features = split_features,
      alpha = alpha,
      audit_v_list = audit_v_list,
      audit_probability_floor = audit_probability_floor,
      cell_limits = local_bind_rows(limit_rows),
      audit_cache = audit_cache,
      audit_features = audit_features
    ),
    class = "local_cams_partition"
  )
}


# ------------------------------------------------
# Apply a frozen partition
# ------------------------------------------------

local_assign_partition <- function(
    partition,
    mdl,
    mdl0,
    data,
    xnames,
    alpha) {

  if (!inherits(partition, "local_cams_partition")) {
    stop("partition must inherit from local_cams_partition.")
  }

  data <- as.data.frame(data)

  features <- local_build_split_features(
    mdl = mdl,
    mdl0 = mdl0,
    data = data,
    xnames = xnames,
    alpha = alpha,
    split_features = partition$split_features
  )

  leaf_id <- rep(
    NA_integer_,
    nrow(data)
  )

  cell <- rep(
    NA_character_,
    nrow(data)
  )

  for (group_value in c(0L, 1L)) {

    idx <- which(
      data$X1 == group_value
    )

    if (length(idx) == 0L) {
      next
    }

    group_model <- partition$group_models[[as.character(group_value)]]

    group_leaf <- rep(
      1L,
      length(idx)
    )

    if (nrow(group_model$rules) > 0L) {

      for (rule_index in seq_len(nrow(group_model$rules))) {

        rule <- group_model$rules[
          rule_index,
          ,
          drop = FALSE
        ]

        active <- which(
          group_leaf == rule$parent_leaf
        )

        if (length(active) == 0L) {
          next
        }

        values <- features[
          idx[active],
          rule$feature
        ]

        go_left <- values <= rule$threshold

        group_leaf[
          active[go_left]
        ] <- rule$left_leaf

        group_leaf[
          active[!go_left]
        ] <- rule$right_leaf
      }
    }

    label_lookup <- setNames(
      group_model$leaf_label_map$cell,
      group_model$leaf_label_map$leaf_id
    )

    leaf_id[idx] <- group_leaf
    cell[idx] <- label_lookup[
      as.character(group_leaf)
    ]
  }

  data.frame(
    group = as.integer(data$X1),
    leaf_id = leaf_id,
    cell = cell,
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------
# Final calibration with the frozen B partition
# ------------------------------------------------

local_calibrate_partition_ipcw <- function(
    mdl,
    mdl0,
    partition,
    final_calibration_data,
    xnames,
    alpha,
    final_v_list = seq(0.001, 0.999, by = 0.001),
    C_0 = 0,
    final_probability_floor = NULL) {

  if (is.null(final_probability_floor)) {
    final_probability_floor <- 1 / log(
      max(nrow(final_calibration_data), 3)
    )
  }

  final_cache <- local_make_ipcw_cache(
    mdl = mdl,
    mdl0 = mdl0,
    data = final_calibration_data,
    xnames = xnames,
    v_list = final_v_list,
    probability_floor = final_probability_floor
  )

  assignment <- local_assign_partition(
    partition = partition,
    mdl = mdl,
    mdl0 = mdl0,
    data = final_calibration_data,
    xnames = xnames,
    alpha = alpha
  )

  selected_maps <- list()
  r_summaries <- list()
  leaf_summaries <- list()
  curve_logs <- list()

  for (group_value in c(0L, 1L)) {

    group_rows <- which(
      final_calibration_data$X1 == group_value
    )

    group_model <- partition$group_models[[as.character(group_value)]]

    expected_leaves <- group_model$final_leaves

    group_assignment <- assignment$leaf_id[
      group_rows
    ]

    if (length(group_rows) == 0L) {

      selected_v <- setNames(
        rep(0, length(expected_leaves)),
        as.character(expected_leaves)
      )

      selected_maps[[as.character(group_value)]] <- list(
        selected_v = selected_v,
        selected_option = setNames(
          rep(1L, length(expected_leaves)),
          as.character(expected_leaves)
        ),
        penalty = NA_real_,
        risk_budget = 0,
        estimated_risk = 0,
        estimated_utility = 0
      )

      next
    }

    cell_curves <- local_compute_cell_curves(
      cache = final_cache,
      row_idx = group_rows,
      leaf_assignment = group_assignment,
      expected_leaves = expected_leaves
    )

    group_weights <- final_cache$weight[
      group_rows
    ]

    n_eff <- (
      length(group_rows)^2
    ) / sum(group_weights^2)

    n_maps <- length(final_v_list)^length(
      expected_leaves
    )

    log_n_maps <- length(expected_leaves) *
      log(length(final_v_list))

    penalty <- C_0 * sqrt(
      (
        log(2) +
          log_n_maps
      ) / max(n_eff, 1)
    )

    risk_budget <- max(
      0,
      alpha - penalty
    )

    local_map <- local_optimize_r_map(
      cell_curves = cell_curves,
      v_list = final_v_list,
      risk_budget = risk_budget
    )

    unsplit_curves <- local_compute_cell_curves(
      cache = final_cache,
      row_idx = group_rows,
      leaf_assignment = rep(1L, length(group_rows)),
      expected_leaves = 1L
    )

    unsplit_map <- local_optimize_r_map(
      cell_curves = unsplit_curves,
      v_list = final_v_list,
      risk_budget = risk_budget
    )

    selected_maps[[as.character(group_value)]] <- list(
      selected_v = local_map$selected_v,
      selected_option = local_map$selected_option,
      penalty = penalty,
      risk_budget = risk_budget,
      estimated_risk = local_map$estimated_risk,
      estimated_utility = local_map$estimated_utility
    )

    r_summaries[[length(r_summaries) + 1L]] <- data.frame(
      group = group_value,
      n_final_calibration = length(group_rows),
      K_selected = length(expected_leaves),
      n_eff = n_eff,
      C_0 = C_0,
      penalty = penalty,
      risk_budget = risk_budget,
      log_number_of_candidate_maps = log_n_maps,
      estimated_R_risk_local = local_map$estimated_risk,
      estimated_R_utility_local = local_map$estimated_utility,
      estimated_R_risk_one_cell = unsplit_map$estimated_risk,
      estimated_R_utility_one_cell = unsplit_map$estimated_utility,
      estimated_absolute_utility_gain = (
        local_map$estimated_utility -
          unsplit_map$estimated_utility
      ),
      estimated_relative_utility_gain = (
        local_map$estimated_utility -
          unsplit_map$estimated_utility
      ) / max(
        abs(unsplit_map$estimated_utility),
        1e-8
      ),
      fraction_G_truncated = mean(
        !is.finite(
          final_cache$G_y_raw[group_rows]
        ) |
          final_cache$G_y_raw[group_rows] <
            final_probability_floor
      ),
      stringsAsFactors = FALSE
    )

    label_lookup <- setNames(
      group_model$leaf_label_map$cell,
      group_model$leaf_label_map$leaf_id
    )

    for (leaf_name in names(cell_curves)) {

      cell_curve <- cell_curves[[leaf_name]]
      selected_option <- local_map$selected_option[
        leaf_name
      ]

      selected_v <- local_map$selected_v[
        leaf_name
      ]

      selected_risk <- if (selected_option <= 1L) {
        0
      } else {
        cell_curve$risk[selected_option - 1L]
      }

      selected_utility <- if (selected_option <= 1L) {
        0
      } else {
        cell_curve$utility[selected_option - 1L]
      }

      leaf_rows <- group_rows[
        group_assignment == as.integer(leaf_name)
      ]

      support <- local_support_summary(
        cache = final_cache,
        row_idx = leaf_rows
      )

      leaf_summaries[[length(leaf_summaries) + 1L]] <- data.frame(
        group = group_value,
        leaf_id = as.integer(leaf_name),
        cell = label_lookup[leaf_name],
        n_final_calibration = support$n,
        final_calibration_fraction_within_R = cell_curve$p,
        observed_events = support$events,
        censoring_rate = support$censoring_rate,
        event_weight_ess = support$event_weight_ess,
        selected_v = selected_v,
        estimated_B_risk_diagnostic_only = selected_risk,
        estimated_B_utility = selected_utility,
        selected_at_upper_grid_boundary = isTRUE(
          selected_v == max(final_v_list)
        ),
        stringsAsFactors = FALSE
      )

      curve_logs[[length(curve_logs) + 1L]] <- data.frame(
        group = group_value,
        leaf_id = as.integer(leaf_name),
        cell = label_lookup[leaf_name],
        v = final_v_list,
        estimated_B_risk_diagnostic_only = cell_curve$risk,
        estimated_B_utility = cell_curve$utility,
        final_calibration_fraction_within_R = cell_curve$p,
        stringsAsFactors = FALSE
      )
    }
  }

  list(
    selected_maps = selected_maps,
    assignment = assignment,
    final_cache = final_cache,
    r_summary = local_bind_rows(r_summaries),
    leaf_summary = local_bind_rows(leaf_summaries),
    cell_curve_log = local_bind_rows(curve_logs)
  )
}


# ------------------------------------------------
# Diagnostics for the audit partition
# ------------------------------------------------

local_audit_partition_diagnostics <- function(
    partition,
    audit_data,
    mdl,
    mdl0,
    xnames,
    alpha) {

  assignment <- local_assign_partition(
    partition = partition,
    mdl = mdl,
    mdl0 = mdl0,
    data = audit_data,
    xnames = xnames,
    alpha = alpha
  )

  audit_cache <- partition$audit_cache

  leaf_rows <- list()

  for (group_value in c(0L, 1L)) {

    group_model <- partition$group_models[[as.character(group_value)]]

    label_lookup <- setNames(
      group_model$leaf_label_map$cell,
      group_model$leaf_label_map$leaf_id
    )

    for (leaf_id in group_model$final_leaves) {

      idx <- which(
        assignment$group == group_value &
          assignment$leaf_id == leaf_id
      )

      support <- local_support_summary(
        cache = audit_cache,
        row_idx = idx
      )

      leaf_rows[[length(leaf_rows) + 1L]] <- data.frame(
        group = group_value,
        leaf_id = leaf_id,
        cell = label_lookup[
          as.character(leaf_id)
        ],
        n_audit = support$n,
        audit_fraction_within_R = if (
          sum(assignment$group == group_value) == 0L
        ) {
          NA_real_
        } else {
          support$n /
            sum(assignment$group == group_value)
        },
        observed_events = support$events,
        censoring_rate = support$censoring_rate,
        event_weight_ess = support$event_weight_ess,
        fraction_G_truncated = support$fraction_G_truncated,
        stringsAsFactors = FALSE
      )
    }
  }

  list(
    assignment = assignment,
    leaf_summary = local_bind_rows(leaf_rows)
  )
}


# ------------------------------------------------
# Logging
# ------------------------------------------------

local_write_mapping_logs <- function(
    mapping_log_dir,
    configuration,
    split_summary,
    partition,
    audit_diagnostics,
    final_calibration_result,
    final_original_row = NULL,
    audit_original_row = NULL) {

  if (is.null(mapping_log_dir)) {
    return(invisible(NULL))
  }

  dir.create(
    mapping_log_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  utils::write.csv(
    configuration,
    file.path(
      mapping_log_dir,
      "mapping_configuration.csv"
    ),
    row.names = FALSE
  )

  utils::write.csv(
    split_summary,
    file.path(
      mapping_log_dir,
      "calibration_split_summary.csv"
    ),
    row.names = FALSE
  )

  utils::write.csv(
    partition$cell_limits,
    file.path(
      mapping_log_dir,
      "data_dependent_cell_limits.csv"
    ),
    row.names = FALSE
  )

  candidate_log <- local_bind_rows(
    lapply(
      partition$group_models,
      `[[`,
      "candidate_log"
    )
  )

  accepted_split_log <- if (
    nrow(candidate_log) > 0L &&
    "accepted" %in% names(candidate_log)
  ) {
    candidate_log[
      candidate_log$accepted %in% TRUE,
      ,
      drop = FALSE
    ]
  } else {
    data.frame()
  }

  partition_rules <- local_bind_rows(
    lapply(
      partition$group_models,
      function(group_model) {
        rules <- group_model$rules
        if (nrow(rules) > 0L) {
          rules$group <- group_model$group
        }
        rules
      }
    )
  )

  accepted_fold_log <- local_bind_rows(
    lapply(
      partition$group_models,
      `[[`,
      "accepted_fold_log"
    )
  )

  utils::write.csv(
    candidate_log,
    file.path(
      mapping_log_dir,
      "candidate_split_log.csv"
    ),
    row.names = FALSE
  )

  utils::write.csv(
    accepted_split_log,
    file.path(
      mapping_log_dir,
      "accepted_split_log.csv"
    ),
    row.names = FALSE
  )

  utils::write.csv(
    partition_rules,
    file.path(
      mapping_log_dir,
      "final_partition_rules.csv"
    ),
    row.names = FALSE
  )

  utils::write.csv(
    accepted_fold_log,
    file.path(
      mapping_log_dir,
      "accepted_split_fold_log.csv"
    ),
    row.names = FALSE
  )

  utils::write.csv(
    audit_diagnostics$leaf_summary,
    file.path(
      mapping_log_dir,
      "audit_leaf_summary.csv"
    ),
    row.names = FALSE
  )

  audit_assignment <- audit_diagnostics$assignment

  if (!is.null(audit_original_row)) {
    audit_assignment$original_calibration_row <- audit_original_row
  }

  utils::write.csv(
    audit_assignment,
    file.path(
      mapping_log_dir,
      "audit_partition_assignments.csv"
    ),
    row.names = FALSE
  )

  utils::write.csv(
    final_calibration_result$leaf_summary,
    file.path(
      mapping_log_dir,
      "final_calibration_leaf_summary.csv"
    ),
    row.names = FALSE
  )

  final_assignment <- final_calibration_result$assignment

  if (!is.null(final_original_row)) {
    final_assignment$original_calibration_row <- final_original_row
  }

  utils::write.csv(
    final_assignment,
    file.path(
      mapping_log_dir,
      "final_calibration_partition_assignments.csv"
    ),
    row.names = FALSE
  )

  utils::write.csv(
    final_calibration_result$r_summary,
    file.path(
      mapping_log_dir,
      "r_level_optimization_summary.csv"
    ),
    row.names = FALSE
  )

  utils::write.csv(
    final_calibration_result$cell_curve_log,
    file.path(
      mapping_log_dir,
      "final_cell_risk_utility_curves.csv"
    ),
    row.names = FALSE
  )

  saveRDS(
    partition,
    file.path(
      mapping_log_dir,
      "local_partition_model.rds"
    )
  )

  invisible(NULL)
}


# ------------------------------------------------
# Main method: local CAMS with IPCW
# ------------------------------------------------

cams_local_ipcw <- function(
    x,
    p,
    len_x,
    xnames,
    data_fit,
    data_calib,
    mdl0,
    alpha,
    use_oracle_sc = FALSE,
    data_audit = NULL,
    data_local_calib = NULL,
    audit_fraction = 0.40,
    calibration_split_seed = 1L,
    audit_v_list = seq(0.01, 0.99, by = 0.01),
    final_v_list = seq(0.001, 0.999, by = 0.001),
    split_features = c(
      "event_lp",
      "logit_G_qalpha",
      "X2",
      "X3",
      "X4"
    ),
    split_quantiles = seq(0.20, 0.80, by = 0.10),
    cv_folds = 4L,
    min_absolute_n = 40L,
    min_group_fraction = 0.15,
    min_events = 30L,
    min_event_ess = 30,
    global_cap = 2,
    min_relative_gain = 0.02,
    min_positive_fold_fraction = 1,
    min_risk_ok_fold_fraction = 1,
    audit_risk_tolerance = 0,
    C_0 = 0,
    mapping_log_dir = NULL) {

  # Keep the original CAMS calling convention.
  # p is not used by the current CAMS implementation either.
  invisible(p)
  invisible(use_oracle_sc)

  newdata <- data.frame(x)
  colnames(newdata) <- xnames

  if (nrow(newdata) != len_x) {
    len_x <- nrow(newdata)
  }

  required_fit_names <- c(
    xnames,
    "censored_T",
    "event"
  )

  missing_fit_names <- setdiff(
    required_fit_names,
    names(data_fit)
  )

  if (length(missing_fit_names) > 0L) {
    stop(
      sprintf(
        "data_fit is missing columns: %s",
        paste(missing_fit_names, collapse = ", ")
      )
    )
  }

  fmla <- stats::as.formula(
    paste(
      "survival::Surv(censored_T, event) ~",
      paste(xnames, collapse = " + ")
    )
  )

  mdl <- survival::survreg(
    fmla,
    data = data_fit,
    dist = "weibull"
  )

  # Either provide both external pieces or let this function split data_calib.
  externally_split <- (
    !is.null(data_audit) ||
      !is.null(data_local_calib)
  )

  if (externally_split) {

    if (
      is.null(data_audit) ||
      is.null(data_local_calib)
    ) {
      stop(
        paste0(
          "Provide both data_audit and data_local_calib, ",
          "or provide neither."
        )
      )
    }

    audit_data <- data_audit
    final_calibration_data <- data_local_calib

    audit_original_row <- NULL
    final_original_row <- NULL

  } else {

    split_result <- local_stratified_calibration_split(
      data_calib = data_calib,
      audit_fraction = audit_fraction,
      split_seed = calibration_split_seed,
      group_name = "X1"
    )

    audit_data <- split_result$audit
    final_calibration_data <- split_result$final_calibration
    audit_original_row <- split_result$audit_original_row
    final_original_row <- split_result$final_original_row
  }

  partition <- local_learn_partition(
    mdl = mdl,
    mdl0 = mdl0,
    audit_data = audit_data,
    xnames = xnames,
    alpha = alpha,
    audit_v_list = audit_v_list,
    split_features = split_features,
    split_quantiles = split_quantiles,
    cv_folds = cv_folds,
    cv_seed = calibration_split_seed,
    min_absolute_n = min_absolute_n,
    min_group_fraction = min_group_fraction,
    min_events = min_events,
    min_event_ess = min_event_ess,
    global_cap = global_cap,
    min_relative_gain = min_relative_gain,
    min_positive_fold_fraction = min_positive_fold_fraction,
    min_risk_ok_fold_fraction = min_risk_ok_fold_fraction,
    audit_risk_tolerance = audit_risk_tolerance
  )

  audit_diagnostics <- local_audit_partition_diagnostics(
    partition = partition,
    audit_data = audit_data,
    mdl = mdl,
    mdl0 = mdl0,
    xnames = xnames,
    alpha = alpha
  )

  final_calibration_result <- local_calibrate_partition_ipcw(
    mdl = mdl,
    mdl0 = mdl0,
    partition = partition,
    final_calibration_data = final_calibration_data,
    xnames = xnames,
    alpha = alpha,
    final_v_list = final_v_list,
    C_0 = C_0
  )

  new_assignment <- local_assign_partition(
    partition = partition,
    mdl = mdl,
    mdl0 = mdl0,
    data = newdata,
    xnames = xnames,
    alpha = alpha
  )

  local_bounds <- rep(
    0,
    nrow(newdata)
  )

  for (group_value in c(0L, 1L)) {

    group_rows <- which(
      newdata$X1 == group_value
    )

    if (length(group_rows) == 0L) {
      next
    }

    selected_map <- final_calibration_result$selected_maps[[as.character(group_value)]]

    if (is.null(selected_map)) {
      next
    }

    group_leaf <- new_assignment$leaf_id[
      group_rows
    ]

    for (leaf_name in names(selected_map$selected_v)) {

      local_rows <- group_rows[
        group_leaf == as.integer(leaf_name)
      ]

      if (length(local_rows) == 0L) {
        next
      }

      selected_v <- selected_map$selected_v[
        leaf_name
      ]

      local_bounds[local_rows] <- local_bound(
        mdl = mdl,
        newdata = newdata[
          local_rows,
          ,
          drop = FALSE
        ],
        v = selected_v
      )
    }
  }

  local_bounds <- pmax(
    local_bounds,
    0
  )

  output <- data.frame(
    "Local-CAMS" = local_bounds,
    check.names = FALSE
  )

  split_summary <- local_bind_rows(
    lapply(
      c(0L, 1L),
      function(group_value) {
        data.frame(
          group = group_value,
          n_total_calibration = if (externally_split) {
            NA_integer_
          } else {
            sum(data_calib$X1 == group_value)
          },
          n_audit = sum(audit_data$X1 == group_value),
          n_final_calibration = sum(
            final_calibration_data$X1 == group_value
          ),
          stringsAsFactors = FALSE
        )
      }
    )
  )

  configuration <- data.frame(
    parameter = c(
      "alpha",
      "audit_fraction",
      "calibration_split_seed",
      "cv_folds",
      "min_absolute_n",
      "min_group_fraction",
      "min_events",
      "min_event_ess",
      "global_cap",
      "min_relative_gain",
      "min_positive_fold_fraction",
      "min_risk_ok_fold_fraction",
      "audit_risk_tolerance",
      "C_0",
      "audit_v_min",
      "audit_v_max",
      "audit_v_count",
      "final_v_min",
      "final_v_max",
      "final_v_count",
      "split_features"
    ),
    value = c(
      alpha,
      audit_fraction,
      calibration_split_seed,
      cv_folds,
      min_absolute_n,
      min_group_fraction,
      min_events,
      min_event_ess,
      global_cap,
      min_relative_gain,
      min_positive_fold_fraction,
      min_risk_ok_fold_fraction,
      audit_risk_tolerance,
      C_0,
      min(audit_v_list),
      max(audit_v_list),
      length(audit_v_list),
      min(final_v_list),
      max(final_v_list),
      length(final_v_list),
      paste(split_features, collapse = ";")
    ),
    stringsAsFactors = FALSE
  )

  local_write_mapping_logs(
    mapping_log_dir = mapping_log_dir,
    configuration = configuration,
    split_summary = split_summary,
    partition = partition,
    audit_diagnostics = audit_diagnostics,
    final_calibration_result = final_calibration_result,
    final_original_row = final_original_row,
    audit_original_row = audit_original_row
  )

  r_diagnostics <- final_calibration_result$r_summary

  if (nrow(r_diagnostics) > 0L) {
    r_diagnostics$method <- "Local-CAMS"
  }

  list(
    output = output,
    times = NA_real_,
    diagnostics = r_diagnostics,
    partition = partition,
    mapping_diagnostics = list(
      calibration_split_summary = split_summary,
      cell_limits = partition$cell_limits,
      candidate_splits = local_bind_rows(
        lapply(
          partition$group_models,
          `[[`,
          "candidate_log"
        )
      ),
      accepted_split_folds = local_bind_rows(
        lapply(
          partition$group_models,
          `[[`,
          "accepted_fold_log"
        )
      ),
      audit_leaf_summary = audit_diagnostics$leaf_summary,
      final_calibration_leaf_summary = final_calibration_result$leaf_summary,
      r_level_optimization_summary = final_calibration_result$r_summary,
      final_cell_risk_utility_curves = final_calibration_result$cell_curve_log,
      audit_assignment = audit_diagnostics$assignment,
      final_calibration_assignment = final_calibration_result$assignment,
      test_assignment = new_assignment
    ),
    event_model = mdl,
    selected_local_levels = lapply(
      final_calibration_result$selected_maps,
      `[[`,
      "selected_v"
    )
  )
}


# Convenient shorter alias
cams_local <- cams_local_ipcw