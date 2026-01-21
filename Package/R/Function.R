#' @keywords internal
"_PACKAGE"

#' MBbioeq - Model-Based Bioequivalence Analysis
#'
#' @name imports
#' @import saemix
#' @import car
#' @import MASS
#' @import VGAM
#' @importFrom stats qnorm qt
#' @importFrom utils read.table
NULL

# =============================================================================
# Model registry (models + metadata embedded)
# =============================================================================

.pk_registry <- list(
  "1cpt" = list(
    model = function(psi, id, xidep) {
      dose <- xidep[, 1]; tim <- xidep[, 2]
      ka <- psi[id, 1]; V <- psi[id, 2]; CL <- psi[id, 3]
      k  <- CL / V
      dose * ka / (V * (ka - k)) * (exp(-k * tim) - exp(-ka * tim))
    },
    par_names = c("ka","V","CL"),
    psi0 = c(1.5, 0.5, 0.04),
    transform = c(1,1,1)
  ),
  "1cpt_tlag" = list(
    model = function(psi, id, xidep) {
      dose <- xidep[, 1]; tim <- xidep[, 2]
      ka <- psi[id, 1]; V <- psi[id, 2]; CL <- psi[id, 3]; tlag <- psi[id, 4]
      k  <- CL / V
      eff_t <- pmax(tim - tlag, 0)
      dose * ka / (V * (ka - k)) * (exp(-k * eff_t) - exp(-ka * eff_t)) * (tim >= tlag)
    },
    par_names = c("ka","V","CL","tlag"),
    psi0 = c(1.5, 0.5, 0.04, 0.5),
    transform = c(1,1,1,1)
  ),
  "2cpt" = list(
    model = function(psi, id, xidep) {
      dose <- xidep[, 1]; tim <- xidep[, 2]
      ka <- psi[id, 1]; Vc <- psi[id, 2]; CL <- psi[id, 3]; Q <- psi[id, 4]; Vp <- psi[id, 5]
      k10 <- CL / Vc; k12 <- Q / Vc; k21 <- Q / Vp
      root  <- sqrt((k10 + k12 + k21)^2 - 4 * k21 * k10)
      alpha <- 0.5 * ((k10 + k12 + k21) + root)
      beta  <- 0.5 * ((k10 + k12 + k21) - root)
      A <- ka * (alpha - k21) / (Vc * (ka - alpha) * (alpha - beta))
      B <- ka * (beta  - k21) / (Vc * (ka - beta)  * (beta  - alpha))
      dose * (A * (exp(-alpha * tim) - exp(-ka * tim)) +
                B * (exp(-beta  * tim) - exp(-ka * tim)))
    },
    par_names = c("ka","Vc","CL","Q","Vp"),
    psi0 = c(1.5, 10, 1, 1, 20),
    transform = c(1,1,1,1,1)
  ),
  "2cpt_tlag" = list(
    model = function(psi, id, xidep) {
      dose <- xidep[, 1]; tim <- xidep[, 2]
      ka <- psi[id, 1]; Vc <- psi[id, 2]; CL <- psi[id, 3]; Q <- psi[id, 4]; Vp <- psi[id, 5]; tlag <- psi[id, 6]
      k10 <- CL / Vc; k12 <- Q / Vc; k21 <- Q / Vp
      eff_t <- pmax(tim - tlag, 0)
      root  <- sqrt((k10 + k12 + k21)^2 - 4 * k21 * k10)
      alpha <- 0.5 * ((k10 + k12 + k21) + root)
      beta  <- 0.5 * ((k10 + k12 + k21) - root)
      A <- ka * (alpha - k21) / (Vc * (ka - alpha) * (alpha - beta))
      B <- ka * (beta  - k21) / (Vc * (ka - beta)  * (beta  - alpha))
      (dose * (A * (exp(-alpha * eff_t) - exp(-ka * eff_t)) +
                 B * (exp(-beta  * eff_t)  - exp(-ka * eff_t)))) * (tim >= tlag)
    },
    par_names = c("ka","Vc","CL","Q","Vp","tlag"),
    psi0 = c(1.5, 10, 1, 1, 20, 0.5),
    transform = c(1,1,1,1,1,1)
  )
)

# =============================================================================
# Flexible builder for omega and error model
# =============================================================================

.build_saemix_model <- function(model,
                                psi0 = NULL,
                                omega_init = NULL,
                                covariance_structure = "diagonal",
                                error_model = "combined",
                                error_init = NULL,
                                verbose = FALSE) {
  stopifnot(model %in% names(.pk_registry))
  spec <- .pk_registry[[model]]
  k <- length(spec$par_names)

  if (is.null(psi0)) psi0 <- spec$psi0
  if (length(psi0) != k) stop("`psi0` must have ", k, " values for model '", model, "'.")
  psi0_mat <- matrix(psi0, ncol = k, byrow = TRUE, dimnames = list(NULL, spec$par_names))

  if (is.null(omega_init)) {
    omega_init <- rep(0.05, k)
  }

  if (is.matrix(omega_init)) {
    if (!all(dim(omega_init) == c(k, k))) {
      stop("`omega_init` must be a ", k, "x", k, " matrix.")
    }
    omega_mat <- omega_init
    cov_model <- (omega_mat != 0) * 1
  } else {
    if (length(omega_init) != k) {
      stop("`omega_init` must have ", k, " values or be a ", k, "x", k, " matrix.")
    }
    omega_mat <- diag(omega_init, k)
    if (is.character(covariance_structure)) {
      cov_model <- switch(covariance_structure,
                          "diagonal" = diag(1, k),
                          "full" = matrix(1, k, k),
                          stop("Unknown covariance_structure: ", covariance_structure))
    } else if (is.matrix(covariance_structure)) {
      cov_model <- covariance_structure
    } else {
      stop("covariance_structure must be 'diagonal', 'full', or a matrix")
    }
  }

  if (is.null(error_init)) {
    error_init <- switch(error_model,
                         "combined" = c(0.1, 0.1),
                         "proportional" = 0.1,
                         "constant" = 0.1,
                         stop("Unknown error_model: ", error_model))
  }

  expected_length <- if (error_model == "combined") 2 else 1
  if (length(error_init) != expected_length) {
    stop("error_init must have ", expected_length, " value(s) for error_model '", error_model, "'")
  }

  saemix::saemixModel(
    model = spec$model,
    description = paste0("PK ", model),
    psi0 = psi0_mat,
    transform.par = spec$transform,
    covariate.model = matrix(1, nrow = 1, ncol = k),
    covariance.model = cov_model,
    omega.init = omega_mat,
    error.model = error_model,
    error.init = error_init,
    verbose = verbose
  )
}

# =============================================================================
# Expression généralisée pour beta_Cmax
# =============================================================================

.get_fe_indices <- function(model) {
  par_names <- .pk_registry[[model]]$par_names
  k <- length(par_names)
  fe_indices <- 1:(2*k)
  fe_names <- character(2*k)
  for (i in 1:k) {
    fe_names[2*i-1] <- par_names[i]
    fe_names[2*i] <- paste0("beta_", par_names[i])
  }
  list(indices = fe_indices, names = fe_names)
}

.beta_cmax_expr <- function(model) {
  switch(model,
         "1cpt" = "-beta_V-(log((ka*V)/CL)+beta_ka+beta_V-beta_CL)*CL*exp(beta_CL)/(ka*V*exp(beta_ka+beta_V)-CL*exp(beta_CL))+(CL/(ka*V-CL))*log(ka*V/CL)",
         "1cpt_tlag" = "-beta_V-(log((ka*V)/CL)+beta_ka+beta_V-beta_CL)*CL*exp(beta_CL)/(ka*V*exp(beta_ka+beta_V)-CL*exp(beta_CL))+(CL/(ka*V-CL))*log(ka*V/CL)",
         "2cpt" = NA_character_,
         "2cpt_tlag" = NA_character_,
         stop("Unsupported model for Cmax delta-method.")
  )
}

# =============================================================================
# Data prep
# =============================================================================

#' Prepare dataset for saemix (parallel-arm T/R)
#' @export
prep_parallel_TR <- function(tab) {
  stopifnot(is.data.frame(tab))
  names(tab) <- tolower(names(tab))

  col_map <- list(
    id = c("id", "subject", "subj", "sid"),
    time = c("time", "t", "hours", "hr"),
    concentration = c("concentration", "conc", "dv", "y", "obs"),
    treatment = c("treatment", "treat", "trt", "tr", "group", "arm"),
    dose = c("dose", "amt")
  )

  rename_one <- function(std) {
    vars <- col_map[[std]]
    hits <- intersect(vars, names(tab))
    if (length(hits) == 0 && std %in% c("id", "time", "concentration", "treatment")) {
      stop("Missing required column: '", std, "'. Allowed variants: ", paste(vars, collapse = ", "), ".")
    }
    if (length(hits) > 1) {
      stop("Ambiguous columns for '", std, "': ", paste(hits, collapse = ", "))
    }
    if (length(hits) == 1) names(tab)[names(tab) == hits] <<- std
  }
  lapply(names(col_map), rename_one)

  req <- c("id", "time", "concentration", "treatment")
  miss <- setdiff(req, names(tab))
  if (length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse = ", "))

  tr_raw <- toupper(as.character(tab$treatment))
  ok <- tr_raw %in% c("T", "R")
  if (!all(ok)) {
    bad <- unique(tr_raw[!ok])
    stop("Invalid treatment values: ", paste(bad, collapse = ", "))
  }
  tab$Treat <- ifelse(tr_raw == "T", 1L, 0L)

  if (!"dose" %in% names(tab)) tab$dose <- NA_real_

  tab <- tab[order(tab$id, tab$time), ]

  ff_subject <- function(x) {
    idx <- cummax(ifelse(!is.na(x), seq_along(x), 0))
    ifelse(idx == 0, NA_real_, x[idx])
  }
  tab$Dose <- ave(tab$dose, tab$id, FUN = ff_subject)

  if (any(is.na(tab$Dose))) {
    bad_ids <- unique(tab$id[is.na(tab$Dose)])
    stop("Could not impute Dose for subject(s): ", paste(bad_ids, collapse = ", "))
  }

  res <- data.frame(
    Id = tab$id,
    Dose = tab$Dose,
    Time = tab$time,
    Concentration = tab$concentration,
    Treat = tab$Treat
  )

  stopifnot(identical(names(res), c("Id","Dose","Time","Concentration","Treat")))
  res[order(res$Id, res$Time), ]
}

# =============================================================================
# Shared fitting function for both TOST and BOT
# =============================================================================

.fit_MB_model <- function(data_path, model, psi0, omega_init,
                          covariance_structure, error_model, error_init,
                          nb_chains, nb_iter, verbose) {

  tab <- utils::read.table(data_path, header = TRUE, sep = " ", dec = ".",
                           na.strings = c(".", "NA", ""))
  tab <- prep_parallel_TR(tab)
  tab <- tab[!is.na(tab$Concentration), ]

  nb_t <- sum(tab$Id == tab$Id[1])
  n <- nrow(tab) / nb_t

  if (!verbose) {
    invisible(capture.output({
      saemix.data <- saemix::saemixData(
        name.data = tab, header = TRUE, sep = " ", na = NA,
        name.group = c("Id"),
        name.predictors = c("Dose", "Time"),
        name.covariates = c("Treat"),
        name.response = c("Concentration"),
        name.X = "Time",
        units = list(x = "hr", y = "mg/L")
      )

      saemix.model <- .build_saemix_model(
        model, psi0, omega_init, covariance_structure, error_model, error_init, verbose = FALSE
      )
    }))
  } else {
    saemix.data <- saemix::saemixData(
      name.data = tab, header = TRUE, sep = " ", na = NA,
      name.group = c("Id"),
      name.predictors = c("Dose", "Time"),
      name.covariates = c("Treat"),
      name.response = c("Concentration"),
      name.X = "Time",
      units = list(x = "hr", y = "mg/L")
    )

    saemix.model <- .build_saemix_model(
      model, psi0, omega_init, covariance_structure, error_model, error_init
    )
  }

  k_params <- length(.pk_registry[[model]]$par_names)

  if (!verbose) {
    invisible(capture.output({
      mod <- suppressMessages(suppressWarnings(
        saemix::saemix(saemix.model, saemix.data,
                       list(nb.chains = nb_chains,
                            nbiter.saemix = nb_iter,
                            print = FALSE,
                            displayProgress = FALSE,
                            save = FALSE,
                            save.graphs = FALSE))
      ))
    }, type = "output"))
  } else {
    mod <- saemix::saemix(saemix.model, saemix.data,
                          list(nb.chains = nb_chains,
                               nbiter.saemix = nb_iter))
  }

  fe <- mod@results@fixed.effects
  sefe <- mod@results@se.fixed
  FIM <- mod@results@fim
  varcov <- tryCatch(solve(FIM), error = function(e) MASS::ginv(FIM))

  par_names <- .pk_registry[[model]]$par_names
  idx_cl <- which(toupper(par_names) == "CL")
  if (length(idx_cl) != 1) stop("Could not locate 'CL' in parameter list")
  beta_cl_pos <- 2 * idx_cl

  beta_auc_est <- -fe[beta_cl_pos]
  beta_auc_se <- sefe[beta_cl_pos]

  cmax_expr <- .beta_cmax_expr(model)
  beta_cmax_est <- NA_real_
  secmax <- NA_real_

  if (!is.na(cmax_expr)) {
    fe_info <- .get_fe_indices(model)
    estimates <- fe[fe_info$indices]
    names(estimates) <- fe_info$names

    var_subset <- varcov[fe_info$indices, fe_info$indices, drop = FALSE]
    colnames(var_subset) <- rownames(var_subset) <- fe_info$names

    SE_deltam <- car::deltaMethod(estimates, cmax_expr, vcov. = var_subset)
    beta_cmax_est <- SE_deltam[, 1]
    secmax <- SE_deltam[, 2]
  }

  list(
    fit = mod,
    model = model,
    n = n,
    nb_t = nb_t,
    k = k_params,
    varcov = varcov,
    beta_auc_est = beta_auc_est,
    beta_auc_se = beta_auc_se,
    beta_cmax_est = beta_cmax_est,
    beta_cmax_se = secmax
  )
}

# =============================================================================
# TOST FUNCTION
# =============================================================================

#' Run Model-Based TOST Bioequivalence Analysis
#'
#' @description Performs Two One-Sided Tests (TOST) for bioequivalence
#' using a model-based approach with SAEMIX. Tests both AUC and Cmax (when available).
#' Provides both asymptotic normal and Gallant-corrected inference.
#'
#' @param data_path Path to dataset file
#' @param model One of: "1cpt", "1cpt_tlag", "2cpt", "2cpt_tlag"
#' @param psi0 Initial parameter values (NULL = defaults)
#' @param omega_init Initial random effects variance (vector or matrix)
#' @param covariance_structure "diagonal" or "full" or custom matrix
#' @param error_model "combined", "proportional", or "constant"
#' @param error_init Initial error parameters (NULL = auto)
#' @param nb_chains Number of SAEM chains
#' @param nb_iter Vector c(iter_saem, iter_mcmc)
#' @param delta Equivalence margin (default log(1.25) for 0.80-1.25 ratio)
#' @param alpha One-sided significance level (default 0.05 for 90 percent CI)
#' @param verbose Print SAEMIX output (default FALSE)
#' @return List containing fit object, estimates, CIs, and TOST decisions for AUC and Cmax
#' @export
#' @examples
#' \dontrun{
#' result <- run_MBTost("data.txt", model = "1cpt")
#' }
run_MBTost <- function(data_path,
                       model = "1cpt",
                       psi0 = NULL,
                       omega_init = NULL,
                       covariance_structure = "diagonal",
                       error_model = "combined",
                       error_init = NULL,
                       nb_chains = 10,
                       nb_iter = c(300, 100),
                       delta = log(1.25),
                       alpha = 0.05,
                       verbose = FALSE) {

  res <- .fit_MB_model(data_path, model, psi0, omega_init,
                       covariance_structure, error_model, error_init,
                       nb_chains, nb_iter, verbose)

  quant <- stats::qnorm(1 - alpha)
  quant_gallant <- stats::qt(1 - alpha, df = res$n - res$k)

  # ========== TOST for beta_AUC ==========
  ic_auc_normal_lower <- res$beta_auc_est - quant * res$beta_auc_se
  ic_auc_normal_upper <- res$beta_auc_est + quant * res$beta_auc_se

  W1_auc <- (res$beta_auc_est + delta) / res$beta_auc_se
  W2_auc <- (res$beta_auc_est - delta) / res$beta_auc_se
  tost_auc_normal <- (W1_auc >= quant) & (W2_auc <= -quant)

  corrected_se_beta_auc <- res$beta_auc_se * sqrt(res$n / quant_gallant)
  ic_auc_gallant_lower <- res$beta_auc_est - quant_gallant * corrected_se_beta_auc
  ic_auc_gallant_upper <- res$beta_auc_est + quant_gallant * corrected_se_beta_auc

  W1_gallant_auc <- (res$beta_auc_est + delta) / corrected_se_beta_auc
  W2_gallant_auc <- (res$beta_auc_est - delta) / corrected_se_beta_auc
  tost_auc_gallant <- (W1_gallant_auc >= quant_gallant) & (W2_gallant_auc <= -quant_gallant)

  # ========== TOST for beta_Cmax ==========
  ic_cmax_normal_lower <- NA_real_; ic_cmax_normal_upper <- NA_real_
  ic_cmax_gallant_lower <- NA_real_; ic_cmax_gallant_upper <- NA_real_
  tost_cmax_normal <- NA; tost_cmax_gallant <- NA
  corrected_se_cmax <- NA_real_

  if (!is.na(res$beta_cmax_est)) {
    ic_cmax_normal_lower <- res$beta_cmax_est - quant * res$beta_cmax_se
    ic_cmax_normal_upper <- res$beta_cmax_est + quant * res$beta_cmax_se

    W1cmax <- (res$beta_cmax_est + delta) / res$beta_cmax_se
    W2cmax <- (res$beta_cmax_est - delta) / res$beta_cmax_se
    tost_cmax_normal <- (W1cmax >= quant) & (W2cmax <= -quant)

    corrected_se_cmax <- res$beta_cmax_se * sqrt(res$n / quant_gallant)
    ic_cmax_gallant_lower <- res$beta_cmax_est - quant_gallant * corrected_se_cmax
    ic_cmax_gallant_upper <- res$beta_cmax_est + quant_gallant * corrected_se_cmax

    W1_gallant_cmax <- (res$beta_cmax_est + delta) / corrected_se_cmax
    W2_gallant_cmax <- (res$beta_cmax_est - delta) / corrected_se_cmax
    tost_cmax_gallant <- (W1_gallant_cmax >= quant_gallant) & (W2_gallant_cmax <= -quant_gallant)
  }

  # ========== Display results ==========
  cat("\n")
  cat("════════════════════════════════════════════════════════════\n")
  cat("        TOST BIOEQUIVALENCE ANALYSIS RESULTS\n")
  cat("════════════════════════════════════════════════════════════\n")
  cat("Model:", model, " | n =", res$n, " | k =", res$k, " | δ =", round(delta, 4), "\n")
  cat("════════════════════════════════════════════════════════════\n\n")

  cat("AUC RATIO (T/R):\n")
  cat("  Estimate:  ", sprintf("%7.4f", exp(res$beta_auc_est)), "\n")
  cat("  90% CI (Normal):  [", sprintf("%7.4f", exp(ic_auc_normal_lower)), ", ",
      sprintf("%7.4f", exp(ic_auc_normal_upper)), "]\n", sep = "")
  cat("  90% CI (Gallant): [", sprintf("%7.4f", exp(ic_auc_gallant_lower)), ", ",
      sprintf("%7.4f", exp(ic_auc_gallant_upper)), "]\n", sep = "")
  cat("\n")
  cat("  Normal:  ", ifelse(tost_auc_normal, "Bioequivalent", "Not bioequivalent"), "\n")
  cat("  Gallant: ", ifelse(tost_auc_gallant, "Bioequivalent", "Not bioequivalent"), "\n")
  cat("\n")

  cat("CMAX RATIO (T/R):\n")
  if (!is.na(res$beta_cmax_est)) {
    cat("  Estimate:  ", sprintf("%7.4f", exp(res$beta_cmax_est)), "\n")
    cat("  90% CI (Normal):  [", sprintf("%7.4f", exp(ic_cmax_normal_lower)), ", ",
        sprintf("%7.4f", exp(ic_cmax_normal_upper)), "]\n", sep = "")
    cat("  90% CI (Gallant): [", sprintf("%7.4f", exp(ic_cmax_gallant_lower)), ", ",
        sprintf("%7.4f", exp(ic_cmax_gallant_upper)), "]\n", sep = "")
    cat("\n")
    cat("  Normal:  ", ifelse(tost_cmax_normal, "Bioequivalent", "Not bioequivalent"), "\n")
    cat("  Gallant: ", ifelse(tost_cmax_gallant, "Bioequivalent", "Not bioequivalent"), "\n")
  } else {
    cat("  Not available for this model\n")
  }
  cat("\n")
  cat("════════════════════════════════════════════════════════════\n\n")

  invisible(list(
    fit = res$fit,
    model = model,
    n = res$n,
    k = res$k,
    varcov = res$varcov,

    auc = list(
      estimate = exp(res$beta_auc_est),
      normal = list(
        ci_lower = exp(ic_auc_normal_lower),
        ci_upper = exp(ic_auc_normal_upper),
        decision = tost_auc_normal
      ),
      gallant = list(
        ci_lower = exp(ic_auc_gallant_lower),
        ci_upper = exp(ic_auc_gallant_upper),
        decision = tost_auc_gallant
      )
    ),

    cmax = list(
      estimate = exp(res$beta_cmax_est),
      normal = list(
        ci_lower = exp(ic_cmax_normal_lower),
        ci_upper = exp(ic_cmax_normal_upper),
        decision = tost_cmax_normal
      ),
      gallant = list(
        ci_lower = exp(ic_cmax_gallant_lower),
        ci_upper = exp(ic_cmax_gallant_upper),
        decision = tost_cmax_gallant
      )
    )
  ))
}

# =============================================================================
# BOT FUNCTION
# =============================================================================

#' Run Model-Based BOT Bioequivalence Analysis
#'
#' @description Performs Bayesian-Optimal Test (BOT) for bioequivalence
#' using a model-based approach with SAEMIX. Tests both AUC and Cmax (when available).
#' Uses the folded normal distribution to determine the critical value.
#'
#' @param data_path Path to dataset file
#' @param model One of: "1cpt", "1cpt_tlag", "2cpt", "2cpt_tlag"
#' @param psi0 Initial parameter values (NULL = defaults)
#' @param omega_init Initial random effects variance (vector or matrix)
#' @param covariance_structure "diagonal" or "full" or custom matrix
#' @param error_model "combined", "proportional", or "constant"
#' @param error_init Initial error parameters (NULL = auto)
#' @param nb_chains Number of SAEM chains
#' @param nb_iter Vector c(iter_saem, iter_mcmc)
#' @param delta Equivalence margin (default log(1.25) for 0.80-1.25 ratio)
#' @param alpha One-sided significance level (default 0.05)
#' @param verbose Print SAEMIX output (default FALSE)
#' @return List containing fit object, estimates, test statistics, critical values, and BOT decisions for AUC and Cmax
#' @export
#' @examples
#' \dontrun{
#' result <- run_MBbot("data.txt", model = "1cpt")
#' }
run_MBbot <- function(data_path,
                      model = "1cpt",
                      psi0 = NULL,
                      omega_init = NULL,
                      covariance_structure = "diagonal",
                      error_model = "combined",
                      error_init = NULL,
                      nb_chains = 10,
                      nb_iter = c(300, 100),
                      delta = log(1.25),
                      alpha = 0.05,
                      verbose = FALSE) {

  res <- .fit_MB_model(data_path, model, psi0, omega_init,
                       covariance_structure, error_model, error_init,
                       nb_chains, nb_iter, verbose)

  # ========== BOT pour AUC ==========
  stat_auc <- abs(res$beta_auc_est)
  limit_auc <- VGAM::qfoldnorm(alpha, delta, res$beta_auc_se)
  bot_auc_decision <- stat_auc <= limit_auc

  # ========== BOT for Cmax ==========
  stat_cmax <- NA_real_
  limit_cmax <- NA_real_
  bot_cmax_decision <- NA

  if (!is.na(res$beta_cmax_est)) {
    stat_cmax <- abs(res$beta_cmax_est)
    limit_cmax <- VGAM::qfoldnorm(alpha, delta, res$beta_cmax_se)
    bot_cmax_decision <- stat_cmax <= limit_cmax
  }

  # ========== Display results ==========
  cat("\n")
  cat("════════════════════════════════════════════════════════════\n")
  cat("         BOT BIOEQUIVALENCE ANALYSIS RESULTS\n")
  cat("════════════════════════════════════════════════════════════\n")
  cat("Model:", model, " | n =", res$n, " | k =", res$k, " | δ =", round(delta, 4), "\n")
  cat("════════════════════════════════════════════════════════════\n\n")

  cat("AUC RATIO (T/R):\n")
  cat("  Estimate:     ", sprintf("%7.4f", exp(res$beta_auc_est)), "\n")
  cat("  Test statistic: ", sprintf("%7.4f", stat_auc), "\n")
  cat("  Critical value: ", sprintf("%7.4f", limit_auc), "\n")
  cat("  Decision:       ", ifelse(bot_auc_decision, "Bioequivalent", "Not bioequivalent"), "\n")
  cat("\n")

  cat("CMAX RATIO (T/R):\n")
  if (!is.na(res$beta_cmax_est)) {
    cat("  Estimate:       ", sprintf("%7.4f", exp(res$beta_cmax_est)), "\n")
    cat("  Test statistic: ", sprintf("%7.4f", stat_cmax), "\n")
    cat("  Critical value: ", sprintf("%7.4f", limit_cmax), "\n")
    cat("  Decision:       ", ifelse(bot_cmax_decision, "Bioequivalent", "Not bioequivalent"), "\n")
  } else {
    cat("  Not available for this model\n")
  }
  cat("\n")
  cat("════════════════════════════════════════════════════════════\n\n")

  invisible(list(
    fit = res$fit,
    model = model,
    n = res$n,
    k = res$k,
    varcov = res$varcov,

    auc = list(
      estimate = exp(res$beta_auc_est),
      statistic = stat_auc,
      critical_value = limit_auc,
      decision = bot_auc_decision
    ),

    cmax = list(
      estimate = exp(res$beta_cmax_est),
      statistic = stat_cmax,
      critical_value = limit_cmax,
      decision = bot_cmax_decision
    )
  ))
}
