#' @keywords internal
"_PACKAGE"

#' 1CPT TOST Toolkit Import Declarations
#'
#' @name imports
#' @import saemix
#' @import car
#' @import MASS
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
# Builders & helpers
# =============================================================================

# Build a saemix model from the registry
#' @noRd
.build_saemix_model <- function(model, psi0 = NULL, omega_init = NULL) {
  stopifnot(model %in% names(.pk_registry))
  spec <- .pk_registry[[model]]
  k <- length(spec$par_names)

  # Initial means
  if (is.null(psi0)) psi0 <- spec$psi0
  if (length(psi0) != k) stop("`psi0` must have ", k, " values for model '", model, "'.")
  psi0_mat <- matrix(psi0, ncol = k, byrow = TRUE, dimnames = list(NULL, spec$par_names))

  # Random-effects covariance (omega)
  if (is.null(omega_init)) omega_init <- rep(0.05, k)
  if (is.matrix(omega_init)) {
    if (!all(dim(omega_init) == c(k, k))) stop("`omega_init` must be a ", k, "x", k, " matrix.")
    omega_mat <- omega_init
  } else {
    if (length(omega_init) != k) stop("`omega_init` must have ", k, " values or be a ", k, "x", k, " matrix.")
    omega_mat <- diag(omega_init, k)
  }

  saemix::saemixModel(
    model = spec$model,
    description = paste0("PK ", model),
    psi0 = psi0_mat,
    transform.par = spec$transform,
    covariate.model = matrix(1, nrow = 1, ncol = k),  # Treat as covariate on all PK params
    covariance.model = diag(1, k),
    omega.init = omega_mat,
    error.model = "combined",
    error.init = c(0.1, 0.1)
  )
}

# Expression for delta-method β_Cmax per model
# Only 1CPT provided; others placeholder (NA)
#' @noRd
.beta_cmax_expr <- function(model) {
  switch(model,
         "1cpt"      = "-beta_V-(log((Ka*V)/Cl)+beta_Ka+beta_V-beta_Cl)*Cl*exp(beta_Cl)/(Ka*V*exp(beta_Ka+beta_V)-Cl*exp(beta_Cl))+(Cl/(Ka*V-Cl))*log(Ka*V/Cl)",
         "1cpt_tlag" = NA_character_,  # TODO
         "2cpt"      = NA_character_,  # TODO
         "2cpt_tlag" = NA_character_,  # TODO
         stop("Unsupported model for Cmax delta-method.")
  )
}


# =============================================================================
# Public: data prep (T/R, column checks, dose forward-fill)
# =============================================================================

#' Prepare dataset for saemix (parallel-arm T/R)
#'
#' Required columns (case-insensitive): Id/Subject, Time/T, Concentration/Conc/DV, Treatment/Tr/Group.
#' Optional column: Dose (if missing or NA, it is forward-filled within subject).
#'
#' Rules:
#' - Validates names, renames to standard, and enforces T/R only for treatment.
#' - Forward-fills missing Dose by subject (value "from above") without explicit loops.
#' - Sorts by Id, then Time; returns columns exactly: Id, Dose, Time, Concentration, Treat.
#' - Treat is 0 for R and 1 for T.
#'
#' @param tab data.frame with the columns above (any order/case).
#' @return data.frame with columns: Id, Dose, Time, Concentration, Treat (0=R, 1=T)
#' @export
prep_parallel_TR <- function(tab) {
  stopifnot(is.data.frame(tab))

  # Normalize names
  names(tab) <- tolower(names(tab))

  # Allowed variants
  col_map <- list(
    id            = c("id", "subject", "subj", "sid"),
    time          = c("time", "t", "hours", "hr"),
    concentration = c("concentration", "conc", "dv", "y", "obs"),
    treatment     = c("treatment", "treat", "trt", "tr", "group", "arm"),
    dose          = c("dose", "amt")
  )

  # Unique match & rename
  rename_one <- function(std) {
    vars <- col_map[[std]]
    hits <- intersect(vars, names(tab))
    if (length(hits) == 0 && std %in% c("id", "time", "concentration", "treatment")) {
      stop("Missing required column: '", std, "'. Allowed variants: ", paste(vars, collapse = ", "), ".")
    }
    if (length(hits) > 1) {
      stop("Ambiguous columns for '", std, "': ", paste(hits, collapse = ", "),
           ". Keep only one of these.")
    }
    if (length(hits) == 1) names(tab)[names(tab) == hits] <<- std
  }
  lapply(names(col_map), rename_one)

  # Final presence check
  req <- c("id", "time", "concentration", "treatment")
  miss <- setdiff(req, names(tab))
  if (length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse = ", "))

  # Enforce T/R only (case-insensitive)
  tr_raw <- toupper(as.character(tab$treatment))
  ok <- tr_raw %in% c("T", "R")
  if (!all(ok)) {
    bad <- unique(tr_raw[!ok])
    stop("Invalid treatment values: ", paste(bad, collapse = ", "),
         ". Only 'T' and 'R' are allowed.")
  }
  tab$Treat <- ifelse(tr_raw == "T", 1L, 0L)

  # Ensure Dose exists
  if (!"dose" %in% names(tab)) tab$dose <- NA_real_

  # Sort for consistent forward-fill
  tab <- tab[order(tab$id, tab$time), ]

  # Forward-fill Dose within subject (no explicit loops)
  ff_subject <- function(x) {
    idx <- cummax(ifelse(!is.na(x), seq_along(x), 0))
    ifelse(idx == 0, NA_real_, x[idx])
  }
  tab$Dose <- ave(tab$dose, tab$id, FUN = ff_subject)

  # If any subject still has all-NA dose, fail explicitly
  if (any(is.na(tab$Dose))) {
    bad_ids <- unique(tab$id[is.na(tab$Dose)])
    stop("Could not impute Dose for subject(s): ",
         paste(bad_ids, collapse = ", "),
         " (no non-missing dose found within those subjects).")
  }

  # Build final result in exact order
  res <- data.frame(
    Id = tab$id,
    Dose = tab$Dose,
    Time = tab$time,
    Concentration = tab$concentration,
    Treat = tab$Treat
  )

  # Sanity checks and return
  stopifnot(identical(names(res), c("Id","Dose","Time","Concentration","Treat")))
  if (!all(res$Treat %in% c(0L, 1L)))
    stop("Treat must contain only 0 (R) or 1 (T).")

  res[order(res$Id, res$Time), ]
}


# =============================================================================
# Public: TOST runner (multi-model)
# =============================================================================

#' Run TOST analysis (1cpt / 1cpt_tlag / 2cpt / 2cpt_tlag)
#'
#' Fits the selected PK model with SAEM, computes delta-method for Cmax when
#' available (1cpt implemented), and performs TOST (Normal/Student/Gallant)
#' for β_AUC via β_CL. Gallant df is dynamic: df = n - k, with k = #PK params.
#'
#' @param data_path Path to dataset file.
#' @param model One of: "1cpt" (default), "1cpt_tlag", "2cpt", "2cpt_tlag".
#' @param psi0 Numeric vector of length k (initial means). If NULL, registry defaults are used.
#' @param omega_init Numeric vector length k (diag variances) or kxk matrix. If NULL, small diag used.
#' @param default_dose (deprecated) ignored; dose is forward-filled instead.
#' @param nb_chains Number of SAEM chains.
#' @param nb_iter Vector c(iter_saem, iter_mcmc).
#' @param delta Equivalence margin (default log(1.25)).
#' @param alpha One-sided alpha (default 0.05).
#' @return A list with fit, n, nb_t, k, varcov, and TOST results for AUC and Cmax.
#' @export
run_tost <- function(data_path,
                     model = "1cpt",
                     psi0 = NULL,
                     omega_init = NULL,
                     default_dose = NULL,
                     nb_chains = 10,
                     nb_iter = c(300, 100),
                     delta = log(1.25),
                     alpha = 0.05) {

  quant <- stats::qnorm(1 - alpha)

  # Load & prepare data
  tab <- utils::read.table(
    data_path, header = TRUE, sep = " ", dec = ".",
    na.strings = c(".", "NA", "")
  )
  tab <- prep_parallel_TR(tab)

  # Remove dosing rows with NA conc (if present)
  tab <- tab[!is.na(tab$Concentration), ]

  # n = subjects, nb_t = timepoints per subject (assumes balanced)
  nb_t <- sum(tab$Id == tab$Id[1])
  n <- nrow(tab) / nb_t

  saemix.data <- saemix::saemixData(
    name.data = tab,
    header = TRUE, sep = " ", na = NA,
    name.group = c("Id"),
    name.predictors = c("Dose", "Time"),
    name.covariates = c("Treat"),
    name.response = c("Concentration"),
    name.X = "Time",
    units = list(x = "hr", y = "mg/L")
  )

  saemix.model <- .build_saemix_model(model, psi0 = psi0, omega_init = omega_init)
  k_params <- length(.pk_registry[[model]]$par_names)

  mod <- saemix::saemix(saemix.model, saemix.data,
                        list(nb.chains = nb_chains, nbiter.saemix = nb_iter))

  fe   <- mod@results@fixed.effects
  sefe <- mod@results@se.fixed
  FIM  <- mod@results@fim
  varcov <- tryCatch(solve(FIM), error = function(e) MASS::ginv(FIM))

  # Locate beta_CL dynamically (Treat effect on CL)
  par_names <- .pk_registry[[model]]$par_names
  idx_cl <- which(toupper(par_names) == "CL")
  if (length(idx_cl) != 1)
    stop("Could not locate 'CL' in parameter list for model '", model, "'.")
  # FE order with 1 covariate Treat is: p1, beta_p1, p2, beta_p2, ...
  beta_cl_pos <- 2 * idx_cl

  beta_cl_est <- -fe[beta_cl_pos]
  beta_cl_se  <-  sefe[beta_cl_pos]

  # Student & Gallant quantiles
  quant_stud    <- stats::qt(1 - alpha, df = (n * nb_t - 5))
  quant_gallant <- stats::qt(1 - alpha, df = n - k_params)

  # -------- TOST on beta_AUC (via beta_CL)
  W1_auc <- ( beta_cl_est + delta) / beta_cl_se
  W2_auc <- ( beta_cl_est - delta) / beta_cl_se
  tost_auc_normal  <- (W1_auc >= quant) & (W2_auc <= -quant)
  tost_auc_student <- (W1_auc >= quant_stud) & (W2_auc <= -quant_stud)

  corrected_se_beta_cl <- beta_cl_se * sqrt(n / quant_gallant)
  W1_gallant_auc <- (beta_cl_est + delta) / corrected_se_beta_cl
  W2_gallant_auc <- (beta_cl_est - delta) / corrected_se_beta_cl
  tost_auc_gallant <- (W1_gallant_auc >= quant_gallant) & (W2_gallant_auc <= -quant_gallant)

  # -------- Cmax delta-method (available for 1cpt only for now)
  cmax_expr <- .beta_cmax_expr(model)
  estimates <- NULL; secmax <- NA_real_; beta_cmax_estim <- NA_real_
  if (!is.na(cmax_expr)) {
    # Expect first six FE: Ka, beta_Ka, V, beta_V, Cl, beta_Cl
    est_idx <- 1:6
    estimates <- c(fe[est_idx])
    names(estimates) <- c("Ka","beta_Ka","V","beta_V","Cl","beta_Cl")
    var6 <- varcov[est_idx, est_idx, drop = FALSE]
    colnames(var6) <- rownames(var6) <- names(estimates)

    SE_deltam <- car::deltaMethod(estimates, cmax_expr, vcov. = var6)
    beta_cmax_estim <- SE_deltam[, 1]
    secmax <- SE_deltam[, 2]

    # TOST Cmax
    W1cmax <- (beta_cmax_estim + delta) / secmax
    W2cmax <- (beta_cmax_estim - delta) / secmax
    tost_cmax_normal  <- (W1cmax >= quant) & (W2cmax <= -quant)
    tost_cmax_student <- (W1cmax >= quant_stud) & (W2cmax <= -quant_stud)

    corrected_se_beta_Cmax <- secmax * sqrt(n / quant_gallant)
    W1_gallant_Cmax <- (beta_cmax_estim + delta) / corrected_se_beta_Cmax
    W2_gallant_Cmax <- (beta_cmax_estim - delta) / corrected_se_beta_Cmax
    tost_cmax_gallant <- (W1_gallant_Cmax >= quant_gallant) & (W2_gallant_Cmax <= -quant_gallant)
  } else {
    tost_cmax_normal <- tost_cmax_student <- tost_cmax_gallant <- NA
    W1cmax <- W2cmax <- corrected_se_beta_Cmax <- NA
  }

  # -------- Print results
  cat("\n========================================\n")
  cat("Model:", model, " | k =", k_params, "\n")
  cat("TOST Results for beta_AUC (via beta_Cl)\n")
  cat("========================================\n")
  cat("Estimate:", beta_cl_est, "\n")
  cat("SE:", beta_cl_se, "\n\n")

  cat("Normal:  ", ifelse(tost_auc_normal, "Reject H0: BE", "Do not reject H0"), "\n")
  cat("  W_plus_delta:", W1_auc, " (crit:", quant, ")\n")
  cat("  W_minus_delta:", W2_auc, " (crit:", -quant, ")\n\n")

  cat("Student: ", ifelse(tost_auc_student, "Reject H0: BE", "Do not reject H0"), "\n")
  cat("  df =", n * nb_t - 5, ", tcrit =", quant_stud, "\n\n")

  cat("Gallant: ", ifelse(tost_auc_gallant, "Reject H0: BE", "Do not reject H0"), "\n")
  cat("  df =", n - k_params, ", tcrit =", quant_gallant, "\n")
  cat("  Corrected SE:", corrected_se_beta_cl, "\n\n")

  cat("========================================\n")
  cat("TOST Results for beta_Cmax\n")
  cat("========================================\n")
  cat("Estimate:", beta_cmax_estim, "\n")
  cat("SE:", secmax, "\n\n")

  cat("Normal:  ", ifelse(is.na(tost_cmax_normal), "ND", ifelse(tost_cmax_normal, "Reject H0: BE", "Do not reject H0")), "\n")
  cat("Student: ", ifelse(is.na(tost_cmax_student), "ND", ifelse(tost_cmax_student, "Reject H0: BE", "Do not reject H0")))
  if (!is.na(tost_cmax_student)) cat("  (df =", n * nb_t - 5, ", tcrit =", quant_stud, ")")
  cat("\n")
  cat("Gallant: ", ifelse(is.na(tost_cmax_gallant), "ND", ifelse(tost_cmax_gallant, "Reject H0: BE", "Do not reject H0")))
  if (!is.na(tost_cmax_gallant)) cat("  (df =", n - k_params, ", tcrit =", quant_gallant, ")")
  if (!is.na(cmax_expr)) cat("\n  Corrected SE:", corrected_se_beta_Cmax)
  cat("\n")

  invisible(list(
    fit = mod,
    n = n,
    nb_t = nb_t,
    k = k_params,
    varcov = varcov,
    beta_auc = list(
      estimate = beta_cl_est,
      se = beta_cl_se,
      normal = list(decision = tost_auc_normal, W1 = W1_auc, W2 = W2_auc, crit = quant),
      student = list(decision = tost_auc_student, df = n * nb_t - 5, crit = quant_stud),
      gallant = list(decision = tost_auc_gallant, df = n - k_params, crit = quant_gallant,
                     se_corrected = corrected_se_beta_cl)
    ),
    beta_cmax = list(
      estimate = beta_cmax_estim,
      se = secmax,
      normal  = list(decision = tost_cmax_normal,  W1 = if (exists("W1cmax")) W1cmax else NA,
                     W2 = if (exists("W2cmax")) W2cmax else NA, crit = quant),
      student = list(decision = tost_cmax_student, df = n * nb_t - 5, crit = quant_stud),
      gallant = list(decision = tost_cmax_gallant, df = n - k_params, crit = quant_gallant,
                     se_corrected = if (exists("corrected_se_beta_Cmax")) corrected_se_beta_Cmax else NA)
    )
  ))
}
