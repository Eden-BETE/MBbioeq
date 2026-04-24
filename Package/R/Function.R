#' MBbioeq: Model-Based Bioequivalence Analysis
#'
#' @description
#' MBbioeq performs **model-based bioequivalence (BE) analysis** for
#' parallel-arm studies using nonlinear mixed-effects pharmacokinetic (PK)
#' modelling via the SAEM algorithm (Stochastic Approximation
#' Expectation-Maximisation, as implemented in the **saemix** package).
#'
#' Two testing procedures are available:
#' \itemize{
#'   \item **TOST** — Two One-Sided Tests with asymptotic normal inference
#'     and a small-sample Gallant correction.
#'   \item **BOT** — Bayesian-Optimal Test based on the folded-normal
#'     distribution.
#' }
#'
#' Both AUC (via the clearance parameter) and Cmax (via the delta method)
#' are tested.  Only one-compartment PK models (with or without absorption
#' lag time) are supported.
#'
#' @section Main functions:
#' \describe{
#'   \item{`prep_parallel_TR()`}{Pre-process a raw dataset for SAEMIX.}
#'   \item{`run_MBTost()`}{Fit the NLMEM and run the TOST procedure.}
#'   \item{`run_MBbot()`}{Fit the NLMEM and run the BOT procedure.}
#'   \item{`plot_spaghetti()`}{Individual PK concentration-time profiles.}
#'   \item{`plot_tost()`}{Forest plot of GMR and 90 % CI from TOST.}
#'   \item{`plot_bot()`}{Folded-normal decision plot from BOT.}
#' }
#'
#' @keywords internal
"_PACKAGE"

# =============================================================================
# Global imports
# =============================================================================

#' @name imports
#' @import saemix
#' @import car
#' @import MASS
#' @import VGAM
#' @importFrom stats qnorm qt dnorm
#' @importFrom utils read.table
#' @importFrom graphics abline legend lines mtext par plot points polygon rect text
#' @importFrom grDevices adjustcolor
NULL

# =============================================================================
# PK model registry — 1-compartment models only
# =============================================================================

#' @keywords internal
.pk_registry <- list(

  "1cpt" = list(
    model = function(psi, id, xidep) {
      dose <- xidep[, 1]; tim <- xidep[, 2]
      ka   <- psi[id, 1]; V  <- psi[id, 2]; CL <- psi[id, 3]
      k    <- CL / V
      dose * ka / (V * (ka - k)) * (exp(-k * tim) - exp(-ka * tim))
    },
    par_names = c("ka", "V", "CL"),
    psi0      = c(1.5, 0.5, 0.04),
    transform = c(1, 1, 1)
  ),

  "1cpt_tlag" = list(
    model = function(psi, id, xidep) {
      dose <- xidep[, 1]; tim  <- xidep[, 2]
      ka   <- psi[id, 1]; V    <- psi[id, 2]
      CL   <- psi[id, 3]; tlag <- psi[id, 4]
      k     <- CL / V
      eff_t <- pmax(tim - tlag, 0)
      dose * ka / (V * (ka - k)) * (exp(-k * eff_t) - exp(-ka * eff_t)) *
        (tim >= tlag)
    },
    par_names = c("ka", "V", "CL", "tlag"),
    psi0      = c(1.5, 0.5, 0.04, 0.5),
    transform = c(1, 1, 1, 1)
  )
)

# =============================================================================
# Internal: auto-detect separator and read a PK dataset file
# =============================================================================

#' @keywords internal
.read_pk_data <- function(data_path) {
  hdr <- readLines(data_path, n = 1L, warn = FALSE)
  sep <- if      (grepl("\t", hdr, fixed = TRUE)) "\t"
         else if (grepl(";",  hdr, fixed = TRUE)) ";"
         else if (grepl(",",  hdr, fixed = TRUE)) ","
         else                                      " "
  utils::read.table(data_path, header = TRUE, sep = sep, dec = ".",
                    na.strings = c(".", "NA", ""))
}

# =============================================================================
# Internal: build a saemixModel object
# =============================================================================

#' @keywords internal
.build_saemix_model <- function(model,
                                psi0                 = NULL,
                                omega_init           = NULL,
                                covariance_structure = "diagonal",
                                error_model          = "combined",
                                error_init           = NULL,
                                verbose              = FALSE) {

  stopifnot(model %in% names(.pk_registry))
  spec <- .pk_registry[[model]]
  k    <- length(spec$par_names)

  # --- Fixed effects initial values ----------------------------------------
  if (is.null(psi0)) psi0 <- spec$psi0
  if (length(psi0) != k)
    stop("`psi0` must have ", k, " values for model '", model, "'.")
  psi0_mat <- matrix(psi0, ncol = k, byrow = TRUE,
                     dimnames = list(NULL, spec$par_names))

  # --- Random effects initial variance ------------------------------------
  if (is.null(omega_init)) omega_init <- rep(0.05, k)

  if (is.matrix(omega_init)) {
    if (!all(dim(omega_init) == c(k, k)))
      stop("`omega_init` must be a ", k, "x", k, " matrix.")
    omega_mat <- omega_init
    cov_model <- (omega_mat != 0) * 1L
  } else {
    if (length(omega_init) != k)
      stop("`omega_init` must have ", k, " values or be a ",
           k, "x", k, " matrix.")
    omega_mat <- diag(omega_init, k)
    if (is.character(covariance_structure)) {
      cov_model <- switch(
        covariance_structure,
        "diagonal" = diag(1L, k),
        "full"     = matrix(1L, k, k),
        stop("Unknown `covariance_structure`: ", covariance_structure,
             ". Use \"diagonal\" or \"full\".")
      )
    } else if (is.matrix(covariance_structure)) {
      cov_model <- covariance_structure
    } else {
      stop("`covariance_structure` must be \"diagonal\", \"full\", or a matrix.")
    }
  }

  # --- Residual error initial values --------------------------------------
  if (is.null(error_init)) {
    error_init <- switch(
      error_model,
      "combined"     = c(0.1, 0.1),
      "proportional" = 0.1,
      "constant"     = 0.1,
      stop("Unknown `error_model`: ", error_model,
           ". Use \"combined\", \"proportional\", or \"constant\".")
    )
  }

  expected_n <- if (error_model == "combined") 2L else 1L
  if (length(error_init) != expected_n)
    stop("`error_init` must have ", expected_n,
         " value(s) for error_model '", error_model, "'.")

  saemix::saemixModel(
    model            = spec$model,
    description      = paste0("PK ", model),
    psi0             = psi0_mat,
    transform.par    = spec$transform,
    covariate.model  = matrix(1L, nrow = 1L, ncol = k),
    covariance.model = cov_model,
    omega.init       = omega_mat,
    error.model      = error_model,
    error.init       = error_init,
    verbose          = verbose
  )
}

# =============================================================================
# Internal: index of fixed effects and analytical Cmax beta expression
# =============================================================================

#' @keywords internal
.get_fe_indices <- function(model) {
  par_names <- .pk_registry[[model]]$par_names
  k         <- length(par_names)
  fe_names  <- character(2L * k)
  for (i in seq_len(k)) {
    fe_names[2L * i - 1L] <- par_names[i]
    fe_names[2L * i]      <- paste0("beta_", par_names[i])
  }
  list(indices = seq_len(2L * k), names = fe_names)
}

#' @keywords internal
.beta_cmax_expr <- function(model) {
  # Analytical log-Cmax treatment effect for 1-compartment models,
  # derived via the delta method.  Expression variables match the SAEM
  # fixed-effects naming convention: population parameter + beta_parameter.
  expr <- paste0(
    "-beta_V",
    " - (log((ka*V)/CL) + beta_ka + beta_V - beta_CL)",
    " * CL*exp(beta_CL)",
    " / (ka*V*exp(beta_ka+beta_V) - CL*exp(beta_CL))",
    " + (CL/(ka*V - CL))*log(ka*V/CL)"
  )
  switch(model,
    "1cpt"      = expr,
    "1cpt_tlag" = expr,   # tlag does not appear in log-Cmax
    stop("Unsupported model '", model, "'.")
  )
}

# =============================================================================
# Data preparation
# =============================================================================

#' Prepare a dataset for SAEMIX (parallel-arm T/R design)
#'
#' @description
#' Reads and standardises a raw pharmacokinetic dataset for use with
#' [run_MBTost()] and [run_MBbot()].  The function accepts flexible column
#' names (see Details), validates treatment coding, forward-fills missing
#' dose values within each subject, and returns a tidy data frame with
#' canonical column names.
#'
#' @param tab A `data.frame` read from the study dataset.
#'
#' @details
#' **Required columns** (case-insensitive; accepted synonyms in parentheses):
#' \describe{
#'   \item{`Id`}{Subject identifier (`subject`, `subj`, `sid`).}
#'   \item{`Time`}{Sampling time in hours (`t`, `hours`, `hr`).}
#'   \item{`Concentration`}{Observed drug concentration (`conc`, `dv`, `y`,
#'     `obs`).}
#'   \item{`Treatment`}{Group indicator: `"T"` for Test, `"R"` for Reference
#'     (`treat`, `trt`, `tr`, `group`, `arm`).}
#' }
#'
#' **Optional column:**
#' \describe{
#'   \item{`Dose`}{Administered dose (`amt`).  If absent or partially
#'     missing, the last observed dose is carried forward within each subject.
#'     At least one dose value per subject must be present.}
#' }
#'
#' Missing concentrations should be coded as `NA`, `"."`, or `""` in the
#' source file; pre-dose rows (time 0) with a missing concentration are
#' retained and filtered inside the fitting functions.
#'
#' @return
#' A `data.frame` with columns `Id` (integer/character), `Dose` (numeric),
#' `Time` (numeric), `Concentration` (numeric), and `Treat` (integer: 1 for
#' Test, 0 for Reference), sorted by `Id` then `Time`.
#'
#' @seealso [run_MBTost()], [run_MBbot()], [plot_spaghetti()]
#'
#' @export
#' @examples
#' \dontrun{
#' raw <- read.table("dataset_1.txt", header = TRUE)
#' tab <- prep_parallel_TR(raw)
#' head(tab)
#' }
prep_parallel_TR <- function(tab) {
  stopifnot(is.data.frame(tab))
  names(tab) <- tolower(names(tab))

  col_map <- list(
    id            = c("id", "subject", "subj", "sid"),
    time          = c("time", "t", "hours", "hr"),
    concentration = c("concentration", "conc", "dv", "y", "obs"),
    treatment     = c("treatment", "treat", "trt", "tr", "group", "arm"),
    dose          = c("dose", "amt")
  )

  rename_one <- function(std) {
    vars <- col_map[[std]]
    hits <- intersect(vars, names(tab))
    if (length(hits) == 0 && std %in% c("id", "time", "concentration", "treatment"))
      stop("Missing required column: '", std,
           "'. Accepted synonyms: ", paste(vars, collapse = ", "), ".")
    if (length(hits) > 1)
      stop("Ambiguous columns for '", std, "': ", paste(hits, collapse = ", "), ".")
    if (length(hits) == 1)
      names(tab)[names(tab) == hits] <<- std
  }
  lapply(names(col_map), rename_one)

  miss <- setdiff(c("id", "time", "concentration", "treatment"), names(tab))
  if (length(miss) > 0)
    stop("Missing required columns: ", paste(miss, collapse = ", "), ".")

  tr_raw <- toupper(as.character(tab$treatment))
  ok     <- tr_raw %in% c("T", "R")
  if (!all(ok))
    stop("Invalid treatment values: ",
         paste(unique(tr_raw[!ok]), collapse = ", "),
         ". Only \"T\" and \"R\" are accepted.")
  tab$Treat <- ifelse(tr_raw == "T", 1L, 0L)

  if (!"dose" %in% names(tab)) tab$dose <- NA_real_
  tab <- tab[order(tab$id, tab$time), ]

  # Forward-fill dose within each subject
  ff <- function(x) {
    idx <- cummax(ifelse(!is.na(x), seq_along(x), 0L))
    ifelse(idx == 0L, NA_real_, x[idx])
  }
  tab$Dose <- ave(tab$dose, tab$id, FUN = ff)

  bad_ids <- unique(tab$id[is.na(tab$Dose)])
  if (length(bad_ids) > 0)
    stop("Could not impute Dose for subject(s): ",
         paste(bad_ids, collapse = ", "),
         ". Ensure at least one non-missing dose per subject.")

  res <- data.frame(
    Id            = tab$id,
    Dose          = tab$Dose,
    Time          = tab$time,
    Concentration = tab$concentration,
    Treat         = tab$Treat
  )
  res[order(res$Id, res$Time), ]
}

# =============================================================================
# Shared fitting core
# =============================================================================

#' @keywords internal
.fit_MB_model <- function(data_path, model, psi0, omega_init,
                          covariance_structure, error_model, error_init,
                          nb_chains, nb_iter, verbose) {

  tab <- .read_pk_data(data_path)
  tab <- prep_parallel_TR(tab)
  tab <- tab[!is.na(tab$Concentration), ]

  n <- length(unique(tab$Id))    # number of subjects

  if (!verbose) {
    invisible(capture.output({
      saemix.data <- saemix::saemixData(
        name.data        = tab,    header = TRUE, sep = " ", na = NA,
        name.group       = "Id",
        name.predictors  = c("Dose", "Time"),
        name.covariates  = "Treat",
        name.response    = "Concentration",
        name.X           = "Time",
        units            = list(x = "hr", y = "mg/L")
      )
      saemix.model <- .build_saemix_model(
        model, psi0, omega_init, covariance_structure,
        error_model, error_init, verbose = FALSE
      )
    }))
  } else {
    saemix.data <- saemix::saemixData(
      name.data        = tab,    header = TRUE, sep = " ", na = NA,
      name.group       = "Id",
      name.predictors  = c("Dose", "Time"),
      name.covariates  = "Treat",
      name.response    = "Concentration",
      name.X           = "Time",
      units            = list(x = "hr", y = "mg/L")
    )
    saemix.model <- .build_saemix_model(
      model, psi0, omega_init, covariance_structure, error_model, error_init
    )
  }

  k_params <- length(.pk_registry[[model]]$par_names)

  saem_opts <- list(
    nb.chains        = nb_chains,
    nbiter.saemix    = nb_iter,
    print            = FALSE,
    displayProgress  = FALSE,
    save             = FALSE,
    save.graphs      = FALSE
  )

  if (!verbose) {
    invisible(capture.output({
      mod <- suppressMessages(suppressWarnings(
        saemix::saemix(saemix.model, saemix.data, saem_opts)
      ))
    }, type = "output"))
  } else {
    mod <- saemix::saemix(saemix.model, saemix.data, saem_opts)
  }

  # --- Extract fixed effects and variance-covariance ----------------------
  fe     <- mod@results@fixed.effects
  sefe   <- mod@results@se.fixed
  FIM    <- mod@results@fim
  varcov <- tryCatch(solve(FIM), error = function(e) MASS::ginv(FIM))

  # --- AUC effect (= -beta_CL) -------------------------------------------
  par_names  <- .pk_registry[[model]]$par_names
  idx_cl     <- which(toupper(par_names) == "CL")
  if (length(idx_cl) != 1L)
    stop("Could not locate 'CL' in the parameter list.")
  beta_cl_pos  <- 2L * idx_cl
  beta_auc_est <- -fe[beta_cl_pos]
  beta_auc_se  <- sefe[beta_cl_pos]

  # --- Cmax effect via delta method ---------------------------------------
  cmax_expr    <- .beta_cmax_expr(model)
  fe_info      <- .get_fe_indices(model)
  estimates    <- fe[fe_info$indices]
  names(estimates) <- fe_info$names

  var_subset <- varcov[fe_info$indices, fe_info$indices, drop = FALSE]
  colnames(var_subset) <- rownames(var_subset) <- fe_info$names

  dm <- car::deltaMethod(estimates, cmax_expr, vcov. = var_subset)
  beta_cmax_est <- dm[, 1L]
  beta_cmax_se  <- dm[, 2L]

  list(
    fit          = mod,
    model        = model,
    n            = n,
    k            = k_params,
    varcov       = varcov,
    beta_auc_est = beta_auc_est,
    beta_auc_se  = beta_auc_se,
    beta_cmax_est = beta_cmax_est,
    beta_cmax_se  = beta_cmax_se
  )
}

# =============================================================================
# TOST
# =============================================================================

#' Run model-based TOST bioequivalence analysis
#'
#' @description
#' Fits a one-compartment nonlinear mixed-effects pharmacokinetic model via
#' the SAEM algorithm and tests bioequivalence for AUC and Cmax using the
#' **Two One-Sided Tests (TOST)** procedure.
#'
#' @section Statistical method:
#' The population PK model is fitted with **saemix** using the SAEM algorithm.
#' Treatment (Test vs. Reference) enters as a covariate on all PK parameters
#' on the log scale, so the treatment log-ratio for AUC is
#' \eqn{\hat\beta_{AUC} = -\hat\beta_{CL}} and for Cmax it is obtained via
#' the delta method on the analytical log-Cmax expression.
#'
#' TOST uses a 90 % confidence interval for the T/R ratio on the original
#' scale:
#' \deqn{CI = \exp\!\bigl(\hat\beta \pm z_{1-\alpha}\,\hat\sigma\bigr)}
#' Bioequivalence is concluded if the CI lies entirely within
#' \eqn{[e^{-\delta},\, e^{\delta}]} (default \eqn{\delta = \log(1.25)},
#' i.e. the 80–125 \% range).
#'
#' Two inference variants are provided:
#' \describe{
#'   \item{**Normal**}{Standard asymptotic normal quantile \eqn{z_{1-\alpha}}.}
#'   \item{**Gallant**}{Small-sample correction replacing \eqn{z_{1-\alpha}}
#'     by the Student-\eqn{t} quantile on \eqn{n-k} degrees of freedom and
#'     inflating the standard error accordingly.  Requires \eqn{n > k};
#'     results are `NA` otherwise with a warning.}
#' }
#'
#' @param data_path Character. Path to the dataset file.  The separator
#'   (space, tab, semicolon, or comma) is detected automatically.
#' @param model Character. PK model: `"1cpt"` (one-compartment) or
#'   `"1cpt_tlag"` (one-compartment with absorption lag time).
#' @param psi0 Named numeric vector of initial fixed-effects values
#'   (`ka`, `V`, `CL`, and optionally `tlag`).  `NULL` uses built-in
#'   defaults.
#' @param omega_init Initial random-effects standard deviations (numeric
#'   vector of length \eqn{k}) or a full \eqn{k \times k} variance matrix.
#'   `NULL` uses `rep(0.05, k)`.
#' @param covariance_structure `"diagonal"` (default), `"full"`, or a
#'   custom \eqn{k \times k} binary matrix indicating which random-effect
#'   covariances are estimated.
#' @param error_model Residual error model: `"combined"` (default),
#'   `"proportional"`, or `"constant"`.
#' @param error_init Initial residual error parameter(s).  `NULL` uses
#'   `c(0.1, 0.1)` for `"combined"` and `0.1` otherwise.
#' @param nb_chains Integer. Number of SAEM Markov chains (default `10`).
#' @param nb_iter Integer vector `c(n_saem, n_mcmc)`. Number of SAEM
#'   iterations followed by MCMC iterations (default `c(300, 100)`).
#' @param delta Numeric. Equivalence margin on the log scale
#'   (default `log(1.25)`, corresponding to the 80–125 \% ratio range).
#' @param alpha Numeric. One-sided significance level for the TOST
#'   (default `0.05`, yielding a 90 \% CI).
#' @param verbose Logical. Print SAEMIX convergence output (default `FALSE`).
#'
#' @return
#' An invisible named list with the following elements:
#' * `fit`: The fitted `SaemixObject` from **saemix**.
#' * `model`: Name of the fitted PK model (character).
#' * `n`: Number of subjects.
#' * `k`: Number of PK parameters.
#' * `varcov`: Variance-covariance matrix of fixed effects.
#' * `auc`: List with elements `estimate` (GMR on ratio scale),
#'   `se` (SE of log-ratio),
#'   `normal` (list: `ci_lower`, `ci_upper`, `decision`),
#'   `gallant` (list: `ci_lower`, `ci_upper`, `decision`).
#' * `cmax`: Same structure as `auc` for the Cmax endpoint.
#'
#' @seealso [plot_tost()] to visualise the result; [run_MBbot()] for the
#'   BOT alternative.
#'
#' @export
#' @examples
#' \dontrun{
#' res <- run_MBTost(
#'   data_path = system.file("extdata", "dataset_1.txt", package = "MBbioeq"),
#'   model     = "1cpt",
#'   psi0      = c(ka = 1.5, V = 0.5, CL = 0.04),
#'   error_model = "combined"
#' )
#' plot_tost(res)
#' }
run_MBTost <- function(data_path,
                       model                = "1cpt",
                       psi0                 = NULL,
                       omega_init           = NULL,
                       covariance_structure = "diagonal",
                       error_model          = "combined",
                       error_init           = NULL,
                       nb_chains            = 10L,
                       nb_iter              = c(300L, 100L),
                       delta                = log(1.25),
                       alpha                = 0.05,
                       verbose              = FALSE) {

  res   <- .fit_MB_model(data_path, model, psi0, omega_init,
                         covariance_structure, error_model, error_init,
                         nb_chains, nb_iter, verbose)
  quant <- stats::qnorm(1 - alpha)

  # --- Gallant correction -------------------------------------------------
  df_gallant <- res$n - res$k
  if (df_gallant > 0L) {
    quant_gallant <- stats::qt(1 - alpha, df = df_gallant)
    if (!is.finite(quant_gallant) || quant_gallant <= 0) {
      warning("Gallant t-quantile is non-finite or non-positive (df = ",
              df_gallant, "). Gallant results set to NA.")
      quant_gallant <- NA_real_
    }
  } else {
    warning("Gallant correction requires n > k (n = ", res$n,
            ", k = ", res$k, "). Gallant results set to NA.")
    quant_gallant <- NA_real_
  }
  gallant_ok <- !is.na(quant_gallant)

  # ========== TOST — AUC ==================================================
  ic_auc_lo <- res$beta_auc_est - quant * res$beta_auc_se
  ic_auc_hi <- res$beta_auc_est + quant * res$beta_auc_se

  W1_auc <- (res$beta_auc_est + delta) / res$beta_auc_se
  W2_auc <- (res$beta_auc_est - delta) / res$beta_auc_se
  tost_auc_normal <- (W1_auc >= quant) & (W2_auc <= -quant)

  if (gallant_ok) {
    se_auc_g  <- res$beta_auc_se * sqrt(res$n / quant_gallant)
    ic_auc_g_lo <- res$beta_auc_est - quant_gallant * se_auc_g
    ic_auc_g_hi <- res$beta_auc_est + quant_gallant * se_auc_g
    W1_g_auc <- (res$beta_auc_est + delta) / se_auc_g
    W2_g_auc <- (res$beta_auc_est - delta) / se_auc_g
    tost_auc_gallant <- (W1_g_auc >= quant_gallant) & (W2_g_auc <= -quant_gallant)
  } else {
    ic_auc_g_lo  <- NA_real_; ic_auc_g_hi   <- NA_real_
    tost_auc_gallant <- NA
  }

  # ========== TOST — Cmax =================================================
  ic_cmax_lo <- res$beta_cmax_est - quant * res$beta_cmax_se
  ic_cmax_hi <- res$beta_cmax_est + quant * res$beta_cmax_se

  W1_cmax <- (res$beta_cmax_est + delta) / res$beta_cmax_se
  W2_cmax <- (res$beta_cmax_est - delta) / res$beta_cmax_se
  tost_cmax_normal <- (W1_cmax >= quant) & (W2_cmax <= -quant)

  if (gallant_ok) {
    se_cmax_g   <- res$beta_cmax_se * sqrt(res$n / quant_gallant)
    ic_cmax_g_lo <- res$beta_cmax_est - quant_gallant * se_cmax_g
    ic_cmax_g_hi <- res$beta_cmax_est + quant_gallant * se_cmax_g
    W1_g_cmax <- (res$beta_cmax_est + delta) / se_cmax_g
    W2_g_cmax <- (res$beta_cmax_est - delta) / se_cmax_g
    tost_cmax_gallant <- (W1_g_cmax >= quant_gallant) & (W2_g_cmax <= -quant_gallant)
  } else {
    ic_cmax_g_lo  <- NA_real_; ic_cmax_g_hi   <- NA_real_
    tost_cmax_gallant <- NA
  }

  # ========== Console output ==============================================
  .be_label <- function(dec) ifelse(isTRUE(dec), "Bioequivalent", "Not bioequivalent")

  cat("\n")
  cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n")
  cat("              TOST BIOEQUIVALENCE RESULTS\n")
  cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n")
  cat("  Model:", model, " | n =", res$n, " | k =", res$k,
      " | \u03b4 =", round(delta, 4L), "\n")
  cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n\n")

  cat("AUC RATIO (T/R):\n")
  cat("  Estimate        : ", sprintf("%.4f", exp(res$beta_auc_est)), "\n", sep = "")
  cat("  90% CI (Normal) : [", sprintf("%.4f", exp(ic_auc_lo)), ", ",
      sprintf("%.4f", exp(ic_auc_hi)), "]  \u2192  ",
      .be_label(tost_auc_normal), "\n", sep = "")
  cat("  90% CI (Gallant): [",
      if (gallant_ok) sprintf("%.4f", exp(ic_auc_g_lo)) else "  NA  ",
      ", ",
      if (gallant_ok) sprintf("%.4f", exp(ic_auc_g_hi)) else "  NA  ",
      "]  \u2192  ", .be_label(tost_auc_gallant), "\n", sep = "")

  cat("\nCmax RATIO (T/R):\n")
  cat("  Estimate        : ", sprintf("%.4f", exp(res$beta_cmax_est)), "\n", sep = "")
  cat("  90% CI (Normal) : [", sprintf("%.4f", exp(ic_cmax_lo)), ", ",
      sprintf("%.4f", exp(ic_cmax_hi)), "]  \u2192  ",
      .be_label(tost_cmax_normal), "\n", sep = "")
  cat("  90% CI (Gallant): [",
      if (gallant_ok) sprintf("%.4f", exp(ic_cmax_g_lo)) else "  NA  ",
      ", ",
      if (gallant_ok) sprintf("%.4f", exp(ic_cmax_g_hi)) else "  NA  ",
      "]  \u2192  ", .be_label(tost_cmax_gallant), "\n", sep = "")

  cat("\n\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n\n")

  # ========== Return value ================================================
  invisible(list(
    fit    = res$fit,
    model  = model,
    n      = res$n,
    k      = res$k,
    varcov = res$varcov,

    auc = list(
      estimate = exp(res$beta_auc_est),
      se       = res$beta_auc_se,
      normal  = list(
        ci_lower = exp(ic_auc_lo),
        ci_upper = exp(ic_auc_hi),
        decision = tost_auc_normal
      ),
      gallant = list(
        ci_lower = exp(ic_auc_g_lo),
        ci_upper = exp(ic_auc_g_hi),
        decision = tost_auc_gallant
      )
    ),

    cmax = list(
      estimate = exp(res$beta_cmax_est),
      se       = res$beta_cmax_se,
      normal  = list(
        ci_lower = exp(ic_cmax_lo),
        ci_upper = exp(ic_cmax_hi),
        decision = tost_cmax_normal
      ),
      gallant = list(
        ci_lower = exp(ic_cmax_g_lo),
        ci_upper = exp(ic_cmax_g_hi),
        decision = tost_cmax_gallant
      )
    )
  ))
}

# =============================================================================
# BOT
# =============================================================================

#' Run model-based BOT bioequivalence analysis
#'
#' @description
#' Fits a one-compartment nonlinear mixed-effects pharmacokinetic model via
#' the SAEM algorithm and tests bioequivalence for AUC and Cmax using the
#' **Bayesian-Optimal Test (BOT)** procedure.
#'
#' @section Statistical method:
#' The population PK model and treatment log-ratios are estimated exactly as
#' in [run_MBTost()].  The BOT test statistic is \eqn{T = |\hat\beta|}.
#' Under the worst-case null hypothesis \eqn{|\beta| = \delta}, \eqn{T}
#' follows a **folded-normal distribution** with location \eqn{\delta} and
#' scale \eqn{\hat\sigma}:
#' \deqn{T \sim FN(\delta,\, \hat\sigma^2)}
#'
#' Bioequivalence is concluded when
#' \deqn{T \;\le\; c^* = q_{FN}(\alpha;\, \delta,\, \hat\sigma),}
#' where \eqn{q_{FN}(\alpha;\cdot)} is the \eqn{\alpha}-quantile of that
#' folded-normal distribution (computed via [VGAM::qfoldnorm()]).  The BOT is
#' uniformly more powerful than TOST in many settings while controlling the
#' type-I error at level \eqn{\alpha}.
#'
#' @inheritParams run_MBTost
#'
#' @return
#' An invisible named list with the following elements:
#' * `fit`: The fitted `SaemixObject` from **saemix**.
#' * `model`: Name of the fitted PK model.
#' * `n`: Number of subjects.
#' * `k`: Number of PK parameters.
#' * `varcov`: Variance-covariance matrix of fixed effects.
#' * `auc`: List with elements `estimate` (GMR), `se` (SE of log-ratio),
#'   `statistic` (\eqn{|\hat\beta_{AUC}|}),
#'   `critical_value` (\eqn{c^*}), `decision` (logical).
#' * `cmax`: Same structure as `auc` for the Cmax endpoint.
#'
#' @seealso [plot_bot()] to visualise the decision; [run_MBTost()] for the
#'   TOST alternative.
#'
#' @export
#' @examples
#' \dontrun{
#' res <- run_MBbot(
#'   data_path = system.file("extdata", "dataset_1.txt", package = "MBbioeq"),
#'   model     = "1cpt",
#'   psi0      = c(ka = 1.5, V = 0.5, CL = 0.04),
#'   error_model = "combined"
#' )
#' plot_bot(res)
#' }
run_MBbot <- function(data_path,
                      model                = "1cpt",
                      psi0                 = NULL,
                      omega_init           = NULL,
                      covariance_structure = "diagonal",
                      error_model          = "combined",
                      error_init           = NULL,
                      nb_chains            = 10L,
                      nb_iter              = c(300L, 100L),
                      delta                = log(1.25),
                      alpha                = 0.05,
                      verbose              = FALSE) {

  res <- .fit_MB_model(data_path, model, psi0, omega_init,
                       covariance_structure, error_model, error_init,
                       nb_chains, nb_iter, verbose)

  # ========== BOT — AUC ===================================================
  stat_auc  <- abs(res$beta_auc_est)
  limit_auc <- VGAM::qfoldnorm(alpha, delta, res$beta_auc_se)
  bot_auc   <- stat_auc <= limit_auc

  # ========== BOT — Cmax ==================================================
  stat_cmax  <- abs(res$beta_cmax_est)
  limit_cmax <- VGAM::qfoldnorm(alpha, delta, res$beta_cmax_se)
  bot_cmax   <- stat_cmax <= limit_cmax

  # ========== Console output ==============================================
  .be_label <- function(dec) ifelse(isTRUE(dec), "Bioequivalent", "Not bioequivalent")

  cat("\n")
  cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n")
  cat("               BOT BIOEQUIVALENCE RESULTS\n")
  cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n")
  cat("  Model:", model, " | n =", res$n, " | k =", res$k,
      " | \u03b4 =", round(delta, 4L), "\n")
  cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n\n")

  cat("AUC RATIO (T/R):\n")
  cat("  Estimate       : ", sprintf("%.4f", exp(res$beta_auc_est)), "\n", sep = "")
  cat("  |beta_hat|     : ", sprintf("%.4f", stat_auc),  "\n", sep = "")
  cat("  Critical value : ", sprintf("%.4f", limit_auc), "\n", sep = "")
  cat("  Decision       : ", .be_label(bot_auc), "\n", sep = "")

  cat("\nCmax RATIO (T/R):\n")
  cat("  Estimate       : ", sprintf("%.4f", exp(res$beta_cmax_est)), "\n", sep = "")
  cat("  |beta_hat|     : ", sprintf("%.4f", stat_cmax),  "\n", sep = "")
  cat("  Critical value : ", sprintf("%.4f", limit_cmax), "\n", sep = "")
  cat("  Decision       : ", .be_label(bot_cmax), "\n", sep = "")

  cat("\n\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n\n")

  # ========== Return value ================================================
  invisible(list(
    fit    = res$fit,
    model  = model,
    n      = res$n,
    k      = res$k,
    varcov = res$varcov,

    auc = list(
      estimate       = exp(res$beta_auc_est),
      se             = res$beta_auc_se,
      statistic      = stat_auc,
      critical_value = limit_auc,
      decision       = bot_auc
    ),

    cmax = list(
      estimate       = exp(res$beta_cmax_est),
      se             = res$beta_cmax_se,
      statistic      = stat_cmax,
      critical_value = limit_cmax,
      decision       = bot_cmax
    )
  ))
}
